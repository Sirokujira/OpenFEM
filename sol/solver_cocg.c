/*
solver_cocg.c

複素対称行列 A = K + jω M に対する COCG 法
(Conjugate Orthogonal Conjugate Gradient、対角前処理付き)。

時間調和渦電流問題 ∫ν∇w・∇Az + jω∫σ w Az = ∫σ w V' の係数行列は
複素**対称** (Hermite ではない) なので、内積に共役を取らない双一次形式
(x, y) = Σ x_i y_i を使う COCG が使える。

複素数は実部・虚部の配列を別に持つ (MSVC の C99 complex 依存を避ける)。
K と M は同じ CRS パターンを共有する。

戻り値 : 反復回数 (収束しなかった場合は負)
*/

#include "fem.h"
#include "fem_prototype.h"

// 複素数の積・商 (スカラー)
static void cmul(double ar, double ai, double br, double bi, double *cr, double *ci)
{
	*cr = (ar * br) - (ai * bi);
	*ci = (ar * bi) + (ai * br);
}

static void cdiv(double ar, double ai, double br, double bi, double *cr, double *ci)
{
	const double d = (br * br) + (bi * bi);
	*cr = ((ar * br) + (ai * bi)) / d;
	*ci = ((ai * br) - (ar * bi)) / d;
}


// 双一次形式 (共役を取らない) の内積
static void cdot(const double *ar, const double *ai, const double *br, const double *bi,
	int n, double *sr, double *si)
{
	double xr = 0, xi = 0;
	int i;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:xr) reduction(+:xi)
#endif
	for (i = 0; i < n; i++) {
		xr += (ar[i] * br[i]) - (ai[i] * bi[i]);
		xi += (ar[i] * bi[i]) + (ai[i] * br[i]);
	}
	*sr = xr;
	*si = xi;
}


// 2 ノルム
static double cnorm(const double *ar, const double *ai, int n)
{
	double s = 0;
	int i;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:s)
#endif
	for (i = 0; i < n; i++) {
		s += (ar[i] * ar[i]) + (ai[i] * ai[i]);
	}

	return sqrt(s);
}


int solver_cocg(const crs_t *K, const crs_t *M, double omega,
	const double *br, const double *bi, double *xr, double *xi,
	const unsigned char *fix, int maxiter, int nout, double converg,
	FILE *fp_log, const char *label)
{
	const int n = (int)K->n;

	// direct = 1 のときは直接解法 (スカイライン LDL^T) に回す。
	// 反復回数の代わりに 0 を返す (実対称系の solver_cg と同じ約束)
	if (Direct) {
		return ((solver_direct_c(K, M, omega, br, bi, xr, xi, fix, fp_log, label) == 0)
			? 0 : -1);
	}

	double *rr = (double *)malloc(n * sizeof(double));
	double *ri = (double *)malloc(n * sizeof(double));
	double *pr = (double *)malloc(n * sizeof(double));
	double *pi = (double *)malloc(n * sizeof(double));
	double *qr = (double *)malloc(n * sizeof(double));
	double *qi = (double *)malloc(n * sizeof(double));
	double *zr = (double *)malloc(n * sizeof(double));
	double *zi = (double *)malloc(n * sizeof(double));
	double *dr = (double *)malloc(n * sizeof(double));
	double *di = (double *)malloc(n * sizeof(double));
	double *w1 = (double *)malloc(n * sizeof(double));
	double *w2 = (double *)malloc(n * sizeof(double));
	double *w3 = (double *)malloc(n * sizeof(double));
	double *w4 = (double *)malloc(n * sizeof(double));

	// 対角前処理 (K + jωM の対角)
	crs_diag(K, dr);
	crs_diag(M, di);
	int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (i = 0; i < n; i++) {
		if (fix[i]) {
			dr[i] = 1;
			di[i] = 0;
		}
		else {
			di[i] *= omega;
			if ((dr[i] == 0) && (di[i] == 0)) dr[i] = 1;
		}
	}

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (i = 0; i < n; i++) {
		xr[i] = xi[i] = 0;
		rr[i] = br[i];
		ri[i] = bi[i];
	}

	const double bnorm = cnorm(br, bi, n);
	if (bnorm <= 0) {
		free(rr); free(ri); free(pr); free(pi); free(qr); free(qi);
		free(zr); free(zi); free(dr); free(di);
		free(w1); free(w2); free(w3); free(w4);
		return 0;
	}

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (i = 0; i < n; i++) {
		cdiv(rr[i], ri[i], dr[i], di[i], &zr[i], &zi[i]);
		pr[i] = zr[i];
		pi[i] = zi[i];
	}
	double rzr, rzi;
	cdot(rr, ri, zr, zi, n, &rzr, &rzi);

	int iter = 0;
	double resid = 1;
	int converged = 0;

	for (iter = 1; iter <= maxiter; iter++) {
		crs_spmv_c(K, M, omega, pr, pi, qr, qi, fix, w1, w2, w3, w4);

		double pqr, pqi;
		cdot(pr, pi, qr, qi, n, &pqr, &pqi);
		if ((pqr == 0) && (pqi == 0)) break;
		double ar, ai;
		cdiv(rzr, rzi, pqr, pqi, &ar, &ai);

#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (i = 0; i < n; i++) {
			double tr, ti;
			cmul(ar, ai, pr[i], pi[i], &tr, &ti);
			xr[i] += tr;
			xi[i] += ti;
			cmul(ar, ai, qr[i], qi[i], &tr, &ti);
			rr[i] -= tr;
			ri[i] -= ti;
		}

		resid = cnorm(rr, ri, n) / bnorm;
		if ((fp_log != NULL) && ((iter % nout) == 0)) {
			fprintf(fp_log, "  %-10s %8d %13.5e\n", label, iter, resid);
			fflush(fp_log);
		}
		if (resid < converg) {
			converged = 1;
			break;
		}

#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (i = 0; i < n; i++) {
			cdiv(rr[i], ri[i], dr[i], di[i], &zr[i], &zi[i]);
		}
		double nr, ni;
		cdot(rr, ri, zr, zi, n, &nr, &ni);
		double betar, betai;
		cdiv(nr, ni, rzr, rzi, &betar, &betai);
		rzr = nr;
		rzi = ni;

#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (i = 0; i < n; i++) {
			double tr, ti;
			cmul(betar, betai, pr[i], pi[i], &tr, &ti);
			pr[i] = zr[i] + tr;
			pi[i] = zi[i] + ti;
		}
	}

	if (fp_log != NULL) {
		fprintf(fp_log, "  %-10s %8d %13.5e %s\n", label, iter, resid,
			(converged ? "converged" : "*** NOT converged"));
		fflush(fp_log);
	}

	free(rr); free(ri); free(pr); free(pi); free(qr); free(qi);
	free(zr); free(zi); free(dr); free(di);
	free(w1); free(w2); free(w3); free(w4);

	return (converged ? iter : -iter);
}
