/*
solver_cg.c

対角スケーリング前処理付き共役勾配法 (Jacobi-PCG)。

Dirichlet 節点 (fix[] != 0) は恒等行として扱い、右辺を 0 にしておくことで
解が自動的に 0 になる (未知量は自由節点のみ)。剛性行列そのものは書き換えない
ので、解いた後に反作用 (電荷・電流) を元の行列から計算できる。

戻り値 : 反復回数 (収束しなかった場合は負)
*/

#include "fem.h"
#include "fem_prototype.h"

static double dotprod(const double *a, const double *b, int n)
{
	double s = 0;
	int i;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:s)
#endif
	for (i = 0; i < n; i++) {
		s += a[i] * b[i];
	}

	return s;
}


int solver_cg(const crs_t *A, const double *b, double *x, const unsigned char *fix,
	int maxiter, int nout, double converg, FILE *fp_log, const char *label)
{
	const int n = (int)A->n;

	double *r = (double *)malloc(n * sizeof(double));
	double *p = (double *)malloc(n * sizeof(double));
	double *q = (double *)malloc(n * sizeof(double));
	double *z = (double *)malloc(n * sizeof(double));
	double *d = (double *)malloc(n * sizeof(double));

	// 前処理 (対角)
	crs_diag(A, d);
	int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (i = 0; i < n; i++) {
		if (fix[i] || (d[i] <= 0)) d[i] = 1;
	}

	// 初期値 x = 0 -> r = b
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (i = 0; i < n; i++) {
		x[i] = 0;
		r[i] = b[i];
	}

	const double bnorm = sqrt(dotprod(b, b, n));
	if (bnorm <= 0) {
		free(r); free(p); free(q); free(z); free(d);
		return 0;
	}

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (i = 0; i < n; i++) {
		z[i] = r[i] / d[i];
		p[i] = z[i];
	}
	double rz = dotprod(r, z, n);

	int iter = 0;
	double resid = 1;
	int converged = 0;

	for (iter = 1; iter <= maxiter; iter++) {
		crs_spmv(A, p, q, fix);

		const double pq = dotprod(p, q, n);
		if (pq <= 0) break;			// 正定値でない (異常)
		const double alpha = rz / pq;

		int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (i = 0; i < n; i++) {
			x[i] += alpha * p[i];
			r[i] -= alpha * q[i];
		}

		resid = sqrt(dotprod(r, r, n)) / bnorm;
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
			z[i] = r[i] / d[i];
		}
		const double rznew = dotprod(r, z, n);
		const double beta = rznew / rz;
		rz = rznew;

#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (i = 0; i < n; i++) {
			p[i] = z[i] + (beta * p[i]);
		}
	}

	if (fp_log != NULL) {
		fprintf(fp_log, "  %-10s %8d %13.5e %s\n", label, iter, resid,
			(converged ? "converged" : "*** NOT converged"));
		fflush(fp_log);
	}

	free(r); free(p); free(q); free(z); free(d);

	return (converged ? iter : -iter);
}
