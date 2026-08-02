/*
solver_direct.c

直接解法 (RCM 並べ替え + スカイライン分解)。

反復解法しか無いと、悪条件な問題 (高コントラストな材料、扁平な要素) で
収束しないときに逃げ道が無い。`direct = 1` でこちらに切り替える。

2 つの系に対応する:
  solver_direct   : 実対称正定値 (C / L / R / M と、非線形反復の内側) -> Cholesky
  solver_direct_c : 複素対称 K + jω M (渦電流 F / A)                  -> LDL^T

複素対称系に Cholesky は使えない (√ が複素になり、そもそも正定値でない) ので
**平方根の要らない LDL^T** で分解する。ピボット選択はしない: 選択すると
対称性が崩れてスカイラインの外に充填が出る。F の系は K が半正定値、
ω M も半正定値で、両者の零空間が交わらないため主小行列式が全て非零になり、
ピボット無しの LDL^T が破綻しないことが示せる。A (A-φ 連成) はそこまで
言えないので、**分解後に真の残差 ‖b − A x‖/‖b‖ を計算してログに出す**
(COCG の残差と同じ意味の数字が並ぶので、破綻すればそこで分かる)。

Dirichlet の扱いは反復解法と完全に同じにする:
固定節点は恒等行 (x_i = b_i = 0) なので、行も列も落として単位対角を入れる。
列を落としてよいのは固定節点の解が 0 だから (右辺は持ち上げ済み)。
**行列そのものは書き換えない** (解いた後に反作用を元の行列から取るため)。

なぜスカイラインか: 並べ替え無しの疎 Cholesky は 3 次元格子で充填が爆発する。
RCM は帯幅 (プロファイル) を減らす古典的な並べ替えで、実装が短く検算しやすい。
AMD + 一般疎分解の方が速いが、コード量と検証コストに見合わない。
プロファイルが大きすぎるときは黙って落ちずに**見積もりを出して断る**。
*/

#include "fem.h"
#include "fem_prototype.h"

// プロファイルの上限 (バイト)。超えたら反復解法を使うよう促す。
// 複素では 1 要素が 2 倍になるので、要素数ではなくバイトで比べる
#define MAXPROFBYTE (320LL * 1024 * 1024)


// 逆 Cuthill-McKee 並べ替え。perm[k] = k 番目に並ぶ元の節点番号
static void rcm_order(const crs_t *A, const unsigned char *fix, int32_t *perm)
{
	const int n = (int)A->n;
	int *deg = (int *)malloc((size_t)n * sizeof(int));
	unsigned char *seen = (unsigned char *)malloc((size_t)n);
	memset(seen, 0, (size_t)n);

	for (int i = 0; i < n; i++) {
		deg[i] = (int)(A->rowptr[i + 1] - A->rowptr[i]);
	}

	int nq = 0;
	for (int pass = 0; pass < n; pass++) {
		// 未訪問のうち次数最小の節点から幅優先で辿る (非連結でも全部拾う)
		int s = -1;
		for (int i = 0; i < n; i++) {
			if (!seen[i] && ((s < 0) || (deg[i] < deg[s]))) s = i;
		}
		if (s < 0) break;

		perm[nq++] = (int32_t)s;
		seen[s] = 1;
		int head = nq - 1;
		while (head < nq) {
			const int32_t r = perm[head++];
			// 隣接節点を次数の昇順に並べて追加する
			const int64_t p0 = A->rowptr[r], p1 = A->rowptr[r + 1];
			const int nq0 = nq;
			for (int64_t p = p0; p < p1; p++) {
				const int32_t c = A->col[p];
				if ((c == r) || seen[c]) continue;
				seen[c] = 1;
				perm[nq++] = c;
			}
			for (int a = nq0; a < nq; a++) {
				for (int b = a + 1; b < nq; b++) {
					if (deg[perm[b]] < deg[perm[a]]) {
						const int32_t t = perm[a]; perm[a] = perm[b]; perm[b] = t;
					}
				}
			}
		}
	}

	// 反転する (Cuthill-McKee -> 逆 Cuthill-McKee)
	for (int i = 0; i < nq / 2; i++) {
		const int32_t t = perm[i];
		perm[i] = perm[nq - 1 - i];
		perm[nq - 1 - i] = t;
	}
	(void)fix;

	free(deg);
	free(seen);
}


/*
並べ替えとスカイラインの形を作る (実・複素で共通)。

perm/ip/first/sp を確保して埋め、プロファイル長を返す。
大きすぎるときは見積もりを報告して -1 を返す (nbyte は 1 要素のバイト数)。
*/
static int64_t skyline_setup(const crs_t *A, const unsigned char *fix, size_t nbyte,
	int32_t **pperm, int32_t **pip, int32_t **pfirst, int64_t **psp,
	FILE *fp_log, const char *label)
{
	const int n = (int)A->n;

	int32_t *perm = (int32_t *)malloc((size_t)n * sizeof(int32_t));
	int32_t *ip = (int32_t *)malloc((size_t)n * sizeof(int32_t));
	for (int i = 0; i < n; i++) perm[i] = (int32_t)i;
	rcm_order(A, fix, perm);
	for (int k = 0; k < n; k++) ip[perm[k]] = (int32_t)k;

	// 行毎のスカイライン先頭列 (並べ替え後)
	int32_t *first = (int32_t *)malloc((size_t)n * sizeof(int32_t));
	int64_t *sp = (int64_t *)malloc(((size_t)n + 1) * sizeof(int64_t));
	for (int k = 0; k < n; k++) {
		const int32_t r = perm[k];
		int32_t f = k;
		if (!(fix != NULL && fix[r])) {
			for (int64_t p = A->rowptr[r]; p < A->rowptr[r + 1]; p++) {
				const int32_t c = A->col[p];
				if (fix != NULL && fix[c]) continue;		// 固定列は落とす
				const int32_t kc = ip[c];
				if ((kc < f) && (kc <= k)) f = kc;
			}
		}
		first[k] = f;
	}
	sp[0] = 0;
	for (int k = 0; k < n; k++) sp[k + 1] = sp[k] + (k - first[k] + 1);
	const int64_t nprof = sp[n];
	const double mb = (double)nprof * nbyte / (1024 * 1024);

	if (fp_log != NULL) {
		fprintf(fp_log, "  %-10s direct : profile = %lld (%.1f MB), mean bandwidth = %.1f\n",
			label, (long long)nprof, mb, (double)nprof / ((n > 0) ? n : 1));
		fflush(fp_log);
	}
	if ((nprof * (int64_t)nbyte) > MAXPROFBYTE) {		// 32bit の size_t で溢れないよう先に広げる
		if (fp_log != NULL) {
			fprintf(fp_log, "*** the direct solver needs %.1f GB for this mesh; "
				"use the iterative solver (remove 'direct = 1')\n", mb / 1024);
		}
		printf("*** direct solver : profile too large (%.1f GB); remove 'direct = 1'\n",
			mb / 1024);
		free(perm); free(ip); free(first); free(sp);
		return -1;
	}

	*pperm = perm;
	*pip = ip;
	*pfirst = first;
	*psp = sp;

	return nprof;
}


int solver_direct(const crs_t *A, const double *b, double *x,
	const unsigned char *fix, FILE *fp_log, const char *label)
{
	const int n = (int)A->n;

	int32_t *perm = NULL, *ip = NULL, *first = NULL;
	int64_t *sp = NULL;
	const int64_t nprof = skyline_setup(A, fix, sizeof(double),
		&perm, &ip, &first, &sp, fp_log, label);
	if (nprof < 0) return -1;

	double *L = (double *)calloc((size_t)nprof, sizeof(double));
	if (L == NULL) {
		printf("%s\n", "*** direct solver : out of memory");
		free(perm); free(ip); free(first); free(sp);
		return -1;
	}

	// 並べ替えた行列をスカイラインに詰める (下三角のみ)
	for (int k = 0; k < n; k++) {
		const int32_t r = perm[k];
		if (fix != NULL && fix[r]) {
			L[sp[k] + (k - first[k])] = 1;			// 恒等行
			continue;
		}
		for (int64_t p = A->rowptr[r]; p < A->rowptr[r + 1]; p++) {
			const int32_t c = A->col[p];
			if (fix != NULL && fix[c]) continue;
			const int32_t kc = ip[c];
			if (kc > k) continue;					// 下三角だけ
			L[sp[k] + (kc - first[k])] = A->val[p];
		}
	}

	// スカイライン Cholesky (L L^T)
	for (int i = 0; i < n; i++) {
		const int32_t fi = first[i];
		for (int j = fi; j < i; j++) {
			const int32_t fj = first[j];
			const int32_t k0 = ((fi > fj) ? fi : fj);
			double s = L[sp[i] + (j - fi)];
			for (int k = k0; k < j; k++) {
				s -= L[sp[i] + (k - fi)] * L[sp[j] + (k - fj)];
			}
			L[sp[i] + (j - fi)] = s / L[sp[j] + (j - fj)];
		}
		double s = L[sp[i] + (i - fi)];
		for (int k = fi; k < i; k++) {
			const double v = L[sp[i] + (k - fi)];
			s -= v * v;
		}
		if (s <= 0) {
			// 正定値でない (材料が負、境界条件が足りない、等)
			printf("*** direct solver : the matrix is not positive definite "
				"(pivot %d = %.4e)\n", i, s);
			if (fp_log != NULL) {
				fprintf(fp_log, "*** direct solver : not positive definite at row %d\n", i);
			}
			free(L); free(perm); free(ip); free(first); free(sp);
			return -1;
		}
		L[sp[i] + (i - fi)] = sqrt(s);
	}

	// 前進代入 L y = P b
	double *y = (double *)malloc((size_t)n * sizeof(double));
	for (int i = 0; i < n; i++) {
		const int32_t fi = first[i];
		double s = b[perm[i]];
		for (int k = fi; k < i; k++) s -= L[sp[i] + (k - fi)] * y[k];
		y[i] = s / L[sp[i] + (i - fi)];
	}
	// 後退代入 L^T z = y (列方向に引くのでスカイラインのまま扱える)
	for (int i = n - 1; i >= 0; i--) {
		const int32_t fi = first[i];
		const double v = y[i] / L[sp[i] + (i - fi)];
		y[i] = v;
		for (int k = fi; k < i; k++) y[k] -= L[sp[i] + (k - fi)] * v;
	}
	for (int i = 0; i < n; i++) x[perm[i]] = y[i];

	free(y);
	free(L);
	free(perm);
	free(ip);
	free(first);
	free(sp);

	return 0;
}


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


/*
L D L^T x = b の三角求解 (並べ替え込み)。

yr/yi は長さ n の作業配列 (並べ替えた順の解)。x に直接書くと並べ替えを
戻すときに未処理の成分を踏むので、作業配列を通して最後に散らす。
*/
static void skyline_solve_c(int n, const double *Lr, const double *Li,
	const double *dr, const double *di, const int32_t *perm,
	const int32_t *first, const int64_t *sp,
	const double *br, const double *bi, double *xr, double *xi,
	double *yr, double *yi)
{
	// 前進代入 L y = P b (L の対角は 1)
	for (int i = 0; i < n; i++) {
		const int32_t fi = first[i];
		double sr = br[perm[i]], si = bi[perm[i]];
		for (int k = fi; k < i; k++) {
			double pr, pi;
			cmul(Lr[sp[i] + (k - fi)], Li[sp[i] + (k - fi)], yr[k], yi[k], &pr, &pi);
			sr -= pr;
			si -= pi;
		}
		yr[i] = sr;
		yi[i] = si;
	}
	// D で割る
	for (int i = 0; i < n; i++) {
		cdiv(yr[i], yi[i], dr[i], di[i], &yr[i], &yi[i]);
	}
	// 後退代入 L^T z = y (列方向に引くのでスカイラインのまま扱える)
	for (int i = n - 1; i >= 0; i--) {
		const int32_t fi = first[i];
		const double vr = yr[i], vi = yi[i];
		for (int k = fi; k < i; k++) {
			double pr, pi;
			cmul(Lr[sp[i] + (k - fi)], Li[sp[i] + (k - fi)], vr, vi, &pr, &pi);
			yr[k] -= pr;
			yi[k] -= pi;
		}
	}
	for (int i = 0; i < n; i++) {
		xr[perm[i]] = yr[i];
		xi[perm[i]] = yi[i];
	}
}


/*
複素対称 A = K + jω M のスカイライン LDL^T。

Cholesky と違い平方根を取らず、単位下三角 L と対角 D に分ける
(A = L D L^T、共役は取らない: 複素**対称**であって Hermite ではない)。
L の非対角は skyline に、D は別配列に持つ。Crout 形で行 i を左から埋め、
途中の t[k] = L_ik D_k を使い回すと D の引き直しが要らない。
*/
int solver_direct_c(const crs_t *K, const crs_t *M, double omega,
	const double *br, const double *bi, double *xr, double *xi,
	const unsigned char *fix, FILE *fp_log, const char *label)
{
	const int n = (int)K->n;

	// K と M は同じ CRS パターンを共有する前提 (虚部を K の格納位置に重ねる)。
	// 崩れていたら黙って誤答を出さずにここで落とす
	if ((M->n != K->n) || (M->nnz != K->nnz)) {
		printf("%s\n", "*** direct solver : K and M must share the CRS pattern");
		return -1;
	}

	int32_t *perm = NULL, *ip = NULL, *first = NULL;
	int64_t *sp = NULL;
	const int64_t nprof = skyline_setup(K, fix, 2 * sizeof(double),
		&perm, &ip, &first, &sp, fp_log, label);
	if (nprof < 0) return -1;

	double *Lr = (double *)calloc((size_t)nprof, sizeof(double));
	double *Li = (double *)calloc((size_t)nprof, sizeof(double));
	double *dr = (double *)malloc((size_t)n * sizeof(double));
	double *di = (double *)malloc((size_t)n * sizeof(double));
	double *tr = (double *)malloc((size_t)n * sizeof(double));
	double *ti = (double *)malloc((size_t)n * sizeof(double));
	if ((Lr == NULL) || (Li == NULL) || (dr == NULL) || (di == NULL)
	 || (tr == NULL) || (ti == NULL)) {
		printf("%s\n", "*** direct solver : out of memory");
		free(Lr); free(Li); free(dr); free(di); free(tr); free(ti);
		free(perm); free(ip); free(first); free(sp);
		return -1;
	}

	// 並べ替えた行列をスカイラインに詰める (下三角のみ)
	for (int k = 0; k < n; k++) {
		const int32_t r = perm[k];
		if (fix != NULL && fix[r]) {
			Lr[sp[k] + (k - first[k])] = 1;			// 恒等行
			continue;
		}
		for (int64_t p = K->rowptr[r]; p < K->rowptr[r + 1]; p++) {
			const int32_t c = K->col[p];
			if (fix != NULL && fix[c]) continue;
			const int32_t kc = ip[c];
			if (kc > k) continue;					// 下三角だけ
			Lr[sp[k] + (kc - first[k])] = K->val[p];
			Li[sp[k] + (kc - first[k])] = omega * M->val[p];
		}
	}

	// スカイライン LDL^T
	for (int i = 0; i < n; i++) {
		const int32_t fi = first[i];
		for (int j = fi; j < i; j++) {
			const int32_t fj = first[j];
			const int32_t k0 = ((fi > fj) ? fi : fj);
			double sr = Lr[sp[i] + (j - fi)], si = Li[sp[i] + (j - fi)];
			for (int k = k0; k < j; k++) {
				double pr, pi;
				cmul(tr[k], ti[k], Lr[sp[j] + (k - fj)], Li[sp[j] + (k - fj)], &pr, &pi);
				sr -= pr;
				si -= pi;
			}
			double lr, li;
			cdiv(sr, si, dr[j], di[j], &lr, &li);
			Lr[sp[i] + (j - fi)] = lr;
			Li[sp[i] + (j - fi)] = li;
			cmul(lr, li, dr[j], di[j], &tr[j], &ti[j]);		// t_j = L_ij D_j
		}
		double sr = Lr[sp[i] + (i - fi)], si = Li[sp[i] + (i - fi)];
		for (int k = fi; k < i; k++) {
			double pr, pi;
			cmul(tr[k], ti[k], Lr[sp[i] + (k - fi)], Li[sp[i] + (k - fi)], &pr, &pi);
			sr -= pr;
			si -= pi;
		}
		if ((sr == 0) && (si == 0)) {
			// 主小行列が特異 (ゲージ固定なしの A-φ 等)
			printf("*** direct solver : the matrix is singular (pivot %d = 0)\n", i);
			if (fp_log != NULL) {
				fprintf(fp_log, "*** direct solver : singular pivot at row %d\n", i);
			}
			free(Lr); free(Li); free(dr); free(di); free(tr); free(ti);
			free(perm); free(ip); free(first); free(sp);
			return -1;
		}
		dr[i] = sr;
		di[i] = si;
	}

	double *zr = (double *)malloc((size_t)n * sizeof(double));
	double *zi = (double *)malloc((size_t)n * sizeof(double));
	double *w1 = (double *)malloc((size_t)n * sizeof(double));
	double *w2 = (double *)malloc((size_t)n * sizeof(double));
	double *w3 = (double *)malloc((size_t)n * sizeof(double));
	double *w4 = (double *)malloc((size_t)n * sizeof(double));

	skyline_solve_c(n, Lr, Li, dr, di, perm, first, sp, br, bi, xr, xi, w1, w2);

	/*
	**真の残差を取る。** ピボット選択をしないので、A (A-φ 連成) のように
	主小行列の非特異性が保証できない系では分解が破綻し得る。反復解法と同じ
	意味の数字 ‖b − A x‖/‖b‖ をログに残せば、壊れたときに黙って通らない。

	機械精度までは落ちない。この系は ν ~ 1/μ0 (1e6) と ωσ が桁で違って
	‖A‖‖x‖ >> ‖b‖ になるため、後退安定な分解でも相対残差はそこで頭打ちに
	なる (実測: plate_line_ac で 5.9e-11、bar_eddy で 4.7e-14。反復改良を
	足しても 3.6e-11 までしか下がらず、条件数で決まっていることが分かる)。
	そのため判定は 1e-6 と緩くしてある。分解が破綻すれば残差は O(1) になる
	ので、この緩さでも検出は効く
	*/
	crs_spmv_c(K, M, omega, xr, xi, zr, zi, fix, w1, w2, w3, w4);
	double bn = 0, rn = 0;
	for (int i = 0; i < n; i++) {
		const double er = br[i] - zr[i], ei = bi[i] - zi[i];
		bn += (br[i] * br[i]) + (bi[i] * bi[i]);
		rn += (er * er) + (ei * ei);
	}
	const double resid = ((bn > 0) ? sqrt(rn / bn) : 0);

	free(Lr); free(Li); free(dr); free(di); free(tr); free(ti);
	free(perm); free(ip); free(first); free(sp);
	free(zr); free(zi); free(w1); free(w2); free(w3); free(w4);

	const int ok = (resid < 1e-6);
	if (fp_log != NULL) {
		fprintf(fp_log, "  %-10s %8d %13.5e %s\n", label, 0, resid,
			(ok ? "direct" : "*** direct solver : inaccurate factorization"));
		fflush(fp_log);
	}
	if (!ok) {
		printf("*** direct solver : the residual is %.3e; use the iterative solver\n", resid);
	}

	return (ok ? 0 : -1);
}
