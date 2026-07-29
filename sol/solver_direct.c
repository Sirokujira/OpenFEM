/*
solver_direct.c

直接解法 (RCM 並べ替え + スカイライン Cholesky)。

反復解法しか無いと、悪条件な問題 (高コントラストな材料、扁平な要素) で
収束しないときに逃げ道が無い。`direct = 1` でこちらに切り替える。

対象は**実対称正定値**の系だけ (C / L / R / M と、非線形反復の内側)。
渦電流 (F / A) の複素対称系は COCG のままで、入力段で弾く。

Dirichlet の扱いは Jacobi-PCG と完全に同じにする:
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

// プロファイルの上限 (要素数)。超えたら反復解法を使うよう促す
#define MAXPROFILE (40000000LL)		// 約 320 MB


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


int solver_direct(const crs_t *A, const double *b, double *x,
	const unsigned char *fix, FILE *fp_log, const char *label)
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

	if (fp_log != NULL) {
		fprintf(fp_log, "  %-10s direct : profile = %lld (%.1f MB), mean bandwidth = %.1f\n",
			label, (long long)nprof, (double)nprof * sizeof(double) / (1024 * 1024),
			(double)nprof / ((n > 0) ? n : 1));
		fflush(fp_log);
	}
	if (nprof > MAXPROFILE) {
		if (fp_log != NULL) {
			fprintf(fp_log, "*** the direct solver needs %.1f GB for this mesh; "
				"use the iterative solver (remove 'direct = 1')\n",
				(double)nprof * sizeof(double) / (1024.0 * 1024 * 1024));
		}
		printf("*** direct solver : profile too large (%.1f GB); remove 'direct = 1'\n",
			(double)nprof * sizeof(double) / (1024.0 * 1024 * 1024));
		free(perm); free(ip); free(first); free(sp);
		return -1;
	}

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
