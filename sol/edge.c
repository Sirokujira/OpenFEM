/*
edge.c

3 次元渦電流 (A-φ) 解析の基盤 : 辺 (edge) の抽出・辺ベース CRS・
Whitney 辺要素 (1 次 Nedelec 要素) の要素行列。

節点要素ではベクトルポテンシャル A の接線連続性を正しく課せないため、
3 次元の渦電流解析には辺要素が要る。自由度が「節点」から「辺」に変わるので、
節点ベースの CRS とは別に辺ベースの連結を作る。

Whitney 基底 (辺 e = (i, j)) :

	W_e     = λ_i ∇λ_j - λ_j ∇λ_i
	∇×W_e   = 2 ∇λ_i × ∇λ_j        (要素内で一定)

要素行列 (1 次四面体、体積 V、g_ab = ∇λ_a・∇λ_b) :

	S_ef = ∫ (∇×W_e)^T ν (∇×W_f) dV = 4 V (c_e)^T ν (c_f),  c = ∇λ_a × ∇λ_b
	T_ef = ∫ σ W_e・W_f dV
	     = σ V/20 [ (1+δ_ac) g_bd - (1+δ_ad) g_bc - (1+δ_bc) g_ad + (1+δ_bd) g_ac ]

	(∫ λ_a λ_b dV = V (1 + δ_ab) / 20 を使うので閉形式で厳密)

辺の向きは全体節点番号の小さい方から大きい方へ揃える。局所の並び (a, b) が
逆向きなら符号 -1 を掛ける。
*/

#include "fem.h"
#include "fem_prototype.h"

// 四面体の局所辺 (節点対)
static const int EDGE_NODE[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};


static int cmp_i32(const void *a, const void *b)
{
	const int32_t x = *(const int32_t *)a;
	const int32_t y = *(const int32_t *)b;

	return ((x < y) ? -1 : ((x > y) ? 1 : 0));
}


// ---- 辺の抽出 ----
//
// 節点 i について「i より大きい隣接節点」を昇順に並べたものを辺とする。
// 辺番号は EdgePtr[i] からの通し番号になり、二分探索で引ける。
void edge_build(void)
{
	const int n = NNode;

	// 節点毎の (i < j) 隣接候補を数える
	int *cnt = (int *)malloc((size_t)n * sizeof(int));
	memset(cnt, 0, (size_t)n * sizeof(int));
	for (int e = 0; e < NTet; e++) {
		const int32_t *nd = &Tet[e * 4];
		for (int k = 0; k < 6; k++) {
			const int32_t a = nd[EDGE_NODE[k][0]];
			const int32_t b = nd[EDGE_NODE[k][1]];
			cnt[(a < b) ? a : b]++;
		}
	}

	int64_t *ptr = (int64_t *)malloc(((size_t)n + 1) * sizeof(int64_t));
	ptr[0] = 0;
	for (int i = 0; i < n; i++) ptr[i + 1] = ptr[i] + cnt[i];
	int32_t *lst = (int32_t *)malloc((size_t)ptr[n] * sizeof(int32_t));
	memset(cnt, 0, (size_t)n * sizeof(int));
	for (int e = 0; e < NTet; e++) {
		const int32_t *nd = &Tet[e * 4];
		for (int k = 0; k < 6; k++) {
			const int32_t a = nd[EDGE_NODE[k][0]];
			const int32_t b = nd[EDGE_NODE[k][1]];
			const int32_t lo = ((a < b) ? a : b);
			const int32_t hi = ((a < b) ? b : a);
			lst[ptr[lo] + cnt[lo]] = hi;
			cnt[lo]++;
		}
	}

	// 整列・重複除去して辺番号を振る
	EdgePtr = (int64_t *)malloc(((size_t)n + 1) * sizeof(int64_t));
	int *un = (int *)malloc((size_t)n * sizeof(int));
	for (int i = 0; i < n; i++) {
		const int64_t p0 = ptr[i], p1 = ptr[i + 1];
		const int m = (int)(p1 - p0);
		qsort(&lst[p0], (size_t)m, sizeof(int32_t), cmp_i32);
		int u = 0;
		for (int q = 0; q < m; q++) {
			if ((q == 0) || (lst[p0 + q] != lst[p0 + q - 1])) u++;
		}
		un[i] = u;
	}
	EdgePtr[0] = 0;
	for (int i = 0; i < n; i++) EdgePtr[i + 1] = EdgePtr[i] + un[i];
	NEdge = (int)EdgePtr[n];

	EdgeTo = (int32_t *)malloc((size_t)NEdge * sizeof(int32_t));
	for (int i = 0; i < n; i++) {
		const int64_t p0 = ptr[i], p1 = ptr[i + 1];
		const int m = (int)(p1 - p0);
		int64_t w = EdgePtr[i];
		for (int q = 0; q < m; q++) {
			if ((q == 0) || (lst[p0 + q] != lst[p0 + q - 1])) EdgeTo[w++] = lst[p0 + q];
		}
	}

	free(cnt);
	free(ptr);
	free(lst);
	free(un);

	// 辺の始点 (節点番号の小さい方)
	EdgeFrom = (int32_t *)malloc((size_t)NEdge * sizeof(int32_t));
	for (int i = 0; i < n; i++) {
		for (int64_t p = EdgePtr[i]; p < EdgePtr[i + 1]; p++) EdgeFrom[p] = (int32_t)i;
	}

	// 四面体毎の辺番号と向き
	TetEdge = (int32_t *)malloc((size_t)NTet * 6 * sizeof(int32_t));
	TetEdgeSgn = (signed char *)malloc((size_t)NTet * 6 * sizeof(signed char));
	for (int e = 0; e < NTet; e++) {
		const int32_t *nd = &Tet[e * 4];
		for (int k = 0; k < 6; k++) {
			const int32_t a = nd[EDGE_NODE[k][0]];
			const int32_t b = nd[EDGE_NODE[k][1]];
			TetEdge[(e * 6) + k] = (int32_t)edge_id(a, b);
			TetEdgeSgn[(e * 6) + k] = (signed char)((a < b) ? +1 : -1);
		}
	}
}


// 節点対 (a, b) の辺番号 (見つからなければ -1)
int64_t edge_id(int32_t a, int32_t b)
{
	const int32_t lo = ((a < b) ? a : b);
	const int32_t hi = ((a < b) ? b : a);

	int64_t p = EdgePtr[lo];
	int64_t q = EdgePtr[lo + 1] - 1;
	while (p <= q) {
		const int64_t m = (p + q) / 2;
		if      (EdgeTo[m] < hi) p = m + 1;
		else if (EdgeTo[m] > hi) q = m - 1;
		else return m;
	}

	return -1;
}


void edge_free(void)
{
	free(EdgePtr);
	free(EdgeTo);
	free(EdgeFrom);
	free(TetEdge);
	free(TetEdgeSgn);
	EdgePtr = NULL;
	EdgeTo = NULL;
	EdgeFrom = NULL;
	TetEdge = NULL;
	TetEdgeSgn = NULL;
	NEdge = 0;
}


// ---- 辺ベース CRS ----

void crs_alloc_edge(crs_t *A)
{
	const int ne = NEdge;

	// 辺毎の四面体リスト
	int *cnt = (int *)malloc((size_t)ne * sizeof(int));
	memset(cnt, 0, (size_t)ne * sizeof(int));
	for (int e = 0; e < NTet; e++) {
		for (int k = 0; k < 6; k++) cnt[TetEdge[(e * 6) + k]]++;
	}
	int64_t *ptr = (int64_t *)malloc(((size_t)ne + 1) * sizeof(int64_t));
	ptr[0] = 0;
	for (int i = 0; i < ne; i++) ptr[i + 1] = ptr[i] + cnt[i];
	int32_t *lst = (int32_t *)malloc((size_t)ptr[ne] * sizeof(int32_t));
	memset(cnt, 0, (size_t)ne * sizeof(int));
	for (int e = 0; e < NTet; e++) {
		for (int k = 0; k < 6; k++) {
			const int32_t ed = TetEdge[(e * 6) + k];
			lst[ptr[ed] + cnt[ed]] = (int32_t)e;
			cnt[ed]++;
		}
	}

	A->n = ne;
	A->rowptr = (int64_t *)malloc(((size_t)ne + 1) * sizeof(int64_t));

	int cap = 64;
	int32_t *work = (int32_t *)malloc((size_t)cap * sizeof(int32_t));
	int *rown = (int *)malloc((size_t)ne * sizeof(int));

	for (int i = 0; i < ne; i++) {
		const int64_t p0 = ptr[i], p1 = ptr[i + 1];
		const int need = (int)(p1 - p0) * 6;
		if (need > cap) {
			cap = need;
			work = (int32_t *)realloc(work, (size_t)cap * sizeof(int32_t));
		}
		int m = 0;
		for (int64_t p = p0; p < p1; p++) {
			const int32_t t = lst[p];
			for (int k = 0; k < 6; k++) work[m++] = TetEdge[(t * 6) + k];
		}
		qsort(work, (size_t)m, sizeof(int32_t), cmp_i32);
		int u = 0;
		for (int q = 0; q < m; q++) {
			if ((q == 0) || (work[q] != work[q - 1])) u++;
		}
		rown[i] = u;
	}

	A->rowptr[0] = 0;
	for (int i = 0; i < ne; i++) A->rowptr[i + 1] = A->rowptr[i] + rown[i];
	A->nnz = A->rowptr[ne];
	A->col = (int32_t *)malloc((size_t)A->nnz * sizeof(int32_t));
	A->val = (double *)malloc((size_t)A->nnz * sizeof(double));

	for (int i = 0; i < ne; i++) {
		const int64_t p0 = ptr[i], p1 = ptr[i + 1];
		int m = 0;
		for (int64_t p = p0; p < p1; p++) {
			const int32_t t = lst[p];
			for (int k = 0; k < 6; k++) work[m++] = TetEdge[(t * 6) + k];
		}
		qsort(work, (size_t)m, sizeof(int32_t), cmp_i32);
		int64_t w = A->rowptr[i];
		for (int q = 0; q < m; q++) {
			if ((q == 0) || (work[q] != work[q - 1])) A->col[w++] = work[q];
		}
	}

	free(cnt);
	free(ptr);
	free(lst);
	free(work);
	free(rown);

	crs_zero(A);
}


static int64_t crs_find_edge(const crs_t *A, int32_t row, int32_t col)
{
	int64_t lo = A->rowptr[row];
	int64_t hi = A->rowptr[row + 1] - 1;

	while (lo <= hi) {
		const int64_t mid = (lo + hi) / 2;
		if      (A->col[mid] < col) lo = mid + 1;
		else if (A->col[mid] > col) hi = mid - 1;
		else return mid;
	}

	return -1;
}


// ---- Whitney 辺要素の要素行列 ----

// 要素 e の回転回転行列 se[6][6] と質量行列 te[6][6]
// nu : 磁気抵抗率テンソル (6 成分)、sig : 導電率 (等方)
void edge_element(int e, const double nu[6], double sig,
	double se[6][6], double te[6][6])
{
	const int32_t *nd = &Tet[e * 4];
	double g[4][3], vol;

	for (int l = 0; l < 6; l++) {
		for (int m = 0; m < 6; m++) {
			se[l][m] = te[l][m] = 0;
		}
	}
	if (tet_grad_pub(nd, g, &vol)) return;

	// ∇λ_a・∇λ_b
	double gg[4][4];
	for (int a = 0; a < 4; a++) {
		for (int b = 0; b < 4; b++) {
			gg[a][b] = (g[a][0] * g[b][0]) + (g[a][1] * g[b][1]) + (g[a][2] * g[b][2]);
		}
	}

	// 各辺の ∇λ_a × ∇λ_b (∇×W_e = 2 * これ)
	double cr[6][3];
	for (int k = 0; k < 6; k++) {
		const int a = EDGE_NODE[k][0], b = EDGE_NODE[k][1];
		cr[k][0] = (g[a][1] * g[b][2]) - (g[a][2] * g[b][1]);
		cr[k][1] = (g[a][2] * g[b][0]) - (g[a][0] * g[b][2]);
		cr[k][2] = (g[a][0] * g[b][1]) - (g[a][1] * g[b][0]);
	}

	for (int k = 0; k < 6; k++) {
		const int a = EDGE_NODE[k][0], b = EDGE_NODE[k][1];
		for (int l = 0; l < 6; l++) {
			const int c = EDGE_NODE[l][0], d = EDGE_NODE[l][1];

			// 回転回転項 : 4 V (cr_k)^T nu (cr_l)
			const double cn = (nu[0] * cr[k][0] * cr[l][0])
			                + (nu[1] * cr[k][1] * cr[l][1])
			                + (nu[2] * cr[k][2] * cr[l][2])
			                + (nu[3] * ((cr[k][0] * cr[l][1]) + (cr[k][1] * cr[l][0])))
			                + (nu[4] * ((cr[k][1] * cr[l][2]) + (cr[k][2] * cr[l][1])))
			                + (nu[5] * ((cr[k][2] * cr[l][0]) + (cr[k][0] * cr[l][2])));
			se[k][l] = 4 * vol * cn;

			// 質量項 : σ V/20 [ (1+δ_ac) g_bd - (1+δ_ad) g_bc
			//                  - (1+δ_bc) g_ad + (1+δ_bd) g_ac ]
			const double t = ((a == c) ? 2.0 : 1.0) * gg[b][d]
			               - ((a == d) ? 2.0 : 1.0) * gg[b][c]
			               - ((b == c) ? 2.0 : 1.0) * gg[a][d]
			               + ((b == d) ? 2.0 : 1.0) * gg[a][c];
			te[k][l] = sig * vol * t / 20;
		}
	}
}


// 全体行列の作成 (辺要素)
//   S : 回転回転行列 (ν)、T : 質量行列 (σ)。どちらも NULL 可
void assemble_edge(crs_t *S, crs_t *T)
{
	if (S != NULL) crs_zero(S);
	if (T != NULL) crs_zero(T);

	for (int e = 0; e < NTet; e++) {
		const int m = TetMat[e];
		double nu[6];
		material_coef_pub(m, 3, nu);			// ν = (μ0 μ~)^-1
		const double sig = Material[m].sigma;

		double se[6][6], te[6][6];
		edge_element(e, nu, sig, se, te);

		const int32_t *ed = &TetEdge[e * 6];
		const signed char *sg = &TetEdgeSgn[e * 6];
		for (int k = 0; k < 6; k++) {
			for (int l = 0; l < 6; l++) {
				const double s = (double)(sg[k] * sg[l]);
				if (S != NULL) {
					const int64_t p = crs_find_edge(S, ed[k], ed[l]);
					if (p >= 0) S->val[p] += s * se[k][l];
				}
				if (T != NULL) {
					const int64_t p = crs_find_edge(T, ed[k], ed[l]);
					if (p >= 0) T->val[p] += s * te[k][l];
				}
			}
		}
	}
}


// ---- 段階 1+2 の検証 (analysis = E) ----
//
// 辺要素の性質を厳密に満たすべき 3 つの量で検査する。
// いずれも解析解が閉形式で分かるので、要素行列・辺の向き・組み立ての
// どこかが間違っていれば必ず落ちる。
//
//  (a) 勾配は回転回転行列の零空間に入る
//      節点値 φ の勾配は辺自由度 u_e = φ_hi - φ_lo で厳密に表せるので u^T S u = 0
//  (b) 一様場 E0 の質量 : u_e = E0・(p_hi - p_lo) として u^T T u = |E0|^2 Σ σ_e V_e
//  (c) 回転場 E = (1/2)(z × r) : curl E = z なので u^T S u = Σ ν_e V_e
//      (この場は最低次 Nedelec 空間 {a + b × r} に含まれるので補間は厳密)
int solve_edge_test(FILE *fp_log)
{
	int ierr = 0;

	edge_build();

	fprintf(fp_log, "\n=== edge element (Nedelec) self test ===\n");
	fprintf(fp_log, "  nodes = %d, tetrahedra = %d, edges = %d\n", NNode, NTet, NEdge);
	fflush(fp_log);

	crs_t S, T;
	crs_alloc_edge(&S);
	crs_alloc_edge(&T);
	assemble_edge(&S, &T);
	fprintf(fp_log, "  edge matrix : %lld nonzeros (%.1f per row)\n",
		(long long)S.nnz, (double)S.nnz / ((NEdge > 0) ? NEdge : 1));

	// 材料の体積積分 (Σ σ V, Σ ν V)
	double svol = 0, nvol = 0, vtot = 0;
	for (int e = 0; e < NTet; e++) {
		double g[4][3], vol;
		if (tet_grad_pub(&Tet[e * 4], g, &vol)) continue;
		const int m = TetMat[e];
		double nu[6];
		material_coef_pub(m, 3, nu);
		svol += Material[m].sigma * vol;
		nvol += nu[0] * vol;
		vtot += vol;
	}
	fprintf(fp_log, "  volume = %.6e [m^3], sum(sigma V) = %.6e, sum(nu V) = %.6e\n",
		vtot, svol, nvol);

	const int ne = NEdge;
	double *u = (double *)malloc((size_t)ne * sizeof(double));
	double *y = (double *)malloc((size_t)ne * sizeof(double));

	// 辺の端点 (lo -> hi)
	int32_t *elo = (int32_t *)malloc((size_t)ne * sizeof(int32_t));
	for (int i = 0; i < NNode; i++) {
		for (int64_t p = EdgePtr[i]; p < EdgePtr[i + 1]; p++) elo[p] = (int32_t)i;
	}

	double dmax = 0;
	for (int i = 0; i < ne; i++) {
		const int64_t p = S.rowptr[i];
		for (int64_t q = p; q < S.rowptr[i + 1]; q++) {
			if ((S.col[q] == i) && (fabs(S.val[q]) > dmax)) dmax = fabs(S.val[q]);
		}
	}

	// (a) 勾配の零空間
	{
		double *phi = (double *)malloc((size_t)NNode * sizeof(double));
		for (int i = 0; i < NNode; i++) {
			phi[i] = fmod(i * 0.61803398875, 1.0) - 0.5;	// 決定的な擬似乱数
		}
		double uu = 0;
		for (int i = 0; i < ne; i++) {
			u[i] = phi[EdgeTo[i]] - phi[elo[i]];
			uu += u[i] * u[i];
		}
		crs_spmv(&S, u, y, NULL);
		double q = 0;
		for (int i = 0; i < ne; i++) q += u[i] * y[i];
		const double rel = fabs(q) / ((dmax * uu > 0) ? (dmax * uu) : 1);
		fprintf(fp_log, "  (a) gradient null space : u^T S u / (max|Sii| u^T u) = %.3e\n", rel);
		if (rel > 1e-10) {
			fprintf(fp_log, "*** the curl-curl matrix does not annihilate gradients\n");
			ierr = 1;
		}
		free(phi);
	}

	// (b) 一様場の質量行列
	{
		const double e0[3] = {0.3, -0.7, 1.1};
		const double e2 = (e0[0] * e0[0]) + (e0[1] * e0[1]) + (e0[2] * e0[2]);
		for (int i = 0; i < ne; i++) {
			const int32_t a = elo[i], b = EdgeTo[i];
			u[i] = (e0[0] * (Xp[b] - Xp[a])) + (e0[1] * (Yp[b] - Yp[a]))
			     + (e0[2] * (Zp[b] - Zp[a]));
		}
		crs_spmv(&T, u, y, NULL);
		double q = 0;
		for (int i = 0; i < ne; i++) q += u[i] * y[i];
		const double exact = e2 * svol;
		const double rel = fabs(q - exact) / ((exact != 0) ? fabs(exact) : 1);
		fprintf(fp_log, "  (b) uniform field mass  : u^T T u = %.10e, exact = %.10e, err = %.3e\n",
			q, exact, rel);
		if (rel > 1e-10) {
			fprintf(fp_log, "*** the edge mass matrix is wrong\n");
			ierr = 1;
		}

		// 一様場は回転が 0 なので回転回転行列でも 0 になる
		crs_spmv(&S, u, y, NULL);
		double q2 = 0, uu = 0;
		for (int i = 0; i < ne; i++) {
			q2 += u[i] * y[i];
			uu += u[i] * u[i];
		}
		const double rel2 = fabs(q2) / ((dmax * uu > 0) ? (dmax * uu) : 1);
		fprintf(fp_log, "  (b') uniform field curl : u^T S u / (max|Sii| u^T u) = %.3e\n", rel2);
		if (rel2 > 1e-10) ierr = 1;
	}

	// (c) 回転場 E = (1/2)(z × r) -> curl E = z
	{
		for (int i = 0; i < ne; i++) {
			const int32_t a = elo[i], b = EdgeTo[i];
			const double xm = (Xp[a] + Xp[b]) / 2;
			const double ym = (Yp[a] + Yp[b]) / 2;
			// E = (1/2)(-y, x, 0)、直線辺上で E は 1 次なので中点則で厳密
			const double ex = -ym / 2, ey = xm / 2;
			u[i] = (ex * (Xp[b] - Xp[a])) + (ey * (Yp[b] - Yp[a]));
		}
		crs_spmv(&S, u, y, NULL);
		double q = 0;
		for (int i = 0; i < ne; i++) q += u[i] * y[i];
		const double rel = fabs(q - nvol) / ((nvol != 0) ? fabs(nvol) : 1);
		fprintf(fp_log, "  (c) rotational field    : u^T S u = %.10e, exact = %.10e, err = %.3e\n",
			q, nvol, rel);
		if (rel > 1e-10) {
			fprintf(fp_log, "*** the curl-curl matrix is wrong\n");
			ierr = 1;
		}
	}

	// (d) 対称性
	{
		double amax = 0, sdif = 0, tdif = 0;
		for (int i = 0; i < ne; i++) {
			for (int64_t p = S.rowptr[i]; p < S.rowptr[i + 1]; p++) {
				const int32_t j = S.col[p];
				const int64_t pj = crs_find_edge(&S, j, (int32_t)i);
				if (pj < 0) continue;
				if (fabs(S.val[p]) > amax) amax = fabs(S.val[p]);
				const double d = fabs(S.val[p] - S.val[pj]);
				if (d > sdif) sdif = d;
				const double dt = fabs(T.val[p] - T.val[pj]);
				if (dt > tdif) tdif = dt;
			}
		}
		const double rel = ((amax > 0) ? (sdif / amax) : 0);
		fprintf(fp_log, "  (d) symmetry            : max|S-S^T|/max|S| = %.3e\n", rel);
		if (rel > 1e-12) ierr = 1;
		(void)tdif;
	}

	// (e) ゲージ固定 (tree-cotree)
	{
		unsigned char *tree = (unsigned char *)malloc((size_t)ne * sizeof(unsigned char));
		const int ncomp = edge_tree(tree);
		int ntree = 0;
		for (int i = 0; i < ne; i++) ntree += tree[i];
		fprintf(fp_log, "  (e) spanning tree       : %d tree edges, %d components "
			"(nodes - components = %d), gauged unknowns = %d\n",
			ntree, ncomp, NNode - ncomp, ne - ntree);
		if (ntree != (NNode - ncomp)) {
			fprintf(fp_log, "*** the spanning tree is not spanning\n");
			ierr = 1;
		}

		// 勾配はゲージ前は S の零空間に入るが、木辺を落とすと外れる
		double *phi = (double *)malloc((size_t)NNode * sizeof(double));
		for (int i = 0; i < NNode; i++) phi[i] = fmod(i * 0.31830988618, 1.0) - 0.5;
		edge_grad(phi, u);
		crs_spmv(&S, u, y, NULL);
		double q0 = 0, uu = 0;
		for (int i = 0; i < ne; i++) {
			q0 += u[i] * y[i];
			uu += u[i] * u[i];
		}
		for (int i = 0; i < ne; i++) {
			if (tree[i]) u[i] = 0;			// ゲージ : 木辺を 0 にする
		}
		crs_spmv(&S, u, y, NULL);
		double q1 = 0, vv = 0;
		for (int i = 0; i < ne; i++) {
			q1 += u[i] * y[i];
			vv += u[i] * u[i];
		}
		const double r0 = fabs(q0) / ((dmax * uu > 0) ? (dmax * uu) : 1);
		const double r1 = fabs(q1) / ((dmax * vv > 0) ? (dmax * vv) : 1);
		fprintf(fp_log, "      gradient energy : before gauge %.3e, after gauge %.3e\n", r0, r1);
		if (r1 < 1e-6) {
			fprintf(fp_log, "*** the gauge did not remove the gradient null space\n");
			ierr = 1;
		}
		free(phi);
		free(tree);
	}

	// (f) 前処理 (Hiptmair vs Jacobi)
	//     A = S + β T (β を小さくすると回転回転支配になり悪条件になる)
	{
		const double beta = 1e-3;
		crs_t A2;
		crs_alloc_edge(&A2);
		for (int64_t p = 0; p < A2.nnz; p++) A2.val[p] = S.val[p] + (beta * T.val[p]);

		crs_t N;
		edge_nodal_aux(&A2, &N);

		unsigned char *nofix = (unsigned char *)malloc((size_t)ne * sizeof(unsigned char));
		memset(nofix, 0, (size_t)ne * sizeof(unsigned char));

		double *xe = (double *)malloc((size_t)ne * sizeof(double));
		double *bb = (double *)malloc((size_t)ne * sizeof(double));
		double *xs = (double *)malloc((size_t)ne * sizeof(double));
		for (int i = 0; i < ne; i++) xe[i] = fmod(i * 0.7548776662, 1.0) - 0.5;
		crs_spmv(&A2, xe, bb, NULL);

		const double tol = 1e-12;
		const int mx = 20000;

		// 誤差は 2 通りで測る。回転回転系の悪条件は「勾配に近い成分」に効くが、
		// その成分は curl A = B に寄与しないので物理量には影響しない。
		//   raw  : 辺自由度そのものの相対誤差 (ゲージ成分を含む)
		//   curl : S ノルムの相対誤差 = 磁束密度の誤差 (物理量)
		double *er = (double *)malloc((size_t)ne * sizeof(double));
		double xn = 0, xs2 = 0;
		for (int i = 0; i < ne; i++) xn += xe[i] * xe[i];
		crs_spmv(&S, xe, y, NULL);
		for (int i = 0; i < ne; i++) xs2 += xe[i] * y[i];

		int it[2];
		double eraw[2], ecurl[2];
		for (int mode = 0; mode < 2; mode++) {
			it[mode] = solver_cg_edge(&A2, ((mode == 1) ? &N : NULL), bb, xs, nofix,
				mode, mx, 0, tol, NULL, ((mode == 1) ? "hiptmair" : "jacobi"));
			double e1 = 0;
			for (int i = 0; i < ne; i++) {
				er[i] = xs[i] - xe[i];
				e1 += er[i] * er[i];
			}
			eraw[mode] = sqrt(e1 / ((xn > 0) ? xn : 1));
			crs_spmv(&S, er, y, NULL);
			double e2 = 0;
			for (int i = 0; i < ne; i++) e2 += er[i] * y[i];
			ecurl[mode] = sqrt(fabs(e2) / ((xs2 > 0) ? xs2 : 1));
		}
		free(er);

		fprintf(fp_log, "  (f) preconditioner      : A = S + %.0e T, tol %.0e\n", beta, tol);
		for (int mode = 0; mode < 2; mode++) {
			fprintf(fp_log, "      %-8s : %6d iterations, raw error %.3e, curl error %.3e\n",
				((mode == 1) ? "Hiptmair" : "Jacobi"),
				(it[mode] < 0 ? -it[mode] : it[mode]), eraw[mode], ecurl[mode]);
		}
		if (it[1] > 0) {
			fprintf(fp_log, "      -> %.1fx fewer iterations, %.0fx smaller raw error\n",
				(double)(it[0] < 0 ? -it[0] : it[0]) / (double)it[1],
				(eraw[1] > 0) ? (eraw[0] / eraw[1]) : 0.0);
		}
		if ((it[0] < 0) || (it[1] < 0)) {
			fprintf(fp_log, "*** the edge solver did not converge\n");
			ierr = 1;
		}
		// 物理量 (curl) は両者とも正しく出ていなければならない
		if ((ecurl[0] > 1e-6) || (ecurl[1] > 1e-6)) {
			fprintf(fp_log, "*** the edge solver did not reproduce the curl of the "
				"manufactured solution\n");
			ierr = 1;
		}
		// 前処理は反復回数と (ゲージ成分を含む) 生の誤差の両方を改善すること
		if ((it[1] > 0) && (it[0] > 0) && (it[1] >= it[0])) {
			fprintf(fp_log, "*** the Hiptmair preconditioner did not reduce the iteration count\n");
			ierr = 1;
		}
		if (eraw[1] > (eraw[0] / 10)) {
			fprintf(fp_log, "*** the Hiptmair preconditioner did not reduce the raw error\n");
			ierr = 1;
		}

		free(nofix);
		free(xe);
		free(bb);
		free(xs);
		crs_free(&A2);
		crs_free(&N);
	}

	fprintf(fp_log, "  %s\n", (ierr ? "*** edge element self test FAILED" : "edge element self test passed"));
	fflush(fp_log);

	free(u);
	free(y);
	free(elo);
	crs_free(&S);
	crs_free(&T);
	edge_free();

	return ierr;
}


// ---- 段階 3 : ゲージ固定 (tree-cotree) ----
//
// 回転回転行列 S の零空間は「節点関数の勾配」全体 (連結成分あたり次元 節点数-1)
// なので、非導電領域を含む系は特異になる。全域木 (spanning tree) の辺で
// A_e = 0 と置くと、この零空間がちょうど消える:
//   勾配 u_e = φ_hi - φ_lo が全ての木辺で 0 <=> φ が木上で一定 <=> φ が一定 <=> u = 0
//
// tree[e] = 1 が木辺。戻り値は連結成分数 (木辺の数は 節点数 - 連結成分数)。
int edge_tree(unsigned char *tree)
{
	const int n = NNode;

	// 双方向の隣接リスト (辺番号つき)
	int *cnt = (int *)malloc((size_t)n * sizeof(int));
	memset(cnt, 0, (size_t)n * sizeof(int));
	for (int e = 0; e < NEdge; e++) {
		cnt[EdgeFrom[e]]++;
		cnt[EdgeTo[e]]++;
	}
	int64_t *ptr = (int64_t *)malloc(((size_t)n + 1) * sizeof(int64_t));
	ptr[0] = 0;
	for (int i = 0; i < n; i++) ptr[i + 1] = ptr[i] + cnt[i];
	int32_t *ato = (int32_t *)malloc((size_t)ptr[n] * sizeof(int32_t));
	int32_t *aed = (int32_t *)malloc((size_t)ptr[n] * sizeof(int32_t));
	memset(cnt, 0, (size_t)n * sizeof(int));
	for (int e = 0; e < NEdge; e++) {
		const int32_t a = EdgeFrom[e], b = EdgeTo[e];
		ato[ptr[a] + cnt[a]] = b;  aed[ptr[a] + cnt[a]] = (int32_t)e;  cnt[a]++;
		ato[ptr[b] + cnt[b]] = a;  aed[ptr[b] + cnt[b]] = (int32_t)e;  cnt[b]++;
	}

	memset(tree, 0, (size_t)NEdge * sizeof(unsigned char));
	unsigned char *vis = (unsigned char *)malloc((size_t)n * sizeof(unsigned char));
	memset(vis, 0, (size_t)n * sizeof(unsigned char));
	int32_t *que = (int32_t *)malloc((size_t)n * sizeof(int32_t));

	int ncomp = 0;
	for (int root = 0; root < n; root++) {
		if (vis[root]) continue;
		ncomp++;
		int head = 0, tail = 0;
		vis[root] = 1;
		que[tail++] = (int32_t)root;
		while (head < tail) {
			const int32_t u = que[head++];
			for (int64_t p = ptr[u]; p < ptr[u + 1]; p++) {
				const int32_t v = ato[p];
				if (vis[v]) continue;
				vis[v] = 1;
				tree[aed[p]] = 1;
				que[tail++] = v;
			}
		}
	}

	free(cnt);
	free(ptr);
	free(ato);
	free(aed);
	free(vis);
	free(que);

	return ncomp;
}


// ---- 段階 4 : 離散勾配と節点補助行列 (Hiptmair 前処理) ----

// G : 節点ベクトル -> 辺ベクトル  (u_e = φ_hi - φ_lo)
void edge_grad(const double *phi, double *u)
{
	for (int e = 0; e < NEdge; e++) {
		u[e] = phi[EdgeTo[e]] - phi[EdgeFrom[e]];
	}
}


// G^T : 辺ベクトル -> 節点ベクトル
void edge_gradT(const double *u, double *c)
{
	memset(c, 0, (size_t)NNode * sizeof(double));
	for (int e = 0; e < NEdge; e++) {
		c[EdgeTo[e]]   += u[e];
		c[EdgeFrom[e]] -= u[e];
	}
}


// 節点補助行列 N = G^T A G (節点の CRS パターンに入る)
//
// 辺 e, f が同じ四面体に属するとき、その端点はすべて同じ四面体の節点なので
// 4 つの寄与先はいずれも節点隣接に含まれる。
void edge_nodal_aux(const crs_t *A, crs_t *N)
{
	crs_alloc_tet(N);			// 節点の隣接パターン

	for (int e = 0; e < NEdge; e++) {
		const int32_t eh = EdgeTo[e], el = EdgeFrom[e];
		for (int64_t p = A->rowptr[e]; p < A->rowptr[e + 1]; p++) {
			const int32_t f = A->col[p];
			const double v = A->val[p];
			if (v == 0) continue;
			const int32_t fh = EdgeTo[f], fl = EdgeFrom[f];
			const int32_t row[4] = {eh, eh, el, el};
			const int32_t col[4] = {fh, fl, fh, fl};
			const double  sgn[4] = {+1, -1, -1, +1};
			for (int q = 0; q < 4; q++) {
				int64_t lo = N->rowptr[row[q]];
				int64_t hi = N->rowptr[row[q] + 1] - 1;
				int64_t pos = -1;
				while (lo <= hi) {
					const int64_t mid = (lo + hi) / 2;
					if      (N->col[mid] < col[q]) lo = mid + 1;
					else if (N->col[mid] > col[q]) hi = mid - 1;
					else { pos = mid; break; }
				}
				if (pos >= 0) N->val[pos] += sgn[q] * v;
			}
		}
	}
}


// ---- 辺系の PCG (Jacobi / Hiptmair 前処理) ----
//
// precond = 0 : 対角 (Jacobi)
//         = 1 : Hiptmair (辺の対角 + 節点空間の補正)
//
//   M^-1 r = D_e^-1 r + G [ N^-1 (G^T r) ]
//
// 回転回転演算子の悪条件は「勾配に近い成分」が原因なので、その成分を節点空間へ
// 移して処理する。節点空間の解は N の Jacobi 多項式 (固定回数) で近似するので、
// 前処理全体は固定の対称正定値作用素になり CG の前提を壊さない。
int solver_cg_edge(const crs_t *A, const crs_t *N, const double *b, double *x,
	const unsigned char *fix, int precond, int maxiter, int nout, double converg,
	FILE *fp_log, const char *label)
{
	const int n = (int)A->n;

	double *r = (double *)malloc((size_t)n * sizeof(double));
	double *rp = (double *)malloc((size_t)n * sizeof(double));
	double *p = (double *)malloc((size_t)n * sizeof(double));
	double *q = (double *)malloc((size_t)n * sizeof(double));
	double *z = (double *)malloc((size_t)n * sizeof(double));
	double *d = (double *)malloc((size_t)n * sizeof(double));
	double *nc = NULL, *ny = NULL, *gy = NULL;
	unsigned char *nfix = NULL;

	crs_diag(A, d);
	for (int i = 0; i < n; i++) {
		if (fix[i] || (d[i] <= 0)) d[i] = 1;
	}
	if (precond == 1) {
		nc = (double *)malloc((size_t)NNode * sizeof(double));
		ny = (double *)malloc((size_t)NNode * sizeof(double));
		gy = (double *)malloc((size_t)n * sizeof(double));
		// 節点補助系は定数 (勾配の定数分) だけ不定なので 1 点を固定する
		nfix = (unsigned char *)malloc((size_t)NNode * sizeof(unsigned char));
		memset(nfix, 0, (size_t)NNode * sizeof(unsigned char));
		nfix[0] = 1;
	}

	for (int i = 0; i < n; i++) {
		x[i] = 0;
		r[i] = (fix[i] ? 0 : b[i]);
		rp[i] = 0;
	}
	double bn = 0;
	for (int i = 0; i < n; i++) bn += r[i] * r[i];
	bn = sqrt(bn);
	if (bn <= 0) {
		free(r); free(rp); free(p); free(q); free(z); free(d);
		free(nc); free(ny); free(gy); free(nfix);
		return 0;
	}

	// 前処理の適用 : z = D_e^-1 r + G [ N^-1 (G^T r) ]
	// 節点空間の解は内側 Jacobi-PCG で近似する。近似度が反復毎に変わるので
	// 外側は flexible CG (Polak-Ribiere 型の β) にして前提を壊さないようにする。
	#define APPLY_M                                                            \
		do {                                                                   \
			for (int i = 0; i < n; i++) z[i] = (fix[i] ? 0 : r[i] / d[i]);     \
			if (precond == 1) {                                                \
				edge_gradT(r, nc);                                             \
				nc[0] = 0;                                                     \
				solver_cg(N, nc, ny, nfix, 2000, 0, 1e-12, NULL, "aux");         \
				edge_grad(ny, gy);                                             \
				for (int i = 0; i < n; i++) {                                  \
					if (!fix[i]) z[i] += gy[i];                                \
				}                                                              \
			}                                                                  \
		} while (0)

	APPLY_M;
	for (int i = 0; i < n; i++) p[i] = z[i];
	double rz = 0;
	for (int i = 0; i < n; i++) rz += r[i] * z[i];

	int iter = 0;
	double resid = 1;
	int converged = 0;
	for (iter = 1; iter <= maxiter; iter++) {
		crs_spmv(A, p, q, fix);
		double pq = 0;
		for (int i = 0; i < n; i++) pq += p[i] * q[i];
		if (pq <= 0) break;
		const double alpha = rz / pq;
		for (int i = 0; i < n; i++) {
			x[i] += alpha * p[i];
			r[i] -= alpha * q[i];
		}
		double rr = 0;
		for (int i = 0; i < n; i++) rr += r[i] * r[i];
		resid = sqrt(rr) / bn;
		if ((fp_log != NULL) && (nout > 0) && ((iter % nout) == 0)) {
			fprintf(fp_log, "  %-12s %8d %13.5e\n", label, iter, resid);
		}
		if (resid < converg) {
			converged = 1;
			break;
		}
		for (int i = 0; i < n; i++) rp[i] = r[i] + (alpha * q[i]);	// 直前の r
		APPLY_M;
		double rzn = 0, ryn = 0;
		for (int i = 0; i < n; i++) {
			rzn += r[i] * z[i];
			ryn += (r[i] - rp[i]) * z[i];		// flexible CG (Polak-Ribiere)
		}
		const double beta = ryn / rz;
		rz = rzn;
		for (int i = 0; i < n; i++) p[i] = z[i] + (beta * p[i]);
	}
	#undef APPLY_M

	free(r); free(rp); free(p); free(q); free(z); free(d);
	free(nc); free(ny); free(gy); free(nfix);

	return (converged ? iter : -iter);
}
