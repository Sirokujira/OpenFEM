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

/*
六面体・角柱の局所辺 (Gmsh の 2 次要素の辺の並びと同じ表。fem.h 参照)。
どの並びでも「小さい節点番号 -> 大きい方」の符号規約が向きを吸収する。

**種別の混在を許す** (四面体 + 角柱 / 六面体 + 角柱)。最低次 Nedelec の
接線トレースは三角形面で三角形の Whitney 空間、四角形面で四角形の最低次
Nedelec 空間になり、どちらも面の 3 (4) 辺の自由度だけで決まるので、
面を共有する要素の種別が違っても接線連続性が保たれる
(solve_edge_test の (g) がこれを直接検査する)。
六面体と四面体は面を共有できないので、その組は間に角柱が要る
(ピラミッドは多項式の Nedelec 基底が無いので setup_unstruct() で弾く)。
*/
static const int EDGE_HEX[12][2] = {
	{0, 1}, {0, 3}, {0, 4}, {1, 2}, {1, 5}, {2, 3},
	{2, 6}, {3, 7}, {4, 5}, {4, 7}, {5, 6}, {6, 7}};
static const int EDGE_PRISM[9][2] = {
	{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 4}, {2, 5}, {3, 4}, {3, 5}, {4, 5}};

// 六面体の頂点の参照座標 (Gmsh の並び。unstruct.c の HEX_SGN と同じ表)
static const signed char EHEX_SGN[8][3] = {
	{-1, -1, -1}, {+1, -1, -1}, {+1, +1, -1}, {-1, +1, -1},
	{-1, -1, +1}, {+1, -1, +1}, {+1, +1, +1}, {-1, +1, +1}};


/*
辺要素で扱う要素の見方。**通し番号は elem3d と同じ「四面体 -> 六面体 ->
角柱」の連番**で、種別ごとの添字は各アクセサの中で引き直す。純粋な格子では
連番が種別内の添字と一致するので、既存の答えはビット単位で変わらない。

(ピラミッドは Nedelec 基底が無く setup_unstruct() で弾かれるので、ここに
来る格子には含まれない)
*/
int edge_elem_count(void)
{
	return (NTet + NHex + NPrism);
}

int edge_elem_kind(int e)
{
	if (e < NTet) return MESHELEM_TET;
	if (e < (NTet + NHex)) return MESHELEM_HEX;

	return MESHELEM_PRISM;
}

int edge_elem_nedge(int e)
{
	const int kind = edge_elem_kind(e);

	if (kind == MESHELEM_HEX) return 12;
	if (kind == MESHELEM_PRISM) return 9;

	return 6;
}

int edge_elem_nen(int e)
{
	const int kind = edge_elem_kind(e);

	if (kind == MESHELEM_HEX) return 8;
	if (kind == MESHELEM_PRISM) return 6;

	return 4;
}

const int32_t *edge_elem_nodes(int e)
{
	if (e < NTet) return &Tet[e * 4];
	if (e < (NTet + NHex)) return &Hex[(e - NTet) * 8];

	return &Prism[(e - NTet - NHex) * 6];
}

const int (*edge_elem_table(int e))[2]
{
	const int kind = edge_elem_kind(e);

	if (kind == MESHELEM_HEX) return EDGE_HEX;
	if (kind == MESHELEM_PRISM) return EDGE_PRISM;

	return EDGE_NODE;
}

int edge_elem_mat(int e)
{
	if (e < NTet) return ((TetMat != NULL) ? TetMat[e] : 0);
	if (e < (NTet + NHex)) return ((HexMat != NULL) ? HexMat[e - NTet] : 0);

	return ((PrismMat != NULL) ? PrismMat[e - NTet - NHex] : 0);
}


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
	const int nelem = edge_elem_count();

	// 要素毎の辺の開始位置 (種別で辺数が違うので通し番号では引けない)
	EdgeOff = (int64_t *)malloc(((size_t)nelem + 1) * sizeof(int64_t));
	EdgeOff[0] = 0;
	for (int e = 0; e < nelem; e++) EdgeOff[e + 1] = EdgeOff[e] + edge_elem_nedge(e);

	// 節点毎の (i < j) 隣接候補を数える
	int *cnt = (int *)malloc((size_t)n * sizeof(int));
	memset(cnt, 0, (size_t)n * sizeof(int));
	for (int e = 0; e < nelem; e++) {
		const int32_t *nd = edge_elem_nodes(e);
		const int nedge = edge_elem_nedge(e);
		const int (*tab)[2] = edge_elem_table(e);
		for (int k = 0; k < nedge; k++) {
			const int32_t a = nd[tab[k][0]];
			const int32_t b = nd[tab[k][1]];
			cnt[(a < b) ? a : b]++;
		}
	}

	int64_t *ptr = (int64_t *)malloc(((size_t)n + 1) * sizeof(int64_t));
	ptr[0] = 0;
	for (int i = 0; i < n; i++) ptr[i + 1] = ptr[i] + cnt[i];
	int32_t *lst = (int32_t *)malloc((size_t)ptr[n] * sizeof(int32_t));
	memset(cnt, 0, (size_t)n * sizeof(int));
	for (int e = 0; e < nelem; e++) {
		const int32_t *nd = edge_elem_nodes(e);
		const int nedge = edge_elem_nedge(e);
		const int (*tab)[2] = edge_elem_table(e);
		for (int k = 0; k < nedge; k++) {
			const int32_t a = nd[tab[k][0]];
			const int32_t b = nd[tab[k][1]];
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

	// 要素毎の辺番号と向き (EdgeOff で引く。純粋な格子では従来と同じ並び)
	TetEdge = (int32_t *)malloc((size_t)EdgeOff[nelem] * sizeof(int32_t));
	TetEdgeSgn = (signed char *)malloc((size_t)EdgeOff[nelem] * sizeof(signed char));
	for (int e = 0; e < nelem; e++) {
		const int32_t *nd = edge_elem_nodes(e);
		const int nedge = edge_elem_nedge(e);
		const int (*tab)[2] = edge_elem_table(e);
		const int64_t off = EdgeOff[e];
		for (int k = 0; k < nedge; k++) {
			const int32_t a = nd[tab[k][0]];
			const int32_t b = nd[tab[k][1]];
			TetEdge[off + k] = (int32_t)edge_id(a, b);
			TetEdgeSgn[off + k] = (signed char)((a < b) ? +1 : -1);
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
	free(EdgeOff);
	EdgePtr = NULL;
	EdgeTo = NULL;
	EdgeFrom = NULL;
	TetEdge = NULL;
	TetEdgeSgn = NULL;
	EdgeOff = NULL;
	NEdge = 0;
}


// ---- 辺ベース CRS ----

void crs_alloc_edge(crs_t *A)
{
	const int ne = NEdge;
	const int nelem = edge_elem_count();

	// 辺毎の要素リスト
	int *cnt = (int *)malloc((size_t)ne * sizeof(int));
	memset(cnt, 0, (size_t)ne * sizeof(int));
	for (int e = 0; e < nelem; e++) {
		for (int64_t k = EdgeOff[e]; k < EdgeOff[e + 1]; k++) cnt[TetEdge[k]]++;
	}
	int64_t *ptr = (int64_t *)malloc(((size_t)ne + 1) * sizeof(int64_t));
	ptr[0] = 0;
	for (int i = 0; i < ne; i++) ptr[i + 1] = ptr[i] + cnt[i];
	int32_t *lst = (int32_t *)malloc((size_t)ptr[ne] * sizeof(int32_t));
	memset(cnt, 0, (size_t)ne * sizeof(int));
	for (int e = 0; e < nelem; e++) {
		for (int64_t k = EdgeOff[e]; k < EdgeOff[e + 1]; k++) {
			const int32_t ed = TetEdge[k];
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
		int need = 0;
		for (int64_t p = p0; p < p1; p++) need += edge_elem_nedge(lst[p]);
		if (need > cap) {
			cap = need;
			work = (int32_t *)realloc(work, (size_t)cap * sizeof(int32_t));
		}
		int m = 0;
		for (int64_t p = p0; p < p1; p++) {
			const int32_t t = lst[p];
			for (int64_t k = EdgeOff[t]; k < EdgeOff[t + 1]; k++) work[m++] = TetEdge[k];
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
			for (int64_t k = EdgeOff[t]; k < EdgeOff[t + 1]; k++) work[m++] = TetEdge[k];
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


/*
六面体・角柱の辺要素 (最低次 Nedelec、等パラメトリック)。

参照要素の辺基底 Ŵ と参照回転 ĉ = ∇̂×Ŵ を作り、共変 Piola 変換で
物理空間に写して数値積分する:

	W = J^{-T} Ŵ,   ∇×W = (J ĉ) / det J
	S_kl = ∫ (∇×W_k)^T ν (∇×W_l) |det J| dξ
	T_kl = ∫ σ W_k・W_l |det J| dξ

**基底の向きは局所辺表の (a -> b)** に合わせる (参照座標で a から b へ)。
全体の向き (小さい節点番号 -> 大きい方) との差は符号 sgn が吸収する —
この規約は四面体と同じで、落とすと隣接要素間で基底が食い違う。

六面体の参照基底 (辺は軸 d に沿い、横方向の座標 s1, s2 が固定):

	Ŵ = dir * ê_d (1 + s1 t1)(1 + s2 t2) / 8    (dir = 参照で a -> b の向き)

辺に沿う線積分が 1 (単位循環)。角柱は三角形の 2 次元 Whitney と 1 次元の
線形関数の積 (水平辺)、および λ_i ê_ζ / 2 (鉛直辺)。

積分は六面体 2x2x2、角柱 三角形 3 点 x ζ 2 点。被積分関数は各方向 2 次まで
なので**アフィンな要素 (平行六面体 / 鉛直に押し出した角柱) では厳密**。
一般の形では近似になる (E の自己検証の恒等式もアフィン要素を前提とする —
{a + b×r} の Nedelec 補間はアフィンでない要素では厳密でないため)。
*/

static double edet3(const double jm[3][3])
{
	return (jm[0][0] * ((jm[1][1] * jm[2][2]) - (jm[1][2] * jm[2][1])))
	     - (jm[0][1] * ((jm[1][0] * jm[2][2]) - (jm[1][2] * jm[2][0])))
	     + (jm[0][2] * ((jm[1][0] * jm[2][1]) - (jm[1][1] * jm[2][0])));
}


static void einv3(const double jm[3][3], double d, double ji[3][3])
{
	ji[0][0] = ((jm[1][1] * jm[2][2]) - (jm[1][2] * jm[2][1])) / d;
	ji[0][1] = ((jm[0][2] * jm[2][1]) - (jm[0][1] * jm[2][2])) / d;
	ji[0][2] = ((jm[0][1] * jm[1][2]) - (jm[0][2] * jm[1][1])) / d;
	ji[1][0] = ((jm[1][2] * jm[2][0]) - (jm[1][0] * jm[2][2])) / d;
	ji[1][1] = ((jm[0][0] * jm[2][2]) - (jm[0][2] * jm[2][0])) / d;
	ji[1][2] = ((jm[0][2] * jm[1][0]) - (jm[0][0] * jm[1][2])) / d;
	ji[2][0] = ((jm[1][0] * jm[2][1]) - (jm[1][1] * jm[2][0])) / d;
	ji[2][1] = ((jm[0][1] * jm[2][0]) - (jm[0][0] * jm[2][1])) / d;
	ji[2][2] = ((jm[0][0] * jm[1][1]) - (jm[0][1] * jm[1][0])) / d;
}


// 六面体 (三重線形) のヤコビアン J[i][j] = ∂x_i/∂ξ_j
static void ehex_jac(const int32_t *nd, double x, double y, double z, double jm[3][3])
{
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) jm[i][j] = 0;
	}
	for (int a = 0; a < 8; a++) {
		const double sx = EHEX_SGN[a][0], sy = EHEX_SGN[a][1], sz = EHEX_SGN[a][2];
		const double dn[3] = {sx * (1 + (sy * y)) * (1 + (sz * z)) / 8,
		                      sy * (1 + (sx * x)) * (1 + (sz * z)) / 8,
		                      sz * (1 + (sx * x)) * (1 + (sy * y)) / 8};
		const int32_t v = nd[a];
		const double p[3] = {Xp[v], Yp[v], Zp[v]};
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) jm[i][j] += p[i] * dn[j];
		}
	}
}


// 角柱 (線形) のヤコビアン ((u, v, ζ)、λ = (1-u-v, u, v)、ζ ∈ [-1,1])
static void eprism_jac(const int32_t *nd, double u, double v, double ze, double jm[3][3])
{
	const double lam[3] = {1 - u - v, u, v};
	static const signed char dlu[3] = {-1, 1, 0};
	static const signed char dlv[3] = {-1, 0, 1};

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) jm[i][j] = 0;
	}
	for (int a = 0; a < 6; a++) {
		const int t = a % 3;
		const double szn = ((a < 3) ? -1.0 : 1.0);
		const double h = (1 + (szn * ze)) / 2;
		const double dn[3] = {dlu[t] * h, dlv[t] * h, lam[t] * szn / 2};
		const int32_t w = nd[a];
		const double p[3] = {Xp[w], Yp[w], Zp[w]};
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) jm[i][j] += p[i] * dn[j];
		}
	}
}


// 六面体の参照辺基底 Ŵ と参照回転 ĉ
static void ehex_basis(double x, double y, double z, double w[12][3], double c[12][3])
{
	const double t[3] = {x, y, z};

	for (int k = 0; k < 12; k++) {
		const int a = EDGE_HEX[k][0], b = EDGE_HEX[k][1];
		// 辺の軸 d (a と b で座標が違う軸) と向き dir
		int d = 0;
		for (int i = 0; i < 3; i++) {
			if (EHEX_SGN[a][i] != EHEX_SGN[b][i]) d = i;
		}
		const double dir = (EHEX_SGN[b][d] - EHEX_SGN[a][d]) / 2.0;
		const int t1 = (d + 1) % 3, t2 = (d + 2) % 3;
		const double s1 = EHEX_SGN[a][t1], s2 = EHEX_SGN[a][t2];

		const double f = dir * (1 + (s1 * t[t1])) * (1 + (s2 * t[t2])) / 8;
		// ∇f (f は t1, t2 だけの関数)
		double gf[3] = {0, 0, 0};
		gf[t1] = dir * s1 * (1 + (s2 * t[t2])) / 8;
		gf[t2] = dir * s2 * (1 + (s1 * t[t1])) / 8;

		for (int i = 0; i < 3; i++) w[k][i] = 0;
		w[k][d] = f;
		// ∇×(f ê_d) = ∇f × ê_d
		double ed[3] = {0, 0, 0};
		ed[d] = 1;
		c[k][0] = (gf[1] * ed[2]) - (gf[2] * ed[1]);
		c[k][1] = (gf[2] * ed[0]) - (gf[0] * ed[2]);
		c[k][2] = (gf[0] * ed[1]) - (gf[1] * ed[0]);
	}
}


// 角柱の参照辺基底 Ŵ と参照回転 ĉ
static void eprism_basis(double u, double v, double ze, double w[9][3], double c[9][3])
{
	const double lam[3] = {1 - u - v, u, v};
	static const signed char dlu[3] = {-1, 1, 0};
	static const signed char dlv[3] = {-1, 0, 1};

	for (int k = 0; k < 9; k++) {
		const int a = EDGE_PRISM[k][0], b = EDGE_PRISM[k][1];
		if ((b - a) == 3) {
			// 鉛直辺 (i, i+3) : Ŵ = (λ_i / 2) ê_ζ (ζ ∈ [-1,1] で循環 1)
			const int i = a;
			w[k][0] = 0;
			w[k][1] = 0;
			w[k][2] = lam[i] / 2;
			// ∇×(0,0,g) = (∂g/∂v, -∂g/∂u, 0)
			c[k][0] = dlv[i] / 2.0;
			c[k][1] = -dlu[i] / 2.0;
			c[k][2] = 0;
		}
		else {
			// 水平辺 : 2 次元 Whitney x 線形。下面 (a,b < 3) は s = -1、上面は +1
			const int i = a % 3, j = b % 3;
			const double szn = ((a < 3) ? -1.0 : 1.0);
			const double h = (1 + (szn * ze)) / 2;
			const double wu = (lam[i] * dlu[j]) - (lam[j] * dlu[i]);
			const double wv = (lam[i] * dlv[j]) - (lam[j] * dlv[i]);
			const double c2 = 2.0 * ((dlu[i] * dlv[j]) - (dlv[i] * dlu[j]));
			w[k][0] = wu * h;
			w[k][1] = wv * h;
			w[k][2] = 0;
			// ∇×(Wu h, Wv h, 0) = (-Wv h', Wu h', curl2d h),  h' = s/2
			c[k][0] = -wv * szn / 2;
			c[k][1] = wu * szn / 2;
			c[k][2] = c2 * h;
		}
	}
}


// 参照点 1 点の寄与を se / te に足す (共変 Piola 変換)
static void edge_piola_point(const int32_t *nd, int nedge, int hex,
	double x, double y, double z, double wq,
	const double nu[6], double sig, double se[12][12], double te[12][12])
{
	double wr[12][3], cr[12][3], jm[3][3], ji[3][3];

	if (hex) {
		ehex_basis(x, y, z, wr, cr);
		ehex_jac(nd, x, y, z, jm);
	}
	else {
		eprism_basis(x, y, z, wr, cr);
		eprism_jac(nd, x, y, z, jm);
	}
	const double det = edet3(jm);
	if (det == 0) return;
	einv3(jm, det, ji);

	// 物理基底 W_i = Σ_j ji[j][i] Ŵ_j、物理回転 C_i = Σ_j jm[i][j] ĉ_j / det
	double wp[12][3], cp[12][3];
	for (int k = 0; k < nedge; k++) {
		for (int i = 0; i < 3; i++) {
			wp[k][i] = (ji[0][i] * wr[k][0]) + (ji[1][i] * wr[k][1]) + (ji[2][i] * wr[k][2]);
			cp[k][i] = ((jm[i][0] * cr[k][0]) + (jm[i][1] * cr[k][1])
			          + (jm[i][2] * cr[k][2])) / det;
		}
	}

	const double dw = wq * fabs(det);
	for (int k = 0; k < nedge; k++) {
		for (int l = 0; l < nedge; l++) {
			const double cn = (nu[0] * cp[k][0] * cp[l][0])
			                + (nu[1] * cp[k][1] * cp[l][1])
			                + (nu[2] * cp[k][2] * cp[l][2])
			                + (nu[3] * ((cp[k][0] * cp[l][1]) + (cp[k][1] * cp[l][0])))
			                + (nu[4] * ((cp[k][1] * cp[l][2]) + (cp[k][2] * cp[l][1])))
			                + (nu[5] * ((cp[k][2] * cp[l][0]) + (cp[k][0] * cp[l][2])));
			se[k][l] += dw * cn;
			te[k][l] += dw * sig * ((wp[k][0] * wp[l][0]) + (wp[k][1] * wp[l][1])
			                      + (wp[k][2] * wp[l][2]));
		}
	}
}


// 要素 e の回転回転行列 se と質量行列 te (種別で分岐。四面体は閉形式のまま)
void edge_elem_matrices(int e, const double nu[6], double sig,
	double se[12][12], double te[12][12])
{
	const int kind = edge_elem_kind(e);

	if (kind == MESHELEM_TET) {
		double s6[6][6], t6[6][6];
		edge_element(e, nu, sig, s6, t6);
		for (int k = 0; k < 6; k++) {
			for (int l = 0; l < 6; l++) {
				se[k][l] = s6[k][l];
				te[k][l] = t6[k][l];
			}
		}
		return;
	}

	for (int k = 0; k < 12; k++) {
		for (int l = 0; l < 12; l++) se[k][l] = te[k][l] = 0;
	}

	if (kind == MESHELEM_HEX) {
		const double gp = 1.0 / sqrt(3.0);
		const int32_t *nd = &Hex[(e - NTet) * 8];
		for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
		for (int k = 0; k < 2; k++) {
			edge_piola_point(nd, 12, 1, (i ? gp : -gp), (j ? gp : -gp),
				(k ? gp : -gp), 1.0, nu, sig, se, te);
		}
		}
		}
		return;
	}

	{
		// 三角形 3 点 (次数 2) x ζ 2 点 Gauss
		static const double tu[3] = {1.0 / 6, 2.0 / 3, 1.0 / 6};
		static const double tv[3] = {1.0 / 6, 1.0 / 6, 2.0 / 3};
		const double gp = 1.0 / sqrt(3.0);
		const int32_t *nd = &Prism[(e - NTet - NHex) * 6];
		for (int q = 0; q < 3; q++) {
			for (int k = 0; k < 2; k++) {
				edge_piola_point(nd, 9, 0, tu[q], tv[q], (k ? gp : -gp),
					1.0 / 6, nu, sig, se, te);
			}
		}
	}
}


// 要素中心での物理基底 W と ∇×W (場の出力用。四面体以外)
void edge_elem_center(int e, double w[12][3], double c[12][3])
{
	double wr[12][3], cr[12][3], jm[3][3], ji[3][3];
	const int hex = (edge_elem_kind(e) == MESHELEM_HEX);
	const int nedge = edge_elem_nedge(e);

	if (hex) {
		ehex_basis(0, 0, 0, wr, cr);
		ehex_jac(&Hex[(e - NTet) * 8], 0, 0, 0, jm);
	}
	else {
		eprism_basis(1.0 / 3, 1.0 / 3, 0, wr, cr);
		eprism_jac(&Prism[(e - NTet - NHex) * 6], 1.0 / 3, 1.0 / 3, 0, jm);
	}
	const double det = edet3(jm);
	if (det == 0) {
		for (int k = 0; k < nedge; k++) {
			for (int i = 0; i < 3; i++) w[k][i] = c[k][i] = 0;
		}
		return;
	}
	einv3(jm, det, ji);
	for (int k = 0; k < nedge; k++) {
		for (int i = 0; i < 3; i++) {
			w[k][i] = (ji[0][i] * wr[k][0]) + (ji[1][i] * wr[k][1]) + (ji[2][i] * wr[k][2]);
			c[k][i] = ((jm[i][0] * cr[k][0]) + (jm[i][1] * cr[k][1])
			         + (jm[i][2] * cr[k][2])) / det;
		}
	}
}


/*
参照点での**物理**基底 W (Piola 変換したもの)。接線連続性の検査 (g) で使う。
ref は種別ごとの参照座標:
  四面体 : (λ1, λ2, λ3)  (λ0 = 1 - λ1 - λ2 - λ3)
  六面体 : (ξ, η, ζ) ∈ [-1,1]^3
  角柱   : (u, v, ζ)     (λ = (1-u-v, u, v)、ζ ∈ [-1,1])
四面体は Whitney 基底が物理座標の閉形式で書けるので Piola を通さない。
*/
static void edge_basis_phys(int e, const double ref[3], double w[12][3])
{
	const int kind = edge_elem_kind(e);
	const int nedge = edge_elem_nedge(e);

	for (int k = 0; k < 12; k++) {
		for (int i = 0; i < 3; i++) w[k][i] = 0;
	}

	if (kind == MESHELEM_TET) {
		double g[4][3], vol;
		if (tet_grad_pub(&Tet[e * 4], g, &vol)) return;
		const double lam[4] = {1 - ref[0] - ref[1] - ref[2], ref[0], ref[1], ref[2]};
		for (int k = 0; k < 6; k++) {
			const int a = EDGE_NODE[k][0], b = EDGE_NODE[k][1];
			for (int i = 0; i < 3; i++) {
				w[k][i] = (lam[a] * g[b][i]) - (lam[b] * g[a][i]);
			}
		}
		return;
	}

	double wr[12][3], jm[3][3], ji[3][3];
	if (kind == MESHELEM_HEX) {
		double cr[12][3];
		ehex_basis(ref[0], ref[1], ref[2], wr, cr);
		ehex_jac(&Hex[(e - NTet) * 8], ref[0], ref[1], ref[2], jm);
	}
	else {
		double wp[9][3], cp[9][3];
		eprism_basis(ref[0], ref[1], ref[2], wp, cp);
		for (int k = 0; k < 9; k++) {
			for (int i = 0; i < 3; i++) wr[k][i] = wp[k][i];
		}
		eprism_jac(&Prism[(e - NTet - NHex) * 6], ref[0], ref[1], ref[2], jm);
	}
	const double det = edet3(jm);
	if (det == 0) return;
	einv3(jm, det, ji);
	for (int k = 0; k < nedge; k++) {
		for (int i = 0; i < 3; i++) {
			w[k][i] = (ji[0][i] * wr[k][0]) + (ji[1][i] * wr[k][1]) + (ji[2][i] * wr[k][2]);
		}
	}
}


// 局所節点の参照座標 (edge_basis_phys の ref と同じ系)
static void edge_ref_node(int kind, int l, double ref[3])
{
	static const signed char TET_REF[4][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
	static const signed char PRISM_REF[6][3] = {
		{0, 0, -1}, {1, 0, -1}, {0, 1, -1}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}};

	for (int i = 0; i < 3; i++) {
		ref[i] = ((kind == MESHELEM_HEX) ? EHEX_SGN[l][i]
		        : (kind == MESHELEM_PRISM) ? PRISM_REF[l][i] : TET_REF[l][i]);
	}
}


// 局所の面 (三角形は 4 番目が -1)。並びは Gmsh の節点順に対応する
static const int FACE_TET[4][4] = {{0, 1, 2, -1}, {0, 1, 3, -1}, {0, 2, 3, -1}, {1, 2, 3, -1}};
static const int FACE_HEX[6][4] = {
	{0, 1, 2, 3}, {4, 5, 6, 7}, {0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}};
static const int FACE_PRISM[5][4] = {
	{0, 1, 2, -1}, {3, 4, 5, -1}, {0, 1, 4, 3}, {1, 2, 5, 4}, {2, 0, 3, 5}};

static int edge_face_table(int kind, const int (**tab)[4])
{
	if (kind == MESHELEM_HEX) { *tab = FACE_HEX; return 6; }
	if (kind == MESHELEM_PRISM) { *tab = FACE_PRISM; return 5; }
	*tab = FACE_TET;

	return 4;
}


// 面の記録 (整列した節点をキーにして共有面を見つける)
typedef struct {
	int32_t key[4];		// 昇順の全体節点番号 (三角形は key[3] = -1)
	int32_t elem;
	int16_t face;
} eface_t;

static int cmp_eface(const void *a, const void *b)
{
	const eface_t *x = (const eface_t *)a;
	const eface_t *y = (const eface_t *)b;

	for (int i = 0; i < 4; i++) {
		if (x->key[i] != y->key[i]) return ((x->key[i] < y->key[i]) ? -1 : 1);
	}

	return 0;
}


/*
(g) 面をまたぐ接線連続性 (Nedelec の適合性)。

**これは (b) (c) (c2) の恒等式では検出できない。** {a + b×r} の場は要素毎に
厳密に補間されるので、隣の要素とトレースが食い違っていてもエネルギーは
要素積分の和として合ってしまう (実測: 角柱の鉛直辺の基底を 1/2 倍しても
(a)〜(f) は全部通り、この検査だけが落ちる)。

任意の辺自由度ベクトル u に対し、面を共有する 2 要素の補間場
  u_h(x) = Σ_k sgn_k u_k W_k(x)
の**接線成分**が面上の各点で一致することを直接見る。同じ点を両側で作るために、
面の全体節点に重み (三角形は面積座標、四角形は双 1 次) を割り当て、
各要素ではその全体節点が入っている局所節点の参照座標を重み付けする。
等パラメトリック写像は面上で (双) 1 次なので、どちらの要素からも同じ物理点になる。

種別をまたぐ面 (四面体 - 角柱の三角形面、六面体 - 角柱の四角形面) が
この検査の本命だが、同じ種別どうしの面も一緒に見る (符号規約の検査になる)。
*/
static int edge_face_conform(FILE *fp_log, const double *u)
{
	const int nelem = edge_elem_count();

	// 全ての面を集めて整列し、2 回現れる面 (内部面) を取り出す
	int64_t nf = 0;
	for (int e = 0; e < nelem; e++) {
		const int (*ft)[4];
		nf += edge_face_table(edge_elem_kind(e), &ft);
	}
	eface_t *fc = (eface_t *)malloc((size_t)nf * sizeof(eface_t));
	int64_t m = 0;
	for (int e = 0; e < nelem; e++) {
		const int kind = edge_elem_kind(e);
		const int32_t *nd = edge_elem_nodes(e);
		const int (*ft)[4];
		const int nfa = edge_face_table(kind, &ft);
		for (int f = 0; f < nfa; f++) {
			int32_t k4[4] = {-1, -1, -1, -1};
			int nv = 0;
			for (int i = 0; i < 4; i++) {
				if (ft[f][i] >= 0) k4[nv++] = nd[ft[f][i]];
			}
			for (int i = 0; i < nv; i++) {			// 小さい方から (nv <= 4)
				for (int j = i + 1; j < nv; j++) {
					if (k4[j] < k4[i]) {
						const int32_t t = k4[i]; k4[i] = k4[j]; k4[j] = t;
					}
				}
			}
			for (int i = 0; i < 4; i++) fc[m].key[i] = k4[i];
			fc[m].elem = (int32_t)e;
			fc[m].face = (int16_t)f;
			m++;
		}
	}
	qsort(fc, (size_t)m, sizeof(eface_t), cmp_eface);

	// 面上の標本点の重み (三角形 : 面積座標、四角形 : 双 1 次の (s, t))
	static const double TW[3][3] = {
		{0.2, 0.3, 0.5}, {0.5, 0.25, 0.25}, {0.1, 0.6, 0.3}};
	static const double QST[3][2] = {{-0.3, 0.4}, {0.6, -0.2}, {0.1, 0.8}};

	double jmax = 0, vmax = 0;
	int64_t nshare = 0, nmix = 0;
	for (int64_t p = 0; p + 1 < m; p++) {
		if (cmp_eface(&fc[p], &fc[p + 1]) != 0) continue;
		const int nv = ((fc[p].key[3] >= 0) ? 4 : 3);
		const int e0 = fc[p].elem, e1 = fc[p + 1].elem;
		nshare++;
		if (edge_elem_kind(e0) != edge_elem_kind(e1)) nmix++;

		// 面の全体節点 (要素 e0 の面の並び) と法線
		const int (*ft0)[4];
		edge_face_table(edge_elem_kind(e0), &ft0);
		const int32_t *nd0 = edge_elem_nodes(e0);
		int32_t gn[4];
		for (int i = 0; i < nv; i++) gn[i] = nd0[ft0[fc[p].face][i]];
		const double a1[3] = {Xp[gn[1]] - Xp[gn[0]], Yp[gn[1]] - Yp[gn[0]], Zp[gn[1]] - Zp[gn[0]]};
		const double a2[3] = {Xp[gn[2]] - Xp[gn[0]], Yp[gn[2]] - Yp[gn[0]], Zp[gn[2]] - Zp[gn[0]]};
		double nv3[3] = {(a1[1] * a2[2]) - (a1[2] * a2[1]),
		                 (a1[2] * a2[0]) - (a1[0] * a2[2]),
		                 (a1[0] * a2[1]) - (a1[1] * a2[0])};
		const double nn = sqrt((nv3[0] * nv3[0]) + (nv3[1] * nv3[1]) + (nv3[2] * nv3[2]));
		if (nn <= 0) continue;
		for (int i = 0; i < 3; i++) nv3[i] /= nn;

		for (int q = 0; q < 3; q++) {
			double wt[4];
			if (nv == 3) {
				for (int i = 0; i < 3; i++) wt[i] = TW[q][i];
			}
			else {
				const double s = QST[q][0], t = QST[q][1];
				wt[0] = (1 - s) * (1 - t) / 4;
				wt[1] = (1 + s) * (1 - t) / 4;
				wt[2] = (1 + s) * (1 + t) / 4;
				wt[3] = (1 - s) * (1 + t) / 4;
			}

			double side[2][3];
			for (int sd = 0; sd < 2; sd++) {
				const int e = ((sd == 0) ? e0 : e1);
				const int kind = edge_elem_kind(e);
				const int nen = edge_elem_nen(e);
				const int nedge = edge_elem_nedge(e);
				const int32_t *nd = edge_elem_nodes(e);
				// 全体節点 -> この要素の局所節点 -> 参照座標を重み付けする
				double ref[3] = {0, 0, 0};
				for (int i = 0; i < nv; i++) {
					int l = -1;
					for (int a = 0; a < nen; a++) {
						if (nd[a] == gn[i]) l = a;
					}
					if (l < 0) break;
					double rp[3];
					edge_ref_node(kind, l, rp);
					for (int c = 0; c < 3; c++) ref[c] += wt[i] * rp[c];
				}
				double wb[12][3];
				edge_basis_phys(e, ref, wb);
				const int32_t *ed = &TetEdge[EdgeOff[e]];
				const signed char *sg = &TetEdgeSgn[EdgeOff[e]];
				for (int c = 0; c < 3; c++) side[sd][c] = 0;
				for (int k = 0; k < nedge; k++) {
					const double uk = sg[k] * u[ed[k]];
					for (int c = 0; c < 3; c++) side[sd][c] += uk * wb[k][c];
				}
			}

			// 接線成分の差 (法線成分は連続でなくてよい)
			double d[3], dn = 0, vn = 0;
			for (int c = 0; c < 3; c++) d[c] = side[0][c] - side[1][c];
			for (int c = 0; c < 3; c++) dn += d[c] * nv3[c];
			for (int c = 0; c < 3; c++) vn += side[0][c] * nv3[c];
			double dt = 0, vt = 0;
			for (int c = 0; c < 3; c++) {
				const double dd = d[c] - (dn * nv3[c]);
				const double vv = side[0][c] - (vn * nv3[c]);
				dt += dd * dd;
				vt += vv * vv;
			}
			if (sqrt(dt) > jmax) jmax = sqrt(dt);
			if (sqrt(vt) > vmax) vmax = sqrt(vt);
		}
	}
	free(fc);

	const double rel = jmax / ((vmax > 0) ? vmax : 1);
	fprintf(fp_log, "  (g) tangential continuity : %lld shared faces (%lld between "
		"different kinds), max jump / max|Et| = %.3e\n",
		(long long)nshare, (long long)nmix, rel);
	if (rel > 1e-10) {
		fprintf(fp_log, "*** the tangential trace is not continuous across faces; "
			"the edge basis is not curl-conforming\n");
		return 1;
	}

	return 0;
}


// 全体行列の作成 (辺要素)
//   S : 回転回転行列 (ν)、T : 質量行列 (σ)。どちらも NULL 可
void assemble_edge(crs_t *S, crs_t *T)
{
	if (S != NULL) crs_zero(S);
	if (T != NULL) crs_zero(T);

	const int nelem = edge_elem_count();

	for (int e = 0; e < nelem; e++) {
		const int m = edge_elem_mat(e);
		const int nedge = edge_elem_nedge(e);
		double nu[6];
		material_coef_pub(m, 3, nu);			// ν = (μ0 μ~)^-1
		const double sig = Material[m].sigma;

		double se[12][12], te[12][12];
		edge_elem_matrices(e, nu, sig, se, te);

		const int32_t *ed = &TetEdge[EdgeOff[e]];
		const signed char *sg = &TetEdgeSgn[EdgeOff[e]];
		for (int k = 0; k < nedge; k++) {
			for (int l = 0; l < nedge; l++) {
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
// 辺要素の性質を厳密に満たすべき量で検査する。
// いずれも解析解が閉形式で分かるので、要素行列・辺の向き・組み立ての
// どこかが間違っていれば必ず落ちる。
//
//  (a) 勾配は回転回転行列の零空間に入る
//      節点値 φ の勾配は辺自由度 u_e = φ_hi - φ_lo で厳密に表せるので u^T S u = 0
//  (b) 一様場 E0 の質量 : u_e = E0・(p_hi - p_lo) として u^T T u = |E0|^2 Σ σ_e V_e
//  (b') 一様場は回転が 0 なので u^T S u = 0
//  (c) 回転場 E = (1/2)(a × r) : curl E = a なので u^T S u = Σ (a^T ν_e a) V_e
//      a は軸に平行でない向き (1,2,3) にとる。軸に平行だと ν の 1 成分しか
//      比較に効かず、他の対角成分や非対角成分の誤りを見逃す
//  (c2) 同じ回転場の質量 : u^T T u = Σ σ_e ∫|E|^2 dV
//      (b) の一様場は Nedelec 補間が要素内で定数になるのでどんな数値積分でも
//      厳密になってしまう。空間変化する場でないと質量行列の階数落ちを検出できない
//  (d) S と T の対称性
//
// (b) (c) (c2) の場はいずれも最低次 Nedelec 空間 {a + b × r} に含まれるので
// 補間は厳密で、比較は機械精度の恒等式になる。
int solve_edge_test(FILE *fp_log)
{
	int ierr = 0;

	edge_build();

	fprintf(fp_log, "\n=== edge element (Nedelec) self test ===\n");
	fprintf(fp_log, "  nodes = %d, elements = %d (%d tet / %d hex / %d prism), edges = %d\n",
		NNode, edge_elem_count(), NTet, NHex, NPrism, NEdge);
	fflush(fp_log);

	crs_t S, T;
	crs_alloc_edge(&S);
	crs_alloc_edge(&T);
	assemble_edge(&S, &T);
	fprintf(fp_log, "  edge matrix : %lld nonzeros (%.1f per row)\n",
		(long long)S.nnz, (double)S.nnz / ((NEdge > 0) ? NEdge : 1));

	// (c) (c2) で使う回転場 E = (1/2)(av × r) の回転ベクトル。
	// 軸に平行にすると ν の 1 成分しか比較に効かないので斜めにとる
	const double av[3] = {1.0, 2.0, 3.0};

	// 材料の体積積分
	//   svol = Σ σ_e V_e                   ... (b) の厳密値 / |E0|^2
	//   nvol = Σ (av^T ν_e av) V_e         ... (c) の厳密値
	//   mvol = Σ σ_e ∫|E|^2 dV             ... (c2) の厳密値
	// ∫|E|^2 は E が 1 次なので節点値から厳密に出せる:
	//   ∫ λ_a λ_b dV = V (1 + δ_ab) / 20 より ∫ p q dV = V[(Σp_a)(Σq_b) + Σ_a p_a q_a] / 20
	/*
	厳密値の積分は**アフィンな要素**の閉形式で組む (組み立てとは独立の計算)。
	{a + b×r} の Nedelec 補間が厳密なのはアフィンな要素だけなので、
	六面体は平行六面体、角柱は「上面が下面の平行移動」であることを検査し、
	外れていれば恒等式の前提が無いものとしてエラーにする
	(四面体は常にアフィン)。

	  四面体       : ∫ p q dV = V [(Σp_a)(Σq_a) + Σ_a p_a q_a] / 20
	  平行六面体   : ∫ p q dV = V [p(c) q(c) + (1/12) Σ_axes (Δ_a p)(Δ_a q)]
	                 (Δ_a = 軸 a に沿う頂点差。アフィンなら定数)
	  角柱         : 3 つの四面体 (0,1,2,5) (0,1,5,4) (0,4,5,3) に割って
	                 四面体の閉形式を使う
	*/
	double svol = 0, nvol = 0, mvol = 0, vtot = 0;
	int naff = 0;
	const int nelem = edge_elem_count();
	for (int e = 0; e < nelem; e++) {
		const int m = edge_elem_mat(e);
		const int kind = edge_elem_kind(e);
		double nu[6];
		material_coef_pub(m, 3, nu);
		const int32_t *nd = edge_elem_nodes(e);
		double vol = 0, ee = 0;

		if (kind == MESHELEM_HEX) {
			// アフィン (平行六面体) か : 非アフィン成分の係数 (ξη, ηζ, ζξ, ξηζ)
			double h = 0;
			for (int a = 0; a < 8; a++) {
				for (int b = a + 1; b < 8; b++) {
					const double dx = Xp[nd[b]] - Xp[nd[a]];
					const double dy = Yp[nd[b]] - Yp[nd[a]];
					const double dz = Zp[nd[b]] - Zp[nd[a]];
					const double d = sqrt((dx * dx) + (dy * dy) + (dz * dz));
					if (d > h) h = d;
				}
			}
			for (int cd = 0; cd < 3; cd++) {
				const double *p = ((cd == 0) ? Xp : (cd == 1) ? Yp : Zp);
				for (int mode = 0; mode < 4; mode++) {
					double sum = 0;
					for (int a = 0; a < 8; a++) {
						const int sx = EHEX_SGN[a][0], sy = EHEX_SGN[a][1], sz = EHEX_SGN[a][2];
						const int w = ((mode == 0) ? (sx * sy) : (mode == 1) ? (sy * sz)
						             : (mode == 2) ? (sz * sx) : (sx * sy * sz));
						sum += w * p[nd[a]];
					}
					if ((h > 0) && ((fabs(sum) / h) > 1e-9)) naff++;
				}
			}
			// 平行六面体の体積と ∫|E|^2
			const double e1[3] = {Xp[nd[1]] - Xp[nd[0]], Yp[nd[1]] - Yp[nd[0]], Zp[nd[1]] - Zp[nd[0]]};
			const double e2[3] = {Xp[nd[3]] - Xp[nd[0]], Yp[nd[3]] - Yp[nd[0]], Zp[nd[3]] - Zp[nd[0]]};
			const double e3[3] = {Xp[nd[4]] - Xp[nd[0]], Yp[nd[4]] - Yp[nd[0]], Zp[nd[4]] - Zp[nd[0]]};
			vol = fabs((e1[0] * ((e2[1] * e3[2]) - (e2[2] * e3[1])))
			         - (e1[1] * ((e2[0] * e3[2]) - (e2[2] * e3[0])))
			         + (e1[2] * ((e2[0] * e3[1]) - (e2[1] * e3[0]))));
			for (int p = 0; p < 3; p++) {
				const int p1 = (p + 1) % 3, p2 = (p + 2) % 3;
				double ea[8];
				for (int a = 0; a < 8; a++) {
					const double r1 = (p1 == 0) ? Xp[nd[a]] : ((p1 == 1) ? Yp[nd[a]] : Zp[nd[a]]);
					const double r2 = (p2 == 0) ? Xp[nd[a]] : ((p2 == 1) ? Yp[nd[a]] : Zp[nd[a]]);
					ea[a] = ((av[p1] * r2) - (av[p2] * r1)) / 2;
				}
				double pc = 0;
				for (int a = 0; a < 8; a++) pc += ea[a] / 8;
				const double d1 = ea[1] - ea[0], d2 = ea[3] - ea[0], d3 = ea[4] - ea[0];
				ee += (pc * pc) + (((d1 * d1) + (d2 * d2) + (d3 * d3)) / 12);
			}
			ee *= vol;
		}
		else if (kind == MESHELEM_PRISM) {
			// アフィン (上面 = 下面の平行移動) か
			double h = 0;
			for (int a = 0; a < 6; a++) {
				for (int b = a + 1; b < 6; b++) {
					const double dx = Xp[nd[b]] - Xp[nd[a]];
					const double dy = Yp[nd[b]] - Yp[nd[a]];
					const double dz = Zp[nd[b]] - Zp[nd[a]];
					const double d = sqrt((dx * dx) + (dy * dy) + (dz * dz));
					if (d > h) h = d;
				}
			}
			for (int cd = 0; cd < 3; cd++) {
				const double *p = ((cd == 0) ? Xp : (cd == 1) ? Yp : Zp);
				for (int l = 1; l < 3; l++) {
					const double d = fabs((p[nd[3 + l]] - p[nd[l]]) - (p[nd[3]] - p[nd[0]]));
					if ((h > 0) && ((d / h) > 1e-9)) naff++;
				}
			}
			// 3 つの四面体に割って四面体の閉形式を使う
			static const int PSPLIT[3][4] = {{0, 1, 2, 5}, {0, 1, 5, 4}, {0, 4, 5, 3}};
			for (int t = 0; t < 3; t++) {
				const int32_t td[4] = {nd[PSPLIT[t][0]], nd[PSPLIT[t][1]],
				                       nd[PSPLIT[t][2]], nd[PSPLIT[t][3]]};
				double g[4][3], tv;
				if (tet_grad_pub(td, g, &tv)) continue;
				vol += tv;
				for (int p = 0; p < 3; p++) {
					const int p1 = (p + 1) % 3, p2 = (p + 2) % 3;
					double sm = 0, s2 = 0;
					for (int a = 0; a < 4; a++) {
						const double r1 = (p1 == 0) ? Xp[td[a]] : ((p1 == 1) ? Yp[td[a]] : Zp[td[a]]);
						const double r2 = (p2 == 0) ? Xp[td[a]] : ((p2 == 1) ? Yp[td[a]] : Zp[td[a]]);
						const double eb = ((av[p1] * r2) - (av[p2] * r1)) / 2;
						sm += eb;
						s2 += eb * eb;
					}
					ee += tv * ((sm * sm) + s2) / 20;
				}
			}
		}
		else {
			double g[4][3];
			if (tet_grad_pub(nd, g, &vol)) continue;
			for (int p = 0; p < 3; p++) {
				const int p1 = (p + 1) % 3, p2 = (p + 2) % 3;
				double sm = 0, s2 = 0;
				for (int a = 0; a < 4; a++) {
					const double r1 = (p1 == 0) ? Xp[nd[a]] : ((p1 == 1) ? Yp[nd[a]] : Zp[nd[a]]);
					const double r2 = (p2 == 0) ? Xp[nd[a]] : ((p2 == 1) ? Yp[nd[a]] : Zp[nd[a]]);
					const double ea = ((av[p1] * r2) - (av[p2] * r1)) / 2;	// E_p at node a
					sm += ea;
					s2 += ea * ea;
				}
				ee += vol * ((sm * sm) + s2) / 20;
			}
		}

		svol += Material[m].sigma * vol;
		nvol += ((nu[0] * av[0] * av[0]) + (nu[1] * av[1] * av[1]) + (nu[2] * av[2] * av[2])
		      + (2 * nu[3] * av[0] * av[1]) + (2 * nu[4] * av[1] * av[2])
		      + (2 * nu[5] * av[2] * av[0])) * vol;
		mvol += Material[m].sigma * ee;
		vtot += vol;
	}
	if (naff > 0) {
		fprintf(fp_log, "*** analysis E needs affine elements (parallelepiped "
			"hexahedra / translated prisms); the exact identities do not hold "
			"on this mesh\n");
		ierr = 1;
	}
	fprintf(fp_log, "  volume = %.6e [m^3], sum(sigma V) = %.6e, sum(a^T nu a V) = %.6e\n",
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

	// (c) 回転場 E = (1/2)(av × r) -> curl E = av
	// (c2) 同じ場の質量。空間変化する場なので質量行列の階数落ちも検出できる
	{
		for (int i = 0; i < ne; i++) {
			const int32_t a = elo[i], b = EdgeTo[i];
			const double xm = (Xp[a] + Xp[b]) / 2;
			const double ym = (Yp[a] + Yp[b]) / 2;
			const double zm = (Zp[a] + Zp[b]) / 2;
			// 直線辺上で E は 1 次なので中点則で厳密
			const double ex = ((av[1] * zm) - (av[2] * ym)) / 2;
			const double ey = ((av[2] * xm) - (av[0] * zm)) / 2;
			const double ez = ((av[0] * ym) - (av[1] * xm)) / 2;
			u[i] = (ex * (Xp[b] - Xp[a])) + (ey * (Yp[b] - Yp[a])) + (ez * (Zp[b] - Zp[a]));
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

		crs_spmv(&T, u, y, NULL);
		double q2 = 0;
		for (int i = 0; i < ne; i++) q2 += u[i] * y[i];
		const double rel2 = fabs(q2 - mvol) / ((mvol != 0) ? fabs(mvol) : 1);
		fprintf(fp_log, "  (c2) rotational mass    : u^T T u = %.10e, exact = %.10e, err = %.3e\n",
			q2, mvol, rel2);
		if (rel2 > 1e-10) {
			fprintf(fp_log, "*** the edge mass matrix is wrong (non-uniform field)\n");
			ierr = 1;
		}
	}

	// (d) 対称性
	{
		double amax = 0, tmax = 0, sdif = 0, tdif = 0;
		for (int i = 0; i < ne; i++) {
			for (int64_t p = S.rowptr[i]; p < S.rowptr[i + 1]; p++) {
				const int32_t j = S.col[p];
				const int64_t pj = crs_find_edge(&S, j, (int32_t)i);
				if (pj < 0) continue;
				if (fabs(S.val[p]) > amax) amax = fabs(S.val[p]);
				if (fabs(T.val[p]) > tmax) tmax = fabs(T.val[p]);
				const double d = fabs(S.val[p] - S.val[pj]);
				if (d > sdif) sdif = d;
				const double dt = fabs(T.val[p] - T.val[pj]);
				if (dt > tdif) tdif = dt;
			}
		}
		const double rel  = ((amax > 0) ? (sdif / amax) : 0);
		const double relt = ((tmax > 0) ? (tdif / tmax) : 0);
		fprintf(fp_log, "  (d) symmetry            : max|S-S^T|/max|S| = %.3e   "
			"max|T-T^T|/max|T| = %.3e\n", rel, relt);
		if ((rel > 1e-12) || (relt > 1e-12)) ierr = 1;
	}

	// (g) 面をまたぐ接線連続性 (種別の混在で本命になる検査)
	{
		double *ur = (double *)malloc((size_t)ne * sizeof(double));
		for (int i = 0; i < ne; i++) {
			ur[i] = fmod((i + 1) * 0.41421356237, 1.0) - 0.5;	// 決定的な擬似乱数
		}
		if (edge_face_conform(fp_log, ur)) ierr = 1;
		free(ur);
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
	// 節点の隣接パターン (同じ要素の辺の端点は同じ要素の節点なので、
	// 種別に依らず G^T A G の寄与は必ずこのパターンに入る)
	crs_alloc_elem3d(N);

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
