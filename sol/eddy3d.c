/*
eddy3d.c

3 次元渦電流解析 (A-φ 定式化、時間調和)。段階 5。

段階 1〜4 で用意した辺 (Nedelec) 要素・ゲージ・前処理の上に、
実際に解く連成系を組む。未知数はベクトルポテンシャル A (辺) と
電気スカラーポテンシャル φ (節点) で、電界は

	E = -jω A - ∇φ,   J = σ E

支配方程式は Ampere 則と電流連続の 2 本 :

	∇×(ν ∇×A) + jωσ A + σ∇φ = 0
	∇・[ σ (jωA + ∇φ) ] = 0

Galerkin 離散化 (辺基底 W、節点基底 N、離散勾配 G) で

	S A + jω T A + T G φ = 0
	jω G^T T A + G^T T G φ = 0

ここで S_ef = ∫(∇×W_e)^T ν (∇×W_f)、T_ef = ∫σ W_e・W_f。
G は節点値を辺自由度に写す離散勾配 (Whitney 恒等式より Σ(Gφ)_e W_e = ∇φ が厳密)
なので、G^T T G はそのまま節点の導電行列 ∫σ∇N_a・∇N_b になる。

**2 本目を jω で割ると係数行列が複素対称になる**:

	[ S + jωT      T G      ] [A]   [0]
	[ G^T T     G^T T G/(jω)] [φ] = [0]

	(1,2) = TG、(2,1) = G^T T = (TG)^T、対角ブロックも対称。
	Hermite ではなく複素**対称**なので、共役を取らない双一次形式の COCG が使える
	(solver_cocg をそのまま流用、A = Kr + j Ki の形に落とす)。

φ を電極で Dirichlet 固定して励振し、電極上の反作用から電流を取り出す:

	Σ_{m∈電極p} ∫∇N_m・σ(jωA + ∇φ) dV = ∮_p J・n dS の符号反転 = 電極 p に
	流れ込む電流 I_p

ω→0 で A が落ちて G^T T G φ = 0 だけが残り、DC の導電解析に厳密に一致する
(検証ケースはこの極限も見る)。

ゲージについて (既定はゲージ固定しない):
	連成系は (A, φ) -> (A + Gψ, φ - jωψ) で不変なので特異だが、右辺を
	b = -A x_D (Dirichlet の持ち上げ) で作る限り b は必ず値域に入るので
	系は無矛盾で、COCG は解に収束する。出力する電流 I は E = -(jωA + ∇φ) から
	決まるゲージ不変量なので、A と φ の非一意性は Z に影響しない。

	実測 (bar_eddy、1e4/1e5 Hz) : ゲージ固定の有無で Z は 7 桁一致し、
	固定しない方が反復回数が約 6 倍少ない (403 対 2546、437 対 2899)。
	そのため既定はゲージ固定なし。ただしその解の A と φ 自体は物理的な値では
	なく、A が物理解の 1000 倍、両端 0V/1V の φ が ±60V まで振れる。
	A や φ そのものが要るとき (場の出力など) は gauge = 1 を指定する。

要素の形について:
	1 次 Nedelec 要素の誤差は「場が変化する方向の刻み」ではなく**要素の最大寸法**
	で決まる。bar_eddy で dz を 1/40 まで細かくしても、場が一様な x 方向の
	dx = 1.7mm を放置すると R が 13% ずれた (dx = 0.42mm で 1%、0.083mm で 0.06%)。
	扁平な要素は避けること。
*/

#include "fem.h"
#include "fem_prototype.h"

// 四面体の局所辺 (節点対) — edge.c と同じ定義
static const int EDGE_NODE3[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};


static int cmp_i32_en(const void *a, const void *b)
{
	const int32_t x = *(const int32_t *)a;
	const int32_t y = *(const int32_t *)b;

	return ((x < y) ? -1 : ((x > y) ? 1 : 0));
}


// 行 row の中で列 col の位置を二分探索で求める (列は昇順)
static int64_t crs_find_en(const crs_t *A, int32_t row, int32_t col)
{
	int64_t lo = A->rowptr[row];
	int64_t hi = A->rowptr[row + 1] - 1;

	while (lo <= hi) {
		const int64_t mid = (lo + hi) / 2;
		if      (A->col[mid] < col) lo = mid + 1;
		else if (A->col[mid] > col) hi = mid - 1;
		else                        return mid;
	}

	return -1;
}


// ---- 連成系 (辺 + 節点) の CRS パターン ----
//
// 自由度は [0, NEdge) が辺、[NEdge, NEdge+NNode) が節点。
// 同じ四面体に属する 10 個の自由度 (辺 6 + 節点 4) は互いに結合する。
void crs_alloc_edge_node(crs_t *A)
{
	const int ne = NEdge;
	const int nn = NNode;
	const int n = ne + nn;
	const int nelem = edge_elem_count();

	// 自由度毎の要素数
	int *cnt = (int *)malloc((size_t)n * sizeof(int));
	memset(cnt, 0, (size_t)n * sizeof(int));
	for (int e = 0; e < nelem; e++) {
		const int32_t *nd = edge_elem_nodes(e);
		const int nen = edge_elem_nen(e);
		for (int64_t k = EdgeOff[e]; k < EdgeOff[e + 1]; k++) cnt[TetEdge[k]]++;
		for (int l = 0; l < nen; l++) cnt[ne + nd[l]]++;
	}
	int64_t *eptr = (int64_t *)malloc(((size_t)n + 1) * sizeof(int64_t));
	eptr[0] = 0;
	for (int i = 0; i < n; i++) eptr[i + 1] = eptr[i] + cnt[i];
	int32_t *elist = (int32_t *)malloc((size_t)eptr[n] * sizeof(int32_t));
	memset(cnt, 0, (size_t)n * sizeof(int));
	for (int e = 0; e < nelem; e++) {
		const int32_t *nd = edge_elem_nodes(e);
		const int nen = edge_elem_nen(e);
		for (int64_t k = EdgeOff[e]; k < EdgeOff[e + 1]; k++) {
			const int32_t d = TetEdge[k];
			elist[eptr[d] + cnt[d]] = (int32_t)e;
			cnt[d]++;
		}
		for (int l = 0; l < nen; l++) {
			const int32_t d = ne + nd[l];
			elist[eptr[d] + cnt[d]] = (int32_t)e;
			cnt[d]++;
		}
	}

	A->n = n;
	A->rowptr = (int64_t *)malloc(((size_t)n + 1) * sizeof(int64_t));

	int cap = 256;
	int32_t *work = (int32_t *)malloc((size_t)cap * sizeof(int32_t));
	int *rown = (int *)malloc((size_t)n * sizeof(int));

	for (int i = 0; i < n; i++) {
		const int64_t p0 = eptr[i], p1 = eptr[i + 1];
		int need = 0;
		for (int64_t p = p0; p < p1; p++) {
			need += edge_elem_nedge(elist[p]) + edge_elem_nen(elist[p]);
		}
		if (need > cap) {
			cap = need;
			work = (int32_t *)realloc(work, (size_t)cap * sizeof(int32_t));
		}
		int m = 0;
		for (int64_t p = p0; p < p1; p++) {
			const int32_t e = elist[p];
			const int32_t *nd = edge_elem_nodes(e);
			const int nen = edge_elem_nen(e);
			for (int64_t k = EdgeOff[e]; k < EdgeOff[e + 1]; k++) work[m++] = TetEdge[k];
			for (int l = 0; l < nen; l++) work[m++] = (int32_t)ne + nd[l];
		}
		qsort(work, (size_t)m, sizeof(int32_t), cmp_i32_en);
		int u = 0;
		for (int q = 0; q < m; q++) {
			if ((q == 0) || (work[q] != work[q - 1])) u++;
		}
		rown[i] = u;
	}

	A->rowptr[0] = 0;
	for (int i = 0; i < n; i++) A->rowptr[i + 1] = A->rowptr[i] + rown[i];
	A->nnz = A->rowptr[n];
	A->col = (int32_t *)malloc((size_t)A->nnz * sizeof(int32_t));
	A->val = (double *)malloc((size_t)A->nnz * sizeof(double));

	for (int i = 0; i < n; i++) {
		const int64_t p0 = eptr[i], p1 = eptr[i + 1];
		int m = 0;
		for (int64_t p = p0; p < p1; p++) {
			const int32_t e = elist[p];
			const int32_t *nd = edge_elem_nodes(e);
			const int nen = edge_elem_nen(e);
			for (int64_t k = EdgeOff[e]; k < EdgeOff[e + 1]; k++) work[m++] = TetEdge[k];
			for (int l = 0; l < nen; l++) work[m++] = (int32_t)ne + nd[l];
		}
		qsort(work, (size_t)m, sizeof(int32_t), cmp_i32_en);
		int64_t w = A->rowptr[i];
		for (int q = 0; q < m; q++) {
			if ((q == 0) || (work[q] != work[q - 1])) A->col[w++] = work[q];
		}
	}

	free(cnt);
	free(eptr);
	free(elist);
	free(work);
	free(rown);

	crs_zero(A);
}


// ---- 連成系の組み立て ----
//
// A = Kr + j Ki (どちらも実対称)。ω は角周波数。
//
//	Kr = [ S     TG ;  (TG)^T   0        ]
//	Ki = [ ωT    0  ;   0      -G^T TG/ω ]
// nuelem が非 NULL なら、B-H 曲線を持つ材料の要素だけ等方性の ν = nuelem[e] を使う
// (非線形の実効磁気抵抗率の反復)。sigelem が非 NULL なら要素ごとの σ = sigelem[e]
// (熱連成の σ(T))。どちらも NULL なら従来と同じ演算
static void assemble_eddy3d_sig(crs_t *Kr, crs_t *Ki, double omega, const double *nuelem,
	const double *sigelem);

void assemble_eddy3d(crs_t *Kr, crs_t *Ki, double omega, const double *nuelem)
{
	assemble_eddy3d_sig(Kr, Ki, omega, nuelem, NULL);
}


static void assemble_eddy3d_sig(crs_t *Kr, crs_t *Ki, double omega, const double *nuelem,
	const double *sigelem)
{
	crs_zero(Kr);
	crs_zero(Ki);

	const int ne = NEdge;
	const int nelem = edge_elem_count();

	for (int e = 0; e < nelem; e++) {
		const int m = edge_elem_mat(e);
		const int nedge = edge_elem_nedge(e);
		const int nen = edge_elem_nen(e);
		const int (*tab)[2] = edge_elem_table(e);
		double nu[6];
		material_coef_pub(m, 3, nu);			// ν = (μ0 μ~)^-1
		if ((nuelem != NULL) && (Material[m].nbh[0] > 0)) {
			nu[0] = nu[1] = nu[2] = nuelem[e];
			nu[3] = nu[4] = nu[5] = 0;
		}
		const double sig = ((sigelem != NULL) ? sigelem[e] : Material[m].sigma);

		double se[12][12], te[12][12];
		edge_elem_matrices(e, nu, sig, se, te);

		// 局所の離散勾配 GL[k][a] : 辺 k = (局所節点 p, q) -> -1 at p, +1 at q
		double gl[12][8];
		for (int k = 0; k < nedge; k++) {
			for (int a = 0; a < nen; a++) gl[k][a] = 0;
			gl[k][tab[k][0]] = -1;
			gl[k][tab[k][1]] = +1;
		}

		// tg = te * GL、gtg = GL^T * te * GL
		double tg[12][8], gtg[8][8];
		for (int k = 0; k < nedge; k++) {
			for (int a = 0; a < nen; a++) {
				double s = 0;
				for (int l = 0; l < nedge; l++) s += te[k][l] * gl[l][a];
				tg[k][a] = s;
			}
		}
		for (int a = 0; a < nen; a++) {
			for (int b = 0; b < nen; b++) {
				double s = 0;
				for (int k = 0; k < nedge; k++) s += gl[k][a] * tg[k][b];
				gtg[a][b] = s;
			}
		}

		const int32_t *ed = &TetEdge[EdgeOff[e]];
		const signed char *sg = &TetEdgeSgn[EdgeOff[e]];
		const int32_t *nd = edge_elem_nodes(e);

		// (1,1) 辺 x 辺
		for (int k = 0; k < nedge; k++) {
			for (int l = 0; l < nedge; l++) {
				const int64_t p = crs_find_en(Kr, ed[k], ed[l]);
				const double w = sg[k] * sg[l];
				Kr->val[p] += w * se[k][l];
				Ki->val[p] += w * omega * te[k][l];
			}
		}

		// (1,2) 辺 x 節点、(2,1) 節点 x 辺 (実部のみ)
		for (int k = 0; k < nedge; k++) {
			for (int a = 0; a < nen; a++) {
				const double w = sg[k] * tg[k][a];
				Kr->val[crs_find_en(Kr, ed[k], (int32_t)ne + nd[a])] += w;
				Kr->val[crs_find_en(Kr, (int32_t)ne + nd[a], ed[k])] += w;
			}
		}

		// (2,2) 節点 x 節点 (虚部のみ、-G^T T G / ω)
		for (int a = 0; a < nen; a++) {
			for (int b = 0; b < nen; b++) {
				const int64_t p = crs_find_en(Ki, (int32_t)ne + nd[a], (int32_t)ne + nd[b]);
				Ki->val[p] -= gtg[a][b] / omega;
			}
		}
	}
}


// ---- A_t = 0 の辺 (awall タグの三角形の 3 辺) ----
static void awall_mark(unsigned char *fix)
{
	for (int t = 0; t < NTri; t++) {
		int hit = 0;
		for (int q = 0; q < NAWall; q++) {
			if (TriTag[t] == AWallTag[q]) hit = 1;
		}
		if (!hit) continue;
		const int32_t *v = &Tri[t * 3];
		for (int l = 0; l < 3; l++) {
			const int64_t id = edge_id(v[l], v[(l + 1) % 3]);	// 向きは edge_id が揃える
			if (id >= 0) fix[id] = 1;
		}
	}
	// 六面体・角柱格子の awall 面は四角形にもなる (4 辺)
	for (int t = 0; t < NQuad; t++) {
		int hit = 0;
		for (int q = 0; q < NAWall; q++) {
			if (QuadTag[t] == AWallTag[q]) hit = 1;
		}
		if (!hit) continue;
		const int32_t *v = &Quad[t * 4];
		for (int l = 0; l < 4; l++) {
			const int64_t id = edge_id(v[l], v[(l + 1) % 4]);
			if (id >= 0) fix[id] = 1;
		}
	}
}


// ---- ゲージ固定 (tree-cotree、awall 優先) ----  gauge = 1 のときだけ使う
//
// 全域木の辺で A_e = 0 と置くと (Gψ)_e = 0 が全木辺で成り立ち ψ = 一定、
// さらに電極の φ 固定から ψ = 0 になって零空間が消える。
//
// ただし A_t = 0 (awall) の辺は物理的な境界条件なので、木がそれを壊さないよう
// **awall の辺を先に木へ入れる** (Kruskal の優先処理)。こうすると木から決まる
// ψ は awall 上で一定になり、awall が連結なら電極上で ψ = 0 にできる
// = 電極間の電位差 V が保たれる。awall の辺だけで電極が連結していないと
// ψ が電極毎に別の定数になって V がずれるので、その場合はエラーで落とす
// (既定の gauge = 0 は常に正しいので、黙って続ける理由が無い)。
//
// Z はゲージ不変なので既定 (gauge = 0) では使わない。A・φ 自体が要るときだけ。
static int uf_find(int *par, int x)
{
	while (par[x] != x) {
		par[x] = par[par[x]];
		x = par[x];
	}

	return x;
}


// pri (awall) の辺を優先して全域木を作る。comp には pri の辺だけで作った
// 連結成分の代表節点を返す。戻り値は木の辺数
static int edge_tree_gauge(const unsigned char *pri, unsigned char *tree, int *comp)
{
	int *par = (int *)malloc((size_t)NNode * sizeof(int));
	for (int i = 0; i < NNode; i++) par[i] = i;
	memset(tree, 0, (size_t)NEdge * sizeof(unsigned char));

	int ntree = 0;
	for (int pass = 0; pass < 2; pass++) {
		for (int e = 0; e < NEdge; e++) {
			const int want = ((pass == 0) ? 1 : 0);
			if ((pri[e] != 0) != want) continue;
			const int a = uf_find(par, EdgeFrom[e]);
			const int b = uf_find(par, EdgeTo[e]);
			if (a == b) continue;
			par[a] = b;
			tree[e] = 1;
			ntree++;
		}
		// 優先辺だけの連結成分を控える
		if (pass == 0) {
			for (int i = 0; i < NNode; i++) comp[i] = uf_find(par, i);
		}
	}

	free(par);

	return ntree;
}


// 行 row の (Kr + j Ki) x の値 (反作用の取り出しに使う)
static void row_dot_c(const crs_t *Kr, const crs_t *Ki, int32_t row,
	const double *xr, const double *xi, double *sr, double *si)
{
	double ar = 0, ai = 0;
	for (int64_t p = Kr->rowptr[row]; p < Kr->rowptr[row + 1]; p++) {
		const int32_t c = Kr->col[p];
		ar += (Kr->val[p] * xr[c]) - (Ki->val[p] * xi[c]);
		ai += (Kr->val[p] * xi[c]) + (Ki->val[p] * xr[c]);
	}
	*sr = ar;
	*si = ai;
}


#define AAM3 (5)		// 非線形の Anderson 加速の深さ

// 要素中心の B = ∇×A (実部・虚部)。四面体は Whitney 基底の閉形式
// (∇×W_e = 2∇λ_a×∇λ_b、要素内で一定)、六面体・角柱は参照基底を Piola 変換した
// 要素中心の値 (場の出力と同じ式)
static void elem_bfield(int e, const double *ur, const double *ui, double bre[3], double bim[3])
{
	const int kind = edge_elem_kind(e);
	const int nedge = edge_elem_nedge(e);
	const int32_t *ed = &TetEdge[EdgeOff[e]];
	const signed char *sg = &TetEdgeSgn[EdgeOff[e]];
	for (int c = 0; c < 3; c++) bre[c] = bim[c] = 0;
	if (kind == MESHELEM_TET) {
		double g[4][3], vol;
		if (tet_grad_pub(&Tet[e * 4], g, &vol)) return;
		for (int k = 0; k < 6; k++) {
			const int a = EDGE_NODE3[k][0], b = EDGE_NODE3[k][1];
			const double cr[3] = {
				(g[a][1] * g[b][2]) - (g[a][2] * g[b][1]),
				(g[a][2] * g[b][0]) - (g[a][0] * g[b][2]),
				(g[a][0] * g[b][1]) - (g[a][1] * g[b][0])};
			const double wr = sg[k] * ur[ed[k]], wi = sg[k] * ui[ed[k]];
			for (int c = 0; c < 3; c++) {
				bre[c] += 2 * wr * cr[c];
				bim[c] += 2 * wi * cr[c];
			}
		}
	}
	else {
		double wc[12][3], cc[12][3];
		edge_elem_center(e, wc, cc);
		for (int k = 0; k < nedge; k++) {
			const double wr = sg[k] * ur[ed[k]], wi = sg[k] * ui[ed[k]];
			for (int c = 0; c < 3; c++) {
				bre[c] += wr * cc[k][c];
				bim[c] += wi * cc[k][c];
			}
		}
	}
}


int solve_eddy3d(FILE *fp_log)
{
	const int np = NPort;
	int ierr = 0;

	fprintf(fp_log, "\n=== 3D eddy current (A-phi, edge elements) ===\n");

	if (!MeshMode) {
		fprintf(fp_log, "*** analysis A requires an unstructured mesh (mesh key)\n");
		return 1;
	}
	if (np < 1) {
		fprintf(fp_log, "*** analysis A requires at least one port\n");
		return 1;
	}
	if (Freq <= 0) {
		fprintf(fp_log, "*** analysis A requires the frequency key\n");
		return 1;
	}

	edge_build();

	const int ne = NEdge;
	const int nn = NNode;
	const int n = ne + nn;
	const double omega = 2 * PI * Freq;

	// 導体の σ μr (最大値) から表皮深さの目安を出す。
	// δ = sqrt(2/(ω μ0 μr σ)) なので μr を落とすと磁性導体で δ を過大評価し、
	// 下の格子分解能の警告が出るべきときに出なくなる。
	// 異方性のときは向きが分からないので μr の最大値 (= δ の最小値) を取る。
	// mu6[0..2] は anisomur 未指定なら mur で埋まっている (input_data.c)
	double sigmax = 0, sigmumax = 0;
	const int nelem = edge_elem_count();
	for (int e = 0; e < nelem; e++) {
		const material_t *mt = &Material[edge_elem_mat(e)];
		const double s = mt->sigma;
		if (s > sigmax) sigmax = s;
		if (s <= 0) continue;
		double mr = mt->mu6[0];
		if (mt->mu6[1] > mr) mr = mt->mu6[1];
		if (mt->mu6[2] > mr) mr = mt->mu6[2];
		if (mr <= 0) mr = 1;
		if ((s * mr) > sigmumax) sigmumax = s * mr;
	}
	if (sigmax <= 0) {
		fprintf(fp_log, "*** analysis A requires a conducting material (sigma > 0)\n");
		edge_free();
		return 1;
	}
	const double delta = sqrt(2 / (omega * MU0 * sigmumax));
	fprintf(fp_log, "  nodes = %d, elements = %d (%d tet / %d hex / %d prism), "
		"edges = %d, unknowns = %d\n", nn, nelem, NTet, NHex, NPrism, ne, n);
	fprintf(fp_log, "  frequency = %.6e [Hz], skin depth = %.4e [m]\n", Freq, delta);

	// 導体要素の最大辺長を表皮深さと並べて報告する。
	// 1 次 Nedelec 要素の誤差は「場が変化する方向の刻み」ではなく
	// **要素の最大寸法**で決まるので、扁平な要素は場が一様な方向にも効く
	// (実測 : 場が x に一様でも dx を 1.7mm -> 0.42mm にすると誤差が 13% -> 1%)
	double hmax = 0;
	for (int e = 0; e < nelem; e++) {
		if (Material[edge_elem_mat(e)].sigma <= 0) continue;
		const int32_t *nd = edge_elem_nodes(e);
		const int nedge = edge_elem_nedge(e);
		const int (*tab)[2] = edge_elem_table(e);
		for (int k = 0; k < nedge; k++) {
			const int32_t a = nd[tab[k][0]], b = nd[tab[k][1]];
			const double dx = Xp[b] - Xp[a], dy = Yp[b] - Yp[a], dz = Zp[b] - Zp[a];
			const double h = sqrt((dx * dx) + (dy * dy) + (dz * dz));
			if (h > hmax) hmax = h;
		}
	}
	fprintf(fp_log, "  max conductor element size = %.4e [m] (%.2f x skin depth)\n",
		hmax, hmax / delta);
	if (hmax > delta) {
		fprintf(fp_log, "*** warning : the largest conductor element exceeds the skin "
			"depth; refine the mesh in every direction, not only across the skin\n");
	}

	crs_t Kr, Ki;
	crs_alloc_edge_node(&Kr);
	crs_alloc_edge_node(&Ki);
	// 非線形 (B-H) の磁性体があれば実効磁気抵抗率法 (solve.c の F と同じ):
	// 要素ごとの ν を波高値 |B| で評価し、H を独立変数にした更新 ν = |H|/B(|H|) に
	// ln ν の Anderson 加速を掛ける。単一ポート・電流振幅 current
	int nlf = 0;
	for (int m = 0; m < NMaterial; m++) {
		if (Material[m].nbh[0] > 0) nlf = 1;
	}
	double *nuelem = NULL, *aay = NULL, *aaf = NULL, *aayp = NULL, *aafp = NULL;
	double *aady = NULL, *aadf = NULL;
	int aan = 0;
	if (nlf) {
		nuelem = (double *)malloc((size_t)nelem * sizeof(double));
		aay  = (double *)malloc((size_t)nelem * sizeof(double));
		aaf  = (double *)malloc((size_t)nelem * sizeof(double));
		aayp = (double *)malloc((size_t)nelem * sizeof(double));
		aafp = (double *)malloc((size_t)nelem * sizeof(double));
		aady = (double *)malloc((size_t)nelem * AAM3 * sizeof(double));
		aadf = (double *)malloc((size_t)nelem * AAM3 * sizeof(double));
		for (int e = 0; e < nelem; e++) {
			const material_t *mt = &Material[edge_elem_mat(e)];
			nuelem[e] = ((mt->nbh[0] > 0) ? bh_nu(mt, 0, 0) : 0);
		}
	}
	assemble_eddy3d(&Kr, &Ki, omega, nuelem);
	fprintf(fp_log, "  coupled matrix : %lld nonzeros (%.1f per row)\n",
		(long long)Kr.nnz, (double)Kr.nnz / ((n > 0) ? n : 1));

	// Dirichlet : awall の辺 (A_t = 0)、電極の節点 (φ = V)、
	//             非導電領域の節点 (φ は定義されない -> 0 に固定)
	unsigned char *fix = (unsigned char *)malloc((size_t)n * sizeof(unsigned char));
	memset(fix, 0, (size_t)n * sizeof(unsigned char));
	unsigned char *wall = (unsigned char *)malloc((size_t)ne * sizeof(unsigned char));
	memset(wall, 0, (size_t)ne * sizeof(unsigned char));
	awall_mark(wall);
	int nwall = 0;
	for (int i = 0; i < ne; i++) {
		fix[i] = wall[i];
		nwall += wall[i];
	}

	// ゲージ固定 (awall 優先の全域木)。gauge = 1 のときだけ
	unsigned char *tree = NULL;
	int *comp = NULL;
	if (GaugeTree) {
		tree = (unsigned char *)malloc((size_t)ne * sizeof(unsigned char));
		comp = (int *)malloc((size_t)nn * sizeof(int));
		const int ntree = edge_tree_gauge(wall, tree, comp);
		int ngauge = 0;
		for (int i = 0; i < ne; i++) {
			if (tree[i] && !fix[i]) ngauge++;
			fix[i] |= tree[i];
		}
		// 木から決まる ψ が電極上で 0 になる条件 : 電極が awall の辺だけで連結
		int bad = 0, ref = -1;
		for (int i = 0; i < nn; i++) {
			if (NodeConductor[i] < 0) continue;
			if (ref < 0) ref = comp[i];
			else if (comp[i] != ref) bad = 1;
		}
		fprintf(fp_log, "  gauge : %d tree edges (%d beyond the wall)\n", ntree, ngauge);
		if (bad) {
			// 木から決まる ψ が電極毎に違う定数になり、端子電圧 V がずれる。
			// 既定 (gauge = 0) なら常に正しいので、黙って続けずに落とす
			fprintf(fp_log, "*** the electrodes are not connected through awall edges, "
				"so the tree gauge would shift the terminal voltage; "
				"use gauge = 0 (the default), or extend awall so that all electrodes "
				"share one A_t = 0 surface\n");
			ierr = 1;
		}
	}

	unsigned char *cond = (unsigned char *)malloc((size_t)nn * sizeof(unsigned char));
	memset(cond, 0, (size_t)nn * sizeof(unsigned char));
	for (int e = 0; e < nelem; e++) {
		if (Material[edge_elem_mat(e)].sigma <= 0) continue;
		const int32_t *nd = edge_elem_nodes(e);
		const int nen2 = edge_elem_nen(e);
		for (int l = 0; l < nen2; l++) cond[nd[l]] = 1;
	}
	int nelec = 0, nfloat = 0;
	for (int i = 0; i < nn; i++) {
		if (!cond[i]) {
			fix[ne + i] = 1;
			nfloat++;
		}
		else if (NodeConductor[i] >= 0) {
			fix[ne + i] = 1;
			nelec++;
		}
	}
	fprintf(fp_log, "  fixed : %d wall edges (A_t = 0), %d electrode nodes, "
		"%d non-conducting nodes\n", nwall, nelec, nfloat);
	if (nwall == 0) {
		fprintf(fp_log, "*** warning : no awall surface; the outer boundary is a "
			"magnetic wall (dA/dn = 0)\n");
	}
	fflush(fp_log);

	double *xdr = (double *)malloc((size_t)n * sizeof(double));
	double *xdi = (double *)malloc((size_t)n * sizeof(double));
	double *br = (double *)malloc((size_t)n * sizeof(double));
	double *bi = (double *)malloc((size_t)n * sizeof(double));
	double *ur = (double *)malloc((size_t)n * sizeof(double));
	double *ui = (double *)malloc((size_t)n * sizeof(double));
	double *w1 = (double *)malloc((size_t)n * sizeof(double));
	double *w2 = (double *)malloc((size_t)n * sizeof(double));
	double *yr = (double *)calloc((size_t)np * np, sizeof(double));
	double *yi = (double *)calloc((size_t)np * np, sizeof(double));

	fprintf(fp_log, "  %-10s %8s %13s\n", "excite", "iter", "residual");
	fflush(fp_log);

	// 線形なら 1 回で終わり。非線形は ν が収束するまで回し、最後にもう 1 回
	// (場の出力つき) 解いて Z を作る
	int final = !nlf;
	if (nlf) {
		fprintf(fp_log, "  nonlinear (B-H) : effective reluctivity at current amplitude "
			"%.6e [A], maxiter = %d, tol = %.1e, relax = %.2f\n", Curr, NlMaxiter, NlTol, NlRelax);
	}
	for (int nlit = 1; !ierr; nlit++) {
	if (nlit > 1) assemble_eddy3d(&Kr, &Ki, omega, nuelem);
	// 基準導体 0 を φ=0、ポート jc を φ=Volt にして順に励振する
	for (int jc = 1; jc <= np; jc++) {
		for (int i = 0; i < n; i++) xdr[i] = xdi[i] = 0;
		for (int i = 0; i < nn; i++) {
			if (fix[ne + i] && (NodeConductor[i] == jc)) xdr[ne + i] = Volt;
		}

		// b = -(A x_D)。x_D は実なので b = -(Kr x_D) - j (Ki x_D)
		crs_spmv(&Kr, xdr, w1, NULL);
		crs_spmv(&Ki, xdr, w2, NULL);
		for (int i = 0; i < n; i++) {
			br[i] = (fix[i] ? 0 : -w1[i]);
			bi[i] = (fix[i] ? 0 : -w2[i]);
		}

		char label[BUFSIZ];
		sprintf(label, "port%d", jc);
		const int iter = solver_cocg(&Kr, &Ki, 1.0, br, bi, ur, ui, fix,
			Solver.maxiter, Solver.nout, Solver.converg, (final ? fp_log : NULL), label);
		if (iter < 0) {
			fprintf(fp_log, "*** solver did not converge (port %d)\n", jc);
			ierr = 1;
		}
		for (int i = 0; i < n; i++) {
			ur[i] += xdr[i];
			ui[i] += xdi[i];
		}

		// 電極上の反作用 -> 流れ込む電流
		// 節点行は 1/(jω) 倍してあるので、元の残差は jω 倍して戻す
		double qr[MAXPORT], qi[MAXPORT];
		for (int p = 0; p < MAXPORT; p++) qr[p] = qi[p] = 0;
		for (int i = 0; i < nn; i++) {
			const int id = NodeConductor[i];
			if ((id < 0) || !cond[i]) continue;
			double rr, ri;
			row_dot_c(&Kr, &Ki, (int32_t)ne + i, ur, ui, &rr, &ri);
			// jω (rr + j ri) = -ω ri + j ω rr
			qr[id] += -omega * ri;
			qi[id] += omega * rr;
		}
		for (int kc = 1; kc <= np; kc++) {
			yr[((kc - 1) * np) + (jc - 1)] = qr[kc] / Volt;
			yi[((kc - 1) * np) + (jc - 1)] = qi[kc] / Volt;
		}
		if (final) {
			fprintf(fp_log, "  port %d : I = %.6e %+.6e j [A] (V = %.6e [V])\n",
				jc, qr[jc], qi[jc], Volt);
		}

		// 非線形 : 端子電流が current になるよう解を s = Curr / I 倍する (Y は倍率に
		// 依らないので上で作った値のまま)。B と場の出力は与えた電流での値になる
		if (nlf) {
			const double ivr = yr[0] * Volt, ivi = yi[0] * Volt;	// V = Volt での端子電流
			const double den = (ivr * ivr) + (ivi * ivi);
			const double svr = Curr * ivr / den, svi = -Curr * ivi / den;
			for (int i = 0; i < n; i++) {
				const double xr = ur[i], xi = ui[i];
				ur[i] = (svr * xr) - (svi * xi);
				ui[i] = (svr * xi) + (svi * xr);
			}
		}

		// 場の出力 : 要素毎に B = ∇×A と J = -σ(jωA + ∇φ) を作る。
		// **要素の種別ごとに分岐する** (混在格子があるため)。四面体は
		// Whitney 基底の閉形式 (∇×W_e = 2∇λ_a×∇λ_b、W_e は 1 次なので重心で
		// 評価)、六面体・角柱は参照基底を Piola 変換した要素中心の値を使う。
		// 純粋な格子では従来と同じ式を同じ順に通るので出力はバイト一致する
		if (FieldOut && (jc == 1) && final) {
			double *bv = (double *)malloc((size_t)nelem * 3 * sizeof(double));
			double *jr = (double *)malloc((size_t)nelem * 3 * sizeof(double));
			double *ji = (double *)malloc((size_t)nelem * 3 * sizeof(double));
			for (int e = 0; e < nelem; e++) {
				const int kind = edge_elem_kind(e);
				const int nedge = edge_elem_nedge(e);
				const int32_t *ed = &TetEdge[EdgeOff[e]];
				const signed char *sg = &TetEdgeSgn[EdgeOff[e]];
				const int32_t *nd = edge_elem_nodes(e);
				double brv[3] = {0, 0, 0}, arv[3] = {0, 0, 0}, aiv[3] = {0, 0, 0};
				double pr[3] = {0, 0, 0}, pi2[3] = {0, 0, 0};	// ∇φ

				if (kind == MESHELEM_TET) {
					double g[4][3], vol;
					for (int c = 0; c < 3; c++) {
						bv[(e * 3) + c] = jr[(e * 3) + c] = ji[(e * 3) + c] = 0;
					}
					if (tet_grad_pub(&Tet[e * 4], g, &vol)) continue;
					for (int k = 0; k < 6; k++) {
						const int a = EDGE_NODE3[k][0], b = EDGE_NODE3[k][1];
						const double cr[3] = {
							(g[a][1] * g[b][2]) - (g[a][2] * g[b][1]),
							(g[a][2] * g[b][0]) - (g[a][0] * g[b][2]),
							(g[a][0] * g[b][1]) - (g[a][1] * g[b][0])};
						const double wr = sg[k] * ur[ed[k]], wi = sg[k] * ui[ed[k]];
						for (int c = 0; c < 3; c++) {
							brv[c] += 2 * wr * cr[c];		// ∇×W_e = 2 cr
							// 重心での W_e = (∇λ_b - ∇λ_a)/4
							arv[c] += wr * (g[b][c] - g[a][c]) / 4;
							aiv[c] += wi * (g[b][c] - g[a][c]) / 4;
						}
					}
					for (int a = 0; a < 4; a++) {
						for (int c = 0; c < 3; c++) {
							pr[c] += ur[ne + nd[a]] * g[a][c];
							pi2[c] += ui[ne + nd[a]] * g[a][c];
						}
					}
				}
				else {
					double wc[12][3], cc[12][3];
					edge_elem_center(e, wc, cc);
					for (int k = 0; k < nedge; k++) {
						const double wr = sg[k] * ur[ed[k]], wi = sg[k] * ui[ed[k]];
						for (int c = 0; c < 3; c++) {
							brv[c] += wr * cc[k][c];
							arv[c] += wr * wc[k][c];
							aiv[c] += wi * wc[k][c];
						}
					}
					if (kind == MESHELEM_HEX) {
						double gn[8][3];
						if (!hex_grad_center(e - NTet, gn)) {
							for (int a = 0; a < 8; a++) {
								for (int c = 0; c < 3; c++) {
									pr[c] += ur[ne + nd[a]] * gn[a][c];
									pi2[c] += ui[ne + nd[a]] * gn[a][c];
								}
							}
						}
					}
					else {
						double gn[6][3];
						if (!prism_grad_center(e - NTet - NHex, gn)) {
							for (int a = 0; a < 6; a++) {
								for (int c = 0; c < 3; c++) {
									pr[c] += ur[ne + nd[a]] * gn[a][c];
									pi2[c] += ui[ne + nd[a]] * gn[a][c];
								}
							}
						}
					}
				}

				const double sg2 = Material[edge_elem_mat(e)].sigma;
				for (int c = 0; c < 3; c++) {
					bv[(e * 3) + c] = brv[c];		// |B| は実部で代表させる
					// J = -σ(jωA + ∇φ)、jωA = jω(ar + j ai) = -ω ai + j ω ar
					jr[(e * 3) + c] = -sg2 * ((-omega * aiv[c]) + pr[c]);
					ji[(e * 3) + c] = -sg2 * ((omega * arv[c]) + pi2[c]);
				}
			}
			field_add_cellvec("B_A_re", bv);
			field_add_cellvec("J_A_re", jr);
			field_add_cellvec("J_A_im", ji);
			free(bv); free(jr); free(ji);
		}

	}
	if (final || ierr) break;

	// ν の更新 (H を独立変数にする) と Anderson 加速 (ln ν)
	double dmax = 0;
	for (int e = 0; e < nelem; e++) {
		const material_t *mt = &Material[edge_elem_mat(e)];
		if (mt->nbh[0] <= 0) {
			aay[e] = aaf[e] = 0;
			continue;
		}
		double bre[3], bim[3];
		elem_bfield(e, ur, ui, bre, bim);
		const double b2 = (bre[0] * bre[0]) + (bre[1] * bre[1]) + (bre[2] * bre[2])
			+ (bim[0] * bim[0]) + (bim[1] * bim[1]) + (bim[2] * bim[2]);
		const double hm = nuelem[e] * sqrt(b2);
		const double nun = ((hm > 0) ? (hm / bh_b_of_h(mt, 0, hm)) : bh_nu(mt, 0, 0));
		const double dn = fabs(nun - nuelem[e]) / nun;
		if (dn > dmax) dmax = dn;
		aay[e] = log(nuelem[e]);
		aaf[e] = log(nun) - aay[e];
	}
	anderson_step(nelem, aay, aaf, aayp, aafp, aady, aadf, AAM3, &aan, NlRelax);
	for (int e = 0; e < nelem; e++) {
		if (Material[edge_elem_mat(e)].nbh[0] > 0) nuelem[e] = exp(aay[e]);
	}
	{
		const double y2 = (yr[0] * yr[0]) + (yi[0] * yi[0]);
		fprintf(fp_log, "  nl-it %4d : R = %.9e, L = %.9e, max dnu/nu = %.3e\n", nlit,
			yr[0] / y2, (-yi[0] / y2) / omega, dmax);
		fflush(fp_log);
	}
	if (dmax < NlTol) {
		final = 1;
	}
	else if (nlit >= NlMaxiter) {
		fprintf(fp_log, "*** nonlinear (B-H) eddy current : the effective reluctivity did not "
			"converge in %d iterations (max dnu/nu = %.3e)\n", NlMaxiter, dmax);
		ierr = 1;
	}
	}
	if (nlf && !ierr) fprintf(fp_log, "  nonlinear (B-H) : converged\n");
	free(nuelem); free(aay); free(aaf); free(aayp); free(aafp); free(aady); free(aadf);

	// Z = inv(Y)
	double *zr = (double *)malloc((size_t)np * np * sizeof(double));
	double *zi = (double *)malloc((size_t)np * np * sizeof(double));
	if (mat_inverse_c(yr, yi, zr, zi, np)) {
		fprintf(fp_log, "*** admittance matrix is singular; R(f)/L(f) is not available\n");
		ierr = 1;
	}
	else {
		const double scale = ((Tline && (TlineLength > 0)) ? TlineLength : 1);
		for (int i = 0; i < np * np; i++) {
			Rfmat[i] = zr[i] / scale;
			Lfmat[i] = zi[i] / (omega * scale);
		}
		HaveF = 1;

		// 受動性の検査。R < 0 や L < 0 は境界条件が物理的でないことを意味する
		// (例 : 外部境界を全部磁気壁にすると囲む電流が 0 に強制され、
		//  端子から電流を流すという条件と矛盾する -> awall を指定していないケース)
		for (int k = 0; k < np; k++) {
			const double rr = zr[(k * np) + k], ll = zi[(k * np) + k];
			if ((rr <= 0) || (ll < 0)) {
				fprintf(fp_log, "*** port %d gives a non-passive impedance "
					"(R = %.6e, L = %.6e); check the boundary conditions "
					"(an all-magnetic-wall boundary forces zero enclosed current, "
					"which contradicts driving a terminal current -- use awall)\n",
					k + 1, rr, ll / omega);
				ierr = 1;
			}
		}
		fprintf(fp_log, "  Z(1,1) = %.6e %+.6e j [ohm], R = %.6e, L = %.6e [H]\n",
			zr[0], zi[0], zr[0], zi[0] / omega);
	}

	free(zr); free(zi);
	free(yr); free(yi);
	free(xdr); free(xdi);
	free(br); free(bi);
	free(ur); free(ui);
	free(w1); free(w2);
	free(fix);
	free(wall);
	free(tree);
	free(comp);
	free(cond);
	crs_free(&Kr);
	crs_free(&Ki);
	// **辺の情報も毎回解放する。** solve() は周波数掃引の点ごとに solve_eddy3d()
	// を呼び直すので、解放しないと点の数だけ漏れる (実測: 4 点の掃引で
	// definitely lost 476,520 バイト)。edge_free() はポインタを NULL に戻すので
	// 次の点で edge_build() をやり直しても安全
	edge_free();

	return ierr;
}



/*
---- 交流の発熱 (thermal + analysis = A) ----

熱連成 (thermal.c) の電気の段を A-φ で解く。σ(T) を要素ごとに与えて port 1 を
V = Volt で解き、端子電流が cur (波高値、実数) になるよう解を複素数倍して、
要素の時間平均の損失 Q_e = ½∫σ|E|² dV を作る。

	e_k = jω a_k + (Gφ)_k   (局所の向き、E = -Σ e_k W_k は Whitney 恒等式から厳密)
	Q_e = ½ eᴴ T_e e,   T_e = σ_e ∫W_k・W_l dV   (組み立てと同じ要素行列)

連成系の A の行 S a + T e = 0 と節点の行 Gᵀ T e = r (電極の反作用 = 端子電流) から
	eᴴ T e = jω aᴴ S a + φᴴ r,   Re(eᴴ T e) = Re(V I*)
なので Σ_e Q_e = ½ Re(V I*) = ½|I|² R が丸めと反復の打ち切りを除いて厳密
(rlc_check.sh で見るエネルギーの恒等式)。節点荷重は f_i = Σ_e Q_e / (要素の節点数)
(要素平均の発熱を一定とみなした 1 次四面体の整合な荷重。和は Σ Q_e のまま)。
ω -> 0 で DC の導電解析に帰着するので、低周波では DC の熱連成 (電流 I/√2) と
一致する (検証の恒等式)。損失はゲージ不変なので gauge は既定 (固定しない) でよい。
*/
static struct {
	crs_t Kr, Ki;
	unsigned char *fix, *cond;
	double *xdr, *br, *bi, *ur, *ui, *w1, *w2;
	double omega;
} EH;


/*
3 次元渦電流の過渡解析 (analysis = A + transient、後退 Euler、電圧駆動)。

時間領域の A-φ 式  ∇×ν∇×A + σ(∂A/∂t + ∇φ) = 0,  ∇・σ(∂A/∂t + ∇φ) = 0  を後退 Euler で
離散化し、φ の代わりに χ = Δt φ を未知数にすると係数行列が**実対称・半正定値**になる:

	[ S + T/Δt     TG/Δt    ] [A^{n+1}]   [ T A^n / Δt   ]
	[ GᵀT/Δt     GᵀTG/Δt    ] [ χ      ] = [ GᵀT A^n / Δt ]

(周波数領域で節点の式を jω で割ったのと同じ形。jω → 1/Δt)。系は (A, χ) -> (A + Gψ, χ - ψ) で
不変なので特異だが、右辺は値域に入るので CG で解ける (周波数領域と同じくゲージ固定しない)。
ポート 1 の電極を χ = Δt V(t^{n+1})、基準導体を 0 に固定し (V(t) = voltage x waveform)、電極の
節点の反作用から端子電流を取る。

**刻みの仕事の恒等式**: z = (A^{n+1} - A^n, χ) とすると zᵀ(Kx - M[A^n;0]) = zᵀ r で、
	V I Δt = (ΔA)ᵀ S A^{n+1} + zᵀ M z
(左辺は端子の仕事、右辺の第 1 項は磁気の仕事、第 2 項は Δt eᵀTe = Joule 損、e = (ΔA + Gχ)/Δt)。
反復の打ち切りを除いて厳密 (ログに出す)。DC の極限 (V 一定で十分長い時間) では ΔA = 0 になり、
節点の式は GᵀTG φ = 0 そのもの (= analysis R の導電行列、Whitney の恒等式) なので I = V / R_dc。
*/
static void assemble_eddy3d_tr(crs_t *K, crs_t *M, double dt, int skipja)
{
	crs_zero(K);
	crs_zero(M);
	const int ne = NEdge;
	const int nelem = edge_elem_count();
	for (int e = 0; e < nelem; e++) {
		const int m = edge_elem_mat(e);
		const int nedge = edge_elem_nedge(e);
		const int nen = edge_elem_nen(e);
		const int (*tab)[2] = edge_elem_table(e);
		double nu[6];
		material_coef_pub(m, 3, nu);
		const double sig = Material[m].sigma;
		double se[12][12], te[12][12];
		edge_elem_matrices(e, nu, sig, se, te);
		double gl[12][8];
		for (int k = 0; k < nedge; k++) {
			for (int a = 0; a < nen; a++) gl[k][a] = 0;
			gl[k][tab[k][0]] = -1;
			gl[k][tab[k][1]] = +1;
		}
		double tg[12][8], gtg[8][8];
		for (int k = 0; k < nedge; k++) {
			for (int a = 0; a < nen; a++) {
				double v = 0;
				for (int l = 0; l < nedge; l++) v += te[k][l] * gl[l][a];
				tg[k][a] = v;
			}
		}
		for (int a = 0; a < nen; a++) {
			for (int b = 0; b < nen; b++) {
				double v = 0;
				for (int k = 0; k < nedge; k++) v += gl[k][a] * tg[k][b];
				gtg[a][b] = v;
			}
		}
		const int32_t *ed = &TetEdge[EdgeOff[e]];
		const signed char *sg = &TetEdgeSgn[EdgeOff[e]];
		const int32_t *nd = edge_elem_nodes(e);
		for (int k = 0; k < nedge; k++) {
			for (int l = 0; l < nedge; l++) {
				const int64_t p = crs_find_en(K, ed[k], ed[l]);
				const double w = sg[k] * sg[l];
				K->val[p] += w * ((skipja && Material[m].ja.on) ? (te[k][l] / dt) : (se[k][l] + (te[k][l] / dt)));
				M->val[p] += w * te[k][l] / dt;
			}
		}
		for (int k = 0; k < nedge; k++) {
			for (int a = 0; a < nen; a++) {
				const double w = sg[k] * tg[k][a] / dt;
				const int64_t p1 = crs_find_en(K, ed[k], (int32_t)ne + nd[a]);
				const int64_t p2 = crs_find_en(K, (int32_t)ne + nd[a], ed[k]);
				K->val[p1] += w;
				K->val[p2] += w;
				M->val[p1] += w;
				M->val[p2] += w;
			}
		}
		for (int a = 0; a < nen; a++) {
			for (int b = 0; b < nen; b++) {
				const int64_t p = crs_find_en(K, (int32_t)ne + nd[a], (int32_t)ne + nd[b]);
				K->val[p] += gtg[a][b] / dt;
				M->val[p] += gtg[a][b] / dt;
			}
		}
	}
}


// 刻みの終わり t の電圧の倍率 (waveform の折れ線の左極限、無ければ 1。静磁場の過渡解析と同じ)
static double ed_wave(double t)
{
	if (WaveN < 1) return 1;
	int k = -1;
	for (int m = 0; m < WaveN; m++) {
		if (WaveT[m] < t) k = m;
	}
	if (k < 0) return WaveF[0];
	if (k == (WaveN - 1)) return WaveF[k];
	const double w = (t - WaveT[k]) / (WaveT[k + 1] - WaveT[k]);
	return WaveF[k] + (w * (WaveF[k + 1] - WaveF[k]));
}


/*
ヒステリシス (ja) の 3 次元過渡解析。J-A の材料の要素は**1 次四面体**に限る (∇×W が要素内で一定なので
B と J-A の状態が要素に 1 つ)。各刻みは静磁場の過渡解析と同じ固定点反復:

	F(x) = K_l x − M [A^n; 0] + Σ_e V_e H_e(B_e)・c_k   (K_l は J-A の要素の S を除いた行列)
	(K_l + S(N_e)) δ = −F,   N_e = dH/dB の対称部分 (テンソル、正定値)

H(B) は確定済みの状態から積分し直す。仕事の恒等式の磁気の仕事は Σ V_e H_e・ΔB_e。
J-A は常に 3 成分のベクトル J-A (ja_eval_vec3。場が 1 方向ならスカラーに厳密に退化する)。
**スカラーの射影 (断面 2 次元の既定) は 3 次元では使わない**: 1 次元の解でも Whitney の辺要素の B は
横方向の成分を持ち、帰路の中央 (B ≈ 0) で向きが乱れる。射影 B・d でも符号付きの |B| でも H の向きが
不連続になり、固定点が 1e-7〜1e-6 の残差で停滞した (実測)。ベクトル J-A は滑らかなので停滞しない。
*/
typedef struct {
	int nj;				// J-A の四面体の数
	int *el;			// [nj] 要素番号
	int32_t *ed;		// [nj*6] 辺 (全体番号)
	double *vol;		// [nj]
	double *cv;			// [nj*18] 符号込みの ∇×W_k (k = 0..5)
	double *b, *h, *m, *d;		// [nj*3] 確定済みの状態 (d は最初に磁化した向き)
	double *bn, *hn, *mn, *dn;	// [nj*3] 反復中
	double *nu6;		// [nj*6] 行列の傾き
} ed_ja_t;


static int ed_ja_open(ed_ja_t *J, FILE *fp_log)
{
	const int nelem = edge_elem_count();
	int nj = 0;
	for (int e = 0; e < nelem; e++) {
		if (!Material[edge_elem_mat(e)].ja.on) continue;
		if (edge_elem_kind(e) != MESHELEM_TET) {
			fprintf(fp_log, "*** hysteresis (ja) with analysis A + transient needs first-order tetrahedra "
				"for the hysteretic material\n");
			return 1;
		}
		nj++;
	}
	J->nj = nj;
	const size_t s1 = (size_t)(nj > 0 ? nj : 1);
	J->el = (int *)malloc(s1 * sizeof(int));
	J->ed = (int32_t *)malloc(s1 * 6 * sizeof(int32_t));
	J->vol = (double *)malloc(s1 * sizeof(double));
	J->cv = (double *)malloc(s1 * 18 * sizeof(double));
	double **v3[8] = {&J->b, &J->h, &J->m, &J->d, &J->bn, &J->hn, &J->mn, &J->dn};
	for (int q = 0; q < 8; q++) *v3[q] = (double *)calloc(s1 * 3, sizeof(double));
	J->nu6 = (double *)calloc(s1 * 6, sizeof(double));
	int j = 0;
	for (int e = 0; e < nelem; e++) {
		if (!Material[edge_elem_mat(e)].ja.on) continue;
		double g[4][3], vol;
		J->el[j] = e;
		J->vol[j] = 0;
		for (int q = 0; q < 18; q++) J->cv[(j * 18) + q] = 0;
		const int32_t *ed = &TetEdge[EdgeOff[e]];
		const signed char *sg = &TetEdgeSgn[EdgeOff[e]];
		for (int k = 0; k < 6; k++) J->ed[(j * 6) + k] = ed[k];
		if (!tet_grad_pub(&Tet[e * 4], g, &vol)) {
			J->vol[j] = vol;
			for (int k = 0; k < 6; k++) {
				const int a = EDGE_NODE3[k][0], c = EDGE_NODE3[k][1];
				double *cv = &J->cv[(j * 18) + (k * 3)];
				cv[0] = 2 * sg[k] * ((g[a][1] * g[c][2]) - (g[a][2] * g[c][1]));
				cv[1] = 2 * sg[k] * ((g[a][2] * g[c][0]) - (g[a][0] * g[c][2]));
				cv[2] = 2 * sg[k] * ((g[a][0] * g[c][1]) - (g[a][1] * g[c][0]));
			}
		}
		j++;
	}
	fprintf(fp_log, "  hysteresis : vector J-A on %d tetrahedra (nlsolver = %d %.3e %.3f)\n",
		nj, NlMaxiter, NlTol, NlRelax);
	// 3 成分のベクトル J-A の自己検証 (向きを固定した経路でスカラーの J-A と一致)
	for (int m = 0; m < NMaterial; m++) {
		if (!Material[m].ja.on) continue;
		fprintf(fp_log, "  J-A model check (material %d) : 3-D vector J-A along a fixed direction == scalar J-A : "
			"max |dH| / max |H| = %.3e\n", m, ja_vec3_selftest(&Material[m].ja));
	}
	return 0;
}


static void ed_ja_close(ed_ja_t *J)
{
	free(J->el); free(J->ed); free(J->vol); free(J->cv);
	free(J->b); free(J->h); free(J->m); free(J->d);
	free(J->bn); free(J->hn); free(J->mn); free(J->dn);
	free(J->nu6);
}


// x の A から各要素の B を作り、確定済みの状態から H (反復中の状態) と行列の傾きを求める
static void ed_ja_eval(ed_ja_t *J, const double *x)
{
	const int nj = J->nj;
	int j;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (j = 0; j < nj; j++) {
		const ja_t *ja = &Material[edge_elem_mat(J->el[j])].ja;
		double bv[3] = {0, 0, 0};
		for (int k = 0; k < 6; k++) {
			const double a = x[J->ed[(j * 6) + k]];
			for (int c = 0; c < 3; c++) bv[c] += a * J->cv[(j * 18) + (k * 3) + c];
		}
		const double bmag = sqrt((bv[0] * bv[0]) + (bv[1] * bv[1]) + (bv[2] * bv[2]));
		const double *dir = &J->d[j * 3];
		double *dn = &J->dn[j * 3];
		const int haved = ((dir[0] != 0) || (dir[1] != 0) || (dir[2] != 0));
		for (int c = 0; c < 3; c++) dn[c] = ((bmag > 1e-12) ? (bv[c] / bmag) : 0);
		double *nu6 = &J->nu6[j * 6];
		double *hv = &J->hn[j * 3];
		double dref[3] = {1, 0, 0};
		if (haved) { dref[0] = dir[0]; dref[1] = dir[1]; dref[2] = dir[2]; }
		else if (bmag > 1e-12) { dref[0] = dn[0]; dref[1] = dn[1]; dref[2] = dn[2]; }
		ja_eval_vec3(ja, &J->b[j * 3], &J->h[j * 3], &J->m[j * 3], bv, dref, hv, &J->mn[j * 3], nu6);
		for (int c = 0; c < 3; c++) J->bn[(j * 3) + c] = bv[c];
	}
}


// 刻みの磁気の仕事 Σ V_e H_e・(B_e − B_e^n) (J-A の要素の分)
static double ed_ja_work(const ed_ja_t *J)
{
	double w = 0;
	for (int j = 0; j < J->nj; j++) {
		double s = 0;
		for (int c = 0; c < 3; c++) s += J->hn[(j * 3) + c] * (J->bn[(j * 3) + c] - J->b[(j * 3) + c]);
		w += J->vol[j] * s;
	}
	return w;
}


static void ed_ja_commit(ed_ja_t *J)
{
	const int nj = J->nj;
	memcpy(J->b, J->bn, (size_t)nj * 3 * sizeof(double));
	memcpy(J->h, J->hn, (size_t)nj * 3 * sizeof(double));
	memcpy(J->m, J->mn, (size_t)nj * 3 * sizeof(double));
	for (int j = 0; j < nj; j++) {
		double *d = &J->d[j * 3];
		if ((d[0] == 0) && (d[1] == 0) && (d[2] == 0)) {
			for (int c = 0; c < 3; c++) d[c] = J->dn[(j * 3) + c];
		}
	}
}


// 残差 r = −F(x) (固定の行は 0) と、その相対ノルム
static double ed_ja_resid(const ed_ja_t *J, const crs_t *Kl, const double *rhs, const unsigned char *fix,
	const double *x, double *r, double bref)
{
	const int n = (int)Kl->n;
	crs_spmv(Kl, x, r, NULL);
	for (int i = 0; i < n; i++) r[i] = rhs[i] - r[i];
	for (int j = 0; j < J->nj; j++) {
		const double *hv = &J->hn[j * 3];
		for (int k = 0; k < 6; k++) {
			const double *cv = &J->cv[(j * 18) + (k * 3)];
			r[J->ed[(j * 6) + k]] -= J->vol[j] * ((hv[0] * cv[0]) + (hv[1] * cv[1]) + (hv[2] * cv[2]));
		}
	}
	double s = 0;
	for (int i = 0; i < n; i++) {
		if (fix[i]) r[i] = 0;
		s += r[i] * r[i];
	}
	return ((bref > 0) ? (sqrt(s) / bref) : sqrt(s));
}


// 1 刻みの固定点反復。x は前の刻みの解で入り、収束した解で出る。戻り値は反復回数 (失敗は -1)
static int ed_ja_step(ed_ja_t *J, const crs_t *Kl, crs_t *Kt, const double *rhs, const unsigned char *fix,
	const double *xd, double *x, double *r, double *del, double *xbest, double bref, int st, FILE *fp_log)
{
	const int n = (int)Kl->n;
	const int64_t nnz = Kl->rowptr[n];
	for (int i = 0; i < n; i++) {
		if (fix[i]) x[i] = xd[i];
	}
	int it, conv = 0, nback = 0;
	double rnorm = 0, rbest = 0, lrel = NlRelax;
	for (it = 1; it <= NlMaxiter; it++) {
		ed_ja_eval(J, x);
		rnorm = ed_ja_resid(J, Kl, rhs, fix, x, r, bref);
		if (rnorm < NlTol) {
			conv = 1;
			break;
		}
		// 発散の検出。H(B) は確定済みの状態で枝が切り替わる折れ点を持ち、反転の近くでは固定点が
		// 2 つの値を往復する (実測: 残差 0.05 と 0.12 の往復が 40 回続いた)。往復は最良の 4 倍程度に
		// 収まるので、静磁場の過渡解析 (mag_transient) の「10 倍」では検出できない。3 次元には
		// 準静的な掃引との一致という制約が無いので、**残差が最良を越えた最初の回**に最良の反復へ
		// 戻して両方の枝の大きい方の傾きにする (以後 5 回程度で収束)。それでも 10 倍を越えたら
		// 更新を半分にする
		if ((it == 1) || (rnorm < rbest)) {
			rbest = rnorm;
			memcpy(xbest, x, (size_t)n * sizeof(double));
		}
		else if ((!JaSafeSlope || (rnorm > (10 * rbest))) && (nback < 30)) {
			memcpy(x, xbest, (size_t)n * sizeof(double));
			if (!JaSafeSlope) JaSafeSlope = 1;
			else lrel /= 2;
			nback++;
			continue;
		}
		memcpy(Kt->val, Kl->val, (size_t)nnz * sizeof(double));
		for (int j = 0; j < J->nj; j++) {
			const double *nu = &J->nu6[j * 6];
			const double *cv = &J->cv[j * 18];
			for (int k = 0; k < 6; k++) {
				const double *ck = &cv[k * 3];
				// ν c_k (対称テンソル xx yy zz xy yz zx)
				const double nc[3] = {
					(nu[0] * ck[0]) + (nu[3] * ck[1]) + (nu[5] * ck[2]),
					(nu[3] * ck[0]) + (nu[1] * ck[1]) + (nu[4] * ck[2]),
					(nu[5] * ck[0]) + (nu[4] * ck[1]) + (nu[2] * ck[2])};
				for (int l = 0; l < 6; l++) {
					const double *cl = &cv[l * 3];
					const int64_t p = crs_find_en(Kt, J->ed[(j * 6) + k], J->ed[(j * 6) + l]);
					Kt->val[p] += J->vol[j] * ((nc[0] * cl[0]) + (nc[1] * cl[1]) + (nc[2] * cl[2]));
				}
			}
		}
		// 内側の許容は刻みの駆動 bref に対する絶対値にする (converg × bref)。反復が進むと r が小さく
		// なり、ゲージの零空間に丸めで漏れる成分 (∇×∇ψ = 0 の丸め) が ‖r‖ に対して無視できなく
		// なるので、‖r‖ に対する相対の許容では CG が止まらない (実測、3 回目で非収束)
		double tol = ((rnorm > 0) ? (Solver.converg / rnorm) : Solver.converg);
		if (tol > 0.1) tol = 0.1;
		if (tol < Solver.converg) tol = Solver.converg;
		for (int i = 0; i < n; i++) del[i] = 0;
		const int li = solver_cg(Kt, r, del, fix, Solver.maxiter, Solver.nout, tol, NULL, "a3j");
		if (li < 0) {
			fprintf(fp_log, "*** linear solver did not converge (step %d)\n", st);
			JaSafeSlope = 0;
			return -1;
		}
		for (int i = 0; i < n; i++) x[i] += lrel * del[i];
	}
	JaSafeSlope = 0;
	if (!conv) {
		fprintf(fp_log, "*** transient : the iteration did not converge at step %d (residual = %.3e)\n", st, rnorm);
		return -1;
	}
	return it;
}


int solve_eddy3d_transient(FILE *fp_log)
{
	int ierr = 0;
	fprintf(fp_log, "\n=== 3D eddy current, transient (A-phi, edge elements, backward Euler) ===\n");
	if (!MeshMode || (NPort != 1)) {
		fprintf(fp_log, "*** analysis A + transient needs an unstructured mesh and exactly one port\n");
		return 1;
	}
	edge_build();
	const int ne = NEdge, nn = NNode, n = ne + nn;
	const int nelem = edge_elem_count();
	const double dt = TransientTime / TransientSteps;
	int hyst = 0;
	for (int m = 0; m < NMaterial; m++) {
		if (Material[m].ja.on) hyst = 1;
	}
	fprintf(fp_log, "  nodes = %d, elements = %d, edges = %d, unknowns = %d\n", nn, nelem, ne, n);
	fprintf(fp_log, "  t_end = %.6e [s], %d steps (dt = %.6e [s]), voltage = %.6e [V]%s\n",
		TransientTime, TransientSteps, dt, Volt, ((WaveN > 0) ? " x waveform" : " (step)"));

	crs_t K, M;
	crs_alloc_edge_node(&K);
	crs_alloc_edge_node(&M);
	assemble_eddy3d_tr(&K, &M, dt, hyst);

	unsigned char *fix = (unsigned char *)calloc((size_t)n, sizeof(unsigned char));
	unsigned char *wall = (unsigned char *)calloc((size_t)ne, sizeof(unsigned char));
	awall_mark(wall);
	int nwall = 0;
	for (int i = 0; i < ne; i++) {
		fix[i] = wall[i];
		nwall += wall[i];
	}
	unsigned char *cond = (unsigned char *)calloc((size_t)nn, sizeof(unsigned char));
	for (int e = 0; e < nelem; e++) {
		if (Material[edge_elem_mat(e)].sigma <= 0) continue;
		const int32_t *nd = edge_elem_nodes(e);
		const int nen2 = edge_elem_nen(e);
		for (int l = 0; l < nen2; l++) cond[nd[l]] = 1;
	}
	int nelec = 0;
	for (int i = 0; i < nn; i++) {
		if (!cond[i]) fix[ne + i] = 1;
		else if (NodeConductor[i] >= 0) {
			fix[ne + i] = 1;
			nelec++;
		}
	}
	fprintf(fp_log, "  fixed : %d wall edges (A_t = 0), %d electrode nodes\n", nwall, nelec);

	double *x = (double *)calloc((size_t)n, sizeof(double));		// [A^{n+1}; χ]
	double *xo = (double *)calloc((size_t)n, sizeof(double));	// [A^n; 0]
	double *xd = (double *)calloc((size_t)n, sizeof(double));
	double *rhs = (double *)malloc((size_t)n * sizeof(double));
	double *b = (double *)malloc((size_t)n * sizeof(double));
	double *w1 = (double *)malloc((size_t)n * sizeof(double));
	double *u = (double *)calloc((size_t)n, sizeof(double));
	double *z = (double *)malloc((size_t)n * sizeof(double));
	ed_ja_t J;
	memset(&J, 0, sizeof(J));
	if (hyst && ed_ja_open(&J, fp_log)) ierr = 1;
	crs_t Kt;
	memset(&Kt, 0, sizeof(Kt));
	if (hyst) crs_alloc_edge_node(&Kt);
	double bref = 0;
	long fptotal = 0;

	FILE *fcsv = fopen("ofe_transient.csv", "w");
	if (fcsv) fprintf(fcsv, "t [s],V [V],I [A],W_in [J],W_mag [J],W_joule [J],iter\n");
	double win = 0, wmag = 0, wjou = 0;
	fprintf(fp_log, "  %6s %13s %13s %13s %6s\n", "step", "t [s]", "V [V]", "I [A]", "iter");
	for (int st = 1; (st <= TransientSteps) && !ierr; st++) {
		const double t = st * dt;
		const double v = Volt * ed_wave(t);
		// 右辺 M [A^n; 0] と Dirichlet の持ち上げ
		for (int i = 0; i < n; i++) xo[i] = ((i < ne) ? x[i] : 0);
		crs_spmv(&M, xo, rhs, NULL);
		for (int i = 0; i < n; i++) xd[i] = 0;
		for (int i = 0; i < nn; i++) {
			if (fix[ne + i] && cond[i] && (NodeConductor[i] == 1)) xd[ne + i] = dt * v;
		}
		crs_spmv(&K, xd, w1, NULL);
		for (int i = 0; i < n; i++) b[i] = (fix[i] ? 0 : (rhs[i] - w1[i]));
		int it;
		if (!hyst) {
			// 初期値は前の刻みの解 (固定の行は 0)
			for (int i = 0; i < n; i++) u[i] = (fix[i] ? 0 : x[i]);
			it = solver_cg(&K, b, u, fix, Solver.maxiter, Solver.nout, Solver.converg, NULL, "a3t");
			if (it < 0) {
				fprintf(fp_log, "*** linear solver did not converge (step %d)\n", st);
				ierr = 1;
				break;
			}
			for (int i = 0; i < n; i++) x[i] = u[i] + xd[i];
		}
		else {
			double bn = 0;
			for (int i = 0; i < n; i++) bn += b[i] * b[i];
			bn = sqrt(bn);
			if (bn > bref) bref = bn;
			it = ed_ja_step(&J, &K, &Kt, rhs, fix, xd, x, b, u, z, bref, st, fp_log);
			if (it < 0) {
				ierr = 1;
				break;
			}
			fptotal += it;
		}
		// 端子電流 : 電極の節点の反作用 (K x - M [A^n; 0])
		crs_spmv(&K, x, w1, NULL);
		double cur = 0;
		for (int i = 0; i < nn; i++) {
			if (fix[ne + i] && cond[i] && (NodeConductor[i] == 1)) cur += w1[ne + i] - rhs[ne + i];
		}
		// 仕事 : 端子 V I Δt、Joule zᵀ M z、磁気 (ΔA)ᵀ S A^{n+1} = zᵀ (K - M) x の A 部分
		for (int i = 0; i < n; i++) z[i] = x[i] - xo[i];
		double zmz = 0, zkx = 0;
		crs_spmv(&M, z, w1, NULL);
		for (int i = 0; i < n; i++) zmz += z[i] * w1[i];
		double zmx = 0;
		crs_spmv(&M, x, w1, NULL);
		for (int i = 0; i < ne; i++) zmx += z[i] * w1[i];
		crs_spmv(&K, x, w1, NULL);
		for (int i = 0; i < ne; i++) zkx += z[i] * w1[i];
		if (hyst) zkx += ed_ja_work(&J);
		win += v * cur * dt;
		wjou += zmz;
		wmag += zkx - zmx;
		if (hyst) ed_ja_commit(&J);
		if (fcsv) fprintf(fcsv, "%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%d\n", t, v, cur, win, wmag, wjou, it);
		if ((st == TransientSteps) || ((st % ((TransientSteps >= 20) ? (TransientSteps / 20) : 1)) == 0)) {
			fprintf(fp_log, "  %6d %13.6e %13.6e %13.6e %6d\n", st, t, v, cur, it);
		}
	}
	if (fcsv) fclose(fcsv);
	if (hyst && !ierr) {
		fprintf(fp_log, "  hysteresis (vector J-A on %d tetrahedra) : fixed-point iterations = %ld (%d steps)\n",
			J.nj, fptotal, TransientSteps);
	}
	fprintf(fp_log, "  energy : terminals = %.12e, magnetic work = %.12e, joule = %.12e [J] "
		"(terminals - magnetic - joule = %.3e)\n", win, wmag, wjou, win - wmag - wjou);
	free(x); free(xo); free(xd); free(rhs); free(b); free(w1); free(u); free(z);
	free(fix); free(wall); free(cond);
	if (hyst) {
		ed_ja_close(&J);
		crs_free(&Kt);
	}
	crs_free(&K);
	crs_free(&M);
	edge_free();
	return ierr;
}


int eddy3d_heat_open(FILE *fp_log)
{
	edge_build();

	const int ne = NEdge;
	const int nn = NNode;
	const int n = ne + nn;
	const int nelem = edge_elem_count();
	EH.omega = 2 * PI * Freq;
	fprintf(fp_log, "  AC heating (analysis = A) : frequency = %.6e [Hz], %d edges, "
		"%d unknowns\n", Freq, ne, n);

	crs_alloc_edge_node(&EH.Kr);
	crs_alloc_edge_node(&EH.Ki);

	// Dirichlet は solve_eddy3d と同じ : awall の辺、電極の節点、非導電領域の節点
	EH.fix = (unsigned char *)calloc((size_t)n, sizeof(unsigned char));
	unsigned char *wall = (unsigned char *)calloc((size_t)ne, sizeof(unsigned char));
	awall_mark(wall);
	int nwall = 0;
	for (int i = 0; i < ne; i++) {
		EH.fix[i] = wall[i];
		nwall += wall[i];
	}
	int ierr = 0;
	if (GaugeTree) {
		unsigned char *tree = (unsigned char *)malloc((size_t)ne * sizeof(unsigned char));
		int *comp = (int *)malloc((size_t)nn * sizeof(int));
		edge_tree_gauge(wall, tree, comp);
		int bad = 0, ref = -1;
		for (int i = 0; i < ne; i++) EH.fix[i] |= tree[i];
		for (int i = 0; i < nn; i++) {
			if (NodeConductor[i] < 0) continue;
			if (ref < 0) ref = comp[i];
			else if (comp[i] != ref) bad = 1;
		}
		if (bad) {
			fprintf(fp_log, "*** the electrodes are not connected through awall edges, "
				"so the tree gauge would shift the terminal voltage; use gauge = 0\n");
			ierr = 1;
		}
		free(tree);
		free(comp);
	}
	free(wall);
	EH.cond = (unsigned char *)calloc((size_t)nn, sizeof(unsigned char));
	for (int e = 0; e < nelem; e++) {
		if (Material[edge_elem_mat(e)].sigma <= 0) continue;
		const int32_t *nd = edge_elem_nodes(e);
		const int nen = edge_elem_nen(e);
		for (int l = 0; l < nen; l++) EH.cond[nd[l]] = 1;
	}
	for (int i = 0; i < nn; i++) {
		if (!EH.cond[i] || (NodeConductor[i] >= 0)) EH.fix[ne + i] = 1;
	}
	if (nwall == 0) {
		fprintf(fp_log, "*** warning : no awall surface; the outer boundary is a "
			"magnetic wall (dA/dn = 0)\n");
	}

	EH.xdr = (double *)malloc((size_t)n * sizeof(double));
	EH.br  = (double *)malloc((size_t)n * sizeof(double));
	EH.bi  = (double *)malloc((size_t)n * sizeof(double));
	EH.ur  = (double *)malloc((size_t)n * sizeof(double));
	EH.ui  = (double *)malloc((size_t)n * sizeof(double));
	EH.w1  = (double *)malloc((size_t)n * sizeof(double));
	EH.w2  = (double *)malloc((size_t)n * sizeof(double));

	return ierr;
}


// σ = sigelem で解き、端子電流の波高値 cur での節点荷重 f[NNode] と Z = zr + j zi を返す
int eddy3d_heat_solve(FILE *fp_log, const double *sigelem, double cur, double *f,
	double *zr, double *zi, int it)
{
	const int ne = NEdge;
	const int nn = NNode;
	const int n = ne + nn;
	const int nelem = edge_elem_count();
	const double omega = EH.omega;

	assemble_eddy3d_sig(&EH.Kr, &EH.Ki, omega, NULL, sigelem);
	for (int i = 0; i < n; i++) EH.xdr[i] = 0;
	for (int i = 0; i < nn; i++) {
		if (EH.fix[ne + i] && (NodeConductor[i] == 1)) EH.xdr[ne + i] = Volt;
	}
	crs_spmv(&EH.Kr, EH.xdr, EH.w1, NULL);
	crs_spmv(&EH.Ki, EH.xdr, EH.w2, NULL);
	for (int i = 0; i < n; i++) {
		EH.br[i] = (EH.fix[i] ? 0 : -EH.w1[i]);
		EH.bi[i] = (EH.fix[i] ? 0 : -EH.w2[i]);
	}
	char label[64];
	sprintf(label, "A%d", it);
	if (solver_cocg(&EH.Kr, &EH.Ki, 1.0, EH.br, EH.bi, EH.ur, EH.ui, EH.fix,
			Solver.maxiter, Solver.nout, Solver.converg, NULL, label) < 0) {
		fprintf(fp_log, "*** thermal : the A-phi solve did not converge (iteration %d)\n", it);
		return 1;
	}
	for (int i = 0; i < n; i++) EH.ur[i] += EH.xdr[i];

	// 端子電流 (節点の行は 1/(jω) 倍してあるので jω 倍して戻す)
	double ir = 0, ii = 0;
	for (int i = 0; i < nn; i++) {
		if ((NodeConductor[i] != 1) || !EH.cond[i]) continue;
		double rr, ri;
		row_dot_c(&EH.Kr, &EH.Ki, (int32_t)ne + i, EH.ur, EH.ui, &rr, &ri);
		ir += -omega * ri;
		ii += omega * rr;
	}
	const double den = (ir * ir) + (ii * ii);
	if (den <= 0) {
		fprintf(fp_log, "*** thermal : no current flows through the port\n");
		return 1;
	}
	*zr = Volt * ir / den;
	*zi = -Volt * ii / den;

	// 端子電流が cur (実数) になるよう s = cur / I_V 倍する
	const double sr = cur * ir / den, si = -cur * ii / den;
	for (int i = 0; i < n; i++) {
		const double xr = EH.ur[i], xi = EH.ui[i];
		EH.ur[i] = (sr * xr) - (si * xi);
		EH.ui[i] = (sr * xi) + (si * xr);
	}

	// 要素の損失 Q_e = ½ eᴴ T_e e を節点に配る
	for (int i = 0; i < nn; i++) f[i] = 0;
	for (int e = 0; e < nelem; e++) {
		const double sg = sigelem[e];
		if (sg <= 0) continue;
		const int m = edge_elem_mat(e);
		const int nedge = edge_elem_nedge(e);
		const int nen = edge_elem_nen(e);
		const int (*tab)[2] = edge_elem_table(e);
		const int32_t *ed = &TetEdge[EdgeOff[e]];
		const signed char *sgn = &TetEdgeSgn[EdgeOff[e]];
		const int32_t *nd = edge_elem_nodes(e);
		double nu[6], se[12][12], te[12][12];
		material_coef_pub(m, 3, nu);
		edge_elem_matrices(e, nu, sg, se, te);
		double er[12], ei[12];
		for (int k = 0; k < nedge; k++) {
			const int32_t p = ne + nd[tab[k][0]], q = ne + nd[tab[k][1]];
			// jω a = -ω a_i + j ω a_r
			er[k] = (-omega * sgn[k] * EH.ui[ed[k]]) + (EH.ur[q] - EH.ur[p]);
			ei[k] = (omega * sgn[k] * EH.ur[ed[k]]) + (EH.ui[q] - EH.ui[p]);
		}
		double qe = 0;
		for (int k = 0; k < nedge; k++) {
			for (int l = 0; l < nedge; l++) {
				qe += te[k][l] * ((er[k] * er[l]) + (ei[k] * ei[l]));
			}
		}
		qe *= 0.5;
		for (int a = 0; a < nen; a++) f[nd[a]] += qe / nen;
	}

	return 0;
}


void eddy3d_heat_close(void)
{
	crs_free(&EH.Kr);
	crs_free(&EH.Ki);
	free(EH.fix); free(EH.cond);
	free(EH.xdr); free(EH.br); free(EH.bi); free(EH.ur); free(EH.ui);
	free(EH.w1); free(EH.w2);
	EH.fix = EH.cond = NULL;
	EH.xdr = EH.br = EH.bi = EH.ur = EH.ui = EH.w1 = EH.w2 = NULL;
	edge_free();
}
