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

	// 自由度毎の四面体数
	int *cnt = (int *)malloc((size_t)n * sizeof(int));
	memset(cnt, 0, (size_t)n * sizeof(int));
	for (int e = 0; e < NTet; e++) {
		for (int k = 0; k < 6; k++) cnt[TetEdge[(e * 6) + k]]++;
		for (int l = 0; l < 4; l++) cnt[ne + Tet[(e * 4) + l]]++;
	}
	int64_t *eptr = (int64_t *)malloc(((size_t)n + 1) * sizeof(int64_t));
	eptr[0] = 0;
	for (int i = 0; i < n; i++) eptr[i + 1] = eptr[i] + cnt[i];
	int32_t *elist = (int32_t *)malloc((size_t)eptr[n] * sizeof(int32_t));
	memset(cnt, 0, (size_t)n * sizeof(int));
	for (int e = 0; e < NTet; e++) {
		for (int k = 0; k < 6; k++) {
			const int32_t d = TetEdge[(e * 6) + k];
			elist[eptr[d] + cnt[d]] = (int32_t)e;
			cnt[d]++;
		}
		for (int l = 0; l < 4; l++) {
			const int32_t d = ne + Tet[(e * 4) + l];
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
		const int need = (int)(p1 - p0) * 10;
		if (need > cap) {
			cap = need;
			work = (int32_t *)realloc(work, (size_t)cap * sizeof(int32_t));
		}
		int m = 0;
		for (int64_t p = p0; p < p1; p++) {
			const int32_t e = elist[p];
			for (int k = 0; k < 6; k++) work[m++] = TetEdge[(e * 6) + k];
			for (int l = 0; l < 4; l++) work[m++] = (int32_t)ne + Tet[(e * 4) + l];
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
			for (int k = 0; k < 6; k++) work[m++] = TetEdge[(e * 6) + k];
			for (int l = 0; l < 4; l++) work[m++] = (int32_t)ne + Tet[(e * 4) + l];
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
void assemble_eddy3d(crs_t *Kr, crs_t *Ki, double omega)
{
	crs_zero(Kr);
	crs_zero(Ki);

	const int ne = NEdge;

	for (int e = 0; e < NTet; e++) {
		const int m = TetMat[e];
		double nu[6];
		material_coef_pub(m, 3, nu);			// ν = (μ0 μ~)^-1
		const double sig = Material[m].sigma;

		double se[6][6], te[6][6];
		edge_element(e, nu, sig, se, te);

		// 局所の離散勾配 GL[k][a] : 辺 k = (局所節点 p, q) -> -1 at p, +1 at q
		double gl[6][4];
		for (int k = 0; k < 6; k++) {
			for (int a = 0; a < 4; a++) gl[k][a] = 0;
			gl[k][EDGE_NODE3[k][0]] = -1;
			gl[k][EDGE_NODE3[k][1]] = +1;
		}

		// tg = te * GL (6x4)、gtg = GL^T * te * GL (4x4)
		double tg[6][4], gtg[4][4];
		for (int k = 0; k < 6; k++) {
			for (int a = 0; a < 4; a++) {
				double s = 0;
				for (int l = 0; l < 6; l++) s += te[k][l] * gl[l][a];
				tg[k][a] = s;
			}
		}
		for (int a = 0; a < 4; a++) {
			for (int b = 0; b < 4; b++) {
				double s = 0;
				for (int k = 0; k < 6; k++) s += gl[k][a] * tg[k][b];
				gtg[a][b] = s;
			}
		}

		const int32_t *ed = &TetEdge[e * 6];
		const signed char *sg = &TetEdgeSgn[e * 6];
		const int32_t *nd = &Tet[e * 4];

		// (1,1) 辺 x 辺
		for (int k = 0; k < 6; k++) {
			for (int l = 0; l < 6; l++) {
				const int64_t p = crs_find_en(Kr, ed[k], ed[l]);
				const double w = sg[k] * sg[l];
				Kr->val[p] += w * se[k][l];
				Ki->val[p] += w * omega * te[k][l];
			}
		}

		// (1,2) 辺 x 節点、(2,1) 節点 x 辺 (実部のみ)
		for (int k = 0; k < 6; k++) {
			for (int a = 0; a < 4; a++) {
				const double w = sg[k] * tg[k][a];
				Kr->val[crs_find_en(Kr, ed[k], (int32_t)ne + nd[a])] += w;
				Kr->val[crs_find_en(Kr, (int32_t)ne + nd[a], ed[k])] += w;
			}
		}

		// (2,2) 節点 x 節点 (虚部のみ、-G^T T G / ω)
		for (int a = 0; a < 4; a++) {
			for (int b = 0; b < 4; b++) {
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
	for (int e = 0; e < NTet; e++) {
		const material_t *mt = &Material[TetMat[e]];
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
		return 1;
	}
	const double delta = sqrt(2 / (omega * MU0 * sigmumax));
	fprintf(fp_log, "  nodes = %d, tetrahedra = %d, edges = %d, unknowns = %d\n",
		nn, NTet, ne, n);
	fprintf(fp_log, "  frequency = %.6e [Hz], skin depth = %.4e [m]\n", Freq, delta);

	// 導体要素の最大辺長を表皮深さと並べて報告する。
	// 1 次 Nedelec 要素の誤差は「場が変化する方向の刻み」ではなく
	// **要素の最大寸法**で決まるので、扁平な要素は場が一様な方向にも効く
	// (実測 : 場が x に一様でも dx を 1.7mm -> 0.42mm にすると誤差が 13% -> 1%)
	double hmax = 0;
	for (int e = 0; e < NTet; e++) {
		if (Material[TetMat[e]].sigma <= 0) continue;
		const int32_t *nd = &Tet[e * 4];
		for (int k = 0; k < 6; k++) {
			const int32_t a = nd[EDGE_NODE3[k][0]], b = nd[EDGE_NODE3[k][1]];
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
	assemble_eddy3d(&Kr, &Ki, omega);
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
	for (int e = 0; e < NTet; e++) {
		if (Material[TetMat[e]].sigma <= 0) continue;
		for (int l = 0; l < 4; l++) cond[Tet[(e * 4) + l]] = 1;
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
			Solver.maxiter, Solver.nout, Solver.converg, fp_log, label);
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
		fprintf(fp_log, "  port %d : I = %.6e %+.6e j [A] (V = %.6e [V])\n",
			jc, qr[jc], qi[jc], Volt);
	}

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

	return ierr;
}
