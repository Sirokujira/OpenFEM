/*
eddy3d_p2.c

3 次元渦電流 (analysis = A、A-φ 連成) の **2 次** の要素。10 節点四面体 (等パラメトリック、曲がって
よい) の格子で、A を 2 次の Nedelec 要素 (第 1 種、20 自由度)、φ を 2 次の節点要素で離散化する。

基底は階層的 (Webb の 1.5 次):

	辺 (a, b)   : W_ab = λa∇λb − λb∇λa (Whitney、全体の向きの符号つき)、∇(λaλb) (符号なし)
	面 (i, j, k) : λk W_ij、λj W_ik    (i < j < k は面の頂点の全体番号の昇順。隣の要素と同じ関数になる)
	φ           : 頂点 λv、辺の泡 λaλb

λ は参照四面体の重心座標を物理空間に引き戻したもの (λ_i = Σ_a λ_i(ξ_a) N_a、2 次の節点基底で厳密)。
**共変 Piola 変換は「引き戻した λ の物理空間の勾配を使う」ことと同じ**なので、曲がった要素でも
∇×W_ab = 2∇λa×∇λb、∇×(λc W_ab) = ∇λc×W_ab + 2λc ∇λa×∇λb が厳密に成り立つ。面の関数の接線成分は
自分の面の外で 0 になる (λk = 0 の面では消え、λj = 0 の面では W_ij が法線を向く) ので、接線連続性は
辺と面の自由度だけで決まる。

**∇φ は A の空間に厳密に入る** (∇λv = Σ_j W_jv、∇(λaλb) は辺の関数そのもの)。そのため結合の行列を
同じ積分点で直接積分しても (A, φ) -> (A + ∇ψ, φ − jωψ) の不変性が厳密に保たれ、右辺は値域に入る
(1 次と同じくゲージは固定しない)。ω -> 0 で節点の式は ∫σ∇N・∇N (2 次の導電行列) そのものになる
(DC の検査)。積分は Duffy + 各方向 4 点 (elem3d_ip、64 点)。

未知数の並び : [Whitney (NEdge) | ∇(λλ) (NEdge) | 面 2 個ずつ (2 NFace) | φ 頂点 (NNode、中間節点の欄は
固定) | φ 泡 (NEdge)]。

対応するのは周波数領域 (線形・B-H 曲線の実効磁気抵抗率法)・過渡解析 (線形)・交流の発熱・場の出力で、
ゲージなし・反復解法だけ (それ以外は入力段で弾く)。
*/

#include "fem.h"
#include "fem_prototype.h"

static const int P2_EDGE[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
static const int P2_FACE[4][3] = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};
// gmsh の tet10 の中間節点 (局所 4..9) の辺 : (0,1) (1,2) (2,0) (3,0) (3,2) (3,1)
#define P2_AAM (5)		// 非線形の Anderson 加速の深さ (1 次の AAM3 と同じ)

static const int P2_MID[6][2] = {{0, 1}, {1, 2}, {2, 0}, {3, 0}, {3, 2}, {3, 1}};

#define P2_NDOF (30)		// 20 (A) + 10 (φ)

typedef struct {
	int32_t v[3];			// 昇順の頂点
	int32_t e, f;			// 要素と局所の面
} p2face_t;

static int cmp_p2face(const void *a, const void *b)
{
	const p2face_t *x = (const p2face_t *)a, *y = (const p2face_t *)b;
	for (int c = 0; c < 3; c++) {
		if (x->v[c] < y->v[c]) return -1;
		if (x->v[c] > y->v[c]) return 1;
	}
	return 0;
}

static void sort3(int32_t v[3])
{
	for (int a = 0; a < 2; a++) {
		for (int b = a + 1; b < 3; b++) {
			if (v[b] < v[a]) { const int32_t t = v[a]; v[a] = v[b]; v[b] = t; }
		}
	}
}

static int cmp_i32p(const void *a, const void *b)
{
	const int32_t x = *(const int32_t *)a, y = *(const int32_t *)b;
	return ((x < y) ? -1 : ((x > y) ? 1 : 0));
}

static void cross3(const double a[3], const double b[3], double c[3])
{
	c[0] = (a[1] * b[2]) - (a[2] * b[1]);
	c[1] = (a[2] * b[0]) - (a[0] * b[2]);
	c[2] = (a[0] * b[1]) - (a[1] * b[0]);
}

static double dot3(const double a[3], const double b[3])
{
	return (a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2]);
}

// ν (対称テンソル xx yy zz xy yz zx) を掛ける
static void nu_mul(const double nu[6], const double c[3], double r[3])
{
	r[0] = (nu[0] * c[0]) + (nu[3] * c[1]) + (nu[5] * c[2]);
	r[1] = (nu[3] * c[0]) + (nu[1] * c[1]) + (nu[4] * c[2]);
	r[2] = (nu[5] * c[0]) + (nu[4] * c[1]) + (nu[2] * c[2]);
}


// 全体の面と、要素ごとの面番号 fe[4 e + f] (3 次以上の eddy3d_pr.c も使う)
int32_t *p2_faces(int *nface, int32_t **fvert)
{
	const int nt = NTet;
	p2face_t *a = (p2face_t *)malloc((size_t)nt * 4 * sizeof(p2face_t));
	for (int e = 0; e < nt; e++) {
		for (int f = 0; f < 4; f++) {
			p2face_t *q = &a[(e * 4) + f];
			for (int c = 0; c < 3; c++) q->v[c] = Tet[(e * 4) + P2_FACE[f][c]];
			sort3(q->v);
			q->e = e;
			q->f = f;
		}
	}
	qsort(a, (size_t)nt * 4, sizeof(p2face_t), cmp_p2face);
	int32_t *fe = (int32_t *)malloc((size_t)nt * 4 * sizeof(int32_t));
	int32_t *fv = (int32_t *)malloc((size_t)nt * 4 * 3 * sizeof(int32_t));
	int nf = -1;
	for (int i = 0; i < (nt * 4); i++) {
		if ((i == 0) || cmp_p2face(&a[i], &a[i - 1])) {
			nf++;
			for (int c = 0; c < 3; c++) fv[(nf * 3) + c] = a[i].v[c];
		}
		fe[(a[i].e * 4) + a[i].f] = nf;
	}
	free(a);
	*nface = nf + 1;
	*fvert = fv;
	return fe;
}


// 面 (昇順の 3 頂点) の番号。無ければ -1
int p2_face_id(const int32_t *fvert, int nface, const int32_t v[3])
{
	int lo = 0, hi = nface - 1;
	while (lo <= hi) {
		const int mid = (lo + hi) / 2;
		const int32_t *w = &fvert[mid * 3];
		int c = 0;
		for (int k = 0; k < 3; k++) {
			if (v[k] != w[k]) { c = ((v[k] < w[k]) ? -1 : 1); break; }
		}
		if (c == 0) return mid;
		if (c < 0) hi = mid - 1;
		else lo = mid + 1;
	}
	return -1;
}


// 要素の 30 自由度の全体番号と符号
static void p2_dofs(int e, const int32_t *fe, int ne, int nf, int nn, int32_t dof[P2_NDOF], double sgn[P2_NDOF])
{
	const int32_t *ed = &TetEdge[EdgeOff[e]];
	const signed char *sg = &TetEdgeSgn[EdgeOff[e]];
	const int32_t *nd = &Tet[e * 4];
	const int oG = ne, oF = 2 * ne, oV = (2 * ne) + (2 * nf), oB = oV + nn;
	for (int k = 0; k < 6; k++) {
		dof[k] = ed[k];
		sgn[k] = sg[k];
		dof[6 + k] = oG + ed[k];
		sgn[6 + k] = 1;
	}
	for (int f = 0; f < 4; f++) {
		dof[12 + (2 * f)] = oF + (2 * fe[(e * 4) + f]);
		dof[13 + (2 * f)] = oF + (2 * fe[(e * 4) + f]) + 1;
		sgn[12 + (2 * f)] = sgn[13 + (2 * f)] = 1;
	}
	for (int a = 0; a < 4; a++) {
		dof[20 + a] = oV + nd[a];
		sgn[20 + a] = 1;
	}
	for (int k = 0; k < 6; k++) {
		dof[24 + k] = oB + ed[k];
		sgn[24 + k] = 1;
	}
}


// 10 節点の形状関数の値 nv と物理空間の勾配 g から、20 個の A の基底 (値と回転、局所の向き) と
// 10 個の φ の基底 (値と勾配) を作る
static void p2_basis_core(int e, const double nv[ELEM_NMAX], const double g[ELEM_NMAX][3],
	double wv[20][3], double wc[20][3], double pn[10], double pg[10][3]);

// 積分点 q で
static int p2_basis(int e, int q, double wv[20][3], double wc[20][3], double pn[10], double pg[10][3], double *wq)
{
	double nv[ELEM_NMAX], g[ELEM_NMAX][3];
	if (elem3d_ip(e, q, nv, g, wq)) return 1;
	p2_basis_core(e, nv, g, wv, wc, pn, pg);
	return 0;
}


// 任意の参照点 (重心座標) で (場の出力)
static int p2_basis_at(int e, const double lam[4], double wv[20][3], double wc[20][3], double pn[10], double pg[10][3])
{
	double nv[ELEM_NMAX], g[ELEM_NMAX][3], det;
	if (tet_ho_point(e, lam, nv, g, &det)) return 1;
	p2_basis_core(e, nv, g, wv, wc, pn, pg);
	return 0;
}


static void p2_basis_core(int e, const double nv[ELEM_NMAX], const double g[ELEM_NMAX][3],
	double wv[20][3], double wc[20][3], double pn[10], double pg[10][3])
{
	// 引き戻した重心座標とその物理空間の勾配 (λ_i = Σ_a λ_i(ξ_a) N_a)
	double lam[4] = {0, 0, 0, 0}, gl[4][3];
	for (int i = 0; i < 4; i++) {
		lam[i] = nv[i];
		for (int c = 0; c < 3; c++) gl[i][c] = g[i][c];
	}
	for (int m = 0; m < 6; m++) {
		const int a = P2_MID[m][0], b = P2_MID[m][1];
		lam[a] += 0.5 * nv[4 + m];
		lam[b] += 0.5 * nv[4 + m];
		for (int c = 0; c < 3; c++) {
			gl[a][c] += 0.5 * g[4 + m][c];
			gl[b][c] += 0.5 * g[4 + m][c];
		}
	}
	const int32_t *nd = &Tet[e * 4];
	for (int k = 0; k < 6; k++) {
		const int a = P2_EDGE[k][0], b = P2_EDGE[k][1];
		double cr[3];
		cross3(gl[a], gl[b], cr);
		for (int c = 0; c < 3; c++) {
			wv[k][c] = (lam[a] * gl[b][c]) - (lam[b] * gl[a][c]);
			wc[k][c] = 2 * cr[c];
			wv[6 + k][c] = (lam[a] * gl[b][c]) + (lam[b] * gl[a][c]);
			wc[6 + k][c] = 0;
		}
	}
	for (int f = 0; f < 4; f++) {
		// 面の頂点を全体番号の昇順 (i, j, k) に
		int lv[3] = {P2_FACE[f][0], P2_FACE[f][1], P2_FACE[f][2]};
		for (int x = 0; x < 2; x++) {
			for (int y = x + 1; y < 3; y++) {
				if (nd[lv[y]] < nd[lv[x]]) { const int t = lv[x]; lv[x] = lv[y]; lv[y] = t; }
			}
		}
		const int i = lv[0], j = lv[1], k = lv[2];
		// F1 = λk W_ij、F2 = λj W_ik
		const int pa[2][3] = {{k, i, j}, {j, i, k}};
		for (int s = 0; s < 2; s++) {
			const int cc = pa[s][0], aa = pa[s][1], bb = pa[s][2];
			double w[3], cr[3], gw[3];
			for (int c = 0; c < 3; c++) w[c] = (lam[aa] * gl[bb][c]) - (lam[bb] * gl[aa][c]);
			cross3(gl[aa], gl[bb], cr);
			cross3(gl[cc], w, gw);
			for (int c = 0; c < 3; c++) {
				wv[12 + (2 * f) + s][c] = lam[cc] * w[c];
				wc[12 + (2 * f) + s][c] = gw[c] + (2 * lam[cc] * cr[c]);
			}
		}
	}
	for (int a = 0; a < 4; a++) {
		pn[a] = lam[a];
		for (int c = 0; c < 3; c++) pg[a][c] = gl[a][c];
	}
	for (int k = 0; k < 6; k++) {
		const int a = P2_EDGE[k][0], b = P2_EDGE[k][1];
		pn[4 + k] = lam[a] * lam[b];
		for (int c = 0; c < 3; c++) pg[4 + k][c] = (lam[a] * gl[b][c]) + (lam[b] * gl[a][c]);
	}
}


static int64_t p2_find(const crs_t *A, int32_t row, int32_t col)
{
	int64_t lo = A->rowptr[row], hi = A->rowptr[row + 1] - 1;
	while (lo <= hi) {
		const int64_t mid = (lo + hi) / 2;
		if (A->col[mid] == col) return mid;
		if (A->col[mid] < col) lo = mid + 1;
		else hi = mid - 1;
	}
	return -1;
}


// 要素の自由度の組から CRS のパターンを作る
static void p2_crs(crs_t *A, int n, const int32_t *fe, int ne, int nf, int nn)
{
	const int nt = NTet;
	int *cnt = (int *)calloc((size_t)n, sizeof(int));
	int32_t dof[P2_NDOF];
	double sgn[P2_NDOF];
	for (int e = 0; e < nt; e++) {
		p2_dofs(e, fe, ne, nf, nn, dof, sgn);
		for (int a = 0; a < P2_NDOF; a++) cnt[dof[a]] += P2_NDOF;
	}
	// 未使用の行 (中間節点の φ) にも対角を置く
	for (int i = 0; i < n; i++) cnt[i] += 1;
	int64_t *ptr = (int64_t *)malloc(((size_t)n + 1) * sizeof(int64_t));
	ptr[0] = 0;
	for (int i = 0; i < n; i++) ptr[i + 1] = ptr[i] + cnt[i];
	int32_t *lst = (int32_t *)malloc((size_t)ptr[n] * sizeof(int32_t));
	memset(cnt, 0, (size_t)n * sizeof(int));
	for (int i = 0; i < n; i++) lst[ptr[i] + (cnt[i]++)] = i;
	for (int e = 0; e < nt; e++) {
		p2_dofs(e, fe, ne, nf, nn, dof, sgn);
		for (int a = 0; a < P2_NDOF; a++) {
			for (int b = 0; b < P2_NDOF; b++) lst[ptr[dof[a]] + (cnt[dof[a]]++)] = dof[b];
		}
	}
	A->n = n;
	A->rowptr = (int64_t *)malloc(((size_t)n + 1) * sizeof(int64_t));
	A->rowptr[0] = 0;
	int64_t nnz = 0;
	for (int i = 0; i < n; i++) {
		int32_t *r = &lst[ptr[i]];
		const int m = cnt[i];
		qsort(r, (size_t)m, sizeof(int32_t), cmp_i32p);
		int u = 0;
		for (int q = 0; q < m; q++) {
			if ((q == 0) || (r[q] != r[q - 1])) u++;
		}
		nnz += u;
		A->rowptr[i + 1] = nnz;
	}
	A->nnz = nnz;
	A->col = (int32_t *)malloc((size_t)nnz * sizeof(int32_t));
	A->val = (double *)calloc((size_t)nnz, sizeof(double));
	for (int i = 0; i < n; i++) {
		const int32_t *r = &lst[ptr[i]];
		int64_t p = A->rowptr[i];
		for (int q = 0; q < cnt[i]; q++) {
			if ((q == 0) || (r[q] != r[q - 1])) A->col[p++] = r[q];
		}
	}
	free(cnt); free(ptr); free(lst);
}


static void p2_row_dot_c(const crs_t *Kr, const crs_t *Ki, int32_t row,
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


// 解く系の一式 (交流の発熱は反復の間これを持ち続ける)
static struct {
	int ne, nf, nn, n, oG, oF, oV, oB;
	int32_t *fe, *fvert;
	crs_t Kr, Ki;
	unsigned char *fix, *cond, *econd, *vert;
	double *xdr, *xdi, *br, *bi, *ur, *ui, *w1, *w2;
	double omega;
	const double *nuq;		// NULL でなければ B-H 曲線の材料の積分点ごとの ν (要素 e の点 q は nuq[e * nq + q])
	int nq;
} P2;


static void p2_close(void);


static int p2_open(FILE *fp_log)
{
	const int np = NPort;
	if (!MeshMode || (TetOrder != 2) || (NHex + NPrism + NPyr > 0) || (np < 1) || (Freq <= 0)) {
		fprintf(fp_log, "*** the second-order A-phi solver needs a pure 10-node tetrahedral mesh, "
			"a port and a frequency\n");
		return 1;
	}
	edge_build();
	const int ne = NEdge, nn = NNode, nt = NTet;
	int nf = 0;
	P2.fe = p2_faces(&nf, &P2.fvert);
	P2.ne = ne; P2.nf = nf; P2.nn = nn;
	P2.oG = ne; P2.oF = 2 * ne; P2.oV = (2 * ne) + (2 * nf); P2.oB = P2.oV + nn;
	const int n = P2.oB + ne;
	P2.n = n;
	P2.omega = 2 * PI * Freq;
	fprintf(fp_log, "  nodes = %d, tetrahedra = %d, edges = %d, faces = %d, unknowns = %d "
		"(A : 2 x %d edges + 2 x %d faces, phi : vertices + %d edge bubbles)\n", nn, nt, ne, nf, n, ne, nf, ne);
	p2_crs(&P2.Kr, n, P2.fe, ne, nf, nn);
	p2_crs(&P2.Ki, n, P2.fe, ne, nf, nn);
	fprintf(fp_log, "  coupled matrix : %lld nonzeros (%.1f per row)\n", (long long)P2.Kr.nnz, (double)P2.Kr.nnz / n);

	// Dirichlet : awall の三角形の辺 (Whitney と勾配) と面 (2 個)、電極の頂点 (φ = V) と電極上の辺の泡 (0)、
	//             非導電領域の φ、使わない欄 (中間節点の φ)
	const int oG = P2.oG, oF = P2.oF, oV = P2.oV, oB = P2.oB;
	unsigned char *fix = P2.fix = (unsigned char *)calloc((size_t)n, sizeof(unsigned char));
	int nwall = 0;
	for (int t = 0; t < NTri; t++) {
		int hit = 0;
		for (int q = 0; q < NAWall; q++) {
			if (TriTag[t] == AWallTag[q]) hit = 1;
		}
		if (!hit) continue;
		const int32_t *v = &Tri[t * 3];
		for (int l = 0; l < 3; l++) {
			const int64_t id = edge_id(v[l], v[(l + 1) % 3]);
			if (id >= 0) {
				fix[id] = 1;
				fix[oG + id] = 1;
			}
		}
		int32_t sv[3] = {v[0], v[1], v[2]};
		sort3(sv);
		const int f = p2_face_id(P2.fvert, nf, sv);
		if (f >= 0) {
			fix[oF + (2 * f)] = fix[oF + (2 * f) + 1] = 1;
			nwall++;
		}
	}
	unsigned char *cond = P2.cond = (unsigned char *)calloc((size_t)nn, sizeof(unsigned char));
	unsigned char *econd = P2.econd = (unsigned char *)calloc((size_t)ne, sizeof(unsigned char));
	unsigned char *vert = P2.vert = (unsigned char *)calloc((size_t)nn, sizeof(unsigned char));
	for (int e = 0; e < nt; e++) {
		const int32_t *nd = &Tet[e * 4];
		for (int a = 0; a < 4; a++) vert[nd[a]] = 1;
		if (Material[elem3d_mat(e)].sigma <= 0) continue;
		for (int a = 0; a < 4; a++) cond[nd[a]] = 1;
		for (int k = 0; k < 6; k++) econd[TetEdge[EdgeOff[e] + k]] = 1;
	}
	int nelec = 0;
	for (int i = 0; i < nn; i++) {
		if (!vert[i] || !cond[i]) fix[oV + i] = 1;
		else if (NodeConductor[i] >= 0) {
			fix[oV + i] = 1;
			nelec++;
		}
	}
	for (int e = 0; e < nt; e++) {
		const int32_t *nd = &Tet[e * 4];
		for (int k = 0; k < 6; k++) {
			const int32_t id = TetEdge[EdgeOff[e] + k];
			const int32_t a = nd[P2_EDGE[k][0]], b = nd[P2_EDGE[k][1]];
			if (!econd[id]) fix[oB + id] = 1;
			else if ((NodeConductor[a] >= 0) && (NodeConductor[a] == NodeConductor[b])) fix[oB + id] = 1;
		}
	}
	fprintf(fp_log, "  fixed : %d wall faces (A_t = 0), %d electrode vertices\n", nwall, nelec);
	// ゲージ固定 (gauge = 1) : 2 次の節点空間の勾配は ∇λv (Whitney の和) と ∇(λaλb) (辺の勾配の自由度
	// そのもの) に分かれるので、**Whitney は 1 次と同じ awall 優先の全域木で、∇(λaλb) は全部 0 に**
	// 固定すると零空間が消える (任意の A は ψ の泡の係数で勾配の自由度を、頂点の値で木の辺を 0 にできる)
	if (GaugeTree) {
		unsigned char *wall = (unsigned char *)malloc((size_t)ne * sizeof(unsigned char));
		unsigned char *tree = (unsigned char *)malloc((size_t)ne * sizeof(unsigned char));
		int *comp = (int *)malloc((size_t)nn * sizeof(int));
		for (int i = 0; i < ne; i++) wall[i] = fix[i];
		const int ntree = edge_tree_gauge(wall, tree, comp);
		// **電極の上の辺 (φ の泡を 0 に固定した導体の辺) の ∇(λaλb) は固定しない**: そこでは ψ の泡も 0 に
		// 縛られていてゲージの自由度が無く、その係数は本物の場 (電極の面の接線の A) になる。固定すると
		// 電極の辺が awall でない格子で Z が黙って変わる (実測 R −5.2%、L +17%)。電極の面を awall にも
		// した格子 (bar_air など) では該当する辺が 0 本なので、そこでは検出できない
		int ngauge = 0, ngrad = 0, nelecg = 0;
		for (int i = 0; i < ne; i++) {
			if (tree[i] && !fix[i]) ngauge++;
			fix[i] |= tree[i];
			if (fix[oG + i]) continue;
			if (econd[i] && fix[oB + i]) {
				nelecg++;
				continue;
			}
			ngrad++;
			fix[oG + i] = 1;
		}
		// 木から決まる ψ が電極上で 0 になる条件 : 電極 (の頂点) が awall の辺だけで連結
		int bad = 0, ref = -1;
		for (int i = 0; i < nn; i++) {
			if ((NodeConductor[i] < 0) || !vert[i]) continue;
			if (ref < 0) ref = comp[i];
			else if (comp[i] != ref) bad = 1;
		}
		fprintf(fp_log, "  gauge : %d tree edges (%d beyond the wall), %d edge gradients "
			"(%d electrode edges left free)\n", ntree, ngauge, ngrad, nelecg);
		free(wall); free(tree); free(comp);
		if (bad) {
			fprintf(fp_log, "*** the electrodes are not connected through awall edges, "
				"so the tree gauge would shift the terminal voltage; "
				"use gauge = 0 (the default), or extend awall so that all electrodes "
				"share one A_t = 0 surface\n");
			p2_close();
			return 1;
		}
	}
	if (nwall == 0) {
		fprintf(fp_log, "*** warning : no awall surface; the outer boundary is a magnetic wall (dA/dn = 0)\n");
	}
	P2.xdr = (double *)calloc((size_t)n, sizeof(double));
	P2.xdi = (double *)calloc((size_t)n, sizeof(double));
	P2.br = (double *)malloc((size_t)n * sizeof(double));
	P2.bi = (double *)malloc((size_t)n * sizeof(double));
	P2.ur = (double *)malloc((size_t)n * sizeof(double));
	P2.ui = (double *)malloc((size_t)n * sizeof(double));
	P2.w1 = (double *)malloc((size_t)n * sizeof(double));
	P2.w2 = (double *)malloc((size_t)n * sizeof(double));
	return 0;
}


static void p2_close(void)
{
	free(P2.xdr); free(P2.xdi); free(P2.br); free(P2.bi); free(P2.ur); free(P2.ui); free(P2.w1); free(P2.w2);
	free(P2.fix); free(P2.cond); free(P2.econd); free(P2.vert);
	free(P2.fe); free(P2.fvert);
	crs_free(&P2.Kr);
	crs_free(&P2.Ki);
	edge_free();
	memset(&P2, 0, sizeof(P2));
}


// 組み立て : [S + jωT, TG; GᵀT, GᵀTG/(jω)] (節点の式を jω で割った複素対称形、1 次と同じ)。
// sigelem が NULL でなければ要素の σ (熱連成の σ(T))
static int p2_assemble(FILE *fp_log, const double *sigelem)
{
	const int ne = P2.ne, nf = P2.nf, nn = P2.nn, nt = NTet;
	const double omega = P2.omega;
	crs_zero(&P2.Kr);
	crs_zero(&P2.Ki);
	int bad = 0;
	for (int e = 0; e < nt; e++) {
		const int m = elem3d_mat(e);
		double nu[6];
		material_coef_pub(m, 3, nu);
		const double sig = ((sigelem != NULL) ? sigelem[e] : Material[m].sigma);
		int32_t dof[P2_NDOF];
		double sgn[P2_NDOF];
		p2_dofs(e, P2.fe, ne, nf, nn, dof, sgn);
		double se[20][20], te[20][20], ce[20][10], pe[10][10];
		memset(se, 0, sizeof(se)); memset(te, 0, sizeof(te));
		memset(ce, 0, sizeof(ce)); memset(pe, 0, sizeof(pe));
		const int nq = elem3d_nip(e);
		const int nlm = ((P2.nuq != NULL) && (Material[m].nbh[0] > 0));
		for (int q = 0; q < nq; q++) {
			double wv[20][3], wc[20][3], pn[10], pg[10][3], w;
			if (p2_basis(e, q, wv, wc, pn, pg, &w)) { bad = 1; break; }
			if (nlm) {
				nu[0] = nu[1] = nu[2] = P2.nuq[((int64_t)e * P2.nq) + q];
				nu[3] = nu[4] = nu[5] = 0;
			}
			double nc[20][3];
			for (int k = 0; k < 20; k++) nu_mul(nu, wc[k], nc[k]);
			for (int k = 0; k < 20; k++) {
				for (int l = 0; l < 20; l++) {
					se[k][l] += w * dot3(nc[k], wc[l]);
					te[k][l] += w * sig * dot3(wv[k], wv[l]);
				}
				for (int a = 0; a < 10; a++) ce[k][a] += w * sig * dot3(wv[k], pg[a]);
			}
			for (int a = 0; a < 10; a++) {
				for (int b = 0; b < 10; b++) pe[a][b] += w * sig * dot3(pg[a], pg[b]);
			}
		}
		for (int k = 0; k < 20; k++) {
			for (int l = 0; l < 20; l++) {
				const int64_t p = p2_find(&P2.Kr, dof[k], dof[l]);
				const double s = sgn[k] * sgn[l];
				P2.Kr.val[p] += s * se[k][l];
				P2.Ki.val[p] += s * omega * te[k][l];
			}
			for (int a = 0; a < 10; a++) {
				const double v = sgn[k] * ce[k][a];
				P2.Kr.val[p2_find(&P2.Kr, dof[k], dof[20 + a])] += v;
				P2.Kr.val[p2_find(&P2.Kr, dof[20 + a], dof[k])] += v;
			}
		}
		for (int a = 0; a < 10; a++) {
			for (int b = 0; b < 10; b++) P2.Ki.val[p2_find(&P2.Ki, dof[20 + a], dof[20 + b])] -= pe[a][b] / omega;
		}
	}
	if (bad) {
		fprintf(fp_log, "*** the second-order A-phi solver : a degenerate element\n");
		return 1;
	}
	return 0;
}


// ポート jc を V = Volt で励振して解き (P2.ur / ui に解)、各導体へ流れ込む電流 qr + j qi を返す
static int p2_solve_port(FILE *fp_log, int jc, FILE *fp_iter, const char *label, double qr[MAXPORT], double qi[MAXPORT])
{
	const int n = P2.n, nn = P2.nn, oV = P2.oV;
	const double omega = P2.omega;
	for (int i = 0; i < n; i++) P2.xdr[i] = P2.xdi[i] = 0;
	for (int i = 0; i < nn; i++) {
		if (P2.fix[oV + i] && P2.cond[i] && P2.vert[i] && (NodeConductor[i] == jc)) P2.xdr[oV + i] = Volt;
	}
	crs_spmv(&P2.Kr, P2.xdr, P2.w1, NULL);
	crs_spmv(&P2.Ki, P2.xdr, P2.w2, NULL);
	for (int i = 0; i < n; i++) {
		P2.br[i] = (P2.fix[i] ? 0 : -P2.w1[i]);
		P2.bi[i] = (P2.fix[i] ? 0 : -P2.w2[i]);
	}
	const int iter = solver_cocg(&P2.Kr, &P2.Ki, 1.0, P2.br, P2.bi, P2.ur, P2.ui, P2.fix,
		Solver.maxiter, Solver.nout, Solver.converg, fp_iter, label);
	if (iter < 0) {
		fprintf(fp_log, "*** the second-order A-phi solve did not converge (%s)\n", label);
		return 1;
	}
	for (int i = 0; i < n; i++) {
		P2.ur[i] += P2.xdr[i];
		P2.ui[i] += P2.xdi[i];
	}
	// 電極の頂点の反作用 (節点の式は 1/(jω) 倍してあるので jω 倍して戻す)。電極の上で 1 になる
	// 試験関数 Σ λv は泡の係数 0 の階層的な 2 次の関数なので、頂点の行だけを足せばよい
	for (int p = 0; p < MAXPORT; p++) qr[p] = qi[p] = 0;
	for (int i = 0; i < nn; i++) {
		const int id = NodeConductor[i];
		if ((id < 0) || !P2.cond[i] || !P2.vert[i]) continue;
		double rr, ri;
		p2_row_dot_c(&P2.Kr, &P2.Ki, (int32_t)(oV + i), P2.ur, P2.ui, &rr, &ri);
		qr[id] += -omega * ri;
		qi[id] += omega * rr;
	}
	return 0;
}


// 場の出力 (ポート 1 の励振) : 要素の重心での B = ∇×A (実部) と J = −σ (jωA + ∇φ) (実部・虚部)。
// どちらもゲージ不変なのでゲージを固定しない解から作ってよい (1 次と同じ名前・同じ符号)
static void p2_fields(void)
{
	const int ne = P2.ne, nf = P2.nf, nn = P2.nn, nt = NTet;
	const double omega = P2.omega;
	const double lc[4] = {0.25, 0.25, 0.25, 0.25};
	double *bv = (double *)calloc((size_t)nt * 3, sizeof(double));
	double *jr = (double *)calloc((size_t)nt * 3, sizeof(double));
	double *ji = (double *)calloc((size_t)nt * 3, sizeof(double));
	for (int e = 0; e < nt; e++) {
		int32_t dof[P2_NDOF];
		double sgn[P2_NDOF];
		p2_dofs(e, P2.fe, ne, nf, nn, dof, sgn);
		double wv[20][3], wc[20][3], pn[10], pg[10][3];
		if (p2_basis_at(e, lc, wv, wc, pn, pg)) continue;
		double br[3] = {0, 0, 0}, er[3] = {0, 0, 0}, ei[3] = {0, 0, 0};
		for (int k = 0; k < 20; k++) {
			const double ar = sgn[k] * P2.ur[dof[k]], ai = sgn[k] * P2.ui[dof[k]];
			for (int c = 0; c < 3; c++) {
				br[c] += ar * wc[k][c];
				er[c] += -omega * ai * wv[k][c];
				ei[c] += omega * ar * wv[k][c];
			}
		}
		for (int a = 0; a < 10; a++) {
			for (int c = 0; c < 3; c++) {
				er[c] += P2.ur[dof[20 + a]] * pg[a][c];
				ei[c] += P2.ui[dof[20 + a]] * pg[a][c];
			}
		}
		const double sg = Material[elem3d_mat(e)].sigma;
		for (int c = 0; c < 3; c++) {
			bv[(e * 3) + c] = br[c];
			jr[(e * 3) + c] = -sg * er[c];
			ji[(e * 3) + c] = -sg * ei[c];
		}
	}
	field_add_cellvec("B_A_re", bv);
	field_add_cellvec("J_A_re", jr);
	field_add_cellvec("J_A_im", ji);
	free(bv); free(jr); free(ji);
}


int solve_eddy3d_p2(FILE *fp_log)
{
	const int np = NPort, nt = NTet;
	int ierr = 0;
	fprintf(fp_log, "\n=== 3D eddy current (A-phi, second-order Nedelec / nodal elements) ===\n");
	if (p2_open(fp_log)) return 1;
	const double omega = P2.omega;

	// 非線形 (B-H) の材料があれば実効磁気抵抗率法 (1 次の eddy3d.c と同じ): ν を波高値 |B| で
	// **積分点ごとに**評価し (2 次の要素では B が要素内で変わる)、H を独立変数にした更新
	// ν = |H|/B(|H|) に ln ν の Anderson 加速を掛ける。単一ポート・電流振幅 current
	int nlf = 0;
	for (int m = 0; m < NMaterial; m++) {
		if (Material[m].nbh[0] > 0) nlf = 1;
	}
	const int nq = elem3d_nip(0);
	const int64_t nqt = (int64_t)nt * nq;
	double *nuq = NULL, *aay = NULL, *aaf = NULL, *aayp = NULL, *aafp = NULL, *aady = NULL, *aadf = NULL;
	int aan = 0;
	if (nlf) {
		nuq  = (double *)malloc((size_t)nqt * sizeof(double));
		aay  = (double *)malloc((size_t)nqt * sizeof(double));
		aaf  = (double *)malloc((size_t)nqt * sizeof(double));
		aayp = (double *)malloc((size_t)nqt * sizeof(double));
		aafp = (double *)malloc((size_t)nqt * sizeof(double));
		aady = (double *)malloc((size_t)nqt * P2_AAM * sizeof(double));
		aadf = (double *)malloc((size_t)nqt * P2_AAM * sizeof(double));
		for (int e = 0; e < nt; e++) {
			const material_t *mt = &Material[elem3d_mat(e)];
			for (int q = 0; q < nq; q++) nuq[((int64_t)e * nq) + q] = ((mt->nbh[0] > 0) ? bh_nu(mt, 0, 0) : 0);
		}
		P2.nuq = nuq;
		P2.nq = nq;
		fprintf(fp_log, "  nonlinear (B-H) : effective reluctivity at %d points per element, current amplitude "
			"%.6e [A], maxiter = %d, tol = %.1e, relax = %.2f\n", nq, Curr, NlMaxiter, NlTol, NlRelax);
	}
	ierr = p2_assemble(fp_log, NULL);
	double *yr = (double *)calloc((size_t)np * np, sizeof(double));
	double *yi = (double *)calloc((size_t)np * np, sizeof(double));
	int final = !nlf;
	for (int nlit = 1; !ierr; nlit++) {
	if (nlit > 1) ierr = p2_assemble(fp_log, NULL);
	if (final) fprintf(fp_log, "  %-10s %8s %13s\n", "excite", "iter", "residual");
	for (int jc = 1; (jc <= np) && !ierr; jc++) {
		char label[64];
		sprintf(label, "port%d", jc);
		double qr[MAXPORT], qi[MAXPORT];
		if (p2_solve_port(fp_log, jc, (final ? fp_log : NULL), label, qr, qi)) {
			ierr = 1;
			break;
		}
		for (int kc = 1; kc <= np; kc++) {
			yr[((kc - 1) * np) + (jc - 1)] = qr[kc] / Volt;
			yi[((kc - 1) * np) + (jc - 1)] = qi[kc] / Volt;
		}
		if (final) fprintf(fp_log, "  port %d : I = %.6e %+.6e j [A] (V = %.6e [V])\n", jc, qr[jc], qi[jc], Volt);
		// 非線形 : 端子電流が current になるよう解を s = Curr / I 倍する (Y は倍率に依らない)
		if (nlf) {
			const double ivr = qr[1], ivi = qi[1];
			const double den = (ivr * ivr) + (ivi * ivi);
			const double svr = Curr * ivr / den, svi = -Curr * ivi / den;
			for (int i = 0; i < P2.n; i++) {
				const double xr = P2.ur[i], xi = P2.ui[i];
				P2.ur[i] = (svr * xr) - (svi * xi);
				P2.ui[i] = (svr * xi) + (svi * xr);
			}
		}
		if (FieldOut && (jc == 1) && final) p2_fields();
	}
	if (final || ierr) break;

	// ν の更新 (H を独立変数にする) と Anderson 加速 (ln ν)。B は組み立てと同じ積分点で
	double dmax = 0;
	for (int e = 0; e < nt; e++) {
		const material_t *mt = &Material[elem3d_mat(e)];
		const int64_t o = (int64_t)e * nq;
		if (mt->nbh[0] <= 0) {
			for (int q = 0; q < nq; q++) aay[o + q] = aaf[o + q] = 0;
			continue;
		}
		int32_t dof[P2_NDOF];
		double sgn[P2_NDOF];
		p2_dofs(e, P2.fe, P2.ne, P2.nf, P2.nn, dof, sgn);
		for (int q = 0; q < nq; q++) {
			double wv[20][3], wc[20][3], pn[10], pg[10][3], w;
			double bre[3] = {0, 0, 0}, bim[3] = {0, 0, 0};
			if (!p2_basis(e, q, wv, wc, pn, pg, &w)) {
				for (int k = 0; k < 20; k++) {
					const double ar = sgn[k] * P2.ur[dof[k]], ai = sgn[k] * P2.ui[dof[k]];
					for (int c = 0; c < 3; c++) {
						bre[c] += ar * wc[k][c];
						bim[c] += ai * wc[k][c];
					}
				}
			}
			const double b2 = dot3(bre, bre) + dot3(bim, bim);
			const double hm = nuq[o + q] * sqrt(b2);
			const double nun = ((hm > 0) ? (hm / bh_b_of_h(mt, 0, hm)) : bh_nu(mt, 0, 0));
			const double dn = fabs(nun - nuq[o + q]) / nun;
			if (dn > dmax) dmax = dn;
			aay[o + q] = log(nuq[o + q]);
			aaf[o + q] = log(nun) - aay[o + q];
		}
	}
	anderson_step(nqt, aay, aaf, aayp, aafp, aady, aadf, P2_AAM, &aan, NlRelax);
	for (int e = 0; e < nt; e++) {
		if (Material[elem3d_mat(e)].nbh[0] <= 0) continue;
		for (int q = 0; q < nq; q++) nuq[((int64_t)e * nq) + q] = exp(aay[((int64_t)e * nq) + q]);
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
	free(nuq); free(aay); free(aaf); free(aayp); free(aafp); free(aady); free(aadf);
	P2.nuq = NULL;

	if (!ierr) {
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
			for (int k = 0; k < np; k++) {
				const double rr = zr[(k * np) + k], ll = zi[(k * np) + k];
				if ((rr <= 0) || (ll < 0)) {
					fprintf(fp_log, "*** port %d gives a non-passive impedance (R = %.6e, L = %.6e); "
						"check the boundary conditions (use awall)\n", k + 1, rr, ll / omega);
					ierr = 1;
				}
			}
			fprintf(fp_log, "  Z(1,1) = %.6e %+.6e j [ohm], R = %.6e, L = %.6e [H]\n",
				zr[0], zi[0], zr[0], zi[0] / omega);
		}
		free(zr); free(zi);
	}
	free(yr); free(yi);
	p2_close();
	return ierr;
}


/*
交流の発熱 (thermal + analysis = A、10 節点四面体)。σ(T) を要素ごとに渡して解き、端子電流 (波高値) cur での
損失 ½σ|e|²、e = jωA + ∇φ を**組み立てと同じ積分点**で節点の形状関数に掛けて節点荷重にする
(f_a = ∫ ½σ|e|² N_a dV)。Σ_a N_a = 1 なので Σ f = ½ eᴴ T e = ½|I|²R が厳密 (1 次の Q_e と同じ恒等式)。
*/
// 過渡解析 (analysis = A + transient、後退 Euler、電圧駆動、線形) の行列。1 次 (eddy3d.c の
// assemble_eddy3d_tr) と同じ形で χ = Δt φ を未知数にする : K = [S + T/Δt, TG/Δt; GᵀT/Δt, GᵀTG/Δt]、
// M は K から S を除いたもの (どちらも実対称)。K は P2.Kr、M は P2.Ki に入れる
static int p2_nonlin(int m)
{
	return (Material[m].ja.on || (Material[m].nbh[0] > 0));
}


// skipnl = 1 : 非線形 (J-A / B-H) の材料の要素の S を入れない (固定点反復が積分点の H から足す)
static int p2_assemble_tr(FILE *fp_log, double dt, int skipnl)
{
	const int ne = P2.ne, nf = P2.nf, nn = P2.nn, nt = NTet;
	crs_zero(&P2.Kr);
	crs_zero(&P2.Ki);
	for (int e = 0; e < nt; e++) {
		const int m = elem3d_mat(e);
		double nu[6];
		material_coef_pub(m, 3, nu);
		const double sig = Material[m].sigma;
		int32_t dof[P2_NDOF];
		double sgn[P2_NDOF];
		p2_dofs(e, P2.fe, ne, nf, nn, dof, sgn);
		double se[20][20], te[20][20], ce[20][10], pe[10][10];
		memset(se, 0, sizeof(se)); memset(te, 0, sizeof(te));
		memset(ce, 0, sizeof(ce)); memset(pe, 0, sizeof(pe));
		const int nq = elem3d_nip(e);
		for (int q = 0; q < nq; q++) {
			double wv[20][3], wc[20][3], pn[10], pg[10][3], w;
			if (p2_basis(e, q, wv, wc, pn, pg, &w)) {
				fprintf(fp_log, "*** the second-order A-phi solver : a degenerate element\n");
				return 1;
			}
			double nc[20][3];
			for (int k = 0; k < 20; k++) nu_mul(nu, wc[k], nc[k]);
			for (int k = 0; k < 20; k++) {
				for (int l = 0; l < 20; l++) {
					se[k][l] += w * dot3(nc[k], wc[l]);
					te[k][l] += w * sig * dot3(wv[k], wv[l]);
				}
				for (int a = 0; a < 10; a++) ce[k][a] += w * sig * dot3(wv[k], pg[a]);
			}
			for (int a = 0; a < 10; a++) {
				for (int b = 0; b < 10; b++) pe[a][b] += w * sig * dot3(pg[a], pg[b]);
			}
		}
		for (int k = 0; k < 20; k++) {
			for (int l = 0; l < 20; l++) {
				const int64_t p = p2_find(&P2.Kr, dof[k], dof[l]);
				const double s = sgn[k] * sgn[l];
				P2.Kr.val[p] += s * (((skipnl && p2_nonlin(m)) ? 0 : se[k][l]) + (te[k][l] / dt));
				P2.Ki.val[p] += s * te[k][l] / dt;
			}
			for (int a = 0; a < 10; a++) {
				const double v = sgn[k] * ce[k][a] / dt;
				const int64_t p1 = p2_find(&P2.Kr, dof[k], dof[20 + a]);
				const int64_t p2 = p2_find(&P2.Kr, dof[20 + a], dof[k]);
				P2.Kr.val[p1] += v;
				P2.Kr.val[p2] += v;
				P2.Ki.val[p1] += v;
				P2.Ki.val[p2] += v;
			}
		}
		for (int a = 0; a < 10; a++) {
			for (int b = 0; b < 10; b++) {
				const int64_t p = p2_find(&P2.Kr, dof[20 + a], dof[20 + b]);
				P2.Kr.val[p] += pe[a][b] / dt;
				P2.Ki.val[p] += pe[a][b] / dt;
			}
		}
	}
	return 0;
}


// 刻みの終わり t の電圧の倍率 (waveform の折れ線の左極限、無ければ 1。1 次と同じ)
static double p2_wave(double t)
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
2 次の過渡解析の非線形 (J-A / B-H)。1 次 (eddy3d.c の ed_ja_*) と同じ固定点反復を**積分点ごとの状態**で解く
(2 次の要素では B = ∇×A が要素内で 1 次に変わる)。各刻み:

	F(x) = K_l x − M [A^n; 0] + Σ_e Σ_q w_q H_q・∇×W_k   (K_l は非線形の要素の S を除いた行列)
	(K_l + Σ w_q (∇×W)ᵀ N_q (∇×W)) δ = −F,   N_q = dH/dB の対称部分 (J-A) / ν I + 2 (dν/dB²) B Bᵀ (B-H)

∇×W は積分点あたり 60 個 (要素あたり 30 KB) なので、合計が P2NL_CACHE_MAX 個以内なら記憶し、越えたら毎回
p2_basis で作り直す (答えは同じ。速さだけが違う)。
*/
typedef struct {
	int nj, nq;
	int *el;			// [nj] 要素番号
	int32_t *dof;		// [nj*20] A の自由度
	double *sg;			// [nj*20] 向きの符号
	double *w;			// [nj*nq] 積分の重み
	double *b, *h, *m, *d;		// [nj*nq*3] 確定済みの状態
	double *bn, *hn, *mn, *dn;	// [nj*nq*3] 反復中
	double *nu6;		// [nj*nq*6] 行列の傾き
	double *loc;		// [nj*400] 要素の作業域
	double *cv;			// [nj*nq*60] 符号込みの ∇×W (記憶するときだけ、しないなら NULL)
} p2nl_t;

#define P2NL_CACHE_MAX ((int64_t)32 * 1024 * 1024)	// ∇×W を記憶する上限 (double の個数、256 MB)


static int p2nl_open(p2nl_t *J, FILE *fp_log)
{
	const int nt = NTet;
	int nj = 0;
	for (int e = 0; e < nt; e++) {
		if (p2_nonlin(elem3d_mat(e))) nj++;
	}
	const int nq = elem3d_nip(0);
	J->nj = nj;
	J->nq = nq;
	const size_t s1 = (size_t)((nj > 0) ? nj : 1), sq = s1 * (size_t)nq;
	J->el = (int *)malloc(s1 * sizeof(int));
	J->dof = (int32_t *)malloc(s1 * 20 * sizeof(int32_t));
	J->sg = (double *)malloc(s1 * 20 * sizeof(double));
	J->w = (double *)calloc(sq, sizeof(double));
	double **v3[8] = {&J->b, &J->h, &J->m, &J->d, &J->bn, &J->hn, &J->mn, &J->dn};
	for (int q = 0; q < 8; q++) *v3[q] = (double *)calloc(sq * 3, sizeof(double));
	J->nu6 = (double *)calloc(sq * 6, sizeof(double));
	J->loc = (double *)calloc(s1 * 400, sizeof(double));
	int j = 0;
	for (int e = 0; e < nt; e++) {
		if (!p2_nonlin(elem3d_mat(e))) continue;
		int32_t dof[P2_NDOF];
		double sgn[P2_NDOF];
		p2_dofs(e, P2.fe, P2.ne, P2.nf, P2.nn, dof, sgn);
		J->el[j] = e;
		for (int k = 0; k < 20; k++) {
			J->dof[(j * 20) + k] = dof[k];
			J->sg[(j * 20) + k] = sgn[k];
		}
		if (j == 0) {
			J->cv = ((((int64_t)nj * nq * 60) <= P2NL_CACHE_MAX) ? (double *)malloc(sq * 60 * sizeof(double)) : NULL);
		}
		for (int q = 0; q < nq; q++) {
			double wv[20][3], wc[20][3], pn[10], pg[10][3], w;
			if (p2_basis(e, q, wv, wc, pn, pg, &w)) {
				fprintf(fp_log, "*** the second-order A-phi solver : a degenerate element\n");
				return 1;
			}
			J->w[((int64_t)j * nq) + q] = w;
			if (J->cv != NULL) {
				double *cv = &J->cv[(((int64_t)j * nq) + q) * 60];
				for (int k = 0; k < 20; k++) {
					for (int c = 0; c < 3; c++) cv[(k * 3) + c] = sgn[k] * wc[k][c];
				}
			}
		}
		j++;
	}
	fprintf(fp_log, "  nonlinear : %d tetrahedra x %d points with J-A (vector) / bh (nlsolver = %d %.3e %.3f)\n",
		nj, nq, NlMaxiter, NlTol, NlRelax);
	for (int m = 0; m < NMaterial; m++) {
		if (!Material[m].ja.on) continue;
		fprintf(fp_log, "  J-A model check (material %d) : 3-D vector J-A along a fixed direction == scalar J-A : "
			"max |dH| / max |H| = %.3e\n", m, ja_vec3_selftest(&Material[m].ja));
	}
	return 0;
}


static void p2nl_close(p2nl_t *J)
{
	free(J->el); free(J->dof); free(J->sg); free(J->w);
	free(J->b); free(J->h); free(J->m); free(J->d);
	free(J->bn); free(J->hn); free(J->mn); free(J->dn);
	free(J->nu6); free(J->loc); free(J->cv);
	memset(J, 0, sizeof(*J));
}


// 要素 j の積分点 q の符号込みの ∇×W_k
static void p2nl_curl(const p2nl_t *J, int j, int q, double wc[20][3])
{
	if (J->cv != NULL) {
		const double *cv = &J->cv[(((int64_t)j * J->nq) + q) * 60];
		for (int k = 0; k < 20; k++) {
			for (int c = 0; c < 3; c++) wc[k][c] = cv[(k * 3) + c];
		}
		return;
	}
	double wv[20][3], pn[10], pg[10][3], w;
	if (p2_basis(J->el[j], q, wv, wc, pn, pg, &w)) {
		memset(wc, 0, 20 * 3 * sizeof(double));
		return;
	}
	for (int k = 0; k < 20; k++) {
		for (int c = 0; c < 3; c++) wc[k][c] *= J->sg[(j * 20) + k];
	}
}


// x の A から各積分点の B を作り、確定済みの状態から H と行列の傾きを求める (1 次の ed_ja_eval と同じ式)
static void p2nl_eval(p2nl_t *J, const double *x)
{
	const int nj = J->nj, nq = J->nq;
	int j;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (j = 0; j < nj; j++) {
		const material_t *mt = &Material[elem3d_mat(J->el[j])];
		const ja_t *ja = &mt->ja;
		for (int q = 0; q < nq; q++) {
			const int64_t p = ((int64_t)j * nq) + q;
			double wc[20][3];
			p2nl_curl(J, j, q, wc);
			double bv[3] = {0, 0, 0};
			for (int k = 0; k < 20; k++) {
				const double a = x[J->dof[(j * 20) + k]];
				for (int c = 0; c < 3; c++) bv[c] += a * wc[k][c];
			}
			const double bmag = sqrt(dot3(bv, bv));
			const double *dir = &J->d[p * 3];
			double *dn = &J->dn[p * 3];
			const int haved = ((dir[0] != 0) || (dir[1] != 0) || (dir[2] != 0));
			for (int c = 0; c < 3; c++) dn[c] = ((bmag > 1e-12) ? (bv[c] / bmag) : 0);
			double *nu6 = &J->nu6[p * 6];
			double *hv = &J->hn[p * 3];
			if (!ja->on) {
				double nu, dnu;
				bh_nu_dnu_pub(mt, bmag, &nu, &dnu);
				for (int c = 0; c < 3; c++) {
					hv[c] = nu * bv[c];
					J->mn[(p * 3) + c] = (bv[c] / MU0) - hv[c];
					J->bn[(p * 3) + c] = bv[c];
				}
				nu6[0] = nu + (2 * dnu * bv[0] * bv[0]);
				nu6[1] = nu + (2 * dnu * bv[1] * bv[1]);
				nu6[2] = nu + (2 * dnu * bv[2] * bv[2]);
				nu6[3] = 2 * dnu * bv[0] * bv[1];
				nu6[4] = 2 * dnu * bv[1] * bv[2];
				nu6[5] = 2 * dnu * bv[2] * bv[0];
				continue;
			}
			double dref[3] = {1, 0, 0};
			if (haved) { dref[0] = dir[0]; dref[1] = dir[1]; dref[2] = dir[2]; }
			else if (bmag > 1e-12) { dref[0] = dn[0]; dref[1] = dn[1]; dref[2] = dn[2]; }
			ja_eval_vec3(ja, &J->b[p * 3], &J->h[p * 3], &J->m[p * 3], bv, dref, hv, &J->mn[p * 3], nu6);
			for (int c = 0; c < 3; c++) J->bn[(p * 3) + c] = bv[c];
		}
	}
}


// 刻みの磁気の仕事 Σ w_q H_q・(B_q − B_q^n)
static double p2nl_work(const p2nl_t *J)
{
	double wsum = 0;
	const int64_t np = (int64_t)J->nj * J->nq;
	for (int64_t p = 0; p < np; p++) {
		double s = 0;
		for (int c = 0; c < 3; c++) s += J->hn[(p * 3) + c] * (J->bn[(p * 3) + c] - J->b[(p * 3) + c]);
		wsum += J->w[p] * s;
	}
	return wsum;
}


static void p2nl_commit(p2nl_t *J)
{
	const int64_t np = (int64_t)J->nj * J->nq;
	memcpy(J->b, J->bn, (size_t)np * 3 * sizeof(double));
	memcpy(J->h, J->hn, (size_t)np * 3 * sizeof(double));
	memcpy(J->m, J->mn, (size_t)np * 3 * sizeof(double));
	for (int64_t p = 0; p < np; p++) {
		double *d = &J->d[p * 3];
		if ((d[0] == 0) && (d[1] == 0) && (d[2] == 0)) {
			for (int c = 0; c < 3; c++) d[c] = J->dn[(p * 3) + c];
		}
	}
}


// 残差 r = −F(x) (固定の行は 0) と、その相対ノルム
static double p2nl_resid(p2nl_t *J, const crs_t *Kl, const double *rhs, const unsigned char *fix,
	const double *x, double *r, double bref)
{
	const int n = (int)Kl->n, nj = J->nj, nq = J->nq;
	crs_spmv(Kl, x, r, NULL);
	for (int i = 0; i < n; i++) r[i] = rhs[i] - r[i];
	int j;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (j = 0; j < nj; j++) {
		double *lv = &J->loc[(int64_t)j * 400];
		for (int k = 0; k < 20; k++) lv[k] = 0;
		for (int q = 0; q < nq; q++) {
			const int64_t p = ((int64_t)j * nq) + q;
			double wc[20][3];
			p2nl_curl(J, j, q, wc);
			for (int k = 0; k < 20; k++) lv[k] += J->w[p] * dot3(&J->hn[p * 3], wc[k]);
		}
	}
	for (j = 0; j < nj; j++) {
		for (int k = 0; k < 20; k++) r[J->dof[(j * 20) + k]] -= J->loc[((int64_t)j * 400) + k];
	}
	double s = 0;
	for (int i = 0; i < n; i++) {
		if (fix[i]) r[i] = 0;
		s += r[i] * r[i];
	}
	return ((bref > 0) ? (sqrt(s) / bref) : sqrt(s));
}


// 接線の行列 Kt = K_l + Σ w_q (∇×W)ᵀ N_q (∇×W)
static void p2nl_tangent(p2nl_t *J, const crs_t *Kl, crs_t *Kt)
{
	const int nj = J->nj, nq = J->nq;
	memcpy(Kt->val, Kl->val, (size_t)Kl->rowptr[Kl->n] * sizeof(double));
	int j;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (j = 0; j < nj; j++) {
		double *ke = &J->loc[(int64_t)j * 400];
		for (int k = 0; k < 400; k++) ke[k] = 0;
		for (int q = 0; q < nq; q++) {
			const int64_t p = ((int64_t)j * nq) + q;
			double wc[20][3], nc[20][3];
			p2nl_curl(J, j, q, wc);
			for (int k = 0; k < 20; k++) nu_mul(&J->nu6[p * 6], wc[k], nc[k]);
			for (int k = 0; k < 20; k++) {
				for (int l = 0; l < 20; l++) ke[(k * 20) + l] += J->w[p] * dot3(nc[k], wc[l]);
			}
		}
	}
	for (j = 0; j < nj; j++) {
		const double *ke = &J->loc[(int64_t)j * 400];
		for (int k = 0; k < 20; k++) {
			for (int l = 0; l < 20; l++) {
				Kt->val[p2_find(Kt, J->dof[(j * 20) + k], J->dof[(j * 20) + l])] += ke[(k * 20) + l];
			}
		}
	}
}


// 1 刻みの固定点反復 (1 次の ed_ja_step と同じ安定化)。x は前の刻みの解で入り、収束した解で出る
static int p2nl_step(p2nl_t *J, const crs_t *Kl, crs_t *Kt, const double *rhs, const unsigned char *fix,
	const double *xd, double *x, double *r, double *del, double *xbest, double bref, int st, FILE *fp_log)
{
	const int n = (int)Kl->n;
	for (int i = 0; i < n; i++) {
		if (fix[i]) x[i] = xd[i];
	}
	int it, conv = 0, nback = 0;
	double rnorm = 0, rbest = 0, lrel = NlRelax;
	for (it = 1; it <= NlMaxiter; it++) {
		p2nl_eval(J, x);
		rnorm = p2nl_resid(J, Kl, rhs, fix, x, r, bref);
		if (rnorm < NlTol) {
			conv = 1;
			break;
		}
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
		p2nl_tangent(J, Kl, Kt);
		double tol = ((rnorm > 0) ? (Solver.converg / rnorm) : Solver.converg);
		if (tol > 0.1) tol = 0.1;
		if (tol < Solver.converg) tol = Solver.converg;
		for (int i = 0; i < n; i++) del[i] = 0;
		const int li = solver_cg(Kt, r, del, fix, Solver.maxiter, Solver.nout, tol, NULL, "a3j2");
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


/*
3 次元渦電流の過渡解析の 2 次の要素 (10 節点四面体)。手順と出力は 1 次 (solve_eddy3d_transient) と同じ :
ポート 1 の電極の頂点を χ = Δt V(t^{n+1}) (泡は 0)、右辺 M [A^n; 0]、端子電流は電極の頂点の行の
反作用 Kx − M[A^n; 0]。仕事の恒等式 V I Δt = (ΔA)ᵀ S A^{n+1} + zᵀMz をログに出す。
*/
int solve_eddy3d_p2_transient(FILE *fp_log)
{
	int ierr = 0;
	fprintf(fp_log, "\n=== 3D eddy current, transient (A-phi, second-order Nedelec / nodal elements, backward Euler) ===\n");
	if (NPort != 1) {
		fprintf(fp_log, "*** analysis A + transient needs exactly one port\n");
		return 1;
	}
	const double fsave = Freq;
	if (Freq <= 0) Freq = 1;		// p2_open の検査を通すだけ (過渡解析は ω を使わない)
	const int oe = p2_open(fp_log);
	Freq = fsave;
	if (oe) return 1;
	const int n = P2.n, nn = P2.nn, oV = P2.oV, oF = P2.oF;
	const double dt = TransientTime / TransientSteps;
	int hyst = 0;
	for (int m = 0; m < NMaterial; m++) {
		if (p2_nonlin(m)) hyst = 1;
	}
	fprintf(fp_log, "  t_end = %.6e [s], %d steps (dt = %.6e [s]), voltage = %.6e [V]%s\n",
		TransientTime, TransientSteps, dt, Volt, ((WaveN > 0) ? " x waveform" : " (step)"));
	if (p2_assemble_tr(fp_log, dt, hyst)) {
		p2_close();
		return 1;
	}
	crs_t *K = &P2.Kr, *M = &P2.Ki;
	double *x = (double *)calloc((size_t)n, sizeof(double));
	double *xo = (double *)calloc((size_t)n, sizeof(double));
	double *xd = (double *)calloc((size_t)n, sizeof(double));
	double *rhs = (double *)malloc((size_t)n * sizeof(double));
	double *b = (double *)malloc((size_t)n * sizeof(double));
	double *w1 = (double *)malloc((size_t)n * sizeof(double));
	double *u = (double *)calloc((size_t)n, sizeof(double));
	double *z = (double *)malloc((size_t)n * sizeof(double));
	p2nl_t J;
	memset(&J, 0, sizeof(J));
	crs_t Kt;
	memset(&Kt, 0, sizeof(Kt));
	if (hyst) {
		if (p2nl_open(&J, fp_log)) ierr = 1;
		p2_crs(&Kt, n, P2.fe, P2.ne, P2.nf, nn);
	}
	double bref = 0;
	long fptotal = 0;
	FILE *fcsv = mpi_fopen_out("ofe_transient.csv", "w");
	if (fcsv) fprintf(fcsv, "t [s],V [V],I [A],W_in [J],W_mag [J],W_joule [J],iter\n");
	double win = 0, wmag = 0, wjou = 0;
	fprintf(fp_log, "  %6s %13s %13s %13s %6s\n", "step", "t [s]", "V [V]", "I [A]", "iter");
	for (int st = 1; (st <= TransientSteps) && !ierr; st++) {
		const double t = st * dt;
		const double v = Volt * p2_wave(t);
		// A の自由度 (φ の欄より前) だけを残した [A^n; 0]
		for (int i = 0; i < n; i++) xo[i] = ((i < oV) ? x[i] : 0);
		crs_spmv(M, xo, rhs, NULL);
		for (int i = 0; i < n; i++) xd[i] = 0;
		for (int i = 0; i < nn; i++) {
			if (P2.fix[oV + i] && P2.cond[i] && P2.vert[i] && (NodeConductor[i] == 1)) xd[oV + i] = dt * v;
		}
		crs_spmv(K, xd, w1, NULL);
		for (int i = 0; i < n; i++) b[i] = (P2.fix[i] ? 0 : (rhs[i] - w1[i]));
		int it;
		if (!hyst) {
			for (int i = 0; i < n; i++) u[i] = (P2.fix[i] ? 0 : x[i]);
			it = solver_cg(K, b, u, P2.fix, Solver.maxiter, Solver.nout, Solver.converg, NULL, "a3t2");
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
			it = p2nl_step(&J, K, &Kt, rhs, P2.fix, xd, x, b, u, z, bref, st, fp_log);
			if (it < 0) {
				ierr = 1;
				break;
			}
			fptotal += it;
		}
		crs_spmv(K, x, w1, NULL);
		double cur = 0;
		for (int i = 0; i < nn; i++) {
			if (P2.fix[oV + i] && P2.cond[i] && P2.vert[i] && (NodeConductor[i] == 1)) cur += w1[oV + i] - rhs[oV + i];
		}
		// 仕事 : 端子 V I Δt、Joule zᵀ M z、磁気 (ΔA)ᵀ S A^{n+1} = zᵀ (K − M) x の A の部分
		for (int i = 0; i < n; i++) z[i] = x[i] - xo[i];
		double zmz = 0, zkx = 0, zmx = 0;
		crs_spmv(M, z, w1, NULL);
		for (int i = 0; i < n; i++) zmz += z[i] * w1[i];
		crs_spmv(M, x, w1, NULL);
		for (int i = 0; i < oV; i++) zmx += z[i] * w1[i];
		crs_spmv(K, x, w1, NULL);
		for (int i = 0; i < oV; i++) zkx += z[i] * w1[i];
		if (hyst) zkx += p2nl_work(&J);
		win += v * cur * dt;
		wjou += zmz;
		wmag += zkx - zmx;
		if (hyst) p2nl_commit(&J);
		if (fcsv) fprintf(fcsv, "%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%d\n", t, v, cur, win, wmag, wjou, it);
		if ((st == TransientSteps) || ((st % ((TransientSteps >= 20) ? (TransientSteps / 20) : 1)) == 0)) {
			fprintf(fp_log, "  %6d %13.6e %13.6e %13.6e %6d\n", st, t, v, cur, it);
		}
	}
	(void)oF;
	if (fcsv) fclose(fcsv);
	if (hyst && !ierr) {
		fprintf(fp_log, "  nonlinear (%d tetrahedra) : fixed-point iterations = %ld (%d steps)\n",
			J.nj, fptotal, TransientSteps);
	}
	if (hyst) {
		p2nl_close(&J);
		crs_free(&Kt);
	}
	fprintf(fp_log, "  energy : terminals = %.12e, magnetic work = %.12e, joule = %.12e [J] "
		"(terminals - magnetic - joule = %.3e)\n", win, wmag, wjou, win - wmag - wjou);
	free(x); free(xo); free(xd); free(rhs); free(b); free(w1); free(u); free(z);
	p2_close();
	return ierr;
}


/*
辺要素の自己検証の 2 次版 (analysis = E、10 節点四面体)。2 次の Nedelec 空間 (第 1 種) は P1³ ⊕ S2
(S2 = 2 次の同次式で p・r = 0) を含むので、そこに入る場を σ で重みをつけた L2 射影すると誤差が 0 になり、
回転のエネルギーも厳密に合う:

	T u = b,   b_k = ∫ σ W_k・E dV,   誤差² = ∫ σ|E|² − uᵀ b = 0,   uᵀ S u = ∫ (∇×E)ᵀ ν (∇×E) dV

(P2a) 一般の 1 次の場 E = a + B r (B は対称部分と反対称部分の両方、軸に平行でない向き) は Whitney と
∇(λaλb) の両方を使い、(P2b) それに S2 の 2 次の場 (yz, −xz, 0) を足すと面の関数が要る (1 次の空間には
入らない)。厳密なのは直線の要素だけ (曲がった要素では共変 Piola で写した多項式になる)。
*/
int solve_edge_test_p2(FILE *fp_log)
{
	fprintf(fp_log, "\n=== edge element (second-order Nedelec) self test ===\n");
	if (!MeshMode || (TetOrder != 2) || (NHex + NPrism + NPyr > 0)) {
		fprintf(fp_log, "*** the second-order edge self test needs a pure 10-node tetrahedral mesh\n");
		return 1;
	}
	// 恒等式が厳密なのは直線の要素だけ (曲がった要素では共変 Piola で写した多項式になり、射影の誤差が
	// 0 にならない)。曲がった格子を誤りと区別できるよう、中間節点が辺の中点から外れていれば先に弾く
	for (int e = 0; e < NTet; e++) {
		for (int m = 0; m < 6; m++) {
			const int32_t a = Tet[(e * 4) + P2_MID[m][0]], b = Tet[(e * 4) + P2_MID[m][1]];
			const int32_t c = Tet2[(e * 6) + m];
			const double d[3] = {Xp[c] - ((Xp[a] + Xp[b]) / 2), Yp[c] - ((Yp[a] + Yp[b]) / 2), Zp[c] - ((Zp[a] + Zp[b]) / 2)};
			const double h[3] = {Xp[b] - Xp[a], Yp[b] - Yp[a], Zp[b] - Zp[a]};
			if (dot3(d, d) > (1e-16 * dot3(h, h))) {
				fprintf(fp_log, "*** the second-order edge self test needs straight-sided elements "
					"(tetrahedron %d has a mid-edge node off the midpoint); the identities are exact "
					"only for affine elements\n", e + 1);
				return 1;
			}
		}
	}
	edge_build();
	const int ne = NEdge, nn = NNode, nt = NTet;
	int nf = 0;
	int32_t *fvert = NULL;
	int32_t *fe = p2_faces(&nf, &fvert);
	const int oV = (2 * ne) + (2 * nf), n = oV + nn + ne;
	fprintf(fp_log, "  nodes = %d, tetrahedra = %d, edges = %d, faces = %d, edge unknowns = %d\n", nn, nt, ne, nf, oV);
	crs_t S, T;
	p2_crs(&S, n, fe, ne, nf, nn);
	p2_crs(&T, n, fe, ne, nf, nn);
	// 中心と大きさ (場の係数を格子の大きさにそろえる)
	double c0[3] = {0, 0, 0}, lmax = 0;
	{
		double mn[3] = {1e300, 1e300, 1e300}, mx[3] = {-1e300, -1e300, -1e300};
		for (int i = 0; i < nn; i++) {
			const double r[3] = {Xp[i], Yp[i], Zp[i]};
			for (int c = 0; c < 3; c++) {
				if (r[c] < mn[c]) mn[c] = r[c];
				if (r[c] > mx[c]) mx[c] = r[c];
			}
		}
		for (int c = 0; c < 3; c++) {
			c0[c] = (mn[c] + mx[c]) / 2;
			if ((mx[c] - mn[c]) > lmax) lmax = mx[c] - mn[c];
		}
	}
	if (!(lmax > 0)) lmax = 1;
	const double av[3] = {1.0, 2.0, 3.0};
	const double bm[3][3] = {{0.3, -1.2, 0.7}, {0.5, 0.2, -0.4}, {1.1, 0.6, -0.9}};
	int ierr = 0;
	for (int e = 0; e < nt; e++) {
		if (!(Material[elem3d_mat(e)].sigma > 0)) {
			fprintf(fp_log, "*** the second-order edge self test needs sigma > 0 in every element\n");
			ierr = 1;
			break;
		}
	}
	double *b1 = (double *)calloc((size_t)n, sizeof(double));
	double *b2 = (double *)calloc((size_t)n, sizeof(double));
	double se1 = 0, se2 = 0, sc1 = 0, sc2 = 0;		// ∫σ|E|²、∫(∇×E)ᵀν(∇×E)
	for (int e = 0; (e < nt) && !ierr; e++) {
		const int m = elem3d_mat(e);
		double nu[6];
		material_coef_pub(m, 3, nu);
		const double sig = Material[m].sigma;
		int32_t dof[P2_NDOF];
		double sgn[P2_NDOF];
		p2_dofs(e, fe, ne, nf, nn, dof, sgn);
		double sk[20][20], tk[20][20], f1[20], f2[20];
		memset(sk, 0, sizeof(sk)); memset(tk, 0, sizeof(tk));
		memset(f1, 0, sizeof(f1)); memset(f2, 0, sizeof(f2));
		const int nq = elem3d_nip(e);
		for (int q = 0; q < nq; q++) {
			double wv[20][3], wc[20][3], pn[10], pg[10][3], w, nv[ELEM_NMAX], g[ELEM_NMAX][3], w2;
			if (p2_basis(e, q, wv, wc, pn, pg, &w) || elem3d_ip(e, q, nv, g, &w2)) { ierr = 1; break; }
			// 積分点の位置 (等パラメトリック : x = Σ x_a N_a)
			int32_t nd[TET_NMAX];
			tet_nodes(e, nd);
			double r[3] = {0, 0, 0};
			for (int a = 0; a < 10; a++) {
				r[0] += nv[a] * (Xp[nd[a]] - c0[0]) / lmax;
				r[1] += nv[a] * (Yp[nd[a]] - c0[1]) / lmax;
				r[2] += nv[a] * (Zp[nd[a]] - c0[2]) / lmax;
			}
			// E1 = a + B r'、E2 = E1 + (y'z', −x'z', 0) (r' = (r − c0)/L)。回転は B の反対称部分 (/L) と
			// 2 次の場の回転 (x', y', −2z')/L
			double e1[3], c1[3], e2[3], c2[3];
			for (int i = 0; i < 3; i++) e1[i] = av[i] + (bm[i][0] * r[0]) + (bm[i][1] * r[1]) + (bm[i][2] * r[2]);
			c1[0] = (bm[2][1] - bm[1][2]) / lmax;
			c1[1] = (bm[0][2] - bm[2][0]) / lmax;
			c1[2] = (bm[1][0] - bm[0][1]) / lmax;
			e2[0] = e1[0] + (r[1] * r[2]); e2[1] = e1[1] - (r[0] * r[2]); e2[2] = e1[2];
			c2[0] = c1[0] + (r[0] / lmax); c2[1] = c1[1] + (r[1] / lmax); c2[2] = c1[2] - (2 * r[2] / lmax);
			double n1[3], n2[3];
			nu_mul(nu, c1, n1);
			nu_mul(nu, c2, n2);
			se1 += w * sig * dot3(e1, e1);
			se2 += w * sig * dot3(e2, e2);
			sc1 += w * dot3(n1, c1);
			sc2 += w * dot3(n2, c2);
			double nc[20][3];
			for (int k = 0; k < 20; k++) nu_mul(nu, wc[k], nc[k]);
			for (int k = 0; k < 20; k++) {
				f1[k] += w * sig * dot3(wv[k], e1);
				f2[k] += w * sig * dot3(wv[k], e2);
				for (int l = 0; l < 20; l++) {
					sk[k][l] += w * dot3(nc[k], wc[l]);
					tk[k][l] += w * sig * dot3(wv[k], wv[l]);
				}
			}
		}
		for (int k = 0; k < 20; k++) {
			b1[dof[k]] += sgn[k] * f1[k];
			b2[dof[k]] += sgn[k] * f2[k];
			for (int l = 0; l < 20; l++) {
				const int64_t p = p2_find(&S, dof[k], dof[l]);
				S.val[p] += sgn[k] * sgn[l] * sk[k][l];
				T.val[p] += sgn[k] * sgn[l] * tk[k][l];
			}
		}
	}
	unsigned char *fix = (unsigned char *)calloc((size_t)n, sizeof(unsigned char));
	for (int i = oV; i < n; i++) fix[i] = 1;
	double *u = (double *)calloc((size_t)n, sizeof(double));
	double *w1 = (double *)malloc((size_t)n * sizeof(double));
	const char *lab[2] = {"(P2a) linear field E = a + B r (general B)", "(P2b) a + B r + (yz, -xz, 0) (needs the face functions)"};
	double *bb[2] = {b1, b2};
	const double sev[2] = {se1, se2}, scv[2] = {sc1, sc2};
	for (int t = 0; (t < 2) && !ierr; t++) {
		for (int i = 0; i < n; i++) u[i] = 0;
		if (solver_cg(&T, bb[t], u, fix, Solver.maxiter, Solver.nout, Solver.converg, NULL, "p2e") < 0) {
			fprintf(fp_log, "*** the L2 projection did not converge\n");
			ierr = 1;
			break;
		}
		double ub = 0, us = 0;
		for (int i = 0; i < oV; i++) ub += u[i] * bb[t][i];
		crs_spmv(&S, u, w1, fix);
		for (int i = 0; i < oV; i++) us += u[i] * w1[i];
		const double ep = fabs(sev[t] - ub) / sev[t];
		const double ec = fabs(us - scv[t]) / scv[t];
		const int ok = ((ep <= 1e-10) && (ec <= 1e-9));
		fprintf(fp_log, "  %s : projection error^2 / int sigma|E|^2 = %.3e, curl energy u^T S u = %.12e (exact %.12e, %.3e) %s\n",
			lab[t], ep, us, scv[t], ec, (ok ? "passed" : "FAILED"));
		if (!ok) ierr = 1;
	}
	fprintf(fp_log, "  second-order edge self test : %s\n", (ierr ? "FAILED" : "passed"));
	free(b1); free(b2); free(u); free(w1); free(fix);
	free(fe); free(fvert);
	crs_free(&S);
	crs_free(&T);
	edge_free();
	return ierr;
}


int eddy3d_p2_heat_open(FILE *fp_log)
{
	fprintf(fp_log, "  AC heating (analysis = A, second-order Nedelec) : frequency = %.6e [Hz]\n", Freq);
	return p2_open(fp_log);
}


int eddy3d_p2_heat_solve(FILE *fp_log, const double *sigelem, double cur, double *f,
	double *zr, double *zi, int it)
{
	const int ne = P2.ne, nf = P2.nf, nn = P2.nn, nt = NTet, n = P2.n;
	const double omega = P2.omega;
	if (p2_assemble(fp_log, sigelem)) return 1;
	char label[64];
	sprintf(label, "A%d", it);
	double qr[MAXPORT], qi[MAXPORT];
	if (p2_solve_port(fp_log, 1, NULL, label, qr, qi)) return 1;
	const double ir = qr[1], ii = qi[1];
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
		const double xr = P2.ur[i], xi = P2.ui[i];
		P2.ur[i] = (sr * xr) - (si * xi);
		P2.ui[i] = (sr * xi) + (si * xr);
	}
	for (int i = 0; i < nn; i++) f[i] = 0;
	for (int e = 0; e < nt; e++) {
		const double sg = sigelem[e];
		if (sg <= 0) continue;
		int32_t dof[P2_NDOF];
		double sgn[P2_NDOF];
		p2_dofs(e, P2.fe, ne, nf, nn, dof, sgn);
		int32_t nd[TET_NMAX];
		tet_nodes(e, nd);
		const int nq = elem3d_nip(e);
		for (int q = 0; q < nq; q++) {
			double wv[20][3], wc[20][3], pn[10], pg[10][3], w;
			double nv[ELEM_NMAX], g[ELEM_NMAX][3], w2;
			if (p2_basis(e, q, wv, wc, pn, pg, &w) || elem3d_ip(e, q, nv, g, &w2)) return 1;
			// e = jω A + ∇φ : 実部 = −ω A_i + ∇φ_r、虚部 = ω A_r + ∇φ_i
			double er[3] = {0, 0, 0}, ei[3] = {0, 0, 0};
			for (int k = 0; k < 20; k++) {
				const double ar = sgn[k] * P2.ur[dof[k]], ai = sgn[k] * P2.ui[dof[k]];
				for (int c = 0; c < 3; c++) {
					er[c] += -omega * ai * wv[k][c];
					ei[c] += omega * ar * wv[k][c];
				}
			}
			for (int a = 0; a < 10; a++) {
				const double pr = P2.ur[dof[20 + a]], pi = P2.ui[dof[20 + a]];
				for (int c = 0; c < 3; c++) {
					er[c] += pr * pg[a][c];
					ei[c] += pi * pg[a][c];
				}
			}
			const double qd = 0.5 * sg * (dot3(er, er) + dot3(ei, ei)) * w;
			for (int a = 0; a < 10; a++) f[nd[a]] += qd * nv[a];
		}
	}
	return 0;
}


void eddy3d_p2_heat_close(void)
{
	p2_close();
}
