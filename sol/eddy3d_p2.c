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

対応するのは周波数領域の線形材料・ゲージなし・場の出力なし (それ以外は入力段で弾く)。
*/

#include "fem.h"
#include "fem_prototype.h"

static const int P2_EDGE[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
static const int P2_FACE[4][3] = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};
// gmsh の tet10 の中間節点 (局所 4..9) の辺 : (0,1) (1,2) (2,0) (3,0) (3,2) (3,1)
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


// 全体の面と、要素ごとの面番号 fe[4 e + f]
static int32_t *p2_faces(int *nface, int32_t **fvert)
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
static int p2_face_id(const int32_t *fvert, int nface, const int32_t v[3])
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


// 積分点での 20 個の A の基底 (値と回転、局所の向き) と 10 個の φ の基底 (値と勾配)
static int p2_basis(int e, int q, double wv[20][3], double wc[20][3], double pn[10], double pg[10][3], double *wq)
{
	double nv[ELEM_NMAX], g[ELEM_NMAX][3];
	if (elem3d_ip(e, q, nv, g, wq)) return 1;
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
	return 0;
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


int solve_eddy3d_p2(FILE *fp_log)
{
	const int np = NPort;
	int ierr = 0;
	fprintf(fp_log, "\n=== 3D eddy current (A-phi, second-order Nedelec / nodal elements) ===\n");
	if (!MeshMode || (TetOrder != 2) || (NHex + NPrism + NPyr > 0) || (np < 1) || (Freq <= 0)) {
		fprintf(fp_log, "*** the second-order A-phi solver needs a pure 10-node tetrahedral mesh, "
			"a port and a frequency\n");
		return 1;
	}
	edge_build();
	const int ne = NEdge, nn = NNode, nt = NTet;
	int nf = 0;
	int32_t *fvert = NULL;
	int32_t *fe = p2_faces(&nf, &fvert);
	const int oG = ne, oF = 2 * ne, oV = (2 * ne) + (2 * nf), oB = oV + nn;
	const int n = oB + ne;
	const double omega = 2 * PI * Freq;
	fprintf(fp_log, "  nodes = %d, tetrahedra = %d, edges = %d, faces = %d, unknowns = %d "
		"(A : 2 x %d edges + 2 x %d faces, phi : vertices + %d edge bubbles)\n", nn, nt, ne, nf, n, ne, nf, ne);

	crs_t Kr, Ki;
	p2_crs(&Kr, n, fe, ne, nf, nn);
	p2_crs(&Ki, n, fe, ne, nf, nn);
	fprintf(fp_log, "  coupled matrix : %lld nonzeros (%.1f per row)\n", (long long)Kr.nnz, (double)Kr.nnz / n);

	// 組み立て : [S + jωT, TG; GᵀT, GᵀTG/(jω)] (節点の式を jω で割った複素対称形、1 次と同じ)
	int bad = 0;
	for (int e = 0; e < nt; e++) {
		const int m = elem3d_mat(e);
		double nu[6];
		material_coef_pub(m, 3, nu);
		const double sig = Material[m].sigma;
		int32_t dof[P2_NDOF];
		double sgn[P2_NDOF];
		p2_dofs(e, fe, ne, nf, nn, dof, sgn);
		double se[20][20], te[20][20], ce[20][10], pe[10][10];
		memset(se, 0, sizeof(se)); memset(te, 0, sizeof(te));
		memset(ce, 0, sizeof(ce)); memset(pe, 0, sizeof(pe));
		const int nq = elem3d_nip(e);
		for (int q = 0; q < nq; q++) {
			double wv[20][3], wc[20][3], pn[10], pg[10][3], w;
			if (p2_basis(e, q, wv, wc, pn, pg, &w)) { bad = 1; break; }
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
				const int64_t p = p2_find(&Kr, dof[k], dof[l]);
				const double s = sgn[k] * sgn[l];
				Kr.val[p] += s * se[k][l];
				Ki.val[p] += s * omega * te[k][l];
			}
			for (int a = 0; a < 10; a++) {
				const double v = sgn[k] * ce[k][a];
				Kr.val[p2_find(&Kr, dof[k], dof[20 + a])] += v;
				Kr.val[p2_find(&Kr, dof[20 + a], dof[k])] += v;
			}
		}
		for (int a = 0; a < 10; a++) {
			for (int b = 0; b < 10; b++) Ki.val[p2_find(&Ki, dof[20 + a], dof[20 + b])] -= pe[a][b] / omega;
		}
	}
	if (bad) {
		fprintf(fp_log, "*** the second-order A-phi solver : a degenerate element\n");
		ierr = 1;
	}

	// Dirichlet : awall の三角形の辺 (Whitney と勾配) と面 (2 個)、電極の頂点 (φ = V) と電極上の辺の泡 (0)、
	//             非導電領域の φ、使わない欄 (中間節点の φ)
	unsigned char *fix = (unsigned char *)calloc((size_t)n, sizeof(unsigned char));
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
		const int f = p2_face_id(fvert, nf, sv);
		if (f >= 0) {
			fix[oF + (2 * f)] = fix[oF + (2 * f) + 1] = 1;
			nwall++;
		}
	}
	unsigned char *cond = (unsigned char *)calloc((size_t)nn, sizeof(unsigned char));
	unsigned char *econd = (unsigned char *)calloc((size_t)ne, sizeof(unsigned char));
	unsigned char *vert = (unsigned char *)calloc((size_t)nn, sizeof(unsigned char));
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
	if (nwall == 0) {
		fprintf(fp_log, "*** warning : no awall surface; the outer boundary is a magnetic wall (dA/dn = 0)\n");
	}
	// 固定の行は恒等行 (crs_spmv の fix) なので、未使用の行の対角は 0 のままでよい

	double *xdr = (double *)calloc((size_t)n, sizeof(double));
	double *xdi = (double *)calloc((size_t)n, sizeof(double));
	double *br = (double *)malloc((size_t)n * sizeof(double));
	double *bi = (double *)malloc((size_t)n * sizeof(double));
	double *ur = (double *)malloc((size_t)n * sizeof(double));
	double *ui = (double *)malloc((size_t)n * sizeof(double));
	double *w1 = (double *)malloc((size_t)n * sizeof(double));
	double *w2 = (double *)malloc((size_t)n * sizeof(double));
	double *yr = (double *)calloc((size_t)np * np, sizeof(double));
	double *yi = (double *)calloc((size_t)np * np, sizeof(double));
	fprintf(fp_log, "  %-10s %8s %13s\n", "excite", "iter", "residual");
	for (int jc = 1; (jc <= np) && !ierr; jc++) {
		for (int i = 0; i < n; i++) xdr[i] = xdi[i] = 0;
		for (int i = 0; i < nn; i++) {
			if (fix[oV + i] && cond[i] && vert[i] && (NodeConductor[i] == jc)) xdr[oV + i] = Volt;
		}
		crs_spmv(&Kr, xdr, w1, NULL);
		crs_spmv(&Ki, xdr, w2, NULL);
		for (int i = 0; i < n; i++) {
			br[i] = (fix[i] ? 0 : -w1[i]);
			bi[i] = (fix[i] ? 0 : -w2[i]);
		}
		char label[64];
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
		// 電極の頂点の反作用 (節点の式は 1/(jω) 倍してあるので jω 倍して戻す)。電極の上で 1 になる
		// 試験関数 Σ λv は泡の係数 0 の階層的な 2 次の関数なので、頂点の行だけを足せばよい
		double qr[MAXPORT], qi[MAXPORT];
		for (int p = 0; p < MAXPORT; p++) qr[p] = qi[p] = 0;
		for (int i = 0; i < nn; i++) {
			const int id = NodeConductor[i];
			if ((id < 0) || !cond[i] || !vert[i]) continue;
			double rr, ri;
			p2_row_dot_c(&Kr, &Ki, (int32_t)(oV + i), ur, ui, &rr, &ri);
			qr[id] += -omega * ri;
			qi[id] += omega * rr;
		}
		for (int kc = 1; kc <= np; kc++) {
			yr[((kc - 1) * np) + (jc - 1)] = qr[kc] / Volt;
			yi[((kc - 1) * np) + (jc - 1)] = qi[kc] / Volt;
		}
		fprintf(fp_log, "  port %d : I = %.6e %+.6e j [A] (V = %.6e [V])\n", jc, qr[jc], qi[jc], Volt);
	}

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
	free(xdr); free(xdi); free(br); free(bi); free(ur); free(ui); free(w1); free(w2);
	free(yr); free(yi); free(fix); free(cond); free(econd); free(vert);
	free(fe); free(fvert);
	crs_free(&Kr);
	crs_free(&Ki);
	edge_free();
	return ierr;
}
