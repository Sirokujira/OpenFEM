/*
eddy3d_pr.c

3 次元渦電流 (analysis = A、A-φ 連成) と辺要素の自己検証 (analysis = E) の **3 次以上** の要素。
次数 r (= 格子の次数 TetOrder、3〜7) の四面体 (等パラメトリック、曲がってよい) の格子で、A を次数 r の
Nedelec 要素 (第 1 種、P_r⁻Λ¹)、φ を次数 r の節点要素 (格子の節点そのもの) で離散化する。

基底は Arnold-Falk-Winther の幾何学的分解 (Geometric decompositions and local bases for spaces of finite
element differential forms, 2009) の P_r⁻Λ¹ の基底:

	λ^α φ_σ,   σ = (a, b) (a < b)、|α| = r − 1、α_i = 0 (i < a)
	φ_σ = λa∇λb − λb∇λa (Whitney)、∇×(λ^α φ_σ) = ∇λ^α × φ_σ + 2 λ^α ∇λa × ∇λb

ここで頂点の「<」は**全体の節点番号の順**。関数は σ ∪ supp α の頂点の組 (辺・面・内部) に属し、個数は
辺 r、面 r(r − 1)、内部 r(r − 1)(r − 2)/2 (r = 3 で 6·3 + 4·6 + 3 = 45)。面 f に属さない関数の接線成分は
f の上で 0 になり (λ_i (i ∉ f) が掛かるか、φ_σ の σ が f の外の頂点を含む)、f に属する関数は f の頂点の
λ と全体の順だけで決まるので、隣の要素と同じ関数になる (接線連続)。向きの符号は要らない (σ も α も全体の
順で決めるので、局所の並びに依らない)。

λ は参照四面体の重心座標を物理空間に引き戻したもの (λ_i = Σ_a λ_i(ξ_a) N_a、節点基底は 1 次の関数を
厳密に補間する)。**共変 Piola 変換は「引き戻した λ の物理空間の勾配を使う」ことと同じ**なので、曲がった
要素でも上の回転の式が厳密に成り立つ。φ は次数 r の等パラメトリックな節点基底で、∇φ (P_{r−1}³ の引き戻し)
は A の空間に厳密に入るので、(A, φ) -> (A + ∇ψ, φ − jωψ) の不変性が保たれ、右辺は値域に入る (2 次と同じく
ゲージは固定しない)。積分は Duffy + 各方向 r + 2 点 (elem3d_ip、質量の 2r 次が直線要素で厳密)。

未知数の並び : [辺 r x NEdge | 面 r(r−1) x NFace | 内部 r(r−1)(r−2)/2 x NTet | φ (NNode)]。

対応するのは周波数領域の線形材料・ゲージなし・反復解法・場の出力・交流の発熱 (thermal + A) (それ以外は setup で弾く)。
*/

#include "fem.h"
#include "fem_prototype.h"

#define PR_RMAX (7)
#define PR_TMAX (128)		// 1 つの実体の関数の数の上限 (内部 r(r−1)(r−2)/2 = 105 @ r = 7)
#define PR_NAMAX (6 * 7 + 4 * 42 + 105)		// 要素の A の自由度 (r = 7 で 315)
#define PR_NPMAX (TET_NMAX)		// 要素の φ の自由度 (r = 7 で 120)
#define PR_NDMAX (PR_NAMAX + PR_NPMAX)

static const int PR_EDGE[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};	// edge.c と同じ
static const int PR_FACE[4][3] = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};		// p2_faces と同じ

// 実体 (頂点 m 個、全体の順の位置 0..m−1) の上の関数 λ^α φ_(a,b)
typedef struct {
	int a, b;
	int al[4];
} prtpl_t;

static prtpl_t TplE[PR_TMAX], TplF[PR_TMAX], TplI[PR_TMAX];
static int NtE, NtF, NtI, TplR;


// 頂点 m 個の実体に属する関数を列挙する (σ ∪ supp α = 全頂点)
static int pr_tpl(int r, int m, prtpl_t *t)
{
	int cnt = 0;
	for (int a = 0; a < m; a++) {
		for (int b = a + 1; b < m; b++) {
			for (int a0 = 0; a0 <= (r - 1); a0++) {
				for (int a1 = 0; a1 <= (r - 1 - a0); a1++) {
					for (int a2 = 0; a2 <= (r - 1 - a0 - a1); a2++) {
						const int al[4] = {a0, a1, a2, r - 1 - a0 - a1 - a2};
						int ok = 1;
						for (int i = 0; i < 4; i++) {
							if ((i >= m) && al[i]) ok = 0;
							if ((i < a) && al[i]) ok = 0;
							if ((i < m) && (i != a) && (i != b) && !al[i]) ok = 0;
						}
						if (!ok) continue;
						if (cnt >= PR_TMAX) return -1;
						t[cnt].a = a;
						t[cnt].b = b;
						for (int i = 0; i < 4; i++) t[cnt].al[i] = al[i];
						cnt++;
					}
				}
			}
		}
	}
	return cnt;
}


static int pr_init(int r)
{
	if ((r < 2) || (r > PR_RMAX)) return 1;
	if (TplR == r) return 0;
	NtE = pr_tpl(r, 2, TplE);
	NtF = pr_tpl(r, 3, TplF);
	NtI = pr_tpl(r, 4, TplI);
	// 個数は AFW の次元の式と一致しなければならない (列挙の条件の誤りの検出)
	if ((NtE != r) || (NtF != (r * (r - 1))) || (NtI != ((r * (r - 1) * (r - 2)) / 2))) return 1;
	TplR = r;
	return 0;
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

static void nu_mul(const double nu[6], const double c[3], double r[3])
{
	r[0] = (nu[0] * c[0]) + (nu[3] * c[1]) + (nu[5] * c[2]);
	r[1] = (nu[3] * c[0]) + (nu[1] * c[1]) + (nu[4] * c[2]);
	r[2] = (nu[5] * c[0]) + (nu[4] * c[1]) + (nu[2] * c[2]);
}


// 局所の頂点の添字 v[0..m-1] を全体の節点番号の昇順に並べる
static void pr_sort(const int32_t *nd, int *v, int m)
{
	for (int x = 0; x < m - 1; x++) {
		for (int y = x + 1; y < m; y++) {
			if (nd[v[y]] < nd[v[x]]) { const int t = v[x]; v[x] = v[y]; v[y] = t; }
		}
	}
}


// 解く系の一式
static struct {
	int r, ne, nf, nn, nt, n, oF, oI, oV, nA, nP;
	int32_t *fe, *fvert;
	crs_t Kr, Ki;
	unsigned char *fix, *cond, *used;
	double *xdr, *xdi, *br, *bi, *ur, *ui, *w1, *w2;
	double omega;
	double *nuq;		// 非線形 (B-H) の積分点ごとの ν [nt * nq] (NULL なら材料の ν)
	int nq;
} PR;

#define PR_AAM (5)		// 非線形の Anderson 加速の深さ (2 次の P2_AAM と同じ)


// 要素の自由度の全体番号 (A : 辺・面・内部、φ : 要素の節点)。戻り値は個数
static int pr_dofs(int e, int32_t dof[PR_NDMAX])
{
	const int32_t *ed = &TetEdge[EdgeOff[e]];
	int k = 0;
	for (int l = 0; l < 6; l++) {
		for (int t = 0; t < NtE; t++) dof[k++] = (ed[l] * NtE) + t;
	}
	for (int f = 0; f < 4; f++) {
		for (int t = 0; t < NtF; t++) dof[k++] = PR.oF + (PR.fe[(e * 4) + f] * NtF) + t;
	}
	for (int t = 0; t < NtI; t++) dof[k++] = PR.oI + (e * NtI) + t;
	int32_t nd[TET_NMAX];
	const int nen = tet_nodes(e, nd);
	for (int a = 0; a < nen; a++) dof[k++] = PR.oV + nd[a];
	return k;
}


// λ^α φ_(a,b) の値と回転 (a, b, α は局所の頂点)
static void pr_fun(const double lam[4], const double gl[4][3], const double pw[4][PR_RMAX + 1],
	int a, int b, const int al[4], double w[3], double c[3])
{
	double f = 1, gf[3] = {0, 0, 0};
	for (int i = 0; i < 4; i++) f *= pw[i][al[i]];
	for (int i = 0; i < 4; i++) {
		if (al[i] == 0) continue;
		double d = al[i] * pw[i][al[i] - 1];
		for (int j = 0; j < 4; j++) {
			if (j != i) d *= pw[j][al[j]];
		}
		for (int q = 0; q < 3; q++) gf[q] += d * gl[i][q];
	}
	double ph[3], cr[3], gp[3];
	for (int q = 0; q < 3; q++) ph[q] = (lam[a] * gl[b][q]) - (lam[b] * gl[a][q]);
	cross3(gl[a], gl[b], cr);
	cross3(gf, ph, gp);
	for (int q = 0; q < 3; q++) {
		w[q] = f * ph[q];
		c[q] = gp[q] + (2 * f * cr[q]);
	}
}


// 形状関数の値 nv と物理空間の勾配 g から、A の基底 (値 wv と回転 wc、pr_dofs の A の並び) を作る
static void pr_basis_core(int e, const double nv[ELEM_NMAX], const double g[ELEM_NMAX][3],
	double (*wv)[3], double (*wc)[3])
{
	const int r = PR.r;
	int32_t nd[TET_NMAX];
	const int nen = tet_nodes(e, nd);
	// 引き戻した重心座標とその勾配 (λ_i = Σ_a (α_a,i / r) N_a)
	double lam[4] = {0, 0, 0, 0}, gl[4][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
	for (int a = 0; a < nen; a++) {
		const int *ix = tet_ho_index(a);
		for (int i = 0; i < 4; i++) {
			if (ix[i] == 0) continue;
			const double s = (double)ix[i] / r;
			lam[i] += s * nv[a];
			for (int q = 0; q < 3; q++) gl[i][q] += s * g[a][q];
		}
	}
	double pw[4][PR_RMAX + 1];
	for (int i = 0; i < 4; i++) {
		pw[i][0] = 1;
		for (int k = 1; k <= r; k++) pw[i][k] = pw[i][k - 1] * lam[i];
	}
	int k = 0;
	for (int l = 0; l < 6; l++) {
		int v[2] = {PR_EDGE[l][0], PR_EDGE[l][1]};
		pr_sort(nd, v, 2);
		for (int t = 0; t < NtE; t++) {
			int al[4] = {0, 0, 0, 0};
			for (int i = 0; i < 2; i++) al[v[i]] = TplE[t].al[i];
			pr_fun(lam, gl, pw, v[TplE[t].a], v[TplE[t].b], al, wv[k], wc[k]);
			k++;
		}
	}
	for (int f = 0; f < 4; f++) {
		int v[3] = {PR_FACE[f][0], PR_FACE[f][1], PR_FACE[f][2]};
		pr_sort(nd, v, 3);
		for (int t = 0; t < NtF; t++) {
			int al[4] = {0, 0, 0, 0};
			for (int i = 0; i < 3; i++) al[v[i]] = TplF[t].al[i];
			pr_fun(lam, gl, pw, v[TplF[t].a], v[TplF[t].b], al, wv[k], wc[k]);
			k++;
		}
	}
	{
		int v[4] = {0, 1, 2, 3};
		pr_sort(nd, v, 4);
		for (int t = 0; t < NtI; t++) {
			int al[4] = {0, 0, 0, 0};
			for (int i = 0; i < 4; i++) al[v[i]] = TplI[t].al[i];
			pr_fun(lam, gl, pw, v[TplI[t].a], v[TplI[t].b], al, wv[k], wc[k]);
			k++;
		}
	}
}


static int64_t pr_find(const crs_t *A, int32_t row, int32_t col)
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


static int cmp_i32(const void *a, const void *b)
{
	const int32_t x = *(const int32_t *)a, y = *(const int32_t *)b;
	return ((x < y) ? -1 : ((x > y) ? 1 : 0));
}


// 要素の自由度の組から CRS のパターンを作る (未使用の行にも対角を置く)
static void pr_crs(crs_t *A)
{
	const int n = PR.n, nt = PR.nt;
	int *cnt = (int *)calloc((size_t)n, sizeof(int));
	int32_t *dof = (int32_t *)malloc(PR_NDMAX * sizeof(int32_t));
	for (int e = 0; e < nt; e++) {
		const int nd = pr_dofs(e, dof);
		for (int a = 0; a < nd; a++) cnt[dof[a]] += nd;
	}
	for (int i = 0; i < n; i++) cnt[i] += 1;
	int64_t *ptr = (int64_t *)malloc(((size_t)n + 1) * sizeof(int64_t));
	ptr[0] = 0;
	for (int i = 0; i < n; i++) ptr[i + 1] = ptr[i] + cnt[i];
	int32_t *lst = (int32_t *)malloc((size_t)ptr[n] * sizeof(int32_t));
	memset(cnt, 0, (size_t)n * sizeof(int));
	for (int i = 0; i < n; i++) lst[ptr[i] + (cnt[i]++)] = i;
	for (int e = 0; e < nt; e++) {
		const int nd = pr_dofs(e, dof);
		for (int a = 0; a < nd; a++) {
			for (int b = 0; b < nd; b++) lst[ptr[dof[a]] + (cnt[dof[a]]++)] = dof[b];
		}
	}
	A->n = n;
	A->rowptr = (int64_t *)malloc(((size_t)n + 1) * sizeof(int64_t));
	A->rowptr[0] = 0;
	int64_t nnz = 0;
	for (int i = 0; i < n; i++) {
		int32_t *rw = &lst[ptr[i]];
		const int m = cnt[i];
		qsort(rw, (size_t)m, sizeof(int32_t), cmp_i32);
		int u = 0;
		for (int q = 0; q < m; q++) {
			if ((q == 0) || (rw[q] != rw[q - 1])) u++;
		}
		nnz += u;
		A->rowptr[i + 1] = nnz;
	}
	A->nnz = nnz;
	A->col = (int32_t *)malloc((size_t)nnz * sizeof(int32_t));
	A->val = (double *)calloc((size_t)nnz, sizeof(double));
	for (int i = 0; i < n; i++) {
		const int32_t *rw = &lst[ptr[i]];
		int64_t p = A->rowptr[i];
		for (int q = 0; q < cnt[i]; q++) {
			if ((q == 0) || (rw[q] != rw[q - 1])) A->col[p++] = rw[q];
		}
	}
	free(cnt); free(ptr); free(lst); free(dof);
}


static void pr_close(void)
{
	free(PR.xdr); free(PR.xdi); free(PR.br); free(PR.bi); free(PR.ur); free(PR.ui); free(PR.w1); free(PR.w2);
	free(PR.fix); free(PR.cond); free(PR.used);
	free(PR.fe); free(PR.fvert);
	crs_free(&PR.Kr);
	crs_free(&PR.Ki);
	edge_free();
	memset(&PR, 0, sizeof(PR));
}


// 実体と自由度の番号付け (A と E で共通)
static int pr_number(FILE *fp_log, const char *what)
{
	if (!MeshMode || (TetOrder < 2) || (NHex + NPrism + NPyr > 0) || pr_init(TetOrder)) {
		fprintf(fp_log, "*** the %s needs a pure tetrahedral mesh of order 2 to %d\n", what, PR_RMAX);
		return 1;
	}
	memset(&PR, 0, sizeof(PR));
	edge_build();
	PR.r = TetOrder;
	PR.ne = NEdge; PR.nn = NNode; PR.nt = NTet;
	PR.fe = p2_faces(&PR.nf, &PR.fvert);
	PR.oF = NtE * PR.ne;
	PR.oI = PR.oF + (NtF * PR.nf);
	PR.oV = PR.oI + (NtI * PR.nt);
	PR.n = PR.oV + PR.nn;
	PR.nA = (6 * NtE) + (4 * NtF) + NtI;
	PR.nP = tet_nen();
	return 0;
}


// 要素の行列 (A-A の回転 se・質量 te、A-φ の結合 ce、φ-φ の pe)。sig は要素の σ
// nuq が NULL でなければ積分点 q の ν を等方性の nuq[q] にする (B-H の実効磁気抵抗率法)
static int pr_element(int e, const double nu0[6], const double *nuq, double sig, double *se, double *te, double *ce,
	double *pe, double (*wv)[3], double (*wc)[3], double (*nc)[3])
{
	const int nA = PR.nA, nP = PR.nP;
	memset(se, 0, (size_t)nA * nA * sizeof(double));
	memset(te, 0, (size_t)nA * nA * sizeof(double));
	memset(ce, 0, (size_t)nA * nP * sizeof(double));
	memset(pe, 0, (size_t)nP * nP * sizeof(double));
	const int nq = elem3d_nip(e);
	for (int q = 0; q < nq; q++) {
		double nv[ELEM_NMAX], g[ELEM_NMAX][3], w;
		if (elem3d_ip(e, q, nv, g, &w)) return 1;
		pr_basis_core(e, nv, g, wv, wc);
		double nu[6];
		for (int c = 0; c < 6; c++) nu[c] = nu0[c];
		if (nuq != NULL) {
			nu[0] = nu[1] = nu[2] = nuq[q];
			nu[3] = nu[4] = nu[5] = 0;
		}
		for (int k = 0; k < nA; k++) nu_mul(nu, wc[k], nc[k]);
		for (int k = 0; k < nA; k++) {
			for (int l = k; l < nA; l++) {
				se[(k * nA) + l] += w * dot3(nc[k], wc[l]);
				te[(k * nA) + l] += w * sig * dot3(wv[k], wv[l]);
			}
			for (int a = 0; a < nP; a++) ce[(k * nP) + a] += w * sig * dot3(wv[k], g[a]);
		}
		for (int a = 0; a < nP; a++) {
			for (int b = a; b < nP; b++) pe[(a * nP) + b] += w * sig * dot3(g[a], g[b]);
		}
	}
	for (int k = 0; k < nA; k++) {
		for (int l = 0; l < k; l++) {
			se[(k * nA) + l] = se[(l * nA) + k];
			te[(k * nA) + l] = te[(l * nA) + k];
		}
	}
	for (int a = 0; a < nP; a++) {
		for (int b = 0; b < a; b++) pe[(a * nP) + b] = pe[(b * nP) + a];
	}
	return 0;
}


static int pr_open(FILE *fp_log)
{
	if ((NPort < 1) || (Freq <= 0)) {
		fprintf(fp_log, "*** the high-order A-phi solver needs a port and a frequency\n");
		return 1;
	}
	if (pr_number(fp_log, "high-order A-phi solver")) return 1;
	const int n = PR.n, ne = PR.ne, nf = PR.nf, nn = PR.nn, nt = PR.nt;
	PR.omega = 2 * PI * Freq;
	fprintf(fp_log, "  order %d : nodes = %d, tetrahedra = %d, edges = %d, faces = %d, unknowns = %d "
		"(A : %d x %d edges + %d x %d faces + %d x %d interiors = %d per element, phi : %d nodes)\n",
		PR.r, nn, nt, ne, nf, n, NtE, ne, NtF, nf, NtI, nt, PR.nA, nn);
	pr_crs(&PR.Kr);
	pr_crs(&PR.Ki);
	// mpisolve = 2 : 実部・虚部とも自分の行ブロックだけを持つ
	crs_distribute_pair(fp_log, &PR.Kr, &PR.Ki);
	fprintf(fp_log, "  coupled matrix : %lld nonzeros (%.1f per row)\n", (long long)PR.Kr.nnz, (double)PR.Kr.nnz / n);

	// Dirichlet : awall の三角形の辺と面の A (接線成分 0)、電極の節点の φ (= V)、非導電領域と未使用の φ
	unsigned char *fix = PR.fix = (unsigned char *)calloc((size_t)n, sizeof(unsigned char));
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
			if (id < 0) continue;
			for (int s = 0; s < NtE; s++) fix[(id * NtE) + s] = 1;
		}
		int32_t sv[3] = {v[0], v[1], v[2]};
		for (int x = 0; x < 2; x++) {
			for (int y = x + 1; y < 3; y++) {
				if (sv[y] < sv[x]) { const int32_t tt = sv[x]; sv[x] = sv[y]; sv[y] = tt; }
			}
		}
		const int f = p2_face_id(PR.fvert, nf, sv);
		if (f >= 0) {
			for (int s = 0; s < NtF; s++) fix[PR.oF + (f * NtF) + s] = 1;
			nwall++;
		}
	}
	unsigned char *cond = PR.cond = (unsigned char *)calloc((size_t)nn, sizeof(unsigned char));
	unsigned char *used = PR.used = (unsigned char *)calloc((size_t)nn, sizeof(unsigned char));
	for (int e = 0; e < nt; e++) {
		int32_t nd[TET_NMAX];
		const int nen = tet_nodes(e, nd);
		const int c = (Material[elem3d_mat(e)].sigma > 0);
		for (int a = 0; a < nen; a++) {
			used[nd[a]] = 1;
			if (c) cond[nd[a]] = 1;
		}
	}
	int nelec = 0;
	for (int i = 0; i < nn; i++) {
		if (!used[i] || !cond[i]) fix[PR.oV + i] = 1;
		else if (NodeConductor[i] >= 0) {
			fix[PR.oV + i] = 1;
			nelec++;
		}
	}
	fprintf(fp_log, "  fixed : %d wall faces (A_t = 0), %d electrode nodes\n", nwall, nelec);
	if (nwall == 0) {
		fprintf(fp_log, "*** warning : no awall surface; the outer boundary is a magnetic wall (dA/dn = 0)\n");
	}
	PR.xdr = (double *)calloc((size_t)n, sizeof(double));
	PR.xdi = (double *)calloc((size_t)n, sizeof(double));
	PR.br = (double *)malloc((size_t)n * sizeof(double));
	PR.bi = (double *)malloc((size_t)n * sizeof(double));
	PR.ur = (double *)malloc((size_t)n * sizeof(double));
	PR.ui = (double *)malloc((size_t)n * sizeof(double));
	PR.w1 = (double *)malloc((size_t)n * sizeof(double));
	PR.w2 = (double *)malloc((size_t)n * sizeof(double));
	return 0;
}


// 組み立て : [S + jωT, TG; GᵀT, GᵀTG/(jω)] (節点の式を jω で割った複素対称形、1 次・2 次と同じ)
static int pr_assemble(FILE *fp_log, const double *sigelem)
{
	const int nt = PR.nt, nA = PR.nA, nP = PR.nP;
	const double omega = PR.omega;
	crs_zero(&PR.Kr);
	crs_zero(&PR.Ki);
	double *se = (double *)malloc((size_t)nA * nA * sizeof(double));
	double *te = (double *)malloc((size_t)nA * nA * sizeof(double));
	double *ce = (double *)malloc((size_t)nA * nP * sizeof(double));
	double *pe = (double *)malloc((size_t)nP * nP * sizeof(double));
	double (*wv)[3] = (double (*)[3])malloc((size_t)nA * 3 * sizeof(double));
	double (*wc)[3] = (double (*)[3])malloc((size_t)nA * 3 * sizeof(double));
	double (*nc)[3] = (double (*)[3])malloc((size_t)nA * 3 * sizeof(double));
	int32_t *dof = (int32_t *)malloc(PR_NDMAX * sizeof(int32_t));
	int bad = 0;
	for (int e = 0; (e < nt) && !bad; e++) {
		const int m = elem3d_mat(e);
		double nu[6];
		material_coef_pub(m, 3, nu);
		const double sig = ((sigelem != NULL) ? sigelem[e] : Material[m].sigma);		// 熱連成なら σ(T)
		pr_dofs(e, dof);
		const double *nuq = (((PR.nuq != NULL) && (Material[m].nbh[0] > 0)) ? &PR.nuq[(int64_t)e * PR.nq] : NULL);
		if (pr_element(e, nu, nuq, sig, se, te, ce, pe, wv, wc, nc)) { bad = 1; break; }
		for (int k = 0; k < nA; k++) {
			for (int l = 0; l < nA; l++) {
				// 行列を分散したとき (mpisolve = 2) は他の rank の行が空で -1 になるので飛ばす
				const int64_t p = pr_find(&PR.Kr, dof[k], dof[l]);
				if (p < 0) continue;
				PR.Kr.val[p] += se[(k * nA) + l];
				PR.Ki.val[p] += omega * te[(k * nA) + l];
			}
			for (int a = 0; a < nP; a++) {
				const double v = ce[(k * nP) + a];
				const int64_t p1 = pr_find(&PR.Kr, dof[k], dof[nA + a]);
				const int64_t p2 = pr_find(&PR.Kr, dof[nA + a], dof[k]);
				if (p1 >= 0) PR.Kr.val[p1] += v;
				if (p2 >= 0) PR.Kr.val[p2] += v;
			}
		}
		for (int a = 0; a < nP; a++) {
			for (int b = 0; b < nP; b++) {
				const int64_t p = pr_find(&PR.Ki, dof[nA + a], dof[nA + b]);
				if (p >= 0) PR.Ki.val[p] -= pe[(a * nP) + b] / omega;
			}
		}
	}
	free(se); free(te); free(ce); free(pe); free(wv); free(wc); free(nc); free(dof);
	if (bad) {
		fprintf(fp_log, "*** the high-order A-phi solver : a degenerate element\n");
		return 1;
	}
	return 0;
}


static void pr_row_dot_c(int32_t row, double *sr, double *si)
{
	double ar = 0, ai = 0;
	for (int64_t p = PR.Kr.rowptr[row]; p < PR.Kr.rowptr[row + 1]; p++) {
		const int32_t c = PR.Kr.col[p];
		ar += (PR.Kr.val[p] * PR.ur[c]) - (PR.Ki.val[p] * PR.ui[c]);
		ai += (PR.Kr.val[p] * PR.ui[c]) + (PR.Ki.val[p] * PR.ur[c]);
	}
	// 行列を分散したとき (mpisolve = 2) は行の持ち主が計算した値を配る (全 rank が同じ順に呼ぶ)
	int64_t d0, d1;
	if (crs_dist_info(&PR.Kr, &d0, &d1)) {
		double v[2] = {ar, ai};
		mpi_bcast_dbl(v, 2, mpi_row_owner((int)PR.Kr.n, (int)row));
		ar = v[0];
		ai = v[1];
	}
	*sr = ar;
	*si = ai;
}


// ポート jc を V = Volt で励振して解き (PR.ur / ui に解)、各導体へ流れ込む電流 qr + j qi を返す
static int pr_solve_port(FILE *fp_log, int jc, FILE *fp_iter, const char *label, double qr[MAXPORT], double qi[MAXPORT])
{
	const int n = PR.n, nn = PR.nn, oV = PR.oV;
	const double omega = PR.omega;
	for (int i = 0; i < n; i++) PR.xdr[i] = PR.xdi[i] = 0;
	for (int i = 0; i < nn; i++) {
		if (PR.fix[oV + i] && PR.cond[i] && PR.used[i] && (NodeConductor[i] == jc)) PR.xdr[oV + i] = Volt;
	}
	crs_spmv(&PR.Kr, PR.xdr, PR.w1, NULL);
	crs_spmv(&PR.Ki, PR.xdr, PR.w2, NULL);
	for (int i = 0; i < n; i++) {
		PR.br[i] = (PR.fix[i] ? 0 : -PR.w1[i]);
		PR.bi[i] = (PR.fix[i] ? 0 : -PR.w2[i]);
	}
	const int iter = solver_cocg(&PR.Kr, &PR.Ki, 1.0, PR.br, PR.bi, PR.ur, PR.ui, PR.fix,
		Solver.maxiter, Solver.nout, Solver.converg, fp_iter, label);
	if (iter < 0) {
		fprintf(fp_log, "*** the high-order A-phi solve did not converge (%s)\n", label);
		return 1;
	}
	for (int i = 0; i < n; i++) {
		PR.ur[i] += PR.xdr[i];
		PR.ui[i] += PR.xdi[i];
	}
	// 電極の節点の反作用 (節点の式は 1/(jω) 倍してあるので jω 倍して戻す)。試験関数 Σ_(電極の節点) N_a は
	// 電極の面で 1 になる (面の節点がすべて電極の節点で、他の節点の関数は面で 0)
	for (int p = 0; p < MAXPORT; p++) qr[p] = qi[p] = 0;
	for (int i = 0; i < nn; i++) {
		const int id = NodeConductor[i];
		if ((id < 0) || !PR.cond[i] || !PR.used[i]) continue;
		double rr, ri;
		pr_row_dot_c((int32_t)(oV + i), &rr, &ri);
		qr[id] += -omega * ri;
		qi[id] += omega * rr;
	}
	return 0;
}


// 場の出力 (ポート 1 の励振) : 要素の重心での B = ∇×A (実部) と J = −σ (jωA + ∇φ) (実部・虚部)
static void pr_fields(void)
{
	const int nt = PR.nt, nA = PR.nA, nP = PR.nP;
	const double omega = PR.omega;
	const double lc[4] = {0.25, 0.25, 0.25, 0.25};
	double *bv = (double *)calloc((size_t)nt * 3, sizeof(double));
	double *jr = (double *)calloc((size_t)nt * 3, sizeof(double));
	double *ji = (double *)calloc((size_t)nt * 3, sizeof(double));
	double (*wv)[3] = (double (*)[3])malloc((size_t)nA * 3 * sizeof(double));
	double (*wc)[3] = (double (*)[3])malloc((size_t)nA * 3 * sizeof(double));
	int32_t *dof = (int32_t *)malloc(PR_NDMAX * sizeof(int32_t));
	for (int e = 0; e < nt; e++) {
		double nv[ELEM_NMAX], g[ELEM_NMAX][3], det;
		if (tet_ho_point(e, lc, nv, g, &det)) continue;
		pr_dofs(e, dof);
		pr_basis_core(e, nv, g, wv, wc);
		double br[3] = {0, 0, 0}, er[3] = {0, 0, 0}, ei[3] = {0, 0, 0};
		for (int k = 0; k < nA; k++) {
			const double ar = PR.ur[dof[k]], ai = PR.ui[dof[k]];
			for (int c = 0; c < 3; c++) {
				br[c] += ar * wc[k][c];
				er[c] += -omega * ai * wv[k][c];
				ei[c] += omega * ar * wv[k][c];
			}
		}
		for (int a = 0; a < nP; a++) {
			for (int c = 0; c < 3; c++) {
				er[c] += PR.ur[dof[nA + a]] * g[a][c];
				ei[c] += PR.ui[dof[nA + a]] * g[a][c];
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
	free(bv); free(jr); free(ji); free(wv); free(wc); free(dof);
}


int solve_eddy3d_pr(FILE *fp_log)
{
	const int np = NPort;
	int ierr = 0;
	fprintf(fp_log, "\n=== 3D eddy current (A-phi, order-%d Nedelec / nodal elements) ===\n", TetOrder);
	if (pr_open(fp_log)) return 1;
	const double omega = PR.omega;
	const int nt = PR.nt;

	// 非線形 (B-H) の材料があれば実効磁気抵抗率法 (2 次の eddy3d_p2.c と同じ): ν を波高値 |B| で
	// **積分点ごとに**評価し、H を独立変数にした更新 ν = |H|/B(|H|) に ln ν の Anderson 加速を掛ける。
	// 単一ポート・電流振幅 current (入力段で検査済み)
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
		aady = (double *)malloc((size_t)nqt * PR_AAM * sizeof(double));
		aadf = (double *)malloc((size_t)nqt * PR_AAM * sizeof(double));
		for (int e = 0; e < nt; e++) {
			const material_t *mt = &Material[elem3d_mat(e)];
			for (int q = 0; q < nq; q++) nuq[((int64_t)e * nq) + q] = ((mt->nbh[0] > 0) ? bh_nu(mt, 0, 0) : 0);
		}
		PR.nuq = nuq;
		PR.nq = nq;
		fprintf(fp_log, "  nonlinear (B-H) : effective reluctivity at %d points per element, current amplitude "
			"%.6e [A], maxiter = %d, tol = %.1e, relax = %.2f\n", nq, Curr, NlMaxiter, NlTol, NlRelax);
	}
	ierr = pr_assemble(fp_log, NULL);
	double *yr = (double *)calloc((size_t)np * np, sizeof(double));
	double *yi = (double *)calloc((size_t)np * np, sizeof(double));
	double (*wv)[3] = (double (*)[3])malloc((size_t)PR.nA * 3 * sizeof(double));
	double (*wc)[3] = (double (*)[3])malloc((size_t)PR.nA * 3 * sizeof(double));
	int32_t *dof = (int32_t *)malloc(PR_NDMAX * sizeof(int32_t));
	int final = !nlf;
	for (int nlit = 1; !ierr; nlit++) {
		if (nlit > 1) ierr = pr_assemble(fp_log, NULL);
		if (final && !ierr) fprintf(fp_log, "  %-10s %8s %13s\n", "excite", "iter", "residual");
		for (int jc = 1; (jc <= np) && !ierr; jc++) {
			char label[64];
			sprintf(label, "port%d", jc);
			double qr[MAXPORT], qi[MAXPORT];
			if (pr_solve_port(fp_log, jc, (final ? fp_log : NULL), label, qr, qi)) {
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
				for (int i = 0; i < PR.n; i++) {
					const double xr = PR.ur[i], xi = PR.ui[i];
					PR.ur[i] = (svr * xr) - (svi * xi);
					PR.ui[i] = (svr * xi) + (svi * xr);
				}
			}
			if (FieldOut && (jc == 1) && final) pr_fields();
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
			pr_dofs(e, dof);
			for (int q = 0; q < nq; q++) {
				double nv[ELEM_NMAX], g[ELEM_NMAX][3], w;
				double bre[3] = {0, 0, 0}, bim[3] = {0, 0, 0};
				if (!elem3d_ip(e, q, nv, g, &w)) {
					pr_basis_core(e, nv, g, wv, wc);
					for (int k = 0; k < PR.nA; k++) {
						const double ar = PR.ur[dof[k]], ai = PR.ui[dof[k]];
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
		anderson_step(nqt, aay, aaf, aayp, aafp, aady, aadf, PR_AAM, &aan, NlRelax);
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
	free(wv); free(wc); free(dof);
	PR.nuq = NULL;
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
	pr_close();
	return ierr;
}


/*
交流の発熱 (thermal + A、3 次以上の四面体)。2 次 (eddy3d_p2_heat_*) と同じ形 : 要素の σ(T) で組み立て直して解き、
端子電流が cur (実数の波高値) になるよう解を倍率し、損失 ½σ|e|² (e = jωA + ∇φ) を**組み立てと同じ積分点で**
節点の形状関数に掛けて荷重 f にする。Σ_a N_a = 1 なので Σf = ∫½σ|e|² = ½|I|²R が厳密 (連成系の行の恒等式)。
*/
int eddy3d_pr_heat_open(FILE *fp_log)
{
	fprintf(fp_log, "  AC heating (analysis = A, order-%d Nedelec) : frequency = %.6e [Hz]\n", TetOrder, Freq);
	return pr_open(fp_log);
}


int eddy3d_pr_heat_solve(FILE *fp_log, const double *sigelem, double cur, double *f,
	double *zr, double *zi, int it)
{
	const int nt = PR.nt, nn = PR.nn, n = PR.n, nA = PR.nA, nP = PR.nP;
	const double omega = PR.omega;
	if (pr_assemble(fp_log, sigelem)) return 1;
	char label[64];
	sprintf(label, "A%d", it);
	double qr[MAXPORT], qi[MAXPORT];
	if (pr_solve_port(fp_log, 1, NULL, label, qr, qi)) return 1;
	const double ir = qr[1], ii = qi[1];
	const double den = (ir * ir) + (ii * ii);
	if (den <= 0) {
		fprintf(fp_log, "*** thermal : no current flows through the port\n");
		return 1;
	}
	*zr = Volt * ir / den;
	*zi = -Volt * ii / den;
	const double sr = cur * ir / den, si = -cur * ii / den;
	for (int i = 0; i < n; i++) {
		const double xr = PR.ur[i], xi = PR.ui[i];
		PR.ur[i] = (sr * xr) - (si * xi);
		PR.ui[i] = (sr * xi) + (si * xr);
	}
	for (int i = 0; i < nn; i++) f[i] = 0;
	double (*wv)[3] = (double (*)[3])malloc((size_t)nA * 3 * sizeof(double));
	double (*wc)[3] = (double (*)[3])malloc((size_t)nA * 3 * sizeof(double));
	int32_t *dof = (int32_t *)malloc(PR_NDMAX * sizeof(int32_t));
	int bad = 0;
	for (int e = 0; (e < nt) && !bad; e++) {
		const double sg = sigelem[e];
		if (sg <= 0) continue;
		pr_dofs(e, dof);
		int32_t nd[TET_NMAX];
		tet_nodes(e, nd);
		const int nq = elem3d_nip(e);
		for (int q = 0; q < nq; q++) {
			double nv[ELEM_NMAX], g[ELEM_NMAX][3], w;
			if (elem3d_ip(e, q, nv, g, &w)) { bad = 1; break; }
			pr_basis_core(e, nv, g, wv, wc);
			// e = jω A + ∇φ : 実部 = −ω A_i + ∇φ_r、虚部 = ω A_r + ∇φ_i
			double er[3] = {0, 0, 0}, ei[3] = {0, 0, 0};
			for (int k = 0; k < nA; k++) {
				const double ar = PR.ur[dof[k]], ai = PR.ui[dof[k]];
				for (int c = 0; c < 3; c++) {
					er[c] += -omega * ai * wv[k][c];
					ei[c] += omega * ar * wv[k][c];
				}
			}
			for (int a = 0; a < nP; a++) {
				const double pr = PR.ur[dof[nA + a]], pi = PR.ui[dof[nA + a]];
				for (int c = 0; c < 3; c++) {
					er[c] += pr * g[a][c];
					ei[c] += pi * g[a][c];
				}
			}
			const double qd = 0.5 * sg * (dot3(er, er) + dot3(ei, ei)) * w;
			for (int a = 0; a < nP; a++) f[nd[a]] += qd * nv[a];
		}
	}
	free(wv); free(wc); free(dof);
	if (bad) {
		fprintf(fp_log, "*** AC heating : a degenerate element\n");
		return 1;
	}
	return 0;
}


void eddy3d_pr_heat_close(void)
{
	pr_close();
}


/*
過渡解析 (analysis = A + transient、後退 Euler、電圧駆動、線形、3 次以上の四面体)。1 次・2 次と同じく
χ = Δt φ を未知数にすると実対称・半正定値になる :

	K = [S + T/Δt, C/Δt; Cᵀ/Δt, P/Δt],   M = K − [S, 0; 0, 0],   K x^{n+1} = M [A^n; 0]

K を PR.Kr、M を PR.Ki に入れる (周波数領域と同じ非零パターン)。ポート 1 の電極の**全節点**を
χ = Δt V(t^{n+1})、端子電流はその行の反作用 Kx − M[A^n; 0] (周波数領域と同じく φ は格子の節点そのもの)。
*/
static int pr_nonlin(int m)
{
	return (Material[m].ja.on || (Material[m].nbh[0] > 0));
}


// skipnl = 1 : 非線形 (J-A / B-H) の材料の要素の S を入れない (固定点反復が積分点の H から足す)
static int pr_assemble_tr(FILE *fp_log, double dt, int skipnl)
{
	const int nt = PR.nt, nA = PR.nA, nP = PR.nP;
	crs_zero(&PR.Kr);
	crs_zero(&PR.Ki);
	double *se = (double *)malloc((size_t)nA * nA * sizeof(double));
	double *te = (double *)malloc((size_t)nA * nA * sizeof(double));
	double *ce = (double *)malloc((size_t)nA * nP * sizeof(double));
	double *pe = (double *)malloc((size_t)nP * nP * sizeof(double));
	double (*wv)[3] = (double (*)[3])malloc((size_t)nA * 3 * sizeof(double));
	double (*wc)[3] = (double (*)[3])malloc((size_t)nA * 3 * sizeof(double));
	double (*nc)[3] = (double (*)[3])malloc((size_t)nA * 3 * sizeof(double));
	int32_t *dof = (int32_t *)malloc(PR_NDMAX * sizeof(int32_t));
	int bad = 0;
	for (int e = 0; (e < nt) && !bad; e++) {
		const int m = elem3d_mat(e);
		double nu[6];
		material_coef_pub(m, 3, nu);
		pr_dofs(e, dof);
		if (pr_element(e, nu, NULL, Material[m].sigma, se, te, ce, pe, wv, wc, nc)) { bad = 1; break; }
		for (int k = 0; k < nA; k++) {
			for (int l = 0; l < nA; l++) {
				const int64_t p = pr_find(&PR.Kr, dof[k], dof[l]);
				if (p < 0) continue;		// 分散した行列 (mpisolve = 2) の他の rank の行
				const double t = te[(k * nA) + l] / dt;
				PR.Kr.val[p] += ((skipnl && pr_nonlin(m)) ? 0 : se[(k * nA) + l]) + t;
				PR.Ki.val[p] += t;
			}
			for (int a = 0; a < nP; a++) {
				const double v = ce[(k * nP) + a] / dt;
				const int64_t p1 = pr_find(&PR.Kr, dof[k], dof[nA + a]);
				const int64_t p2 = pr_find(&PR.Kr, dof[nA + a], dof[k]);
				if (p1 >= 0) {
					PR.Kr.val[p1] += v;
					PR.Ki.val[p1] += v;
				}
				if (p2 >= 0) {
					PR.Kr.val[p2] += v;
					PR.Ki.val[p2] += v;
				}
			}
		}
		for (int a = 0; a < nP; a++) {
			for (int b = 0; b < nP; b++) {
				const int64_t p = pr_find(&PR.Kr, dof[nA + a], dof[nA + b]);
				if (p < 0) continue;
				PR.Kr.val[p] += pe[(a * nP) + b] / dt;
				PR.Ki.val[p] += pe[(a * nP) + b] / dt;
			}
		}
	}
	free(se); free(te); free(ce); free(pe); free(wv); free(wc); free(nc); free(dof);
	if (bad) {
		fprintf(fp_log, "*** the high-order A-phi solver : a degenerate element\n");
		return 1;
	}
	return 0;
}


/*
過渡解析の非線形 (J-A / B-H)。2 次 (eddy3d_p2.c の p2nl_*) と同じ固定点反復を**積分点ごとの状態**で解く:

	F(x) = K_l x − M [A^n; 0] + Σ_e Σ_q w_q H_q・∇×W_k   (K_l は非線形の要素の S を除いた行列)
	(K_l + Σ w_q (∇×W)ᵀ N_q (∇×W)) δ = −F,   N_q = dH/dB の対称部分 (J-A) / ν I + 2 (dν/dB²) B Bᵀ (B-H)

∇×W は積分点あたり 3 nA 個。合計が PRNL_CACHE_MAX 個以内なら記憶し、越えたら毎回作り直す (答えは同じ)。
接線の要素行列 (nA² 個) は要素の塊ごとに並列に作って直列に足す (塊の大きさは PRNL_LOC_MAX で決める)。
*/
typedef struct {
	int nj, nq, nA;
	int *el;			// [nj] 要素番号
	int32_t *dof;		// [nj*nA] A の自由度
	double *w;			// [nj*nq] 積分の重み
	double *b, *h, *m, *d;		// [nj*nq*3] 確定済みの状態
	double *bn, *hn, *mn, *dn;	// [nj*nq*3] 反復中
	double *nu6;		// [nj*nq*6] 行列の傾き
	double *lv;			// [nj*nA] 残差の要素の作業域
	double *loc;		// [nb*nA*nA] 接線の要素行列の作業域 (塊 nb 要素分)
	int nb;
	double *cv;			// [nj*nq*nA*3] ∇×W (記憶するときだけ、しないなら NULL)
} prnl_t;

#define PRNL_CACHE_MAX ((int64_t)32 * 1024 * 1024)	// ∇×W を記憶する上限 (double の個数、256 MB)
#define PRNL_LOC_MAX ((int64_t)4 * 1024 * 1024)	// 接線の要素行列の塊の上限 (double の個数、32 MB)


// 要素 j の積分点 q の ∇×W_k (wc は nA x 3、wv は作り直すときの作業域)
static void prnl_curl(const prnl_t *J, int j, int q, double (*wc)[3], double (*wv)[3])
{
	const int nA = J->nA;
	if (J->cv != NULL) {
		const double *cv = &J->cv[(((int64_t)j * J->nq) + q) * nA * 3];
		for (int k = 0; k < nA; k++) {
			for (int c = 0; c < 3; c++) wc[k][c] = cv[(k * 3) + c];
		}
		return;
	}
	double nv[ELEM_NMAX], g[ELEM_NMAX][3], w;
	if (elem3d_ip(J->el[j], q, nv, g, &w)) {
		memset(wc, 0, (size_t)nA * 3 * sizeof(double));
		return;
	}
	pr_basis_core(J->el[j], nv, g, wv, wc);
}


static void prnl_close(prnl_t *J)
{
	free(J->el); free(J->dof); free(J->w);
	free(J->b); free(J->h); free(J->m); free(J->d);
	free(J->bn); free(J->hn); free(J->mn); free(J->dn);
	free(J->nu6); free(J->lv); free(J->loc); free(J->cv);
	memset(J, 0, sizeof(*J));
}


static int prnl_open(prnl_t *J, FILE *fp_log)
{
	const int nt = PR.nt, nA = PR.nA;
	int nj = 0;
	for (int e = 0; e < nt; e++) {
		if (pr_nonlin(elem3d_mat(e))) nj++;
	}
	const int nq = elem3d_nip(0);
	J->nj = nj;
	J->nq = nq;
	J->nA = nA;
	const size_t s1 = (size_t)((nj > 0) ? nj : 1), sq = s1 * (size_t)nq;
	J->el = (int *)malloc(s1 * sizeof(int));
	J->dof = (int32_t *)malloc(s1 * nA * sizeof(int32_t));
	J->w = (double *)calloc(sq, sizeof(double));
	double **v3[8] = {&J->b, &J->h, &J->m, &J->d, &J->bn, &J->hn, &J->mn, &J->dn};
	for (int q = 0; q < 8; q++) *v3[q] = (double *)calloc(sq * 3, sizeof(double));
	J->nu6 = (double *)calloc(sq * 6, sizeof(double));
	J->lv = (double *)calloc(s1 * nA, sizeof(double));
	int64_t nb = PRNL_LOC_MAX / ((int64_t)nA * nA);
	if (nb < 1) nb = 1;
	if (nb > nj) nb = ((nj > 0) ? nj : 1);
	J->nb = (int)nb;
	J->loc = (double *)calloc((size_t)nb * nA * nA, sizeof(double));
	J->cv = ((((int64_t)nj * nq * nA * 3) <= PRNL_CACHE_MAX) ? (double *)malloc(sq * nA * 3 * sizeof(double)) : NULL);
	double (*wv)[3] = (double (*)[3])malloc((size_t)nA * 3 * sizeof(double));
	double (*wc)[3] = (double (*)[3])malloc((size_t)nA * 3 * sizeof(double));
	int32_t *dof = (int32_t *)malloc(PR_NDMAX * sizeof(int32_t));
	int j = 0, bad = 0;
	for (int e = 0; (e < nt) && !bad; e++) {
		if (!pr_nonlin(elem3d_mat(e))) continue;
		pr_dofs(e, dof);
		J->el[j] = e;
		for (int k = 0; k < nA; k++) J->dof[((int64_t)j * nA) + k] = dof[k];
		for (int q = 0; q < nq; q++) {
			double nv[ELEM_NMAX], g[ELEM_NMAX][3], w;
			if (elem3d_ip(e, q, nv, g, &w)) { bad = 1; break; }
			J->w[((int64_t)j * nq) + q] = w;
			if (J->cv != NULL) {
				pr_basis_core(e, nv, g, wv, wc);
				double *cv = &J->cv[(((int64_t)j * nq) + q) * nA * 3];
				for (int k = 0; k < nA; k++) {
					for (int c = 0; c < 3; c++) cv[(k * 3) + c] = wc[k][c];
				}
			}
		}
		j++;
	}
	free(wv); free(wc); free(dof);
	if (bad) {
		fprintf(fp_log, "*** the high-order A-phi solver : a degenerate element\n");
		return 1;
	}
	fprintf(fp_log, "  nonlinear : %d tetrahedra x %d points with J-A (vector) / bh (nlsolver = %d %.3e %.3f)%s\n",
		nj, nq, NlMaxiter, NlTol, NlRelax, ((J->cv != NULL) ? "" : ", curl W recomputed"));
	for (int m = 0; m < NMaterial; m++) {
		if (!Material[m].ja.on) continue;
		fprintf(fp_log, "  J-A model check (material %d) : 3-D vector J-A along a fixed direction == scalar J-A : "
			"max |dH| / max |H| = %.3e\n", m, ja_vec3_selftest(&Material[m].ja));
	}
	return 0;
}


// x の A から各積分点の B を作り、確定済みの状態から H と行列の傾きを求める (2 次の p2nl_eval と同じ式)
static void prnl_eval(prnl_t *J, const double *x)
{
	const int nj = J->nj, nq = J->nq, nA = J->nA;
	int j;
#ifdef _OPENMP
#pragma omp parallel
#endif
	{
	double (*wc)[3] = (double (*)[3])malloc((size_t)nA * 3 * sizeof(double));
	double (*wv)[3] = (double (*)[3])malloc((size_t)nA * 3 * sizeof(double));
#ifdef _OPENMP
#pragma omp for
#endif
	for (j = 0; j < nj; j++) {
		const material_t *mt = &Material[elem3d_mat(J->el[j])];
		const ja_t *ja = &mt->ja;
		for (int q = 0; q < nq; q++) {
			const int64_t p = ((int64_t)j * nq) + q;
			prnl_curl(J, j, q, wc, wv);
			double bv[3] = {0, 0, 0};
			for (int k = 0; k < nA; k++) {
				const double a = x[J->dof[((int64_t)j * nA) + k]];
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
	free(wc); free(wv);
	}
}


// 刻みの磁気の仕事 Σ w_q H_q・(B_q − B_q^n)
static double prnl_work(const prnl_t *J)
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


static void prnl_commit(prnl_t *J)
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
static double prnl_resid(prnl_t *J, const crs_t *Kl, const double *rhs, const unsigned char *fix,
	const double *x, double *r, double bref)
{
	const int n = (int)Kl->n, nj = J->nj, nq = J->nq, nA = J->nA;
	crs_spmv(Kl, x, r, NULL);
	for (int i = 0; i < n; i++) r[i] = rhs[i] - r[i];
	int j;
#ifdef _OPENMP
#pragma omp parallel
#endif
	{
	double (*wc)[3] = (double (*)[3])malloc((size_t)nA * 3 * sizeof(double));
	double (*wv)[3] = (double (*)[3])malloc((size_t)nA * 3 * sizeof(double));
#ifdef _OPENMP
#pragma omp for
#endif
	for (j = 0; j < nj; j++) {
		double *lv = &J->lv[(int64_t)j * nA];
		for (int k = 0; k < nA; k++) lv[k] = 0;
		for (int q = 0; q < nq; q++) {
			const int64_t p = ((int64_t)j * nq) + q;
			prnl_curl(J, j, q, wc, wv);
			for (int k = 0; k < nA; k++) lv[k] += J->w[p] * dot3(&J->hn[p * 3], wc[k]);
		}
	}
	free(wc); free(wv);
	}
	for (j = 0; j < nj; j++) {
		for (int k = 0; k < nA; k++) r[J->dof[((int64_t)j * nA) + k]] -= J->lv[((int64_t)j * nA) + k];
	}
	double s = 0;
	for (int i = 0; i < n; i++) {
		if (fix[i]) r[i] = 0;
		s += r[i] * r[i];
	}
	return ((bref > 0) ? (sqrt(s) / bref) : sqrt(s));
}


// 接線の行列 Kt = K_l + Σ w_q (∇×W)ᵀ N_q (∇×W)。要素の塊ごとに並列に作って直列に足す
static void prnl_tangent(prnl_t *J, const crs_t *Kl, crs_t *Kt)
{
	const int nj = J->nj, nq = J->nq, nA = J->nA, nb = J->nb;
	memcpy(Kt->val, Kl->val, (size_t)Kl->rowptr[Kl->n] * sizeof(double));
	for (int j0 = 0; j0 < nj; j0 += nb) {
		const int j1 = (((j0 + nb) < nj) ? (j0 + nb) : nj);
		int j;
#ifdef _OPENMP
#pragma omp parallel
#endif
		{
		double (*wc)[3] = (double (*)[3])malloc((size_t)nA * 3 * sizeof(double));
		double (*wv)[3] = (double (*)[3])malloc((size_t)nA * 3 * sizeof(double));
		double (*nc)[3] = (double (*)[3])malloc((size_t)nA * 3 * sizeof(double));
#ifdef _OPENMP
#pragma omp for
#endif
		for (j = j0; j < j1; j++) {
			double *ke = &J->loc[(int64_t)(j - j0) * nA * nA];
			for (int k = 0; k < nA * nA; k++) ke[k] = 0;
			for (int q = 0; q < nq; q++) {
				const int64_t p = ((int64_t)j * nq) + q;
				prnl_curl(J, j, q, wc, wv);
				for (int k = 0; k < nA; k++) nu_mul(&J->nu6[p * 6], wc[k], nc[k]);
				for (int k = 0; k < nA; k++) {
					for (int l = 0; l < nA; l++) ke[(k * nA) + l] += J->w[p] * dot3(nc[k], wc[l]);
				}
			}
		}
		free(wc); free(wv); free(nc);
		}
		for (j = j0; j < j1; j++) {
			const double *ke = &J->loc[(int64_t)(j - j0) * nA * nA];
			const int32_t *dj = &J->dof[(int64_t)j * nA];
			for (int k = 0; k < nA; k++) {
				for (int l = 0; l < nA; l++) {
					const int64_t p = pr_find(Kt, dj[k], dj[l]);
					if (p >= 0) Kt->val[p] += ke[(k * nA) + l];		// 分散した行列の他の rank の行は飛ばす
				}
			}
		}
	}
}


// 1 刻みの固定点反復 (2 次の p2nl_step と同じ安定化)。x は前の刻みの解で入り、収束した解で出る
static int prnl_step(prnl_t *J, const crs_t *Kl, crs_t *Kt, const double *rhs, const unsigned char *fix,
	const double *xd, double *x, double *r, double *del, double *xbest, double bref, int st, FILE *fp_log)
{
	const int n = (int)Kl->n;
	for (int i = 0; i < n; i++) {
		if (fix[i]) x[i] = xd[i];
	}
	int it, conv = 0, nback = 0;
	double rnorm = 0, rbest = 0, lrel = NlRelax;
	for (it = 1; it <= NlMaxiter; it++) {
		prnl_eval(J, x);
		rnorm = prnl_resid(J, Kl, rhs, fix, x, r, bref);
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
		prnl_tangent(J, Kl, Kt);
		double tol = ((rnorm > 0) ? (Solver.converg / rnorm) : Solver.converg);
		if (tol > 0.1) tol = 0.1;
		if (tol < Solver.converg) tol = Solver.converg;
		for (int i = 0; i < n; i++) del[i] = 0;
		const int li = solver_cg(Kt, r, del, fix, Solver.maxiter, Solver.nout, tol, NULL, "a3jr");
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


// 刻みの終わり t の電圧の倍率 (waveform の折れ線の左極限、無ければ 1。1 次・2 次と同じ)
static double pr_wave(double t)
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


// 電極の節点 (ポート 1) : χ を固定して励振し、端子電流をその行から取る
static int pr_port1_node(int i)
{
	return (PR.fix[PR.oV + i] && PR.cond[i] && PR.used[i] && (NodeConductor[i] == 1));
}


int solve_eddy3d_pr_transient(FILE *fp_log)
{
	int ierr = 0;
	fprintf(fp_log, "\n=== 3D eddy current, transient (A-phi, order-%d Nedelec / nodal elements, backward Euler) ===\n", TetOrder);
	if (NPort != 1) {
		fprintf(fp_log, "*** analysis A + transient needs exactly one port\n");
		return 1;
	}
	const double fsave = Freq;
	if (Freq <= 0) Freq = 1;		// pr_open の検査を通すだけ (過渡解析は ω を使わない)
	const int oe = pr_open(fp_log);
	Freq = fsave;
	if (oe) return 1;
	const int n = PR.n, nn = PR.nn, oV = PR.oV;
	const double dt = TransientTime / TransientSteps;
	int hyst = 0;
	for (int m = 0; m < NMaterial; m++) {
		if (pr_nonlin(m)) hyst = 1;
	}
	fprintf(fp_log, "  t_end = %.6e [s], %d steps (dt = %.6e [s]), voltage = %.6e [V]%s\n",
		TransientTime, TransientSteps, dt, Volt, ((WaveN > 0) ? " x waveform" : " (step)"));
	if (pr_assemble_tr(fp_log, dt, hyst)) {
		pr_close();
		return 1;
	}
	const crs_t *K = &PR.Kr, *M = &PR.Ki;
	prnl_t J;
	memset(&J, 0, sizeof(J));
	crs_t Kt;
	memset(&Kt, 0, sizeof(Kt));
	if (hyst) {
		if (prnl_open(&J, fp_log)) ierr = 1;
		pr_crs(&Kt);
		crs_distribute(&Kt);		// mpisolve = 2 : K と同じ分け方 (Kt = K の値の写し + 接線)
	}
	double bref = 0;
	const unsigned char *fix = PR.fix;
	double *x = (double *)calloc((size_t)n, sizeof(double));
	double *xo = (double *)calloc((size_t)n, sizeof(double));
	double *xd = (double *)calloc((size_t)n, sizeof(double));
	double *rhs = (double *)malloc((size_t)n * sizeof(double));
	double *b = (double *)malloc((size_t)n * sizeof(double));
	double *w1 = (double *)malloc((size_t)n * sizeof(double));
	double *u = (double *)calloc((size_t)n, sizeof(double));
	double *z = (double *)malloc((size_t)n * sizeof(double));
	FILE *fcsv = mpi_fopen_out("ofe_transient.csv", "w");
	if (fcsv) fprintf(fcsv, "t [s],V [V],I [A],W_in [J],W_mag [J],W_joule [J],iter\n");
	double win = 0, wmag = 0, wjou = 0;
	long itotal = 0;
	fprintf(fp_log, "  %6s %13s %13s %13s %6s\n", "step", "t [s]", "V [V]", "I [A]", "iter");
	for (int st = 1; (st <= TransientSteps) && !ierr; st++) {
		const double t = st * dt;
		const double v = Volt * pr_wave(t);
		// A の自由度 (φ の欄より前) だけを残した [A^n; 0]
		for (int i = 0; i < n; i++) xo[i] = ((i < oV) ? x[i] : 0);
		crs_spmv(M, xo, rhs, NULL);
		for (int i = 0; i < n; i++) xd[i] = 0;
		for (int i = 0; i < nn; i++) {
			if (pr_port1_node(i)) xd[oV + i] = dt * v;
		}
		crs_spmv(K, xd, w1, NULL);
		for (int i = 0; i < n; i++) b[i] = (fix[i] ? 0 : (rhs[i] - w1[i]));
		int it;
		if (!hyst) {
			for (int i = 0; i < n; i++) u[i] = (fix[i] ? 0 : x[i]);
			it = solver_cg(K, b, u, fix, Solver.maxiter, Solver.nout, Solver.converg, NULL, "a3tr");
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
			it = prnl_step(&J, K, &Kt, rhs, fix, xd, x, b, u, z, bref, st, fp_log);
			if (it < 0) {
				ierr = 1;
				break;
			}
		}
		itotal += it;
		crs_spmv(K, x, w1, NULL);
		double cur = 0;
		for (int i = 0; i < nn; i++) {
			if (pr_port1_node(i)) cur += w1[oV + i] - rhs[oV + i];
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
		if (hyst) zkx += prnl_work(&J);
		win += v * cur * dt;
		wjou += zmz;
		wmag += zkx - zmx;
		if (hyst) prnl_commit(&J);
		if (fcsv) fprintf(fcsv, "%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%d\n", t, v, cur, win, wmag, wjou, it);
		if ((st == TransientSteps) || ((st % ((TransientSteps >= 20) ? (TransientSteps / 20) : 1)) == 0)) {
			fprintf(fp_log, "  %6d %13.6e %13.6e %13.6e %6d\n", st, t, v, cur, it);
		}
	}
	if (fcsv) fclose(fcsv);
	if (!ierr && !hyst) fprintf(fp_log, "  CG iterations = %ld (%d steps)\n", itotal, TransientSteps);
	if (!ierr && hyst) {
		fprintf(fp_log, "  nonlinear (%d tetrahedra) : fixed-point iterations = %ld (%d steps)\n",
			J.nj, itotal, TransientSteps);
	}
	if (hyst) {
		prnl_close(&J);
		crs_free(&Kt);
	}
	fprintf(fp_log, "  energy : terminals = %.12e, magnetic work = %.12e, joule = %.12e [J] "
		"(terminals - magnetic - joule = %.3e)\n", win, wmag, wjou, win - wmag - wjou);
	free(x); free(xo); free(xd); free(rhs); free(b); free(w1); free(u); free(z);
	pr_close();
	return ierr;
}


/*
辺要素の自己検証の高次版 (analysis = E、次数 r の四面体)。次数 r の Nedelec 空間 (第 1 種) は
P_{r−1}³ ⊕ r × P̃_{r−1}³ (P̃ は同次式) なので、そこに入る場を σ で重みをつけた L2 射影すると誤差が 0 になり、
回転のエネルギーも厳密に合う:

	T u = b,   b_k = ∫ σ W_k・E dV,   誤差² = ∫ σ|E|² − uᵀ b = 0,   uᵀ S u = ∫ (∇×E)ᵀ ν (∇×E) dV

(Pra) 次数 r − 1 の一般の場 (1 次の a + B r に r − 1 次と r − 2 次の単項式を足したもの)、(Prb) それに
r' × q (q は r − 1 次の同次式) を足した次数 r の場 (P_{r−1}³ に入らない部分。面・内部の関数が要る)。
場は多項式の単項式の表で持ち、回転も表から作る (式の手書きの誤りを避ける)。厳密なのは直線の要素だけ。
*/
#define PR_MMAX (32)
typedef struct { double c; int e[3]; } prmono_t;
typedef struct { int n; prmono_t m[PR_MMAX]; } prpoly_t;

static void poly_add(prpoly_t *p, double c, int i, int j, int k)
{
	if ((c == 0) || (i < 0) || (j < 0) || (k < 0) || (p->n >= PR_MMAX)) return;
	p->m[p->n].c = c;
	p->m[p->n].e[0] = i; p->m[p->n].e[1] = j; p->m[p->n].e[2] = k;
	p->n++;
}

static double poly_eval(const prpoly_t *p, const double r[3])
{
	double s = 0;
	for (int q = 0; q < p->n; q++) {
		double v = p->m[q].c;
		for (int c = 0; c < 3; c++) {
			for (int k = 0; k < p->m[q].e[c]; k++) v *= r[c];
		}
		s += v;
	}
	return s;
}

// d (∂/∂x_d) に係数 s を掛けたものを o に足す
static void poly_dadd(prpoly_t *o, const prpoly_t *p, int d, double s)
{
	for (int q = 0; q < p->n; q++) {
		const prmono_t *m = &p->m[q];
		if (m->e[d] == 0) continue;
		int e[3] = {m->e[0], m->e[1], m->e[2]};
		e[d]--;
		poly_add(o, s * m->c * m->e[d], e[0], e[1], e[2]);
	}
}

// 回転 (r' = (r − c0)/L の多項式なので 1/L を掛ける)
static void poly_curl(const prpoly_t f[3], prpoly_t c[3], double il)
{
	for (int i = 0; i < 3; i++) c[i].n = 0;
	poly_dadd(&c[0], &f[2], 1, il); poly_dadd(&c[0], &f[1], 2, -il);
	poly_dadd(&c[1], &f[0], 2, il); poly_dadd(&c[1], &f[2], 0, -il);
	poly_dadd(&c[2], &f[1], 0, il); poly_dadd(&c[2], &f[0], 1, -il);
}

// o += r' × q (q の各成分は単項式の表)
static void poly_rcross(prpoly_t o[3], const prpoly_t q[3])
{
	for (int i = 0; i < 3; i++) {
		const int j = (i + 1) % 3, k = (i + 2) % 3;
		// (r × q)_i = r_j q_k − r_k q_j
		for (int t = 0; t < q[k].n; t++) {
			const prmono_t *m = &q[k].m[t];
			int e[3] = {m->e[0], m->e[1], m->e[2]};
			e[j]++;
			poly_add(&o[i], m->c, e[0], e[1], e[2]);
		}
		for (int t = 0; t < q[j].n; t++) {
			const prmono_t *m = &q[j].m[t];
			int e[3] = {m->e[0], m->e[1], m->e[2]};
			e[k]++;
			poly_add(&o[i], -m->c, e[0], e[1], e[2]);
		}
	}
}


int solve_edge_test_pr(FILE *fp_log)
{
	fprintf(fp_log, "\n=== edge element (order-%d Nedelec) self test ===\n", TetOrder);
	if (pr_number(fp_log, "high-order edge self test")) return 1;
	const int r = PR.r, n = PR.n, nt = PR.nt, nn = PR.nn, nA = PR.nA, oV = PR.oV;
	// 直線の要素だけ (節点が頂点の 1 次の写像の位置にあること)
	for (int e = 0; e < nt; e++) {
		int32_t nd[TET_NMAX];
		const int nen = tet_nodes(e, nd);
		const double h[3] = {Xp[nd[1]] - Xp[nd[0]], Yp[nd[1]] - Yp[nd[0]], Zp[nd[1]] - Zp[nd[0]]};
		for (int a = 4; a < nen; a++) {
			const int *ix = tet_ho_index(a);
			double x[3] = {0, 0, 0};
			for (int i = 0; i < 4; i++) {
				x[0] += ix[i] * Xp[nd[i]] / r;
				x[1] += ix[i] * Yp[nd[i]] / r;
				x[2] += ix[i] * Zp[nd[i]] / r;
			}
			const double d[3] = {Xp[nd[a]] - x[0], Yp[nd[a]] - x[1], Zp[nd[a]] - x[2]};
			if (dot3(d, d) > (1e-16 * dot3(h, h))) {
				fprintf(fp_log, "*** the high-order edge self test needs straight-sided elements "
					"(tetrahedron %d has a node off its affine position); the identities are exact "
					"only for affine elements\n", e + 1);
				pr_close();
				return 1;
			}
		}
	}
	fprintf(fp_log, "  order %d : nodes = %d, tetrahedra = %d, edges = %d, faces = %d, edge unknowns = %d "
		"(%d per element)\n", r, nn, nt, PR.ne, PR.nf, oV, nA);
	crs_t S, T;
	pr_crs(&S);
	pr_crs(&T);
	double c0[3] = {0, 0, 0}, lmax = 0;
	{
		double mn[3] = {1e300, 1e300, 1e300}, mx[3] = {-1e300, -1e300, -1e300};
		for (int i = 0; i < nn; i++) {
			const double x[3] = {Xp[i], Yp[i], Zp[i]};
			for (int c = 0; c < 3; c++) {
				if (x[c] < mn[c]) mn[c] = x[c];
				if (x[c] > mx[c]) mx[c] = x[c];
			}
		}
		for (int c = 0; c < 3; c++) {
			c0[c] = (mn[c] + mx[c]) / 2;
			if ((mx[c] - mn[c]) > lmax) lmax = mx[c] - mn[c];
		}
	}
	if (!(lmax > 0)) lmax = 1;
	// 場 (r' = (r − c0)/L の多項式)
	prpoly_t f1[3], f2[3], q[3], g1[3], g2[3];
	memset(f1, 0, sizeof(f1)); memset(q, 0, sizeof(q));
	{
		const double av[3] = {1.0, 2.0, 3.0};
		const double bm[3][3] = {{0.3, -1.2, 0.7}, {0.5, 0.2, -0.4}, {1.1, 0.6, -0.9}};
		for (int i = 0; i < 3; i++) {
			poly_add(&f1[i], av[i], 0, 0, 0);
			poly_add(&f1[i], bm[i][0], 1, 0, 0);
			poly_add(&f1[i], bm[i][1], 0, 1, 0);
			poly_add(&f1[i], bm[i][2], 0, 0, 1);
		}
		const int d = r - 1;
		// 次数 r − 1 の項 (各方向の最高次と混ざった項) と r − 2 次の項
		poly_add(&f1[0], 1.0, 0, d, 0);
		poly_add(&f1[0], 0.4, d - 2, 1, 1);
		poly_add(&f1[0], -0.7, 1, 0, d - 1);
		poly_add(&f1[1], -0.8, 0, 0, d);
		poly_add(&f1[1], 0.6, 1, d - 2, 1);
		poly_add(&f1[1], 0.5, d - 1, 1, 0);
		poly_add(&f1[2], 0.9, d, 0, 0);
		poly_add(&f1[2], -0.5, 1, d - 1, 0);
		poly_add(&f1[2], 0.3, 0, 1, d - 2);
		// q : 次数 r − 1 の同次式
		poly_add(&q[0], 1.0, 0, d, 0);
		poly_add(&q[0], -0.6, 1, 0, d - 1);
		poly_add(&q[1], 0.7, d - 1, 0, 1);
		poly_add(&q[1], 0.2, 0, 0, d);
		poly_add(&q[2], -0.9, d, 0, 0);
		poly_add(&q[2], 0.5, 0, d - 1, 1);
	}
	memcpy(f2, f1, sizeof(f1));
	poly_rcross(f2, q);
	poly_curl(f1, g1, 1 / lmax);
	poly_curl(f2, g2, 1 / lmax);

	int ierr = 0;
	for (int e = 0; e < nt; e++) {
		if (!(Material[elem3d_mat(e)].sigma > 0)) {
			fprintf(fp_log, "*** the high-order edge self test needs sigma > 0 in every element\n");
			ierr = 1;
			break;
		}
	}
	double *b1 = (double *)calloc((size_t)n, sizeof(double));
	double *b2 = (double *)calloc((size_t)n, sizeof(double));
	double se1 = 0, se2 = 0, sc1 = 0, sc2 = 0;
	double *sk = (double *)malloc((size_t)nA * nA * sizeof(double));
	double *tk = (double *)malloc((size_t)nA * nA * sizeof(double));
	double *h1 = (double *)malloc((size_t)nA * sizeof(double));
	double *h2 = (double *)malloc((size_t)nA * sizeof(double));
	double (*wv)[3] = (double (*)[3])malloc((size_t)nA * 3 * sizeof(double));
	double (*wc)[3] = (double (*)[3])malloc((size_t)nA * 3 * sizeof(double));
	double (*nc)[3] = (double (*)[3])malloc((size_t)nA * 3 * sizeof(double));
	int32_t *dof = (int32_t *)malloc(PR_NDMAX * sizeof(int32_t));
	for (int e = 0; (e < nt) && !ierr; e++) {
		const int m = elem3d_mat(e);
		double nu[6];
		material_coef_pub(m, 3, nu);
		const double sig = Material[m].sigma;
		pr_dofs(e, dof);
		memset(sk, 0, (size_t)nA * nA * sizeof(double));
		memset(tk, 0, (size_t)nA * nA * sizeof(double));
		memset(h1, 0, (size_t)nA * sizeof(double));
		memset(h2, 0, (size_t)nA * sizeof(double));
		int32_t nd[TET_NMAX];
		const int nen = tet_nodes(e, nd);
		const int nq = elem3d_nip(e);
		for (int qq = 0; qq < nq; qq++) {
			double nv[ELEM_NMAX], g[ELEM_NMAX][3], w;
			if (elem3d_ip(e, qq, nv, g, &w)) { ierr = 1; break; }
			pr_basis_core(e, nv, g, wv, wc);
			double x[3] = {0, 0, 0};
			for (int a = 0; a < nen; a++) {
				x[0] += nv[a] * (Xp[nd[a]] - c0[0]) / lmax;
				x[1] += nv[a] * (Yp[nd[a]] - c0[1]) / lmax;
				x[2] += nv[a] * (Zp[nd[a]] - c0[2]) / lmax;
			}
			double e1[3], e2[3], k1[3], k2[3], n1[3], n2[3];
			for (int i = 0; i < 3; i++) {
				e1[i] = poly_eval(&f1[i], x);
				e2[i] = poly_eval(&f2[i], x);
				k1[i] = poly_eval(&g1[i], x);
				k2[i] = poly_eval(&g2[i], x);
			}
			nu_mul(nu, k1, n1);
			nu_mul(nu, k2, n2);
			se1 += w * sig * dot3(e1, e1);
			se2 += w * sig * dot3(e2, e2);
			sc1 += w * dot3(n1, k1);
			sc2 += w * dot3(n2, k2);
			for (int k = 0; k < nA; k++) nu_mul(nu, wc[k], nc[k]);
			for (int k = 0; k < nA; k++) {
				h1[k] += w * sig * dot3(wv[k], e1);
				h2[k] += w * sig * dot3(wv[k], e2);
				for (int l = k; l < nA; l++) {
					sk[(k * nA) + l] += w * dot3(nc[k], wc[l]);
					tk[(k * nA) + l] += w * sig * dot3(wv[k], wv[l]);
				}
			}
		}
		for (int k = 0; k < nA; k++) {
			for (int l = 0; l < k; l++) {
				sk[(k * nA) + l] = sk[(l * nA) + k];
				tk[(k * nA) + l] = tk[(l * nA) + k];
			}
		}
		for (int k = 0; k < nA; k++) {
			b1[dof[k]] += h1[k];
			b2[dof[k]] += h2[k];
			for (int l = 0; l < nA; l++) {
				const int64_t p = pr_find(&S, dof[k], dof[l]);
				S.val[p] += sk[(k * nA) + l];
				T.val[p] += tk[(k * nA) + l];
			}
		}
	}
	free(sk); free(tk); free(h1); free(h2); free(wv); free(wc); free(nc); free(dof);
	unsigned char *fix = (unsigned char *)calloc((size_t)n, sizeof(unsigned char));
	for (int i = oV; i < n; i++) fix[i] = 1;
	double *u = (double *)calloc((size_t)n, sizeof(double));
	double *w1 = (double *)malloc((size_t)n * sizeof(double));
	const char *lab[2] = {"(Pra) general field of degree r-1", "(Prb) + r x q (q homogeneous of degree r-1, needs the face / interior functions)"};
	double *bb[2] = {b1, b2};
	const double sev[2] = {se1, se2}, scv[2] = {sc1, sc2};
	for (int t = 0; (t < 2) && !ierr; t++) {
		for (int i = 0; i < n; i++) u[i] = 0;
		const int it = solver_cg(&T, bb[t], u, fix, Solver.maxiter, Solver.nout, Solver.converg, NULL, "pre");
		if (it < 0) {
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
		fprintf(fp_log, "  %s : projection error^2 / int sigma|E|^2 = %.3e, curl energy u^T S u = %.12e (exact %.12e, %.3e), %d CG iterations %s\n",
			lab[t], ep, us, scv[t], ec, it, (ok ? "passed" : "FAILED"));
		if (!ok) ierr = 1;
	}
	fprintf(fp_log, "  order-%d edge self test : %s\n", r, (ierr ? "FAILED" : "passed"));
	free(b1); free(b2); free(u); free(w1); free(fix);
	crs_free(&S);
	crs_free(&T);
	pr_close();
	return ierr;
}
