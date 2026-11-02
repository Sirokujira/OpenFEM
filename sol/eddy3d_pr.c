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
#define PR_NAMAX (3 * 5 * 6 * 6)		// 要素の A の自由度 (四面体 r = 7 で 315、六面体 p = 5 で 3p(p+1)² = 540、角柱 p = 5 で 315)
#define PR_NPMAX (ELEM_NMAX)		// 要素の φ の自由度 (四面体 r = 7 で 120、六面体 p = 5 で 216)
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
	int hex;			// 高次の六面体の格子なら次数 p (A は部分辺の自由度、下の「高次の六面体」)
	int mix;			// 種別の混在した高次の格子なら次数 p (下の「種別の混在」。部分辺・三角形の面・四面体の内部)
	int dummy;			// 種別の混在 : 要素の自由度を最大の数に揃える詰め物の未知数 (固定、基底は 0)
	int32_t *hse;		// [nt * nA] 要素の局所の部分辺 -> 全体の部分辺
	signed char *hsg;	// [nt * nA] 向きの符号 (全体の向きは小さい節点番号 -> 大きい方)
	int32_t *hsa, *hsb;	// [ne] 部分辺の両端の節点 (a < b)
	int prism;			// 高次の角柱の格子なら次数 p (下の「高次の角柱」)
	int nz;				// 角柱 : ζ 方向の部分辺の数 (未知数の並びは [部分辺 | 層の辺 x NtE | 層の面 x NtF | φ])
	int32_t *pe, *pf;	// [nt * 3(p+1)] / [nt * (p+1)] 要素の層の辺 / 層の面 -> 全体の番号
} PR;

/*
---- 高次の六面体 (27 / 64 / 125 / 216 節点、完全 Lagrange、等パラメトリック) ----

A を次数 p の Nedelec 要素 (第 1 種、x 成分が Q_{p−1,p,p}、y・z も同様) で、φ を六面体の節点そのもので離散化する。
自由度は**格子の隣り合う節点を結ぶ部分辺** (要素あたり 3p(p+1)²) の線積分で、基底は 1 次元の histopolation の関数

	e_i(ξ) = −Σ_{m ≤ i} ℓ_m'(ξ)     (ℓ_m は等間隔の p 次 Lagrange。∫_{ξ_m}^{ξ_{m+1}} e_i dξ = δ_im)

と Lagrange の積 : 方向 d の部分辺 (始点の格子位置 i) の関数は W = e_{i_d}(ξ_d) Π_{c≠d} ℓ_{i_c}(ξ_c) ∇ξ_d。参照座標 ξ_c は
節点の値の補間 Σ_a ξ_c(a) N_a で厳密に表せるので、∇ξ_c = J⁻ᵀ e_c (共変 Piola) が形状関数の勾配から作れ、
∇×W = ∇f × ∇ξ_d が曲がった要素でも厳密。ℓ_a' = e_{a−1} − e_a なので**節点の関数の勾配は部分辺の関数の接続行列
の和に厳密に一致**し (∇N の係数は端点の ±1)、ゲージの零空間は 1 次と同じく部分辺のグラフの全域木で消える。
隣の要素と共有する面の部分辺は節点の組で同じ自由度になり、面の上のトレース (e_i ℓ_j) は節点だけで決まるので
接線連続 (向きの反転は e_i(−ξ) = −e_{p−1−i} で、全体の向きの符号に入る)。
*/
static int HXp, HXnen;
static int HXI[ELEM_NMAX][3];			// 局所節点の格子位置
static int HXinv[ELEM_NMAX];			// 格子位置 (i + (p+1)(j + (p+1)k)) -> 局所節点
static int HXd[PR_NAMAX], HXa[PR_NAMAX], HXb[PR_NAMAX];		// 局所の部分辺の方向・始点・終点 (局所節点)


static int pr_hex_init(void)
{
	int ix[3];
	const int p = hex_lattice(0, ix);
	if (p < 2) return 1;
	HXp = p;
	HXnen = HexNen;
	const int q = p + 1;
	for (int i = 0; i < q * q * q; i++) HXinv[i] = -1;
	for (int a = 0; a < HexNen; a++) {
		if (hex_lattice(a, HXI[a]) != p) return 1;
		const int id = HXI[a][0] + (q * (HXI[a][1] + (q * HXI[a][2])));
		if ((id < 0) || (id >= q * q * q) || (HXinv[id] >= 0)) return 1;
		HXinv[id] = a;
	}
	int k = 0;
	for (int d = 0; d < 3; d++) {
		for (int a = 0; a < HexNen; a++) {
			if (HXI[a][d] >= p) continue;
			int jx[3] = {HXI[a][0], HXI[a][1], HXI[a][2]};
			jx[d]++;
			HXd[k] = d;
			HXa[k] = a;
			HXb[k] = HXinv[jx[0] + (q * (jx[1] + (q * jx[2])))];
			k++;
		}
	}
	return (k != (3 * p * q * q));
}


// 等間隔の p 次 Lagrange ℓ_m の値・1 階・2 階微分 (節点 −1 + 2m/p)
static void pr_lag(int p, int m, double x, double *v, double *d1, double *d2)
{
	double val = 1, a1 = 0, a2 = 0;
	const double xm = -1 + ((2.0 * m) / p);
	for (int j = 0; j <= p; j++) {
		if (j == m) continue;
		const double xj = -1 + ((2.0 * j) / p);
		const double f = (x - xj) / (xm - xj), fd = 1 / (xm - xj);
		a2 = (a2 * f) + (2 * a1 * fd);
		a1 = (a1 * f) + (val * fd);
		val *= f;
	}
	*v = val;
	*d1 = a1;
	*d2 = a2;
}


// 高次の六面体の A の基底 (値 wv と回転 wc、pr_dofs の並び)。nv / g は積分点の形状関数の値と物理空間の勾配
static void pr_hex_basis_s(const signed char *sg, const double nv[ELEM_NMAX], const double g[ELEM_NMAX][3], double (*wv)[3],
	double (*wc)[3])
{
	const int p = HXp, nen = HXnen, nA = 3 * HXp * (HXp + 1) * (HXp + 1);
	double xi[3] = {0, 0, 0}, gx[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
	for (int a = 0; a < nen; a++) {
		for (int c = 0; c < 3; c++) {
			const double s = -1 + ((2.0 * HXI[a][c]) / p);
			xi[c] += s * nv[a];
			for (int q = 0; q < 3; q++) gx[c][q] += s * g[a][q];
		}
	}
	double lv[3][PR_RMAX + 1], ld[3][PR_RMAX + 1], ev[3][PR_RMAX + 1], ed[3][PR_RMAX + 1];
	for (int c = 0; c < 3; c++) {
		double sv = 0, sd = 0;
		for (int m = 0; m <= p; m++) {
			double d2;
			pr_lag(p, m, xi[c], &lv[c][m], &ld[c][m], &d2);
			sv -= ld[c][m];
			sd -= d2;
			ev[c][m] = sv;
			ed[c][m] = sd;
		}
	}
	for (int k = 0; k < nA; k++) {
		const int d = HXd[k], *ix = HXI[HXa[k]];
		double f = sg[k], gf[3] = {0, 0, 0};
		double fc[3], dc[3];
		for (int c = 0; c < 3; c++) {
			fc[c] = ((c == d) ? ev[c][ix[c]] : lv[c][ix[c]]);
			dc[c] = ((c == d) ? ed[c][ix[c]] : ld[c][ix[c]]);
		}
		f *= fc[0] * fc[1] * fc[2];
		for (int c = 0; c < 3; c++) {
			double t = sg[k] * dc[c];
			for (int o = 0; o < 3; o++) {
				if (o != c) t *= fc[o];
			}
			for (int q = 0; q < 3; q++) gf[q] += t * gx[c][q];
		}
		for (int q = 0; q < 3; q++) wv[k][q] = f * gx[d][q];
		cross3(gf, gx[d], wc[k]);
	}
}


static void pr_hex_basis(int e, const double nv[ELEM_NMAX], const double g[ELEM_NMAX][3], double (*wv)[3], double (*wc)[3])
{
	pr_hex_basis_s(&PR.hsg[(int64_t)e * PR.nA], nv, g, wv, wc);
}


static int cmp_i64(const void *a, const void *b)
{
	const int64_t x = *(const int64_t *)a, y = *(const int64_t *)b;
	return ((x < y) ? -1 : ((x > y) ? 1 : 0));
}


// 高次の六面体の部分辺の番号付け (全要素の部分辺を節点の組で重複なく並べる)
static int pr_hex_number(void)
{
	const int nt = PR.nt, nA = PR.nA;
	const int64_t m = (int64_t)nt * nA, nn = NNode;
	int64_t *key = (int64_t *)malloc((size_t)m * sizeof(int64_t));
	PR.hse = (int32_t *)malloc((size_t)m * sizeof(int32_t));
	PR.hsg = (signed char *)malloc((size_t)m * sizeof(signed char));
	for (int e = 0; e < nt; e++) {
		int32_t nd[ELEM_NMAX];
		hex_nodes(e, nd);
		for (int k = 0; k < nA; k++) {
			const int64_t a = nd[HXa[k]], b = nd[HXb[k]];
			key[((int64_t)e * nA) + k] = ((a < b) ? ((a * nn) + b) : ((b * nn) + a));
			PR.hsg[((int64_t)e * nA) + k] = (signed char)((a < b) ? 1 : -1);
		}
	}
	int64_t *srt = (int64_t *)malloc((size_t)m * sizeof(int64_t));
	memcpy(srt, key, (size_t)m * sizeof(int64_t));
	qsort(srt, (size_t)m, sizeof(int64_t), cmp_i64);
	int64_t u = 0;
	for (int64_t i = 0; i < m; i++) {
		if ((i == 0) || (srt[i] != srt[i - 1])) srt[u++] = srt[i];
	}
	if (u >= INT32_MAX) {
		free(key); free(srt);
		return 1;
	}
	PR.ne = (int)u;
	PR.hsa = (int32_t *)malloc((size_t)u * sizeof(int32_t));
	PR.hsb = (int32_t *)malloc((size_t)u * sizeof(int32_t));
	for (int64_t i = 0; i < u; i++) {
		PR.hsa[i] = (int32_t)(srt[i] / nn);
		PR.hsb[i] = (int32_t)(srt[i] % nn);
	}
	for (int64_t i = 0; i < m; i++) {
		const int64_t *f = (const int64_t *)bsearch(&key[i], srt, (size_t)u, sizeof(int64_t), cmp_i64);
		PR.hse[i] = (int32_t)(f - srt);
	}
	free(key); free(srt);
	return 0;
}

/*
---- 高次の角柱 (18 / 40 / 75 / 126 節点、完全 Lagrange、等パラメトリック) ----

角柱の第 1 種 Nedelec 要素 (次数 p) は「三角形の次数 p の Nedelec x ζ の p 次」と「三角形の p 次 Lagrange x ζ の p−1 次 (∇ζ 向き)」
の和 (次元 3p(p+1)(p+2)/2)。前者は格子の層 k (ζ = ζ_k) ごとに三角形の AFW の基底 (四面体の面と同じ λ^α φ_σ、辺 p 個・面 p(p−1) 個)
に ℓ_k(ζ) を掛け、後者は六面体と同じ histopolation の e_k(ζ) に三角形の Lagrange L_T(λ) を掛けた ζ 方向の部分辺の関数にする:

	W = ℓ_k(ζ) λ^α φ_σ,    ∇×W = ℓ_k ∇×(λ^α φ_σ) + ℓ_k' ∇ζ × λ^α φ_σ
	W = e_k(ζ) L_T(λ) ∇ζ,  ∇×W = e_k ∇L_T × ∇ζ

層の辺 (頂点 a, b) の関数は層 k の 2 つの「頂点の列」の節点で、層の面は 3 つで番号を付け、σ・α の順は**その節点の全体番号**で
決める (四面体と同じく向きの符号は要らない)。λ・ζ は節点の値の補間で厳密に引き戻す。∇N は空間に厳密に入る (ζ 方向は部分辺の
接続行列、面内は三角形の Lagrange の勾配) が面内は ±1 の接続行列でないので、ゲージ固定と直接解法は弾く (gauge = 0 の反復解法)。
*/
static void pr_fun(const double lam[4], const double gl[4][3], const double pw[4][PR_RMAX + 1],
	int a, int b, const int al[4], double w[3], double c[3]);
static int PXp, PXnen, PXnt;			// 次数、節点数、三角形の格子点の数 (p+1)(p+2)/2
static int PXA[ELEM_NMAX][4];			// 局所節点 -> (α0, α1, α2, k)
static int PXinv[64 * 6];				// (三角形の格子点 T, 層 k) -> 局所節点 (T + PXnt k)
static int PXT[64][3];					// 三角形の格子点 T -> α
static const int PX_TE[3][2] = {{0, 1}, {1, 2}, {0, 2}};	// 三角形の辺 (PR_FACE と同じく頂点の昇順)


static int pr_px_tidx(const int al[3])
{
	for (int t = 0; t < PXnt; t++) {
		if ((PXT[t][0] == al[0]) && (PXT[t][1] == al[1]) && (PXT[t][2] == al[2])) return t;
	}
	return -1;
}


static int pr_prism_init(void)
{
	int al[4];
	const int p = prism_lattice(0, al);
	if (p < 2) return 1;
	PXp = p;
	PXnen = PrismNen;
	PXnt = ((p + 1) * (p + 2)) / 2;
	int t = 0;
	for (int a0 = p; a0 >= 0; a0--) {
		for (int a1 = p - a0; a1 >= 0; a1--) {
			PXT[t][0] = a0;
			PXT[t][1] = a1;
			PXT[t][2] = p - a0 - a1;
			t++;
		}
	}
	for (int i = 0; i < PXnt * (p + 1); i++) PXinv[i] = -1;
	for (int a = 0; a < PrismNen; a++) {
		if (prism_lattice(a, PXA[a]) != p) return 1;
		const int ti = pr_px_tidx(PXA[a]);
		if ((ti < 0) || (PXA[a][3] < 0) || (PXA[a][3] > p)) return 1;
		const int id = ti + (PXnt * PXA[a][3]);
		if (PXinv[id] >= 0) return 1;
		PXinv[id] = a;
	}
	return 0;
}


// 三角形の格子点 T の頂点 v (α = p e_v) の局所の添字
static int pr_px_vert(int v)
{
	const int al[3] = {(v == 0) ? PXp : 0, (v == 1) ? PXp : 0, (v == 2) ? PXp : 0};
	return pr_px_tidx(al);
}


// 三角形の p 次 Lagrange L_α(λ) の値と勾配 (λ の勾配 gl から)
static double pr_tri_lag(int p, const int al[3], const double lam[3], const double gl[4][3], double gv[3])
{
	double g[3], gd[3];
	for (int m = 0; m < 3; m++) {
		double val = 1, der = 0;
		for (int s = 0; s < al[m]; s++) {
			const double f = ((p * lam[m]) - s) / (s + 1);
			der = (der * f) + (val * p / (s + 1));
			val *= f;
		}
		g[m] = val;
		gd[m] = der;
	}
	for (int q = 0; q < 3; q++) {
		gv[q] = (gd[0] * g[1] * g[2] * gl[0][q]) + (g[0] * gd[1] * g[2] * gl[1][q]) + (g[0] * g[1] * gd[2] * gl[2][q]);
	}
	return g[0] * g[1] * g[2];
}


// 角柱の部分辺・層の辺・層の面の番号付け
typedef struct { int32_t v[3]; } pr_key3;

static int cmp_key3(const void *a, const void *b)
{
	const pr_key3 *x = (const pr_key3 *)a, *y = (const pr_key3 *)b;
	for (int i = 0; i < 3; i++) {
		if (x->v[i] != y->v[i]) return ((x->v[i] < y->v[i]) ? -1 : 1);
	}
	return 0;
}

// 3 つ (辺なら 3 つ目を −1) の鍵の表を作り、要素ごとの全体番号 idx に詰める。戻り値は表の大きさ (失敗は −1)
static int pr_key_table(pr_key3 *key, int64_t m, int32_t *idx)
{
	pr_key3 *srt = (pr_key3 *)malloc((size_t)m * sizeof(pr_key3));
	memcpy(srt, key, (size_t)m * sizeof(pr_key3));
	qsort(srt, (size_t)m, sizeof(pr_key3), cmp_key3);
	int64_t u = 0;
	for (int64_t i = 0; i < m; i++) {
		if ((i == 0) || cmp_key3(&srt[i], &srt[i - 1])) srt[u++] = srt[i];
	}
	for (int64_t i = 0; i < m; i++) {
		const pr_key3 *f = (const pr_key3 *)bsearch(&key[i], srt, (size_t)u, sizeof(pr_key3), cmp_key3);
		idx[i] = (int32_t)(f - srt);
	}
	free(srt);
	return ((u < INT32_MAX) ? (int)u : -1);
}

static void sort3(int32_t *v, int m)
{
	for (int x = 0; x < m - 1; x++) {
		for (int y = x + 1; y < m; y++) {
			if (v[y] < v[x]) { const int32_t t = v[x]; v[x] = v[y]; v[y] = t; }
		}
	}
}

static int pr_prism_number(void)
{
	const int nt = PR.nt, p = PXp, nzl = PXnt * p, ne3 = 3 * (p + 1), nf1 = p + 1;
	const int64_t mz = (int64_t)nt * nzl, me = (int64_t)nt * ne3, mf = (int64_t)nt * nf1, nn = NNode;
	int64_t *zk = (int64_t *)malloc((size_t)mz * sizeof(int64_t));
	pr_key3 *ek = (pr_key3 *)malloc((size_t)me * sizeof(pr_key3));
	pr_key3 *fk = (pr_key3 *)malloc((size_t)mf * sizeof(pr_key3));
	PR.hse = (int32_t *)malloc((size_t)mz * sizeof(int32_t));
	PR.hsg = (signed char *)malloc((size_t)mz * sizeof(signed char));
	PR.pe = (int32_t *)malloc((size_t)me * sizeof(int32_t));
	PR.pf = (int32_t *)malloc((size_t)mf * sizeof(int32_t));
	for (int e = 0; e < nt; e++) {
		int32_t nd[ELEM_NMAX];
		prism_nodes(e, nd);
		for (int t = 0; t < PXnt; t++) {
			for (int k = 0; k < p; k++) {
				const int64_t a = nd[PXinv[t + (PXnt * k)]], b = nd[PXinv[t + (PXnt * (k + 1))]];
				const int64_t o = ((int64_t)e * nzl) + (t * p) + k;
				zk[o] = ((a < b) ? ((a * nn) + b) : ((b * nn) + a));
				PR.hsg[o] = (signed char)((a < b) ? 1 : -1);
			}
		}
		for (int k = 0; k <= p; k++) {
			int32_t lv[3];
			for (int v = 0; v < 3; v++) lv[v] = nd[PXinv[pr_px_vert(v) + (PXnt * k)]];
			for (int l = 0; l < 3; l++) {
				pr_key3 *q = &ek[((int64_t)e * ne3) + (k * 3) + l];
				q->v[0] = lv[PX_TE[l][0]];
				q->v[1] = lv[PX_TE[l][1]];
				q->v[2] = -1;
				sort3(q->v, 2);
			}
			pr_key3 *q = &fk[((int64_t)e * nf1) + k];
			for (int v = 0; v < 3; v++) q->v[v] = lv[v];
			sort3(q->v, 3);
		}
	}
	// ζ 方向の部分辺 (六面体と同じ int64 の鍵)
	int64_t *srt = (int64_t *)malloc((size_t)mz * sizeof(int64_t));
	memcpy(srt, zk, (size_t)mz * sizeof(int64_t));
	qsort(srt, (size_t)mz, sizeof(int64_t), cmp_i64);
	int64_t u = 0;
	for (int64_t i = 0; i < mz; i++) {
		if ((i == 0) || (srt[i] != srt[i - 1])) srt[u++] = srt[i];
	}
	PR.nz = (int)u;
	PR.hsa = (int32_t *)malloc((size_t)u * sizeof(int32_t));
	PR.hsb = (int32_t *)malloc((size_t)u * sizeof(int32_t));
	for (int64_t i = 0; i < u; i++) {
		PR.hsa[i] = (int32_t)(srt[i] / nn);
		PR.hsb[i] = (int32_t)(srt[i] % nn);
	}
	for (int64_t i = 0; i < mz; i++) {
		const int64_t *f = (const int64_t *)bsearch(&zk[i], srt, (size_t)u, sizeof(int64_t), cmp_i64);
		PR.hse[i] = (int32_t)(f - srt);
	}
	PR.ne = pr_key_table(ek, me, PR.pe);
	PR.nf = pr_key_table(fk, mf, PR.pf);
	free(zk); free(ek); free(fk); free(srt);
	return ((u >= INT32_MAX) || (PR.ne < 0) || (PR.nf < 0));
}


// 角柱の A の基底 (pr_dofs の並び : ζ 方向の部分辺 [T][k]、層の辺 [k][l][t]、層の面 [k][t])
static void pr_prism_basis_s(int e, const signed char *sg, const double nv[ELEM_NMAX], const double g[ELEM_NMAX][3],
	double (*wv)[3], double (*wc)[3])
{
	const int p = PXp, nen = PXnen;
	double lam[4] = {0, 0, 0, 0}, gl[4][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, ze = 0, gz[3] = {0, 0, 0};
	for (int a = 0; a < nen; a++) {
		for (int i = 0; i < 3; i++) {
			if (PXA[a][i] == 0) continue;
			const double s = (double)PXA[a][i] / p;
			lam[i] += s * nv[a];
			for (int q = 0; q < 3; q++) gl[i][q] += s * g[a][q];
		}
		const double sz = -1 + ((2.0 * PXA[a][3]) / p);
		ze += sz * nv[a];
		for (int q = 0; q < 3; q++) gz[q] += sz * g[a][q];
	}
	double lv[PR_RMAX + 1], ld[PR_RMAX + 1], ev[PR_RMAX + 1];
	{
		double sv = 0;
		for (int m = 0; m <= p; m++) {
			double d2;
			pr_lag(p, m, ze, &lv[m], &ld[m], &d2);
			sv -= ld[m];
			ev[m] = sv;
		}
	}
	double pw[4][PR_RMAX + 1];
	for (int i = 0; i < 4; i++) {
		pw[i][0] = 1;
		for (int k = 1; k <= p; k++) pw[i][k] = pw[i][k - 1] * lam[i];
	}
	int32_t nd[ELEM_NMAX];
	prism_nodes(e, nd);
	int k = 0;
	for (int t = 0; t < PXnt; t++) {
		double gL[3];
		const double L = pr_tri_lag(p, PXT[t], lam, gl, gL);
		double cr[3];
		cross3(gL, gz, cr);
		for (int z = 0; z < p; z++) {
			const double s = sg[(t * p) + z];
			for (int q = 0; q < 3; q++) {
				wv[k][q] = s * ev[z] * L * gz[q];
				wc[k][q] = s * ev[z] * cr[q];
			}
			k++;
		}
	}
	for (int z = 0; z <= p; z++) {
		int32_t lnd[4];
		for (int v = 0; v < 3; v++) lnd[v] = nd[PXinv[pr_px_vert(v) + (PXnt * z)]];
		lnd[3] = INT32_MAX;
		for (int pass = 0; pass < 2; pass++) {
			const int nent = ((pass == 0) ? 3 : 1);
			for (int l = 0; l < nent; l++) {
				int v[3] = {0, 1, 2};
				int m = 3;
				if (pass == 0) {
					v[0] = PX_TE[l][0];
					v[1] = PX_TE[l][1];
					m = 2;
				}
				pr_sort(lnd, v, m);
				const prtpl_t *tp = ((pass == 0) ? TplE : TplF);
				const int nt0 = ((pass == 0) ? NtE : NtF);
				for (int t = 0; t < nt0; t++) {
					int al[4] = {0, 0, 0, 0};
					for (int i = 0; i < m; i++) al[v[i]] = tp[t].al[i];
					double w[3], c[3], gw[3];
					pr_fun(lam, gl, pw, v[tp[t].a], v[tp[t].b], al, w, c);
					cross3(gz, w, gw);
					for (int q = 0; q < 3; q++) {
						wv[k][q] = lv[z] * w[q];
						wc[k][q] = (lv[z] * c[q]) + (ld[z] * gw[q]);
					}
					k++;
				}
			}
		}
	}
}

#define PR_AAM (5)		// 非線形の Anderson 加速の深さ (2 次の P2_AAM と同じ)


static void pr_prism_basis(int e, const double nv[ELEM_NMAX], const double g[ELEM_NMAX][3], double (*wv)[3], double (*wc)[3])
{
	pr_prism_basis_s(e, &PR.hsg[(int64_t)e * PXnt * PXp], nv, g, wv, wc);
}

static double PXC[PR_RMAX][PR_RMAX];	// 辺の AFW の関数 -> 部分線分の双対の基底 (列 m が部分線分 m)
static int32_t *MXof;					// [nt + 1] 要素の A の自由度の位置
static int32_t *MXdof;
static signed char *MXsg;				// 部分辺の向きの符号 (面・内部の関数は +1)
static int MXnA[3], MXnen[3];			// 種別 (四面体 / 六面体 / 角柱) ごとの A の自由度の数と節点の数
static int TXinv[(PR_RMAX + 1) * (PR_RMAX + 1) * (PR_RMAX + 1)];	// 四面体の多重指数 -> 局所節点
static int pr_mix_init(int p);
static int pr_mix_number(FILE *fp_log);

// 積分点 (種別の混在では節点の数・積分点の数を最大に揃えて 0 で詰める : 足りない積分点は重み 0)
static int pr_ip(int e, int q, double nv[ELEM_NMAX], double g[ELEM_NMAX][3], double *w)
{
	if (!PR.mix) return elem3d_ip(e, q, nv, g, w);
	for (int a = 0; a < PR.nP; a++) {
		nv[a] = 0;
		g[a][0] = g[a][1] = g[a][2] = 0;
	}
	if (q >= elem3d_nip(e)) {
		*w = 0;
		return 0;
	}
	return elem3d_ip(e, q, nv, g, w);
}

// 積分点の数 (種別の混在では最大)
static int pr_nq(void)
{
	if (!PR.mix) return elem3d_nip(0);
	int m = 0;
	for (int e = 0; e < PR.nt; e++) {
		const int q = elem3d_nip(e);
		if (q > m) m = q;
	}
	return m;
}

static int pr_mix_dofs(int e, int32_t dof[PR_NDMAX]);
static void pr_mix_basis(int e, const double nv[ELEM_NMAX], const double g[ELEM_NMAX][3], double (*wv)[3], double (*wc)[3]);


// 要素の自由度の全体番号 (A : 辺・面・内部、φ : 要素の節点)。戻り値は個数
static int pr_dofs(int e, int32_t dof[PR_NDMAX])
{
	if (PR.mix) return pr_mix_dofs(e, dof);
	if (PR.prism) {
		const int p = PXp, nzl = PXnt * p;
		int k = 0;
		for (int i = 0; i < nzl; i++) dof[k++] = PR.hse[((int64_t)e * nzl) + i];
		for (int z = 0; z <= p; z++) {
			for (int l = 0; l < 3; l++) {
				const int32_t id = PR.pe[((int64_t)e * 3 * (p + 1)) + (z * 3) + l];
				for (int t = 0; t < NtE; t++) dof[k++] = PR.nz + (id * NtE) + t;
			}
			const int32_t f = PR.pf[((int64_t)e * (p + 1)) + z];
			for (int t = 0; t < NtF; t++) dof[k++] = PR.oF + (f * NtF) + t;
		}
		int32_t nd[ELEM_NMAX];
		const int nen = prism_nodes(e, nd);
		for (int a = 0; a < nen; a++) dof[k++] = PR.oV + nd[a];
		return k;
	}
	if (PR.hex) {
		const int nA = PR.nA;
		for (int k = 0; k < nA; k++) dof[k] = PR.hse[((int64_t)e * nA) + k];
		int32_t nd[ELEM_NMAX];
		const int nen = hex_nodes(e, nd);
		for (int a = 0; a < nen; a++) dof[nA + a] = PR.oV + nd[a];
		return nA + nen;
	}
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
static void pr_tet_basis(int e, const double nv[ELEM_NMAX], const double g[ELEM_NMAX][3], double (*wv)[3], double (*wc)[3]);

static void pr_basis_core(int e, const double nv[ELEM_NMAX], const double g[ELEM_NMAX][3],
	double (*wv)[3], double (*wc)[3])
{
	if (PR.mix) {
		pr_mix_basis(e, nv, g, wv, wc);
		return;
	}
	if (PR.hex) {
		pr_hex_basis(e, nv, g, wv, wc);
		return;
	}
	if (PR.prism) {
		pr_prism_basis(e, nv, g, wv, wc);
		return;
	}
	pr_tet_basis(e, nv, g, wv, wc);
}


// 四面体の AFW の基底 (辺 r x 6・面 r(r−1) x 4・内部の順)
static void pr_tet_basis(int e, const double nv[ELEM_NMAX], const double g[ELEM_NMAX][3], double (*wv)[3], double (*wc)[3])
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
	free(PR.hse); free(PR.hsg); free(PR.hsa); free(PR.hsb); free(PR.pe); free(PR.pf);
	free(MXof); free(MXdof); free(MXsg);
	MXof = MXdof = NULL;
	MXsg = NULL;
	crs_free(&PR.Kr);
	crs_free(&PR.Ki);
	edge_free();
	memset(&PR, 0, sizeof(PR));
}


// 実体と自由度の番号付け (A と E で共通)
static int pr_number(FILE *fp_log, const char *what)
{
	// 種別の混在した高次の格子 (四面体・六面体・角柱、同じ次数、完全 Lagrange)
	if (MeshMode && (TetOrder >= 2) && (NPyr == 0) && ((((NTet > 0) + (NHex > 0) + (NPrism > 0)) >= 2))) {
		memset(&PR, 0, sizeof(PR));
		const int p = TetOrder;
		if ((p > 5) || pr_init(p) || ((NHex > 0) && pr_hex_init()) || ((NPrism > 0) && pr_prism_init())
			|| ((NHex > 0) && (HXp != p)) || ((NPrism > 0) && (PXp != p)) || pr_mix_init(p)) {
			fprintf(fp_log, "*** the %s on a mixed mesh needs elements of one order (2 to 5) with the full Lagrange nodes "
				"(not the 20-node / 15-node serendipity elements)\n", what);
			return 1;
		}
		PR.mix = PR.r = p;
		PR.nn = NNode;
		PR.nt = NTet + NHex + NPrism;
		MXnA[0] = (6 * p) + (4 * NtF) + NtI;
		MXnA[1] = 3 * p * (p + 1) * (p + 1);
		MXnA[2] = (((p + 1) * (p + 2)) / 2 * p) + ((p + 1) * ((3 * p) + NtF));
		MXnen[0] = ((NTet > 0) ? tet_nen() : 0);
		MXnen[1] = ((NHex > 0) ? HexNen : 0);
		MXnen[2] = ((NPrism > 0) ? PrismNen : 0);
		PR.nA = PR.nP = 0;
		for (int K = 0; K < 3; K++) {
			const int cnt = ((K == 0) ? NTet : (K == 1) ? NHex : NPrism);
			if (cnt == 0) continue;
			if (MXnA[K] > PR.nA) PR.nA = MXnA[K];
			if (MXnen[K] > PR.nP) PR.nP = MXnen[K];
		}
		if (pr_mix_number(fp_log)) {
			fprintf(fp_log, "*** the %s : too many unknowns\n", what);
			return 1;
		}
		return 0;
	}
	// 高次の六面体だけの格子 (完全 Lagrange) : A は部分辺の自由度
	if (MeshMode && (NHex > 0) && ((NTet + NPrism + NPyr) == 0) && (TetOrder >= 2)) {
		memset(&PR, 0, sizeof(PR));
		if (pr_hex_init()) {
			fprintf(fp_log, "*** the %s needs hexahedra with the full Lagrange nodes (27 / 64 / 125 / 216, not the "
				"20-node serendipity element)\n", what);
			return 1;
		}
		PR.hex = PR.r = HXp;
		PR.nn = NNode;
		PR.nt = NHex;
		PR.nA = 3 * HXp * (HXp + 1) * (HXp + 1);
		PR.nP = HexNen;
		if (pr_hex_number()) {
			fprintf(fp_log, "*** the %s : too many sub-edges\n", what);
			return 1;
		}
		PR.oF = PR.oI = PR.oV = PR.ne;
		PR.n = PR.oV + PR.nn;
		return 0;
	}
	// 高次の角柱だけの格子 (完全 Lagrange) : A は ζ 方向の部分辺と層ごとの三角形の AFW の関数
	if (MeshMode && (NPrism > 0) && ((NTet + NHex + NPyr) == 0) && (TetOrder >= 2)) {
		memset(&PR, 0, sizeof(PR));
		if (pr_prism_init() || pr_init(PXp)) {
			fprintf(fp_log, "*** the %s needs prisms with the full Lagrange nodes (18 / 40 / 75 / 126, not the "
				"15-node serendipity element)\n", what);
			return 1;
		}
		PR.prism = PR.r = PXp;
		PR.nn = NNode;
		PR.nt = NPrism;
		PR.nA = (PXnt * PXp) + ((PXp + 1) * ((3 * NtE) + NtF));
		PR.nP = PrismNen;
		if (pr_prism_number()) {
			fprintf(fp_log, "*** the %s : too many sub-edges\n", what);
			return 1;
		}
		PR.oF = PR.nz + (NtE * PR.ne);
		PR.oI = PR.oV = PR.oF + (NtF * PR.nf);
		PR.n = PR.oV + PR.nn;
		return 0;
	}
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
		if (pr_ip(e, q, nv, g, &w)) return 1;
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

/*
ゲージ固定 (gauge = 1)。A-φ 系は (A, φ) → (A + ∇ψ, φ − jωψ) で不変で、ψ は次数 r の節点空間 (電極の節点で 0、awall の上で一定)。
ψ を階層的に「頂点の λ_v + 辺の泡 + 面の泡 + 内部の泡」に分けると、∇ψ の AFW の係数は実体の次元について**ブロック三角**になる
(泡の勾配は低い次元の実体の接線のトレースが 0 なので、その実体の関数の係数も 0)。そこで
  ・頂点 : 1 次・2 次と同じ awall 優先の全域木。木の辺は辺の r 個の自由度を全部 0 (定数のトレース + 泡の勾配)
  ・辺 / 面 / 内部 : その実体の自由度のうち、泡の勾配の「対角ブロック」D (実体の関数 x 泡) の行で正則な小行列になる g 個を 0
    (g = 泡の数 : 辺 r − 1、面 (r − 1)(r − 2)/2、内部 (r − 1)(r − 2)(r − 3)/6)
とすれば ψ が一意に決まり零空間が消える (実体ごとに下から順に解ける)。D は重心座標の代数だけで決まり要素の形に依らないので、
参照の四面体 (頂点の順 = 全体の順) で 1 度だけ最小二乗で作り (残差と、他の実体の係数が 0 になることを確かめる)、行は
ピボット付きの Gram-Schmidt で選ぶ。**電極の上の実体 (泡の節点の φ を固定した辺・面) は ψ の泡も 0 なのでゲージの自由度でなく、
固定しない** (2 次の eddy3d_p2.c と同じ理由)。
*/
static int GTe[PR_TMAX], GTf[PR_TMAX], GTi[PR_TMAX];
static int NGe, NGf, NGi, GTr;


// Householder の QR で最小二乗 min |A x − b| を nb 本の右辺について解く (A は m x n の列優先、上書きする)。戻り値は最大の相対残差
static double pr_lsq(double *A, int m, int n, double *B, int nb, double *X)
{
	double *v = (double *)malloc((size_t)m * sizeof(double));
	for (int k = 0; k < n; k++) {
		double nrm = 0;
		for (int i = k; i < m; i++) nrm += A[((int64_t)k * m) + i] * A[((int64_t)k * m) + i];
		nrm = sqrt(nrm);
		const double akk = A[((int64_t)k * m) + k];
		const double alpha = ((akk > 0) ? -nrm : nrm);
		for (int i = 0; i < m; i++) v[i] = ((i < k) ? 0 : A[((int64_t)k * m) + i]);
		v[k] -= alpha;
		double vn = 0;
		for (int i = k; i < m; i++) vn += v[i] * v[i];
		if (vn <= 0) continue;
		for (int j = k; j < n; j++) {
			double d = 0;
			for (int i = k; i < m; i++) d += v[i] * A[((int64_t)j * m) + i];
			d = 2 * d / vn;
			for (int i = k; i < m; i++) A[((int64_t)j * m) + i] -= d * v[i];
		}
		for (int j = 0; j < nb; j++) {
			double d = 0;
			for (int i = k; i < m; i++) d += v[i] * B[((int64_t)j * m) + i];
			d = 2 * d / vn;
			for (int i = k; i < m; i++) B[((int64_t)j * m) + i] -= d * v[i];
		}
	}
	double rmax = 0;
	for (int j = 0; j < nb; j++) {
		double rr = 0, bb = 0;
		for (int i = n; i < m; i++) rr += B[((int64_t)j * m) + i] * B[((int64_t)j * m) + i];
		for (int i = n - 1; i >= 0; i--) {
			double t = B[((int64_t)j * m) + i];
			for (int c = i + 1; c < n; c++) t -= A[((int64_t)c * m) + i] * X[((int64_t)j * n) + c];
			X[((int64_t)j * n) + i] = t / A[((int64_t)i * m) + i];
			bb += t * t;
		}
		const double rel = sqrt(rr / ((bb + rr) > 0 ? (bb + rr) : 1));
		if (rel > rmax) rmax = rel;
	}
	free(v);
	return rmax;
}


// 参照の四面体での λ^α φ_(a,b) と、泡 Π λ^m の勾配
static void pr_ref_point(const double x[3], double lam[4], double gl[4][3], double pw[4][PR_RMAX + 1], int r)
{
	lam[0] = 1 - x[0] - x[1] - x[2];
	lam[1] = x[0]; lam[2] = x[1]; lam[3] = x[2];
	const double g0[4][3] = {{-1, -1, -1}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
	for (int i = 0; i < 4; i++) {
		for (int q = 0; q < 3; q++) gl[i][q] = g0[i][q];
		pw[i][0] = 1;
		for (int k = 1; k <= r; k++) pw[i][k] = pw[i][k - 1] * lam[i];
	}
}


static void pr_bubble_grad(const double lam[4], const double gl[4][3], const double pw[4][PR_RMAX + 1], const int m[4], double gb[3])
{
	for (int q = 0; q < 3; q++) gb[q] = 0;
	for (int i = 0; i < 4; i++) {
		if (m[i] == 0) continue;
		double d = m[i] * pw[i][m[i] - 1];
		for (int j = 0; j < 4; j++) {
			if (j != i) d *= pw[j][m[j]];
		}
		for (int q = 0; q < 3; q++) gb[q] += d * gl[i][q];
	}
	(void)lam;
}


// 実体 (頂点 nv 個 : 2 辺 / 3 面 / 4 内部) の泡の勾配の対角ブロックを作り、固定する自由度 (テンプレートの番号) を選ぶ
static int pr_gauge_kind(FILE *fp_log, int nv, int *sel, int *nsel)
{
	const int r = PR.r;
	// 列 : 実体の閉包に属する AFW の関数 (辺 / 面 / 内部、参照の頂点 0..nv-1 の上)
	typedef struct { int a, b, al[4], own; } fn_t;
	fn_t *fn = (fn_t *)malloc((size_t)PR_NAMAX * sizeof(fn_t));
	int nf = 0;
	for (int u = 0; u < nv; u++) {
		for (int w = u + 1; w < nv; w++) {
			for (int t = 0; t < NtE; t++) {
				fn_t z = {0, 0, {0, 0, 0, 0}, (nv == 2)};
				const int v2[2] = {u, w};
				z.a = v2[TplE[t].a]; z.b = v2[TplE[t].b];
				z.al[u] = TplE[t].al[0]; z.al[w] = TplE[t].al[1];
				fn[nf++] = z;
			}
		}
	}
	if (nv >= 3) {
		for (int u = 0; u < nv; u++) {
			for (int w = u + 1; w < nv; w++) {
				for (int y = w + 1; y < nv; y++) {
					const int v3[3] = {u, w, y};
					for (int t = 0; t < NtF; t++) {
						fn_t z = {0, 0, {0, 0, 0, 0}, (nv == 3)};
						z.a = v3[TplF[t].a]; z.b = v3[TplF[t].b];
						for (int i = 0; i < 3; i++) z.al[v3[i]] = TplF[t].al[i];
						fn[nf++] = z;
					}
				}
			}
		}
	}
	if (nv == 4) {
		for (int t = 0; t < NtI; t++) {
			fn_t z = {TplI[t].a, TplI[t].b, {TplI[t].al[0], TplI[t].al[1], TplI[t].al[2], TplI[t].al[3]}, 1};
			fn[nf++] = z;
		}
	}
	// 泡 : λ_0 … λ_(nv-1) x 次数 r − nv の単項式
	int (*bm)[4] = (int (*)[4])malloc((size_t)PR_TMAX * sizeof(int[4]));
	int nb = 0;
	const int dr = r - nv;
	for (int a0 = 0; a0 <= dr; a0++) {
		for (int a1 = 0; a1 <= dr - a0; a1++) {
			for (int a2 = 0; a2 <= ((nv >= 3) ? (dr - a0 - a1) : 0); a2++) {
				const int a3 = ((nv == 4) ? (dr - a0 - a1 - a2) : 0);
				if ((nv == 2) && ((a0 + a1) != dr)) continue;
				if ((nv == 3) && ((a0 + a1 + a2) != dr)) continue;
				int m[4] = {a0 + 1, a1 + 1, ((nv >= 3) ? (a2 + 1) : 0), ((nv == 4) ? (a3 + 1) : 0)};
				for (int i = 0; i < 4; i++) bm[nb][i] = m[i];
				nb++;
			}
		}
	}
	// 標本点 (実体の上) と接線の成分
	const int K = (2 * r) + 3;
	int np = 0;
	double (*pt)[3] = (double (*)[3])malloc((size_t)K * K * K * sizeof(double[3]));
	for (int i = 0; i < K; i++) {
		for (int j = 0; j < ((nv >= 3) ? K : 1); j++) {
			for (int k = 0; k < ((nv == 4) ? K : 1); k++) {
				const double x = (i + 0.37) / K, y = ((nv >= 3) ? (j + 0.21) / K : 0), z = ((nv == 4) ? (k + 0.13) / K : 0);
				if ((x + y + z) >= 1) continue;
				pt[np][0] = x; pt[np][1] = y; pt[np][2] = z;
				np++;
			}
		}
	}
	const int nc = nv - 1;		// 接線の成分の数 (辺 1 : x、面 2 : x y、内部 3)
	const int m = np * nc;
	double *A = (double *)calloc((size_t)m * nf, sizeof(double));
	double *B = (double *)calloc((size_t)m * nb, sizeof(double));
	double *X = (double *)calloc((size_t)nf * nb, sizeof(double));
	for (int p = 0; p < np; p++) {
		double lam[4], gl[4][3], pw[4][PR_RMAX + 1];
		pr_ref_point(pt[p], lam, gl, pw, r);
		for (int c = 0; c < nf; c++) {
			double w[3], cu[3];
			pr_fun(lam, gl, pw, fn[c].a, fn[c].b, fn[c].al, w, cu);
			for (int q = 0; q < nc; q++) A[((int64_t)c * m) + (p * nc) + q] = w[q];
		}
		for (int c = 0; c < nb; c++) {
			double gb[3];
			pr_bubble_grad(lam, gl, pw, bm[c], gb);
			for (int q = 0; q < nc; q++) B[((int64_t)c * m) + (p * nc) + q] = gb[q];
		}
	}
	const double res = pr_lsq(A, m, nf, B, nb, X);
	// 対角ブロック (実体自身の関数の行) と、他の実体の係数の大きさ
	int nown = 0, own[PR_TMAX];
	double offmax = 0, onmax = 0;
	for (int c = 0; c < nf; c++) {
		for (int j = 0; j < nb; j++) {
			const double v = fabs(X[((int64_t)j * nf) + c]);
			if (fn[c].own) { if (v > onmax) onmax = v; }
			else if (v > offmax) offmax = v;
		}
		if (fn[c].own) own[nown++] = c;
	}
	// 行の選択 : ピボット付きの Gram-Schmidt (残りの成分が最大の行を順に)
	double *Q = (double *)calloc((size_t)nb * nb, sizeof(double));
	unsigned char *used = (unsigned char *)calloc((size_t)nown, sizeof(unsigned char));
	double minpiv = 1e300;
	*nsel = 0;
	for (int s = 0; s < nb; s++) {
		int best = -1;
		double bn = -1, bv[PR_TMAX];
		for (int o = 0; o < nown; o++) {
			if (used[o]) continue;
			double v[PR_TMAX];
			for (int j = 0; j < nb; j++) v[j] = X[((int64_t)j * nf) + own[o]];
			for (int t = 0; t < s; t++) {
				double d = 0;
				for (int j = 0; j < nb; j++) d += v[j] * Q[(t * nb) + j];
				for (int j = 0; j < nb; j++) v[j] -= d * Q[(t * nb) + j];
			}
			double nn = 0;
			for (int j = 0; j < nb; j++) nn += v[j] * v[j];
			if (nn > bn) {
				bn = nn;
				best = o;
				for (int j = 0; j < nb; j++) bv[j] = v[j];
			}
		}
		if ((best < 0) || (bn <= 0)) break;
		used[best] = 1;
		bn = sqrt(bn);
		if ((bn / ((onmax > 0) ? onmax : 1)) < minpiv) minpiv = bn / ((onmax > 0) ? onmax : 1);
		for (int j = 0; j < nb; j++) Q[(s * nb) + j] = bv[j] / bn;
		sel[(*nsel)++] = own[best] - (nf - nown);		// 実体の関数の中での番号 (テンプレートの番号)
	}
	// テンプレートの番号の昇順に並べる
	for (int i = 0; i < *nsel; i++) {
		for (int j = i + 1; j < *nsel; j++) {
			if (sel[j] < sel[i]) { const int t = sel[i]; sel[i] = sel[j]; sel[j] = t; }
		}
	}
	const char *nm = ((nv == 2) ? "edge" : (nv == 3) ? "face" : "interior");
	fprintf(fp_log, "  gauge (%s) : %d bubble gradients in %d %s functions, least-squares residual %.1e, other entities %.1e, "
		"smallest pivot %.1e\n", nm, nb, nown, nm, res, offmax / ((onmax > 0) ? onmax : 1), minpiv);
	int bad = ((*nsel != nb) || !(res < 1e-9) || !(offmax <= (1e-8 * onmax)) || !(minpiv > 1e-6));
	free(fn); free(bm); free(pt); free(A); free(B); free(X); free(Q); free(used);
	if (bad) {
		fprintf(fp_log, "*** the gauge of the order-%d Nedelec space : the %s bubble gradients are not resolved\n", r, nm);
		return 1;
	}
	return 0;
}


static int pr_gauge(FILE *fp_log, unsigned char *fix)
{
	const int r = PR.r, ne = PR.ne, nf = PR.nf, nt = PR.nt, nn = PR.nn;
	NGe = NGf = NGi = 0;
	if (pr_gauge_kind(fp_log, 2, GTe, &NGe)) return 1;
	if ((r >= 3) && pr_gauge_kind(fp_log, 3, GTf, &NGf)) return 1;
	if ((r >= 4) && pr_gauge_kind(fp_log, 4, GTi, &NGi)) return 1;
	GTr = r;
	// 電極の上の実体 (泡の節点が電極の節点) と awall の実体
	unsigned char *eel = (unsigned char *)calloc((size_t)ne, sizeof(unsigned char));
	unsigned char *fel = (unsigned char *)calloc((size_t)nf, sizeof(unsigned char));
	for (int e = 0; e < nt; e++) {
		int32_t nd[TET_NMAX];
		const int nen = tet_nodes(e, nd);
		for (int a = 4; a < nen; a++) {
			const int i = nd[a];
			if (!((NodeConductor[i] >= 0) && PR.cond[i] && PR.used[i])) continue;
			const int *ix = tet_ho_index(a);
			int v[4], k = 0;
			for (int q = 0; q < 4; q++) {
				if (ix[q] > 0) v[k++] = nd[q];
			}
			if (k == 2) {
				const int64_t id = edge_id(v[0], v[1]);
				if (id >= 0) eel[id] = 1;
			}
			else if (k == 3) {
				int32_t sv[3] = {v[0], v[1], v[2]};
				for (int x = 0; x < 2; x++) {
					for (int y = x + 1; y < 3; y++) {
						if (sv[y] < sv[x]) { const int32_t tt = sv[x]; sv[x] = sv[y]; sv[y] = tt; }
					}
				}
				const int f = p2_face_id(PR.fvert, nf, sv);
				if (f >= 0) fel[f] = 1;
			}
		}
	}
	unsigned char *wall = (unsigned char *)malloc((size_t)ne * sizeof(unsigned char));
	unsigned char *tree = (unsigned char *)malloc((size_t)ne * sizeof(unsigned char));
	int *comp = (int *)malloc((size_t)nn * sizeof(int));
	for (int i = 0; i < ne; i++) wall[i] = fix[(int64_t)i * NtE];
	const int ntree = edge_tree_gauge(wall, tree, comp);
	int nte = 0, ngr = 0, nfe = 0, nfg = 0, nig = 0, nee = 0;
	for (int i = 0; i < ne; i++) {
		if (wall[i]) continue;
		if (tree[i]) {
			for (int t = 0; t < NtE; t++) fix[((int64_t)i * NtE) + t] = 1;
			nte++;
		}
		else if (eel[i]) nee++;
		else {
			for (int t = 0; t < NGe; t++) fix[((int64_t)i * NtE) + GTe[t]] = 1;
			ngr++;
		}
	}
	for (int f = 0; (f < nf) && (NGf > 0); f++) {
		if (fix[PR.oF + ((int64_t)f * NtF)]) continue;
		if (fel[f]) { nfe++; continue; }
		for (int t = 0; t < NGf; t++) fix[PR.oF + ((int64_t)f * NtF) + GTf[t]] = 1;
		nfg++;
	}
	for (int e = 0; (e < nt) && (NGi > 0); e++) {
		for (int t = 0; t < NGi; t++) fix[PR.oI + ((int64_t)e * NtI) + GTi[t]] = 1;
		nig++;
	}
	// 木から決まる ψ が電極上で 0 になる条件 : 電極 (の頂点) が awall の辺だけで連結
	int bad = 0, ref = -1;
	for (int e = 0; e < nt; e++) {
		int32_t nd[TET_NMAX];
		tet_nodes(e, nd);
		for (int a = 0; a < 4; a++) {
			const int i = nd[a];
			if ((NodeConductor[i] < 0) || !PR.cond[i]) continue;
			if (ref < 0) ref = comp[i];
			else if (comp[i] != ref) bad = 1;
		}
	}
	fprintf(fp_log, "  gauge : %d tree edges (%d beyond the wall, all %d functions), edge gradients on %d edges (%d of %d functions; "
		"%d electrode edges left free), face gradients on %d faces (%d of %d; %d electrode faces free), interior gradients on %d (%d of %d)\n",
		ntree, nte, NtE, ngr, NGe, NtE, nee, nfg, NGf, NtF, nfe, nig, NGi, NtI);
	free(eel); free(fel); free(wall); free(tree); free(comp);
	if (bad) {
		fprintf(fp_log, "*** the electrodes are not connected through awall edges, "
			"so the tree gauge would shift the terminal voltage; "
			"use gauge = 0 (the default), or extend awall so that all electrodes "
			"share one A_t = 0 surface\n");
		return 1;
	}
	return 0;
}


// 高次の六面体の awall : 要素の面の 4 隅が awall の四角形に一致したら、その面に乗る部分辺 (両端が面の上) を固定する。
// 戻り値は固定した面の数 (内部の awall は両側から数える)
static int cmp_quad4(const void *a, const void *b)
{
	const int32_t *x = (const int32_t *)a, *y = (const int32_t *)b;
	for (int i = 0; i < 4; i++) {
		if (x[i] != y[i]) return ((x[i] < y[i]) ? -1 : 1);
	}
	return 0;
}

static void sort4(int32_t *v)
{
	for (int x = 0; x < 3; x++) {
		for (int y = x + 1; y < 4; y++) {
			if (v[y] < v[x]) { const int32_t t = v[x]; v[x] = v[y]; v[y] = t; }
		}
	}
}

static int uf_find_pr(int *par, int x)
{
	while (par[x] != x) {
		par[x] = par[par[x]];
		x = par[x];
	}
	return x;
}

static int pr_hex_wall(unsigned char *fix)
{
	const int nt = PR.nt, nA = PR.nA, p = HXp, q = p + 1;
	int32_t *key = (int32_t *)malloc(((size_t)NQuad + 1) * 4 * sizeof(int32_t));
	int nk = 0;
	for (int t = 0; t < NQuad; t++) {
		int hit = 0;
		for (int w = 0; w < NAWall; w++) {
			if (QuadTag[t] == AWallTag[w]) hit = 1;
		}
		if (!hit) continue;
		for (int c = 0; c < 4; c++) key[(nk * 4) + c] = Quad[(t * 4) + c];
		sort4(&key[nk * 4]);
		nk++;
	}
	qsort(key, (size_t)nk, 4 * sizeof(int32_t), cmp_quad4);
	int nwall = 0;
	for (int e = 0; (e < nt) && (nk > 0); e++) {
		int32_t nd[ELEM_NMAX];
		hex_nodes(e, nd);
		for (int d = 0; d < 3; d++) {
			for (int side = 0; side <= p; side += p) {
				int32_t cv[4];
				int m = 0;
				for (int s1 = 0; s1 <= p; s1 += p) {
					for (int s2 = 0; s2 <= p; s2 += p) {
						int ix[3];
						ix[d] = side;
						ix[(d + 1) % 3] = s1;
						ix[(d + 2) % 3] = s2;
						cv[m++] = nd[HXinv[ix[0] + (q * (ix[1] + (q * ix[2])))]];
					}
				}
				sort4(cv);
				if (bsearch(cv, key, (size_t)nk, 4 * sizeof(int32_t), cmp_quad4) == NULL) continue;
				nwall++;
				for (int k = 0; k < nA; k++) {
					if ((HXd[k] != d) && (HXI[HXa[k]][d] == side)) fix[PR.hse[((int64_t)e * nA) + k]] = 1;
				}
			}
		}
	}
	free(key);
	return nwall;
}


// 高次の角柱の awall : 上下の三角形の面 (層 0 / p) は 3 隅を awall の三角形と、側面 (三角形の辺 x ζ) は 4 隅を awall の四角形と
// 照合し、面に乗る自由度 (三角形の面なら層の辺 3 本と層の面、側面なら全層のその辺と ζ 方向の部分辺のうち辺の上のもの) を固定する
static int pr_prism_wall(unsigned char *fix)
{
	const int nt = PR.nt, p = PXp, nzl = PXnt * p;
	int32_t *qk = (int32_t *)malloc(((size_t)NQuad + 1) * 4 * sizeof(int32_t));
	pr_key3 *tk = (pr_key3 *)malloc(((size_t)NTri + 1) * sizeof(pr_key3));
	int nq = 0, ntr = 0;
	for (int t = 0; t < NQuad; t++) {
		int hit = 0;
		for (int w = 0; w < NAWall; w++) {
			if (QuadTag[t] == AWallTag[w]) hit = 1;
		}
		if (!hit) continue;
		for (int c = 0; c < 4; c++) qk[(nq * 4) + c] = Quad[(t * 4) + c];
		sort4(&qk[nq * 4]);
		nq++;
	}
	for (int t = 0; t < NTri; t++) {
		int hit = 0;
		for (int w = 0; w < NAWall; w++) {
			if (TriTag[t] == AWallTag[w]) hit = 1;
		}
		if (!hit) continue;
		for (int c = 0; c < 3; c++) tk[ntr].v[c] = Tri[(t * 3) + c];
		sort3(tk[ntr].v, 3);
		ntr++;
	}
	qsort(qk, (size_t)nq, 4 * sizeof(int32_t), cmp_quad4);
	qsort(tk, (size_t)ntr, sizeof(pr_key3), cmp_key3);
	int nwall = 0;
	for (int e = 0; e < nt; e++) {
		int32_t nd[ELEM_NMAX];
		prism_nodes(e, nd);
		for (int z = 0; (z <= p) && (ntr > 0); z += p) {
			pr_key3 c;
			for (int v = 0; v < 3; v++) c.v[v] = nd[PXinv[pr_px_vert(v) + (PXnt * z)]];
			sort3(c.v, 3);
			if (bsearch(&c, tk, (size_t)ntr, sizeof(pr_key3), cmp_key3) == NULL) continue;
			nwall++;
			for (int l = 0; l < 3; l++) {
				const int32_t id = PR.pe[((int64_t)e * 3 * (p + 1)) + (z * 3) + l];
				for (int t = 0; t < NtE; t++) fix[PR.nz + (id * NtE) + t] = 1;
			}
			const int32_t f = PR.pf[((int64_t)e * (p + 1)) + z];
			for (int t = 0; t < NtF; t++) fix[PR.oF + (f * NtF) + t] = 1;
		}
		for (int l = 0; (l < 3) && (nq > 0); l++) {
			const int va = pr_px_vert(PX_TE[l][0]), vb = pr_px_vert(PX_TE[l][1]);
			int32_t c[4] = {nd[PXinv[va]], nd[PXinv[vb]], nd[PXinv[va + (PXnt * p)]], nd[PXinv[vb + (PXnt * p)]]};
			sort4(c);
			if (bsearch(c, qk, (size_t)nq, 4 * sizeof(int32_t), cmp_quad4) == NULL) continue;
			nwall++;
			for (int z = 0; z <= p; z++) {
				const int32_t id = PR.pe[((int64_t)e * 3 * (p + 1)) + (z * 3) + l];
				for (int t = 0; t < NtE; t++) fix[PR.nz + (id * NtE) + t] = 1;
			}
			// 辺の上の格子点 (辺の外の頂点の α が 0) の ζ 方向の部分辺
			const int other = 3 - PX_TE[l][0] - PX_TE[l][1];
			for (int t = 0; t < PXnt; t++) {
				if (PXT[t][other] != 0) continue;
				for (int z = 0; z < p; z++) fix[PR.hse[((int64_t)e * nzl) + (t * p) + z]] = 1;
			}
		}
	}
	free(qk); free(tk);
	return nwall;
}


// 高次の六面体のゲージ固定 (gauge = 1)。∇ψ は部分辺の接続行列そのものなので、部分辺のグラフの全域木の辺で A = 0 と置く。
// ψ は電極の節点で 0 (φ を固定した節点。電極ごとの別の値を許さない) で awall の上で一定なので、**電極の節点を全部
// 1 つの根につなぎ、awall の部分辺を先に入れた** Kruskal の木を作る。木から決まる ψ は根で 0・awall の連結成分で一定に
// なって零空間が消え、電極の上の部分辺 (根の中の辺) は木に入らないので固定されない (本物の接線の A)。根に awall の
// 成分が入らなくても成分ごとの定数を木が決めるので、1 次と違い電極が awall でつながっている必要は無い
static int pr_hex_gauge(FILE *fp_log, unsigned char *fix)
{
	const int ne = PR.ne, nn = PR.nn;
	int *par = (int *)malloc((size_t)nn * sizeof(int));
	for (int i = 0; i < nn; i++) par[i] = i;
	int root = -1, nel = 0;
	for (int i = 0; i < nn; i++) {
		if ((NodeConductor[i] < 0) || !PR.cond[i] || !PR.used[i]) continue;
		if (root < 0) root = i;
		else par[i] = root;
		nel++;
	}
	int ntree = 0, nwall = 0;
	for (int pass = 0; pass < 2; pass++) {
		for (int i = 0; i < ne; i++) {
			if ((fix[i] != 0) != (pass == 0)) continue;
			const int a = uf_find_pr(par, PR.hsa[i]), b = uf_find_pr(par, PR.hsb[i]);
			if (a == b) continue;
			par[a] = b;
			if (pass == 0) nwall++;
			else {
				fix[i] = 1;
				ntree++;
			}
		}
	}
	int ncomp = 0;
	for (int i = 0; i < nn; i++) {
		if (PR.used[i] && (uf_find_pr(par, i) == i)) ncomp++;
	}
	fprintf(fp_log, "  gauge : %d tree sub-edges beyond %d awall forest sub-edges, %d electrode nodes as the root, "
		"%d connected component(s)\n", ntree, nwall, nel, ncomp);
	free(par);
	return 0;
}


/*
高次の角柱のゲージ固定 (gauge = 1)。ψ = Σ c_(T,k) L_T(λ) ℓ_k(ζ) の勾配は、ζ 方向の部分辺 (列 T の c の差、接続行列) と、層 k の
三角形の AFW の関数 (層の実体の閉包の c_(·,k) だけに依る、四面体の面と同じブロック三角) に分かれる。そこで節点を三角形の
格子点の種類 (頂点 / 辺の泡 / 面の泡) に分け、電極の節点を 1 つの根・awall の自由度を先に入れた union-find で、下の次元から
順に「新しい c を 1 つずつ決める」自由度だけを 0 にする:
  1. 頂点の列 : 頂点どうしの ζ の部分辺、次に層の辺 (木の辺は p 個を全部 0 にし、その層の辺の泡も同時に決まる)
  2. 辺の泡の列 : ζ の部分辺 (列の中の木)、次に列がまだ決まっていない層の辺の泡の行 (四面体と同じ GTe の行)
  3. 面の泡の列 : ζ の部分辺、次に層の面の泡の行 (GTf)
電極の上の実体は根の中なので固定されない (六面体と同じく電極が awall でつながっている必要は無い)。
*/
static int pr_prism_gauge(FILE *fp_log, unsigned char *fix)
{
	const int nt = PR.nt, nn = PR.nn, p = PXp, nzl = PXnt * p, ne3 = 3 * (p + 1);
	NGe = NGf = NGi = 0;
	if (pr_gauge_kind(fp_log, 2, GTe, &NGe)) return 1;
	if ((p >= 3) && pr_gauge_kind(fp_log, 3, GTf, &NGf)) return 1;
	int *par = (int *)malloc((size_t)nn * sizeof(int));
	signed char *typ = (signed char *)malloc((size_t)nn * sizeof(signed char));
	unsigned char *seenE = (unsigned char *)calloc((size_t)PR.ne + 1, sizeof(unsigned char));
	unsigned char *seenF = (unsigned char *)calloc((size_t)PR.nf + 1, sizeof(unsigned char));
	for (int i = 0; i < nn; i++) { par[i] = i; typ[i] = -1; }
	int root = -1, nel = 0;
	for (int i = 0; i < nn; i++) {
		if ((NodeConductor[i] < 0) || !PR.cond[i] || !PR.used[i]) continue;
		if (root < 0) root = i;
		else par[i] = root;
		nel++;
	}
	// 節点の種類 (三角形の格子点の 0 でない重心座標の数 − 1 : 頂点 0、辺の泡 1、面の泡 2)
	for (int e = 0; e < nt; e++) {
		int32_t nd[ELEM_NMAX];
		prism_nodes(e, nd);
		for (int a = 0; a < PXnen; a++) typ[nd[a]] = (signed char)((PXA[a][0] > 0) + (PXA[a][1] > 0) + (PXA[a][2] > 0) - 1);
	}
	// awall : 固定済みの自由度でつながる節点をまとめる (ζ の部分辺、全部固定の層の辺・層の面)
	int nwz = 0;
	for (int i = 0; i < PR.nz; i++) {
		if (!fix[i]) continue;
		const int a = uf_find_pr(par, PR.hsa[i]), b = uf_find_pr(par, PR.hsb[i]);
		if (a != b) { par[a] = b; nwz++; }
	}
	int32_t lst[ELEM_NMAX];
	for (int e = 0; e < nt; e++) {
		int32_t nd[ELEM_NMAX];
		prism_nodes(e, nd);
		for (int z = 0; z <= p; z++) {
			for (int l = 0; l < 3; l++) {
				const int32_t id = PR.pe[((int64_t)e * ne3) + (z * 3) + l];
				if (!fix[PR.nz + (id * NtE)]) continue;
				const int other = 3 - PX_TE[l][0] - PX_TE[l][1];
				int m = 0;
				for (int t = 0; t < PXnt; t++) {
					if (PXT[t][other] == 0) lst[m++] = nd[PXinv[t + (PXnt * z)]];
				}
				for (int q = 1; q < m; q++) {
					const int a = uf_find_pr(par, lst[0]), b = uf_find_pr(par, lst[q]);
					if (a != b) par[a] = b;
				}
			}
			const int32_t f = PR.pf[((int64_t)e * (p + 1)) + z];
			if ((NtF > 0) && fix[PR.oF + (f * NtF)]) {
				for (int t = 1; t < PXnt; t++) {
					const int a = uf_find_pr(par, nd[PXinv[PXnt * z]]), b = uf_find_pr(par, nd[PXinv[t + (PXnt * z)]]);
					if (a != b) par[a] = b;
				}
			}
		}
	}
	int nzt[3] = {0, 0, 0}, nle = 0, nge = 0, ngf = 0;
	for (int stage = 0; stage < 3; stage++) {
		// ζ の部分辺 (両端がこの種類の節点)
		for (int i = 0; i < PR.nz; i++) {
			if (fix[i] || (typ[PR.hsa[i]] != stage)) continue;
			const int a = uf_find_pr(par, PR.hsa[i]), b = uf_find_pr(par, PR.hsb[i]);
			if (a == b) continue;
			par[a] = b;
			fix[i] = 1;
			nzt[stage]++;
		}
		// 層の実体
		for (int e = 0; e < nt; e++) {
			int32_t nd[ELEM_NMAX];
			prism_nodes(e, nd);
			for (int z = 0; z <= p; z++) {
				const int32_t lv0 = nd[PXinv[pr_px_vert(0) + (PXnt * z)]];
				if (stage < 2) {
					for (int l = 0; l < 3; l++) {
						const int32_t id = PR.pe[((int64_t)e * ne3) + (z * 3) + l];
						if (seenE[id] || fix[PR.nz + (id * NtE)]) continue;
						const int other = 3 - PX_TE[l][0] - PX_TE[l][1];
						const int32_t na = nd[PXinv[pr_px_vert(PX_TE[l][0]) + (PXnt * z)]];
						const int32_t nb = nd[PXinv[pr_px_vert(PX_TE[l][1]) + (PXnt * z)]];
						int m = 0, nfree = 0;
						for (int t = 0; t < PXnt; t++) {
							if ((PXT[t][other] == 0) && (PXT[t][PX_TE[l][0]] > 0) && (PXT[t][PX_TE[l][1]] > 0)) {
								lst[m] = nd[PXinv[t + (PXnt * z)]];
								if (uf_find_pr(par, lst[m]) != uf_find_pr(par, na)) nfree++;
								m++;
							}
						}
						if (stage == 0) {
							// 木の辺 : 両端が別の成分で、泡がどれもまだ決まっていない (自分だけの成分)
							if (uf_find_pr(par, na) == uf_find_pr(par, nb)) continue;
							int ok = 1;
							for (int q = 0; q < m; q++) {
								if (uf_find_pr(par, lst[q]) != lst[q]) ok = 0;
							}
							if (!ok) continue;
							seenE[id] = 1;
							for (int t = 0; t < NtE; t++) fix[PR.nz + (id * NtE) + t] = 1;
							par[uf_find_pr(par, na)] = uf_find_pr(par, nb);
							for (int q = 0; q < m; q++) par[uf_find_pr(par, lst[q])] = uf_find_pr(par, nb);
							nle++;
						}
						else if ((nfree == m) && (m > 0) && (uf_find_pr(par, na) == uf_find_pr(par, nb))) {
							// 泡の行 : この層で列がまとめて決まる (どの泡もまだ頂点の成分に入っていないときだけ)
							int distinct = 1;
							for (int x = 0; x < m; x++) {
								for (int y = x + 1; y < m; y++) {
									if (uf_find_pr(par, lst[x]) == uf_find_pr(par, lst[y])) distinct = 0;
								}
							}
							if (!distinct) continue;
							seenE[id] = 1;
							for (int t = 0; t < NGe; t++) fix[PR.nz + (id * NtE) + GTe[t]] = 1;
							for (int q = 0; q < m; q++) par[uf_find_pr(par, lst[q])] = uf_find_pr(par, na);
							nge++;
						}
					}
				}
				else if (NGf > 0) {
					const int32_t f = PR.pf[((int64_t)e * (p + 1)) + z];
					if (seenF[f] || fix[PR.oF + (f * NtF)]) continue;
					int m = 0, nfree = 0;
					for (int t = 0; t < PXnt; t++) {
						if ((PXT[t][0] > 0) && (PXT[t][1] > 0) && (PXT[t][2] > 0)) {
							lst[m] = nd[PXinv[t + (PXnt * z)]];
							if (uf_find_pr(par, lst[m]) != uf_find_pr(par, lv0)) nfree++;
							m++;
						}
					}
					if ((nfree != m) || (m == 0)) continue;
					int distinct = 1;
					for (int x = 0; x < m; x++) {
						for (int y = x + 1; y < m; y++) {
							if (uf_find_pr(par, lst[x]) == uf_find_pr(par, lst[y])) distinct = 0;
						}
					}
					if (!distinct) continue;
					seenF[f] = 1;
					for (int t = 0; t < NGf; t++) fix[PR.oF + (f * NtF) + GTf[t]] = 1;
					for (int q = 0; q < m; q++) par[uf_find_pr(par, lst[q])] = uf_find_pr(par, lv0);
					ngf++;
				}
			}
		}
	}
	int ncomp = 0;
	for (int i = 0; i < nn; i++) {
		if (PR.used[i] && (uf_find_pr(par, i) == i)) ncomp++;
	}
	fprintf(fp_log, "  gauge : zeta sub-edges %d / %d / %d (vertex / edge / face columns), %d layer-edge tree edges (all %d functions), "
		"edge bubbles on %d layer edges (%d of %d), face bubbles on %d layer faces (%d of %d), %d electrode nodes as the root, "
		"%d awall zeta sub-edges, %d connected component(s)\n", nzt[0], nzt[1], nzt[2], nle, NtE, nge, NGe, NtE, ngf, NGf, NtF,
		nel, nwz, ncomp);
	free(par); free(typ); free(seenE); free(seenF);
	return 0;
}



/*
---- 種別の混在した高次の格子 (四面体・六面体・角柱、同じ次数 p、完全 Lagrange) ----

辺の自由度を**すべて部分辺** (格子の隣り合う節点を結ぶ線分の線積分) にそろえる。六面体はもともと部分辺の histopolation の
基底で、四面体の辺の AFW の関数 λ_A^a λ_B^b φ_AB (a + b = p − 1) と角柱の層の辺の関数は、辺の上で P_{p−1} を張るので、辺の
p 本の部分線分の積分 M_ij = ∫_{i/p}^{(i+1)/p} (1−s)^a_j s^b_j ds (s = λ_B、φ_AB・dx/ds = 1) の逆で双対の基底に付け替える
(面・内部の関数は辺の接線のトレースが 0 なので付け替え不要)。こうすると四角形の面 (六面体 / 角柱の側面) の接線トレースは
どちらも histopolation x Lagrange で同じ基底、三角形の面 (四面体 / 角柱の上下) は AFW の面の関数 (全体の頂点の順) と
付け替えた辺の関数で同じ基底になり、自由度は節点の組 (部分辺)・頂点の 3 つ組 (三角形の面)・要素 (四面体の内部) で共有される。
部分辺の自由度は線積分なので ∇N の係数は接続行列 (±1) で、ゲージは部分辺の木 + 面・内部の泡の行。要素ごとの自由度の数が
違うので、最大の数に 0 の基底と固定した詰め物の未知数で揃える (組み立て・過渡解析・発熱・非線形の経路を変えずに通す)。
*/


static int pr_mix_kind(int e)
{
	return ((e < NTet) ? 0 : (e < (NTet + NHex)) ? 1 : 2);
}


// 小さな正方行列の逆 (Gauss-Jordan、部分ピボット)。戻り値 1 は特異
static int pr_inv_small(int n, double a[PR_RMAX][PR_RMAX], double x[PR_RMAX][PR_RMAX])
{
	double m[PR_RMAX][2 * PR_RMAX];
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			m[i][j] = a[i][j];
			m[i][n + j] = (i == j);
		}
	}
	for (int c = 0; c < n; c++) {
		int piv = c;
		for (int i = c + 1; i < n; i++) {
			if (fabs(m[i][c]) > fabs(m[piv][c])) piv = i;
		}
		if (fabs(m[piv][c]) < 1e-14) return 1;
		for (int j = 0; j < 2 * n; j++) { const double t = m[c][j]; m[c][j] = m[piv][j]; m[piv][j] = t; }
		const double d = m[c][c];
		for (int j = 0; j < 2 * n; j++) m[c][j] /= d;
		for (int i = 0; i < n; i++) {
			if (i == c) continue;
			const double f = m[i][c];
			for (int j = 0; j < 2 * n; j++) m[i][j] -= f * m[c][j];
		}
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) x[i][j] = m[i][n + j];
	}
	return 0;
}


static int pr_mix_init(int p)
{
	// 辺の付け替え : M_ij = ∫_{i/p}^{(i+1)/p} (1−s)^a_j s^b_j ds (8 点 Gauss-Legendre、被積分関数は p − 1 ≤ 6 次)
	static const double gx[8] = {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
		0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363};
	static const double gw[8] = {0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
		0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763};
	double M[PR_RMAX][PR_RMAX];
	for (int i = 0; i < p; i++) {
		for (int j = 0; j < p; j++) {
			double v = 0;
			for (int q = 0; q < 8; q++) {
				const double sq = (i + 0.5 * (1 + gx[q])) / p;
				v += 0.5 * gw[q] / p * pow(1 - sq, TplE[j].al[0]) * pow(sq, TplE[j].al[1]);
			}
			M[i][j] = v;
		}
	}
	if (pr_inv_small(p, M, PXC)) return 1;
	// 四面体の多重指数の表
	if (NTet > 0) {
		const int q = p + 1;
		for (int i = 0; i < q * q * q; i++) TXinv[i] = -1;
		for (int a = 0; a < tet_nen(); a++) {
			const int *ix = tet_ho_index(a);
			TXinv[ix[0] + (q * (ix[1] + (q * ix[2])))] = a;
		}
	}
	return 0;
}


// 四面体の局所の頂点 A から B への辺の上の格子点 i (λ_B = i/p) の局所節点
static int pr_tet_edge_node(int A, int B, int i, int p)
{
	int ix[4] = {0, 0, 0, 0};
	ix[A] = p - i;
	ix[B] = i;
	const int q = p + 1;
	return TXinv[ix[0] + (q * (ix[1] + (q * ix[2])))];
}


// 角柱の層 z の三角形の頂点 A から B への辺の上の格子点 i の局所節点
static int pr_px_edge_node(int A, int B, int i, int z)
{
	int al[3] = {0, 0, 0};
	al[A] = PXp - i;
	al[B] = i;
	return PXinv[pr_px_tidx(al) + (PXnt * z)];
}


/*
要素の A の自由度を局所の順に列挙する。kind[k] : 0 部分辺 (key に節点の組、sg に向き)、1 三角形の面 (k3 に頂点の 3 つ組、
t に関数の番号)、2 四面体の内部 (t)。戻り値は個数
*/
static int pr_mix_local(int e, signed char *kind, int64_t *key, pr_key3 *k3, int *tt, signed char *sg)
{
	const int p = PR.r, K = pr_mix_kind(e);
	const int64_t nn = NNode;
	int32_t nd[ELEM_NMAX];
	elem3d_nodes(e, nd);
	int k = 0;
	#define PR_SUB(na, nb) do { const int64_t a_ = (na), b_ = (nb); kind[k] = 0; \
		key[k] = ((a_ < b_) ? ((a_ * nn) + b_) : ((b_ * nn) + a_)); sg[k] = (signed char)((a_ < b_) ? 1 : -1); k++; } while (0)
	if (K == 0) {
		for (int l = 0; l < 6; l++) {
			int v[2] = {PR_EDGE[l][0], PR_EDGE[l][1]};
			pr_sort(nd, v, 2);
			for (int i = 0; i < p; i++) PR_SUB(nd[pr_tet_edge_node(v[0], v[1], i, p)], nd[pr_tet_edge_node(v[0], v[1], i + 1, p)]);
		}
		for (int f = 0; f < 4; f++) {
			pr_key3 c = {{nd[PR_FACE[f][0]], nd[PR_FACE[f][1]], nd[PR_FACE[f][2]]}};
			sort3(c.v, 3);
			for (int t = 0; t < NtF; t++) { kind[k] = 1; k3[k] = c; tt[k] = t; sg[k] = 1; k++; }
		}
		for (int t = 0; t < NtI; t++) { kind[k] = 2; tt[k] = t; sg[k] = 1; k++; }
	}
	else if (K == 1) {
		const int nAh = 3 * HXp * (HXp + 1) * (HXp + 1);
		for (int q = 0; q < nAh; q++) PR_SUB(nd[HXa[q]], nd[HXb[q]]);
	}
	else {
		for (int t = 0; t < PXnt; t++) {
			for (int z = 0; z < p; z++) PR_SUB(nd[PXinv[t + (PXnt * z)]], nd[PXinv[t + (PXnt * (z + 1))]]);
		}
		for (int z = 0; z <= p; z++) {
			int32_t lnd[4];
			for (int v = 0; v < 3; v++) lnd[v] = nd[PXinv[pr_px_vert(v) + (PXnt * z)]];
			lnd[3] = INT32_MAX;
			for (int l = 0; l < 3; l++) {
				int v[2] = {PX_TE[l][0], PX_TE[l][1]};
				pr_sort(lnd, v, 2);
				for (int i = 0; i < p; i++) PR_SUB(nd[pr_px_edge_node(v[0], v[1], i, z)], nd[pr_px_edge_node(v[0], v[1], i + 1, z)]);
			}
			pr_key3 c = {{lnd[0], lnd[1], lnd[2]}};
			sort3(c.v, 3);
			for (int t = 0; t < NtF; t++) { kind[k] = 1; k3[k] = c; tt[k] = t; sg[k] = 1; k++; }
		}
	}
	#undef PR_SUB
	return k;
}


static int pr_mix_number(FILE *fp_log)
{
	const int nt = PR.nt;
	signed char kind[PR_NAMAX], sg[PR_NAMAX];
	int64_t key[PR_NAMAX];
	pr_key3 k3[PR_NAMAX];
	int tt[PR_NAMAX];
	MXof = (int32_t *)malloc(((size_t)nt + 1) * sizeof(int32_t));
	MXof[0] = 0;
	int64_t ns = 0, nfk = 0;
	for (int e = 0; e < nt; e++) {
		const int m = pr_mix_local(e, kind, key, k3, tt, sg);
		MXof[e + 1] = MXof[e] + m;
		for (int q = 0; q < m; q++) {
			if (kind[q] == 0) ns++;
			if ((kind[q] == 1) && (tt[q] == 0)) nfk++;
		}
	}
	int64_t *sk = (int64_t *)malloc(((size_t)ns + 1) * sizeof(int64_t));
	pr_key3 *fk = (pr_key3 *)malloc(((size_t)nfk + 1) * sizeof(pr_key3));
	ns = nfk = 0;
	for (int e = 0; e < nt; e++) {
		const int m = pr_mix_local(e, kind, key, k3, tt, sg);
		for (int q = 0; q < m; q++) {
			if (kind[q] == 0) sk[ns++] = key[q];
			if ((kind[q] == 1) && (tt[q] == 0)) fk[nfk++] = k3[q];
		}
	}
	qsort(sk, (size_t)ns, sizeof(int64_t), cmp_i64);
	int64_t u = 0;
	for (int64_t i = 0; i < ns; i++) {
		if ((i == 0) || (sk[i] != sk[i - 1])) sk[u++] = sk[i];
	}
	qsort(fk, (size_t)nfk, sizeof(pr_key3), cmp_key3);
	int64_t uf = 0;
	for (int64_t i = 0; i < nfk; i++) {
		if ((i == 0) || cmp_key3(&fk[i], &fk[i - 1])) fk[uf++] = fk[i];
	}
	if ((u + (uf * NtF) + ((int64_t)NTet * NtI) + NNode + 1) >= INT32_MAX) {
		free(sk); free(fk);
		return 1;
	}
	PR.ne = (int)u;
	PR.nf = (int)uf;
	PR.oF = PR.ne;
	PR.oI = PR.oF + (NtF * PR.nf);
	PR.oV = PR.oI + (NtI * NTet);
	PR.dummy = PR.oV + NNode;
	PR.n = PR.dummy + 1;
	PR.hsa = (int32_t *)malloc(((size_t)u + 1) * sizeof(int32_t));
	PR.hsb = (int32_t *)malloc(((size_t)u + 1) * sizeof(int32_t));
	for (int64_t i = 0; i < u; i++) {
		PR.hsa[i] = (int32_t)(sk[i] / NNode);
		PR.hsb[i] = (int32_t)(sk[i] % NNode);
	}
	MXdof = (int32_t *)malloc(((size_t)MXof[nt] + 1) * sizeof(int32_t));
	MXsg = (signed char *)malloc(((size_t)MXof[nt] + 1) * sizeof(signed char));
	for (int e = 0; e < nt; e++) {
		const int m = pr_mix_local(e, kind, key, k3, tt, sg);
		for (int q = 0; q < m; q++) {
			int32_t d = 0;
			if (kind[q] == 0) d = (int32_t)((const int64_t *)bsearch(&key[q], sk, (size_t)u, sizeof(int64_t), cmp_i64) - sk);
			else if (kind[q] == 1) {
				const pr_key3 *f = (const pr_key3 *)bsearch(&k3[q], fk, (size_t)uf, sizeof(pr_key3), cmp_key3);
				d = PR.oF + ((int32_t)(f - fk) * NtF) + tt[q];
			}
			else d = PR.oI + (e * NtI) + tt[q];
			MXdof[MXof[e] + q] = d;
			MXsg[MXof[e] + q] = sg[q];
		}
	}
	free(sk); free(fk);
	fprintf(fp_log, "  order %d mixed mesh : %d tetrahedra, %d hexahedra, %d prisms, nodes = %d, unknowns = %d (A : %d sub-edges + "
		"%d x %d triangle faces + %d x %d tetrahedron interiors, phi : %d nodes, 1 padding)\n", PR.r, NTet, NHex, NPrism, NNode,
		PR.n, PR.ne, NtF, PR.nf, NtI, NTet, NNode);
	return 0;
}


static int pr_mix_dofs(int e, int32_t dof[PR_NDMAX])
{
	const int nA = PR.nA, nP = PR.nP;
	const int m = MXof[e + 1] - MXof[e];
	for (int k = 0; k < nA; k++) dof[k] = ((k < m) ? MXdof[MXof[e] + k] : PR.dummy);
	int32_t nd[ELEM_NMAX];
	const int nen = elem3d_nodes(e, nd);
	for (int a = 0; a < nP; a++) dof[nA + a] = ((a < nen) ? (PR.oV + nd[a]) : PR.dummy);
	return nA + nP;
}


// 辺の AFW の関数 (局所の位置 k0 から p 個) を部分線分の双対の基底に付け替えて向きの符号を掛ける
static void pr_mix_edge(int k0, const signed char *sg, double (*wv)[3], double (*wc)[3])
{
	const int p = PR.r;
	double av[PR_RMAX][3], ac[PR_RMAX][3];
	for (int j = 0; j < p; j++) {
		for (int c = 0; c < 3; c++) { av[j][c] = wv[k0 + j][c]; ac[j][c] = wc[k0 + j][c]; }
	}
	for (int m = 0; m < p; m++) {
		for (int c = 0; c < 3; c++) {
			double v = 0, w = 0;
			for (int j = 0; j < p; j++) { v += av[j][c] * PXC[j][m]; w += ac[j][c] * PXC[j][m]; }
			wv[k0 + m][c] = sg[k0 + m] * v;
			wc[k0 + m][c] = sg[k0 + m] * w;
		}
	}
}


static void pr_mix_basis(int e, const double nv[ELEM_NMAX], const double g[ELEM_NMAX][3], double (*wv)[3], double (*wc)[3])
{
	const int p = PR.r, K = pr_mix_kind(e), m = MXof[e + 1] - MXof[e];
	const signed char *sg = &MXsg[MXof[e]];
	if (K == 0) {
		pr_tet_basis(e, nv, g, wv, wc);
		for (int l = 0; l < 6; l++) pr_mix_edge(l * p, sg, wv, wc);
	}
	else if (K == 1) pr_hex_basis_s(sg, nv, g, wv, wc);
	else {
		pr_prism_basis_s(e - NTet - NHex, sg, nv, g, wv, wc);
		for (int z = 0; z <= p; z++) {
			for (int l = 0; l < 3; l++) pr_mix_edge((PXnt * p) + (z * ((3 * p) + NtF)) + (l * p), sg, wv, wc);
		}
	}
	for (int k = m; k < PR.nA; k++) wv[k][0] = wv[k][1] = wv[k][2] = wc[k][0] = wc[k][1] = wc[k][2] = 0;
}


// 種別の混在の awall : 要素の面 (三角形 / 四角形) の頂点を awall の面と照合し、面に乗る自由度を固定する
static int pr_mix_wall(unsigned char *fix)
{
	const int nt = PR.nt, p = PR.r;
	int32_t *qk = (int32_t *)malloc(((size_t)NQuad + 1) * 4 * sizeof(int32_t));
	pr_key3 *tk = (pr_key3 *)malloc(((size_t)NTri + 1) * sizeof(pr_key3));
	int nq = 0, ntr = 0;
	for (int t = 0; t < NQuad; t++) {
		int hit = 0;
		for (int w = 0; w < NAWall; w++) hit |= (QuadTag[t] == AWallTag[w]);
		if (!hit) continue;
		for (int c = 0; c < 4; c++) qk[(nq * 4) + c] = Quad[(t * 4) + c];
		sort4(&qk[nq * 4]);
		nq++;
	}
	for (int t = 0; t < NTri; t++) {
		int hit = 0;
		for (int w = 0; w < NAWall; w++) hit |= (TriTag[t] == AWallTag[w]);
		if (!hit) continue;
		for (int c = 0; c < 3; c++) tk[ntr].v[c] = Tri[(t * 3) + c];
		sort3(tk[ntr].v, 3);
		ntr++;
	}
	qsort(qk, (size_t)nq, 4 * sizeof(int32_t), cmp_quad4);
	qsort(tk, (size_t)ntr, sizeof(pr_key3), cmp_key3);
	int nwall = 0;
	for (int e = 0; e < nt; e++) {
		const int K = pr_mix_kind(e);
		const int32_t *dof = &MXdof[MXof[e]];
		int32_t nd[ELEM_NMAX];
		elem3d_nodes(e, nd);
		if ((K == 0) && (ntr > 0)) {
			for (int f = 0; f < 4; f++) {
				pr_key3 c = {{nd[PR_FACE[f][0]], nd[PR_FACE[f][1]], nd[PR_FACE[f][2]]}};
				sort3(c.v, 3);
				if (bsearch(&c, tk, (size_t)ntr, sizeof(pr_key3), cmp_key3) == NULL) continue;
				nwall++;
				for (int l = 0; l < 6; l++) {
					int in = 0;
					for (int v = 0; v < 3; v++) in += ((PR_FACE[f][v] == PR_EDGE[l][0]) || (PR_FACE[f][v] == PR_EDGE[l][1]));
					if (in < 2) continue;
					for (int i = 0; i < p; i++) fix[dof[(l * p) + i]] = 1;
				}
				for (int t = 0; t < NtF; t++) fix[dof[(6 * p) + (f * NtF) + t]] = 1;
			}
		}
		else if (K == 1) {
			const int q = HXp + 1, nAh = 3 * HXp * q * q;
			for (int d = 0; (d < 3) && (nq > 0); d++) {
				for (int side = 0; side <= HXp; side += HXp) {
					int32_t cv[4];
					int mm = 0;
					for (int s1 = 0; s1 <= HXp; s1 += HXp) {
						for (int s2 = 0; s2 <= HXp; s2 += HXp) {
							int ix[3];
							ix[d] = side;
							ix[(d + 1) % 3] = s1;
							ix[(d + 2) % 3] = s2;
							cv[mm++] = nd[HXinv[ix[0] + (q * (ix[1] + (q * ix[2])))]];
						}
					}
					sort4(cv);
					if (bsearch(cv, qk, (size_t)nq, 4 * sizeof(int32_t), cmp_quad4) == NULL) continue;
					nwall++;
					for (int k = 0; k < nAh; k++) {
						if ((HXd[k] != d) && (HXI[HXa[k]][d] == side)) fix[dof[k]] = 1;
					}
				}
			}
		}
		else if (K == 2) {
			const int nzl = PXnt * p, blk = (3 * p) + NtF;
			for (int z = 0; (z <= p) && (ntr > 0); z += p) {
				pr_key3 c;
				for (int v = 0; v < 3; v++) c.v[v] = nd[PXinv[pr_px_vert(v) + (PXnt * z)]];
				sort3(c.v, 3);
				if (bsearch(&c, tk, (size_t)ntr, sizeof(pr_key3), cmp_key3) == NULL) continue;
				nwall++;
				for (int k = 0; k < blk; k++) fix[dof[nzl + (z * blk) + k]] = 1;
			}
			for (int l = 0; (l < 3) && (nq > 0); l++) {
				const int va = pr_px_vert(PX_TE[l][0]), vb = pr_px_vert(PX_TE[l][1]);
				int32_t c[4] = {nd[PXinv[va]], nd[PXinv[vb]], nd[PXinv[va + (PXnt * p)]], nd[PXinv[vb + (PXnt * p)]]};
				sort4(c);
				if (bsearch(c, qk, (size_t)nq, 4 * sizeof(int32_t), cmp_quad4) == NULL) continue;
				nwall++;
				for (int z = 0; z <= p; z++) {
					for (int i = 0; i < p; i++) fix[dof[nzl + (z * blk) + (l * p) + i]] = 1;
				}
				const int other = 3 - PX_TE[l][0] - PX_TE[l][1];
				for (int t = 0; t < PXnt; t++) {
					if (PXT[t][other] != 0) continue;
					for (int z = 0; z < p; z++) fix[dof[(t * p) + z]] = 1;
				}
			}
		}
	}
	free(qk); free(tk);
	return nwall;
}


/*
種別の混在のゲージ固定 : 部分辺の自由度は線積分なので ∇N の係数は接続行列。電極の節点を根に、awall の自由度を先にまとめた
union-find で部分辺の全域木を 0 にし、まだ決まっていない三角形の面の泡 (GTf の行) と四面体の内部の泡 (GTi の行) を 0 にする
(角柱の面の泡の列は ζ の部分辺で列ごとにまとまっているので、層の面のどれか 1 つの行で決まる)
*/
// 種別の混在のゲージ : 要素 e の三角形の実体 (pass 0 : 四面体の面 f / 角柱の層 f、pass 1 : 四面体の内部) の泡の節点と、
// その関数の局所の位置 fo・基準の頂点 ref。戻り値は泡の数 (実体が無ければ −1)
static int pr_mix_bubbles(int e, int pass, int f, int *fo, int32_t *ref, int32_t *bub)
{
	const int p = PR.r, K = pr_mix_kind(e);
	if ((K == 1) || ((pass == 1) && (K != 0))) return -1;
	if ((K == 0) && (f >= ((pass == 0) ? 4 : 1))) return -1;
	if ((K == 2) && (f > p)) return -1;
	int32_t nd[ELEM_NMAX];
	const int nen = elem3d_nodes(e, nd);
	int nb = 0;
	if (K == 0) {
		*fo = (pass == 0) ? ((6 * p) + (f * NtF)) : ((6 * p) + (4 * NtF));
		*ref = nd[(pass == 0) ? PR_FACE[f][0] : 0];
		for (int a = 4; a < nen; a++) {
			const int *ix = tet_ho_index(a);
			int nz = 0, inface = 1;
			for (int i = 0; i < 4; i++) nz += (ix[i] > 0);
			if (pass == 0) {
				if (nz != 3) continue;
				for (int i = 0; i < 4; i++) {
					const int onf = ((PR_FACE[f][0] == i) || (PR_FACE[f][1] == i) || (PR_FACE[f][2] == i));
					if ((ix[i] > 0) && !onf) inface = 0;
				}
				if (!inface) continue;
			}
			else if (nz != 4) continue;
			bub[nb++] = nd[a];
		}
	}
	else {
		*fo = (PXnt * p) + (f * ((3 * p) + NtF)) + (3 * p);
		*ref = nd[PXinv[pr_px_vert(0) + (PXnt * f)]];
		for (int t = 0; t < PXnt; t++) {
			if ((PXT[t][0] > 0) && (PXT[t][1] > 0) && (PXT[t][2] > 0)) bub[nb++] = nd[PXinv[t + (PXnt * f)]];
		}
	}
	return nb;
}


static int pr_mix_gauge(FILE *fp_log, unsigned char *fix)
{
	const int nt = PR.nt, nn = PR.nn, p = PR.r;
	NGe = NGf = NGi = 0;
	if ((p >= 3) && pr_gauge_kind(fp_log, 3, GTf, &NGf)) return 1;
	if ((p >= 4) && (NTet > 0) && pr_gauge_kind(fp_log, 4, GTi, &NGi)) return 1;
	int *par = (int *)malloc((size_t)nn * sizeof(int));
	for (int i = 0; i < nn; i++) par[i] = i;
	int root = -1, nel = 0;
	for (int i = 0; i < nn; i++) {
		if ((NodeConductor[i] < 0) || !PR.cond[i] || !PR.used[i]) continue;
		if (root < 0) root = i;
		else par[i] = root;
		nel++;
	}
	int32_t bub[ELEM_NMAX], ref;
	int fo;
	// awall : 固定済みの部分辺の両端と、関数を固定した三角形の面 (awall の面) の泡の節点をまとめる
	for (int i = 0; i < PR.ne; i++) {
		if (!fix[i]) continue;
		const int a = uf_find_pr(par, PR.hsa[i]), b = uf_find_pr(par, PR.hsb[i]);
		if (a != b) par[a] = b;
	}
	for (int e = 0; (e < nt) && (NtF > 0); e++) {
		for (int f = 0; f <= p; f++) {
			const int nb = pr_mix_bubbles(e, 0, f, &fo, &ref, bub);
			if ((nb <= 0) || !fix[MXdof[MXof[e] + fo]]) continue;
			for (int x = 0; x < nb; x++) {
				const int a = uf_find_pr(par, bub[x]), b = uf_find_pr(par, ref);
				if (a != b) par[a] = b;
			}
		}
	}
	int ntree = 0, nfr = 0, nir = 0;
	for (int i = 0; i < PR.ne; i++) {
		if (fix[i]) continue;
		const int a = uf_find_pr(par, PR.hsa[i]), b = uf_find_pr(par, PR.hsb[i]);
		if (a == b) continue;
		par[a] = b;
		fix[i] = 1;
		ntree++;
	}
	// 三角形の面の泡 (四面体の面・角柱の層の面) と四面体の内部の泡 : まだ決まっていない泡だけを行で決める
	for (int pass = 0; pass < 2; pass++) {
		const int ng = ((pass == 0) ? NGf : NGi);
		for (int e = 0; (e < nt) && (ng > 0); e++) {
			for (int f = 0; f <= p; f++) {
				const int nb = pr_mix_bubbles(e, pass, f, &fo, &ref, bub);
				if ((nb <= 0) || fix[MXdof[MXof[e] + fo]]) continue;
				int ok = 1;
				for (int x = 0; (x < nb) && ok; x++) {
					if (uf_find_pr(par, bub[x]) == uf_find_pr(par, ref)) ok = 0;
					for (int y = x + 1; (y < nb) && ok; y++) {
						if (uf_find_pr(par, bub[x]) == uf_find_pr(par, bub[y])) ok = 0;
					}
				}
				if (!ok) continue;
				for (int t = 0; t < ng; t++) fix[MXdof[MXof[e] + fo + ((pass == 0) ? GTf[t] : GTi[t])]] = 1;
				for (int x = 0; x < nb; x++) par[uf_find_pr(par, bub[x])] = uf_find_pr(par, ref);
				if (pass == 0) nfr++;
				else nir++;
			}
		}
	}
	int ncomp = 0;
	for (int i = 0; i < nn; i++) {
		if (PR.used[i] && (uf_find_pr(par, i) == i)) ncomp++;
	}
	fprintf(fp_log, "  gauge : %d tree sub-edges, face bubbles on %d triangle faces (%d of %d), interior bubbles in %d tetrahedra "
		"(%d of %d), %d electrode nodes as the root, %d connected component(s)\n", ntree, nfr, NGf, NtF, nir, NGi, NtI, nel, ncomp);
	free(par);
	return 0;
}


// 要素の中心 (四面体 λ = 1/4、六面体 ξ = 0、角柱 λ = 1/3, ζ = 0) の形状関数と勾配 (種別の混在の場の出力)
static int pr_mix_center(int e, double nv[ELEM_NMAX], double g[ELEM_NMAX][3])
{
	const int K = pr_mix_kind(e);
	for (int a = 0; a < PR.nP; a++) { nv[a] = 0; g[a][0] = g[a][1] = g[a][2] = 0; }
	if (K == 0) {
		const double lc[4] = {0.25, 0.25, 0.25, 0.25};
		double det;
		return tet_ho_point(e, lc, nv, g, &det);
	}
	int nen;
	if (K == 1) {
		if (hex2_grad_center(e - NTet, g, &nen)) return 1;
		for (int a = 0; a < nen; a++) {
			nv[a] = 1;
			for (int c = 0; c < 3; c++) {
				double v, d1, d2;
				pr_lag(HXp, HXI[a][c], 0, &v, &d1, &d2);
				nv[a] *= v;
			}
		}
		return 0;
	}
	if (prism2_grad_center(e - NTet - NHex, g, &nen)) return 1;
	const double l3[3] = {1.0 / 3, 1.0 / 3, 1.0 / 3};
	const double g0[4][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
	for (int a = 0; a < nen; a++) {
		double gv[3], v, d1, d2;
		pr_lag(PXp, PXA[a][3], 0, &v, &d1, &d2);
		nv[a] = pr_tri_lag(PXp, PXA[a], l3, g0, gv) * v;
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
	if (PR.mix) {
		// (pr_mix_number が書いた)
	}
	else if (PR.hex) {
		fprintf(fp_log, "  order %d hexahedra : nodes = %d, hexahedra = %d, unknowns = %d (A : %d sub-edges, %d per element, "
			"phi : %d nodes)\n", PR.r, nn, nt, n, ne, PR.nA, nn);
	}
	else if (PR.prism) {
		fprintf(fp_log, "  order %d prisms : nodes = %d, prisms = %d, unknowns = %d (A : %d zeta sub-edges + %d x %d layer edges + "
			"%d x %d layer faces = %d per element, phi : %d nodes)\n", PR.r, nn, nt, n, PR.nz, NtE, ne, NtF, nf, PR.nA, nn);
	}
	else {
		fprintf(fp_log, "  order %d : nodes = %d, tetrahedra = %d, edges = %d, faces = %d, unknowns = %d "
			"(A : %d x %d edges + %d x %d faces + %d x %d interiors = %d per element, phi : %d nodes)\n",
			PR.r, nn, nt, ne, nf, n, NtE, ne, NtF, nf, NtI, nt, PR.nA, nn);
	}
	pr_crs(&PR.Kr);
	pr_crs(&PR.Ki);
	// mpisolve = 2 : 実部・虚部とも自分の行ブロックだけを持つ
	crs_distribute_pair(fp_log, &PR.Kr, &PR.Ki);
	fprintf(fp_log, "  coupled matrix : %lld nonzeros (%.1f per row)\n", (long long)PR.Kr.nnz, (double)PR.Kr.nnz / n);

	// Dirichlet : awall の三角形の辺と面の A (接線成分 0)、電極の節点の φ (= V)、非導電領域と未使用の φ
	unsigned char *fix = PR.fix = (unsigned char *)calloc((size_t)n, sizeof(unsigned char));
	int nwall = 0;
	if (PR.hex) nwall = pr_hex_wall(fix);
	if (PR.prism) nwall = pr_prism_wall(fix);
	if (PR.mix) {
		nwall = pr_mix_wall(fix);
		fix[PR.dummy] = 1;
	}
	for (int t = 0; (t < NTri) && !PR.hex && !PR.prism && !PR.mix; t++) {
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
		int32_t nd[ELEM_NMAX];
		const int nen = elem3d_nodes(e, nd);
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
	if (GaugeTree && (PR.mix ? pr_mix_gauge(fp_log, fix) : PR.prism ? pr_prism_gauge(fp_log, fix) : PR.hex ? pr_hex_gauge(fp_log, fix)
		: pr_gauge(fp_log, fix))) return 1;
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
		if (PR.mix) {
			if (pr_mix_center(e, nv, g)) continue;
		}
		else if (PR.prism) {
			// 角柱の中心 (λ = 1/3, ζ = 0)
			int nen;
			if (prism2_grad_center(e, g, &nen)) continue;
			const double l3[3] = {1.0 / 3, 1.0 / 3, 1.0 / 3};
			const double g0[4][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
			for (int a = 0; a < nen; a++) {
				double gv[3], v, d1, d2;
				pr_lag(PXp, PXA[a][3], 0, &v, &d1, &d2);
				nv[a] = pr_tri_lag(PXp, PXA[a], l3, g0, gv) * v;
			}
		}
		else if (PR.hex) {
			// 六面体の中心 (ξ = 0) : 形状関数は 1 次元の Lagrange の積
			int nen;
			if (hex2_grad_center(e, g, &nen)) continue;
			for (int a = 0; a < nen; a++) {
				nv[a] = 1;
				for (int c = 0; c < 3; c++) {
					double v, d1, d2;
					pr_lag(HXp, HXI[a][c], 0, &v, &d1, &d2);
					nv[a] *= v;
				}
			}
		}
		else if (tet_ho_point(e, lc, nv, g, &det)) continue;
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


// ν の更新の 1 段 : いまの解 (PR.ur / ui) の B から H を独立変数にした ν = |H|/B(|H|) を積分点ごとに作り、
// Anderson 加速の入力 (ln ν と その変化) を aay / aaf に入れる。戻り値は max |Δν| / ν (周波数領域と交流の発熱で共有)
static double pr_nu_update(const double *nuq, double *aay, double *aaf, double (*wv)[3], double (*wc)[3], int32_t *dof)
{
	const int nt = PR.nt, nq = PR.nq;
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
			if (!pr_ip(e, q, nv, g, &w)) {
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
	return dmax;
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
	const int nq = pr_nq();
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
		const double dmax = pr_nu_update(nuq, aay, aaf, wv, wc, dof);
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
// 交流の発熱の非線形 (B-H) : 積分点ごとの ν と Anderson の作業域 (熱連成の回をまたいで持ち越す。無ければ NULL)
static double *HNuq, *HAay, *HAaf, *HAayp, *HAafp, *HAady, *HAadf;
static int HNlTotal;


int eddy3d_pr_heat_open(FILE *fp_log)
{
	fprintf(fp_log, "  AC heating (analysis = A, order-%d Nedelec) : frequency = %.6e [Hz]\n", TetOrder, Freq);
	if (pr_open(fp_log)) return 1;
	HNuq = NULL;
	HNlTotal = 0;
	int nlf = 0;
	for (int m = 0; m < NMaterial; m++) {
		if (Material[m].nbh[0] > 0) nlf = 1;
	}
	if (nlf) {
		const int nq = pr_nq();
		const size_t nqt = (size_t)PR.nt * (size_t)nq;
		HNuq  = (double *)malloc(nqt * sizeof(double));
		HAay  = (double *)malloc(nqt * sizeof(double));
		HAaf  = (double *)malloc(nqt * sizeof(double));
		HAayp = (double *)malloc(nqt * sizeof(double));
		HAafp = (double *)malloc(nqt * sizeof(double));
		HAady = (double *)malloc(nqt * PR_AAM * sizeof(double));
		HAadf = (double *)malloc(nqt * PR_AAM * sizeof(double));
		for (int e = 0; e < PR.nt; e++) {
			const material_t *mt = &Material[elem3d_mat(e)];
			for (int q = 0; q < nq; q++) HNuq[((int64_t)e * nq) + q] = ((mt->nbh[0] > 0) ? bh_nu(mt, 0, 0) : 0);
		}
		PR.nuq = HNuq;
		PR.nq = nq;
		fprintf(fp_log, "  nonlinear (B-H) : effective reluctivity at %d points per element, solved to convergence in "
			"every thermal iteration (maxiter = %d, tol = %.1e); the thermal tolerance is at least 10 x that\n",
			nq, NlMaxiter, NlTol);
	}
	return 0;
}


int eddy3d_pr_heat_solve(FILE *fp_log, const double *sigelem, double cur, double *f,
	double *zr, double *zi, int it)
{
	const int nt = PR.nt, nn = PR.nn, n = PR.n, nA = PR.nA, nP = PR.nP;
	const double omega = PR.omega;
	double (*wv)[3] = (double (*)[3])malloc((size_t)nA * 3 * sizeof(double));
	double (*wc)[3] = (double (*)[3])malloc((size_t)nA * 3 * sizeof(double));
	int32_t *dof = (int32_t *)malloc(PR_NDMAX * sizeof(int32_t));
	double ir = 0, ii = 0, den = 0;
	int aan = 0, ierr = 0;		// Anderson の履歴は σ(T) が変わるたびに捨てる
	// 非線形 (B-H) なら 1 次の eddy3d_heat_solve と同じく、いまの σ(T) で ν が動かなくなるまで解き直す
	for (int nlit = 1; ; nlit++) {
		if (pr_assemble(fp_log, sigelem)) { ierr = 1; break; }
		char label[64];
		sprintf(label, "A%d", it);
		double qr[MAXPORT], qi[MAXPORT];
		if (pr_solve_port(fp_log, 1, NULL, label, qr, qi)) { ierr = 1; break; }
		ir = qr[1];
		ii = qi[1];
		den = (ir * ir) + (ii * ii);
		if (den <= 0) {
			fprintf(fp_log, "*** thermal : no current flows through the port\n");
			ierr = 1;
			break;
		}
		const double sr = cur * ir / den, si = -cur * ii / den;
		for (int i = 0; i < n; i++) {
			const double xr = PR.ur[i], xi = PR.ui[i];
			PR.ur[i] = (sr * xr) - (si * xi);
			PR.ui[i] = (sr * xi) + (si * xr);
		}
		if (HNuq == NULL) break;
		const double dmax = pr_nu_update(HNuq, HAay, HAaf, wv, wc, dof);
		HNlTotal++;
		if (dmax < NlTol) break;
		if (nlit >= NlMaxiter) {
			fprintf(fp_log, "*** thermal : the effective reluctivity did not converge in %d iterations "
				"(thermal iteration %d, max dnu/nu = %.3e)\n", NlMaxiter, it, dmax);
			ierr = 1;
			break;
		}
		const int64_t nqt = (int64_t)nt * PR.nq;
		anderson_step(nqt, HAay, HAaf, HAayp, HAafp, HAady, HAadf, PR_AAM, &aan, NlRelax);
		for (int e = 0; e < nt; e++) {
			if (Material[elem3d_mat(e)].nbh[0] <= 0) continue;
			for (int q = 0; q < PR.nq; q++) HNuq[((int64_t)e * PR.nq) + q] = exp(HAay[((int64_t)e * PR.nq) + q]);
		}
	}
	if (ierr) {
		free(wv); free(wc); free(dof);
		return 1;
	}
	*zr = Volt * ir / den;
	*zi = -Volt * ii / den;
	if (HNuq != NULL) {
		fprintf(fp_log, "  thermal iteration %d : effective reluctivity converged (%d B-H iterations so far)\n",
			it, HNlTotal);
	}
	for (int i = 0; i < nn; i++) f[i] = 0;
	int bad = 0;
	for (int e = 0; (e < nt) && !bad; e++) {
		const double sg = sigelem[e];
		if (sg <= 0) continue;
		pr_dofs(e, dof);
		int32_t nd[ELEM_NMAX];
		const int nen = elem3d_nodes(e, nd);
		const int nq = elem3d_nip(e);
		for (int q = 0; q < nq; q++) {
			double nv[ELEM_NMAX], g[ELEM_NMAX][3], w;
			if (pr_ip(e, q, nv, g, &w)) { bad = 1; break; }
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
			for (int a = 0; a < nen; a++) f[nd[a]] += qd * nv[a];
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
	free(HNuq); free(HAay); free(HAaf); free(HAayp); free(HAafp); free(HAady); free(HAadf);
	HNuq = HAay = HAaf = HAayp = HAafp = HAady = HAadf = NULL;
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
	if (pr_ip(J->el[j], q, nv, g, &w)) {
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
	const int nq = pr_nq();
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
			if (pr_ip(e, q, nv, g, &w)) { bad = 1; break; }
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
	// 直線の要素だけ (節点が要素のアフィンな写像の位置にあること。四面体は 4 頂点、六面体は角 0 と 3 本の辺 (平行六面体)、
	// 角柱は下の三角形と鉛直の辺 1 本 (上の三角形が平行移動) で決まる写像)
	for (int e = 0; e < nt; e++) {
		int32_t nd[ELEM_NMAX];
		const int nen = elem3d_nodes(e, nd);
		int32_t o = nd[0], v1 = nd[1], v2 = nd[2], v3 = nd[3];
		const int K = (PR.mix ? pr_mix_kind(e) : PR.hex ? 1 : PR.prism ? 2 : 0);
		if (K == 1) {
			const int q = HXp + 1;
			o = nd[HXinv[0]];
			v1 = nd[HXinv[HXp]];
			v2 = nd[HXinv[q * HXp]];
			v3 = nd[HXinv[q * q * HXp]];
		}
		else if (K == 2) {
			o = nd[PXinv[pr_px_vert(0)]];
			v1 = nd[PXinv[pr_px_vert(1)]];
			v2 = nd[PXinv[pr_px_vert(2)]];
			v3 = nd[PXinv[pr_px_vert(0) + (PXnt * PXp)]];
		}
		const double h[3] = {Xp[v1] - Xp[o], Yp[v1] - Yp[o], Zp[v1] - Zp[o]};
		for (int a = 0; a < nen; a++) {
			double c[3];		// アフィンな写像の座標 (o からの v1 / v2 / v3 の向きの重み)
			if (K == 1) {
				for (int k = 0; k < 3; k++) c[k] = (double)HXI[a][k] / HXp;
			}
			else if (K == 2) {
				c[0] = (double)PXA[a][1] / PXp;
				c[1] = (double)PXA[a][2] / PXp;
				c[2] = (double)PXA[a][3] / PXp;
			}
			else {
				if (a < 4) continue;
				const int *ix = tet_ho_index(a);
				c[0] = (double)ix[1] / r;
				c[1] = (double)ix[2] / r;
				c[2] = (double)ix[3] / r;
			}
			const double x[3] = {
				Xp[o] + (c[0] * (Xp[v1] - Xp[o])) + (c[1] * (Xp[v2] - Xp[o])) + (c[2] * (Xp[v3] - Xp[o])),
				Yp[o] + (c[0] * (Yp[v1] - Yp[o])) + (c[1] * (Yp[v2] - Yp[o])) + (c[2] * (Yp[v3] - Yp[o])),
				Zp[o] + (c[0] * (Zp[v1] - Zp[o])) + (c[1] * (Zp[v2] - Zp[o])) + (c[2] * (Zp[v3] - Zp[o]))};
			const double d[3] = {Xp[nd[a]] - x[0], Yp[nd[a]] - x[1], Zp[nd[a]] - x[2]};
			if (dot3(d, d) > (1e-16 * dot3(h, h))) {
				fprintf(fp_log, "*** the high-order edge self test needs straight-sided elements "
					"(element %d has a node off its affine position); the identities are exact "
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
		int32_t nd[ELEM_NMAX];
		const int nen = elem3d_nodes(e, nd);
		const int nq = elem3d_nip(e);
		for (int qq = 0; qq < nq; qq++) {
			double nv[ELEM_NMAX], g[ELEM_NMAX][3], w;
			if (pr_ip(e, qq, nv, g, &w)) { ierr = 1; break; }
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
