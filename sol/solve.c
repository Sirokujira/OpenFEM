/*
solve.c

準静的 FEM の求解と回路パラメータの抽出。

ポート k を Volt [V]、他の導体と基準導体 (id=0) を 0 [V] にして解き、
各導体上の反作用 (Σ (K φ)_i) から電荷 Q_j (静電界) または電流 I_j
(定常電流界) を得る。これを全ポートについて繰り返して行列を作る。

	静電界   : C[k][j] = Q_j / V   短絡容量行列 (Maxwell 行列)
	真空静電 : C0 -> L = μ0 ε0 inv(C0)   (TEM 仮定)
	電流界   : G[k][j] = I_j / V,  R = inv(G)
*/

#include "fem.h"
#include "fem_prototype.h"

// ---- 要素の走査 (M / F の断面 2 次元定式化を格子の種類で分けずに書く) ----
//
// M / F は「導体要素を走査して節点に重みを配る」形の処理が 4 か所あり
// (∫σdV、右辺 ∫σN V'、電流 ∫σAz dV、場の出力)、構造格子と 2 次元三角形で
// 書き分けると必ずどこかが食い違う。要素の情報だけ 1 か所に集める。
//
// **「体積」は 2 次元格子では面積**。単位長あたりで扱う (TlineLength = 1) ので、
// 出てくる R [Ω/m]・L [H/m] は構造格子と同じ意味になる。
typedef struct {
	int nen;					// 節点数 (構造格子 8、1 次三角形 3、2 次三角形 6、3 次三角形 10)
	int64_t nd[TRI_NMAX];		// 節点番号 (TRI_NMAX = 10 >= 8)
	double vol;					// 体積 [m^3] または面積 [m^2]
	// **節点あたりの重み ∫N_i dV。均等とは限らない。**
	// 1 次要素は vol/nen で均等だが、2 次三角形は ∫N dA が
	// 頂点で 0、辺の中間節点で A/3 になる。「vol/nen を配る」と書くと
	// 2 次要素で右辺も端子電流も壊れるので、重みを要素側に持たせる
	double w[TRI_NMAX];
	// 要素中心での形状関数の値 N_i(重心)。場の出力で「中心の値」を作るのに使う。
	// 2 次三角形では頂点が -1/9、中間節点が 4/9 で、単純平均ではない
	// (3 次三角形では頂点・辺が 0、面心が 1)
	double wc[TRI_NMAX];
	int cond;					// 導体番号 (-1 = 導体でない)
	int mat;					// 材料番号
} elem_t;


static int64_t elem_count(void)
{
	return ((MeshMode && (MeshDim == 2)) ? (int64_t)NTri : ((int64_t)Nx * Ny * Nz));
}


static void elem_get(int64_t e, elem_t *el)
{
	if (MeshMode && (MeshDim == 2)) {
		int32_t nd[TRI_NMAX];
		el->nen = tri_nodes((int)e, nd);
		for (int l = 0; l < el->nen; l++) el->nd[l] = nd[l];
		el->vol = TriArea[e];
		if (el->nen >= 10) {
			// 3 次三角形 : ∫N dA を等パラメトリック写像で積分する
			// (直線要素では頂点 A/30、辺 3A/40、面心 9A/20)
			tri_weights((int)e, el->w, el->wc);
		}
		else if (el->nen == 6) {
			// 2 次三角形 : ∫N dA は頂点で 0、辺の中間節点で A/3
			for (int l = 0; l < 3; l++) { el->w[l] = 0;                el->wc[l] = -1.0 / 9; }
			for (int l = 3; l < 6; l++) { el->w[l] = TriArea[e] / 3;   el->wc[l] =  4.0 / 9; }
		}
		else {
			for (int l = 0; l < 3; l++) { el->w[l] = TriArea[e] / 3;   el->wc[l] = 1.0 / 3; }
		}
		el->cond = TriCond[e];
		el->mat = TriMat[e];
		return;
	}

	const int k = (int)(e % Nz);
	const int j = (int)((e / Nz) % Ny);
	const int i = (int)(e / ((int64_t)Ny * Nz));
	el->nen = 8;
	for (int l = 0; l < 8; l++) {
		el->nd[l] = node_index(i + ((l >> 2) & 1), j + ((l >> 1) & 1), k + (l & 1));
	}
	el->vol = (Xn[i + 1] - Xn[i]) * (Yn[j + 1] - Yn[j]) * (Zn[k + 1] - Zn[k]);
	for (int l = 0; l < 8; l++) {
		el->w[l] = el->vol / 8;			// 3 重線形は均等
		el->wc[l] = 1.0 / 8;
	}
	el->cond = CellConductor[e];
	el->mat = CellMaterial[e];
}


// 1 ポート励振の求解 -> 各導体の反作用 q[0..NPort]
int solve_port(const crs_t *A, const unsigned char *fix, int kport,
	double *phi, double *work, double *q, double *energy,
	FILE *fp_log, const char *label)
{
	const int n = (int)A->n;

	// Dirichlet 値 (phi は φ_D として使い、後で u を足す)
	int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (i = 0; i < n; i++) {
		const int id = NodeConductor[i];
		phi[i] = ((id == kport) ? Volt : 0);
	}

	// 右辺 b = -(A φ_D) (自由節点のみ、固定節点は 0)
	crs_spmv(A, phi, work, NULL);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (i = 0; i < n; i++) {
		work[i] = (fix[i] ? 0 : -work[i]);
	}

	// 求解 (u は work2 として phi とは別に持つ必要がある)
	double *u = (double *)malloc(n * sizeof(double));
	const int iter = solver_cg(A, work, u, fix,
		Solver.maxiter, Solver.nout, Solver.converg, fp_log, label);

	// φ = φ_D + u
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (i = 0; i < n; i++) {
		phi[i] += u[i];
	}
	free(u);

	// 反作用 (電荷または電流) と静電エネルギー
	for (int p = 0; p <= NPort; p++) {
		q[p] = 0;
	}
	double w2 = 0;
	for (int i = 0; i < n; i++) {
		const int id = NodeConductor[i];
		if (id < 0) continue;
		const double f = crs_row_dot(A, i, phi);
		q[id] += f;
		w2 += f * phi[i];
	}
	*energy = w2 / 2;

	return iter;
}


// 行列の対称化 (離散化・収束による非対称のみのはず) と非対称度の報告
static void symmetrize(double *m, int np, FILE *fp_log, const char *name)
{
	double amax = 0, dmax = 0;
	for (int i = 0; i < np; i++) {
		for (int j = 0; j < np; j++) {
			const double a = fabs(m[(i * np) + j]);
			if (a > amax) amax = a;
		}
	}
	for (int i = 0; i < np; i++) {
		for (int j = i + 1; j < np; j++) {
			const double d = fabs(m[(i * np) + j] - m[(j * np) + i]);
			if (d > dmax) dmax = d;
			const double a = (m[(i * np) + j] + m[(j * np) + i]) / 2;
			m[(i * np) + j] = m[(j * np) + i] = a;
		}
	}
	if ((fp_log != NULL) && (np > 1) && (amax > 0)) {
		fprintf(fp_log, "  %s matrix asymmetry = %.3e (relative)\n", name, dmax / amax);
	}
}


// 鉄損 (Bertotti の損失分離)
//
// 静磁場解析で得た要素毎の |B| を**正弦波励磁の波高値**とみなし、
// 単位体積あたりの損失を 3 項に分けて積分する:
//
//   P = kh f B^alpha  +  (pi^2 sigma d^2 / 6) f^2 B^2  +  ke (f B)^1.5   [W/m^3]
//        ヒステリシス       古典渦電流 (積層厚 d)          異常 (過剰)
//
// 古典項の係数 pi^2 sigma d^2/6 は、厚さ d の薄板を一様な正弦波磁束が貫くときの
// 渦電流損の閉形式 (d << 表皮深さ)。異常項の指数 1.5 は Bertotti の統計理論による。
//
// **静磁場の解を波高値に使う**ので、渦電流による磁束の押し出し (表皮効果) は
// 入らない。d << 表皮深さ という積層鉄心の前提で使うこと。
//
// 戻り値は単位長あたりの損失 [W/m] (Tline があるとき) または全損失 [W]。
static double iron_loss(const double *az, FILE *fp_log, int kport)
{
	const int64_t ncell = elem_count();
	double *bv = (double *)malloc((size_t)ncell * 3 * sizeof(double));
	field_cell_grad(az, 1, bv);			// B = ∇×(Az ê_t)

	double ph = 0, pc = 0, pe = 0, vsum = 0, bmax = 0;
	for (int64_t c = 0; c < ncell; c++) {
		elem_t el;
		elem_get(c, &el);
		const material_t *mt = &Material[el.mat];
		if ((mt->be_kh <= 0) && (mt->be_ke <= 0) && (mt->be_d <= 0)) continue;

		const double bx = bv[(c * 3)], by = bv[(c * 3) + 1], bz = bv[(c * 3) + 2];
		const double b = sqrt((bx * bx) + (by * by) + (bz * bz));
		if (b > bmax) bmax = b;
		vsum += el.vol;

		if (mt->be_kh > 0) ph += mt->be_kh * Freq * pow(b, mt->be_alpha) * el.vol;
		if ((mt->be_d > 0) && (mt->sigma > 0)) {
			pc += (PI * PI * mt->sigma * mt->be_d * mt->be_d / 6)
			    * Freq * Freq * b * b * el.vol;
		}
		if (mt->be_ke > 0) pe += mt->be_ke * pow(Freq * b, 1.5) * el.vol;
	}
	free(bv);

	const double scale = ((TlineLength > 0) ? (1 / TlineLength) : 1);
	ph *= scale;
	pc *= scale;
	pe *= scale;

	// max|B| は検証で閉形式と突き合わせるので有効桁を落とさない
	fprintf(fp_log, "  iron loss (port %d) : volume = %.6e, max|B| = %.10e [T], "
		"f = %.6e [Hz]\n", kport, vsum * scale, bmax, Freq);
	fprintf(fp_log, "    hysteresis = %.6e, classical = %.6e, excess = %.6e, "
		"total = %.6e [W/m]\n", ph, pc, pe, ph + pc + pe);

	return ph + pc + pe;
}


// 静磁場解析 (断面 2 次元、ベクトルポテンシャル Az 定式化)
//
//   ∇・(ν ∇Az) = -Jz,  ν = 1/(μ0 μr)
//
// ポート k に +I、基準導体 (id=0) に -I を一様電流密度で流して解く。
// 正味電流が 0 なので純 Neumann 系は可解 (定数分の不定性のみ) であり、
// 1 節点を固定して解いた解は真の解を定数だけずらしたものになる。
// 相互エネルギー ∫J^(j)・A^(k) dV = L_kj I_k I_j は定数分に依らない
// (∫J dV = 0) ので、そのまま L 行列が得られる。
// 導体内部の電流分布を含むため内部インダクタンスが自動的に入り、μr にも対応する。
// 材料に B-H 曲線 (bh キー) があるときは ν = H(|B|)/|B| を反復更新する非線形解析に
// なる (単一ポートのみ)。得られる L は与えた電流での磁束鎖交 (割線) インダクタンス。
/*
静磁場の過渡解析 (analysis = M + transient、断面 2 次元の構造格子)。

	-∇・(ν ∇A) + σ (∂A/∂t - c_k(t)) = J_s(t)

J_s は駆動導体 (conductor) の一様な電流密度 I(t)/S (撚り線とみなし、導体の中の渦電流は
入れない)。渦電流は導体でない材料の σ (Material[].sigma) で流れる。断面 2 次元では
**孤立した導電領域ごとに正味の電流が 0** (無限に長い板の端は開いている) なので、
c_k = (1/S_k)∫_k σ ∂A/∂t を消去した射影質量
	M_eff y = M y - Σ_k m_k (m_k・y) / S_k,   m_k = ∫_k σ N,  S_k = ∫_k σ
を使う (M_eff は半正定値で、定数 (ゲージ) を零空間に持つので 1 節点固定のまま解ける)。

時間は後退 Euler。各刻みで
	R(A) = b(t_{n+1}) - r_ja(A) - M_eff (A - A^n) / Δt = 0
を、行列 K(dH/dB) + M_eff/Δt の固定点反復で解く (ヒステリシスの静磁場と同じ構成。
行列は常に正の微分磁気抵抗率、残差だけ正確な H(B))。収束したら J-A の状態を確定させる。

**離散エネルギーの恒等式**: 刻みの式に w = ΔA = A^{n+1} - A^n を掛けると
	b(t_{n+1})・ΔA = ΔA・r_ja(A^{n+1}) + ΔA・M_eff ΔA / Δt
左辺は端子から入れた仕事 I ΔΛ、右辺第 1 項は磁気の仕事 ∫H・ΔB (残差と同じ Gauss 点)、
第 2 項は渦電流損 Δt ∫σ|E|² (E = -(ΔA/Δt - c))。刻みを足し合わせた
	入力 = 磁気の仕事 + 渦電流損
が反復の打ち切りを除いて厳密に成り立つ (ログに出し、rlc_check.sh で見る)。
σ = 0 なら刻みの式は静磁場そのもので、waveform の点と刻みを揃えれば準静的な
電流掃引 (currentsweep) と同じ解の列になる (検証の恒等式)。
*/

// 刻みの終わり t の電流の倍率 f(t^-) (waveform の折れ線の左極限。無ければ 1)
static double mt_wave(double t)
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


static int mt_find(int *par, int x)
{
	while (par[x] != x) {
		par[x] = par[par[x]];
		x = par[x];
	}
	return x;
}


// M_eff y = M y - Σ_k m_k (m_k・y)/S_k   (out = M_eff y、comp[i] < 0 は導電領域の外)
static void mt_meff(const crs_t *M, const double *mv, const int *comp, const double *sk,
	double *dk, int ncomp, const double *y, double *out)
{
	const int n = (int)M->n;
	crs_spmv(M, y, out, NULL);
	for (int k = 0; k < ncomp; k++) dk[k] = 0;
	for (int i = 0; i < n; i++) {
		if (comp[i] >= 0) dk[comp[i]] += mv[i] * y[i];
	}
	for (int i = 0; i < n; i++) {
		if (comp[i] >= 0) out[i] -= mv[i] * dk[comp[i]] / sk[comp[i]];
	}
}


// (A - sc Σ_k m_k m_kᵀ/S_k) x = b の Jacobi-PCG (A は CRS、固定節点は恒等行で x = 0)
static int mt_pcg(const crs_t *A, const double *mv, const int *comp, const double *sk,
	int ncomp, double sc, const double *b, double *x, const unsigned char *fix,
	int maxiter, double tol)
{
	const int n = (int)A->n;
	double *r = (double *)malloc((size_t)n * sizeof(double));
	double *z = (double *)malloc((size_t)n * sizeof(double));
	double *p = (double *)malloc((size_t)n * sizeof(double));
	double *q = (double *)malloc((size_t)n * sizeof(double));
	double *dg = (double *)malloc((size_t)n * sizeof(double));
	double *dk = (double *)malloc(((size_t)ncomp + 1) * sizeof(double));
	crs_diag(A, dg);
	double bn = 0;
	for (int i = 0; i < n; i++) {
		x[i] = 0;
		r[i] = (fix[i] ? 0 : b[i]);
		bn += r[i] * r[i];
	}
	bn = sqrt(bn);
	int it = 0, ok = (bn == 0);
	double rz = 0;
	for (int i = 0; i < n; i++) {
		z[i] = (fix[i] ? 0 : r[i] / dg[i]);
		p[i] = z[i];
		rz += r[i] * z[i];
	}
	for (it = 1; !ok && (it <= maxiter); it++) {
		crs_spmv(A, p, q, fix);
		for (int k = 0; k < ncomp; k++) dk[k] = 0;
		for (int i = 0; i < n; i++) {
			if (comp[i] >= 0) dk[comp[i]] += mv[i] * p[i];
		}
		double pq = 0;
		for (int i = 0; i < n; i++) {
			if (fix[i]) { q[i] = 0; continue; }
			if (comp[i] >= 0) q[i] -= sc * mv[i] * dk[comp[i]] / sk[comp[i]];
			pq += p[i] * q[i];
		}
		const double al = rz / pq;
		double rn = 0, rz2 = 0;
		for (int i = 0; i < n; i++) {
			x[i] += al * p[i];
			r[i] -= al * q[i];
			rn += r[i] * r[i];
		}
		if (sqrt(rn) <= (tol * bn)) { ok = 1; break; }
		for (int i = 0; i < n; i++) {
			z[i] = (fix[i] ? 0 : r[i] / dg[i]);
			rz2 += r[i] * z[i];
		}
		const double be = rz2 / rz;
		rz = rz2;
		for (int i = 0; i < n; i++) p[i] = z[i] + (be * p[i]);
	}
	free(r); free(z); free(p); free(q); free(dg); free(dk);
	return (ok ? it : -1);
}


static int mag_transient(FILE *fp_log, crs_t *A, const double *b1, unsigned char *fix,
	double *nucell, double *az)
{
	const int n = (int)num_node();
	const int64_t ncell = (int64_t)Nx * Ny * Nz;
	const double dt = TransientTime / TransientSteps;
	const double scale = 1 / TlineLength;
	int ierr = 0;

	// 渦電流の σ (導体でない材料のセル)
	double *sig = (double *)malloc((size_t)ncell * sizeof(double));
	int nsigc = 0;
	for (int64_t c = 0; c < ncell; c++) {
		const double s = Material[CellMaterial[c]].sigma;
		sig[c] = (((CellConductor[c] < 0) && (s > 0)) ? s : 0);
		if (sig[c] > 0) nsigc++;
	}
	crs_t Ms;
	crs_alloc(&Ms);
	crs_distribute(&Ms);		// mpisolve = 2 : 剛性と同じ分け方 (A += Ms/Δt を値の配列で足す)
	assemble_mass_cell(&Ms, sig);

	// 孤立した導電領域 (節点を共有するセルどうしを 1 つにまとめる)
	int *par = (int *)malloc((size_t)n * sizeof(int));
	int *comp = (int *)malloc((size_t)n * sizeof(int));
	for (int i = 0; i < n; i++) par[i] = i;
	for (int i = 0; i < Nx; i++) {
	for (int j = 0; j < Ny; j++) {
	for (int k = 0; k < Nz; k++) {
		if (sig[((int64_t)i * Ny + j) * Nz + k] <= 0) continue;
		const int a0 = mt_find(par, (int)node_index(i, j, k));
		for (int l = 1; l < 8; l++) {
			const int b = mt_find(par, (int)node_index(i + ((l >> 2) & 1), j + ((l >> 1) & 1), k + (l & 1)));
			if (b != mt_find(par, a0)) par[b] = mt_find(par, a0);
		}
	}
	}
	}
	double *mv = (double *)malloc((size_t)n * sizeof(double));
	double *one = (double *)malloc((size_t)n * sizeof(double));
	for (int i = 0; i < n; i++) one[i] = 1;
	crs_spmv(&Ms, one, mv, NULL);			// m = M 1 = ∫σ N_i
	int ncomp = 0;
	int *cid = (int *)malloc((size_t)n * sizeof(int));
	for (int i = 0; i < n; i++) cid[i] = -1;
	for (int i = 0; i < n; i++) {
		comp[i] = -1;
		if (mv[i] <= 0) continue;
		const int r = mt_find(par, i);
		if (cid[r] < 0) cid[r] = ncomp++;
		comp[i] = cid[r];
	}
	double *sk = (double *)calloc((size_t)ncomp + 1, sizeof(double));
	double *dk = (double *)malloc(((size_t)ncomp + 1) * sizeof(double));
	for (int i = 0; i < n; i++) {
		if (comp[i] >= 0) sk[comp[i]] += mv[i];
	}
	free(par); free(cid); free(one);

	fprintf(fp_log, "  transient (eddy currents%s) : t_end = %.6e [s], %d steps (dt = %.6e [s]), "
		"%d conducting cells in %d isolated regions (zero net current each)\n",
		((JaB != NULL) ? " + hysteresis" : ""), TransientTime, TransientSteps, dt, nsigc, ncomp);
	FILE *fcsv = mpi_fopen_out(FN_transient, "w");
	if (fcsv == NULL) {
		fprintf(fp_log, "*** cannot open %s\n", FN_transient);
		ierr = 1;
	}
	else {
		fprintf(fcsv, "t [s],I [A],Lambda [Wb/m],H_mean [A/m],B_mean [T],iter\n");
	}
	fprintf(fp_log, "  %6s %13s %13s %13s %13s %13s %5s\n", "step", "t [s]", "I [A]",
		"Lambda [Wb/m]", "H_mean [A/m]", "B_mean [T]", "iter");

	double *res = (double *)malloc((size_t)n * sizeof(double));
	double *del = (double *)malloc((size_t)n * sizeof(double));
	double *r = (double *)malloc((size_t)n * sizeof(double));
	double *aold = (double *)malloc((size_t)n * sizeof(double));
	double *da = (double *)malloc((size_t)n * sizeof(double));
	double *w = (double *)malloc((size_t)n * sizeof(double));
	double *abest = (double *)malloc((size_t)n * sizeof(double));
	memset(az, 0, (size_t)n * sizeof(double));
	memset(aold, 0, (size_t)n * sizeof(double));

	// 残差の基準 : 波形の最大電流での右辺のノルム
	double fmax = 1;
	if (WaveN > 0) {
		fmax = 0;
		for (int m = 0; m < WaveN; m++) {
			if (fabs(WaveF[m]) > fmax) fmax = fabs(WaveF[m]);
		}
	}
	double bref = 0;
	for (int i = 0; i < n; i++) {
		if (!fix[i]) bref += b1[i] * b1[i];
	}
	bref = sqrt(bref) * ((fmax > 0) ? fmax : 1);
	if (bref <= 0) bref = 1;

	double win = 0, wmag = 0, wed = 0, lam = 0, cur = 0;
	for (int step = 1; !ierr && (step <= TransientSteps); step++) {
		const double t = step * dt;
		cur = Curr * mt_wave(t);
		const double fac = cur / Curr;
		field_free();
		int it, conv = 0, nback = 0;
		double rnorm = 0, rbest = 0, lrel = NlRelax;
		for (it = 1; it <= NlMaxiter; it++) {
			assemble_ja(az, res, nucell);
			for (int i = 0; i < n; i++) da[i] = az[i] - aold[i];
			if (nsigc > 0) mt_meff(&Ms, mv, comp, sk, dk, ncomp, da, w);
			rnorm = 0;
			for (int i = 0; i < n; i++) {
				r[i] = (fix[i] ? 0 : (((fac * b1[i]) - res[i]) - ((nsigc > 0) ? (w[i] / dt) : 0)));
				rnorm += r[i] * r[i];
			}
			rnorm = sqrt(rnorm) / bref;
			if (rnorm < NlTol) {
				conv = 1;
				break;
			}
			// 残差が増えたら直前の更新を半分に戻してやり直す (反転の近くでは微分磁気
			// 抵抗率の行列が真の接線より小さく、固定点が不安定になることがある : 実測)。
			// 単調に減っている間は何もしないので、σ = 0 の電流掃引との一致は崩れない
			// 発散の検出 : 残差がそれまでの最良の 10 倍を越えたら最良の反復に戻り、
			// この刻みの残りは行列の傾きを両方の枝の大きい方にする (単調な H(B) に対して
			// 縮小写像になる)。それでも越えるときは更新を半分にする。固定点の反復は
			// 単調に減るとは限らない (実測で途中 4 倍に増えてから収束する) ので、
			// 1 回増えただけでは何もしない (σ = 0 の電流掃引との一致を崩さないため)
			if ((it == 1) || (rnorm < rbest)) {
				rbest = rnorm;
				memcpy(abest, az, (size_t)n * sizeof(double));
			}
			else if ((rnorm > (10 * rbest)) && (nback < 30)) {
				memcpy(az, abest, (size_t)n * sizeof(double));
				if (!JaSafeSlope) JaSafeSlope = 1;
				else lrel /= 2;
				nback++;
				continue;
			}
			if (JaNu6 != NULL) assemble_nu6(A, JaNu6);
			else assemble_nu(A, nucell);
			int li;
			if (nsigc > 0) {
				const int64_t nnz = A->rowptr[n];
				for (int64_t e = 0; e < nnz; e++) A->val[e] += Ms.val[e] / dt;
				li = mt_pcg(A, mv, comp, sk, ncomp, 1 / dt, r, del, fix, Solver.maxiter, Solver.converg);
			}
			else {
				li = solver_cg(A, r, del, fix, Solver.maxiter, Solver.nout, Solver.converg, NULL, "ja");
			}
			if (li < 0) {
				fprintf(fp_log, "*** linear solver did not converge (step %d)\n", step);
				ierr = 1;
				break;
			}
			for (int i = 0; i < n; i++) az[i] += lrel * del[i];
		}
		JaSafeSlope = 0;
		if (ierr) break;
		if (!conv) {
			fprintf(fp_log, "*** transient : the iteration did not converge at step %d "
				"(residual = %.3e)\n", step, rnorm);
			ierr = 1;
			break;
		}

		// J-A の状態を確定させる (静磁場の電流掃引と同じ)
		ja_commit(ncell);

		// エネルギー : 入力 b・ΔA、磁気の仕事 ΔA・r_ja、渦電流損 ΔA・M_eff ΔA / Δt
		double din = 0, dmag = 0, ded = 0, lb = 0;
		for (int i = 0; i < n; i++) {
			if (fix[i]) continue;
			din += fac * b1[i] * da[i];
			dmag += da[i] * res[i];
			lb += b1[i] * az[i];
		}
		if (nsigc > 0) {
			mt_meff(&Ms, mv, comp, sk, dk, ncomp, da, w);
			for (int i = 0; i < n; i++) {
				if (!fix[i]) ded += da[i] * w[i] / dt;
			}
		}
		win += din * scale;
		wmag += dmag * scale;
		wed += ded * scale;
		lam = lb / Curr * scale;

		// 平均の H, B : ヒステリシス材料は J-A の状態 (符号つき)、それ以外の導電領域は
		// セル中心の |B|
		double bsum = 0, hsum = 0, vsum = 0;
		for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
		for (int k = 0; k < Nz; k++) {
			const int64_t c = ((int64_t)i * Ny + j) * Nz + k;
			const material_t *mt = &Material[CellMaterial[c]];
			const double dx = Xn[i + 1] - Xn[i], dy = Yn[j + 1] - Yn[j], dz = Zn[k + 1] - Zn[k];
			const double v = dx * dy * dz;
			if (mt->ja.on && (JaB != NULL)) {
				// 静磁場の電流掃引と同じ式・同じ順に足す (σ = 0 で解の列が一致する恒等式のため)
				const double wg = (Xn[i + 1] - Xn[i]) * (Yn[j + 1] - Yn[j]) * (Zn[k + 1] - Zn[k]) / 8;
				for (int g = 0; g < 8; g++) {
					bsum += JaB[(c * 8) + g] * wg;
					hsum += JaH[(c * 8) + g] * wg;
					vsum += wg;
				}
			}
			else if (sig[c] > 0) {
				double gr[3] = {0, 0, 0};
				for (int l = 0; l < 8; l++) {
					const double al = az[node_index(i + ((l >> 2) & 1), j + ((l >> 1) & 1), k + (l & 1))];
					gr[0] += al * (((l >> 2) & 1) ? 1 : -1) / (4 * dx);
					gr[1] += al * (((l >> 1) & 1) ? 1 : -1) / (4 * dy);
					gr[2] += al * ((l & 1) ? 1 : -1) / (4 * dz);
				}
				const double bm = sqrt((gr[0] * gr[0]) + (gr[1] * gr[1]) + (gr[2] * gr[2]));
				bsum += bm * v;
				hsum += bm * v / (MU0 * mt->mu6[0] * temp_factor(mt, 1));
				vsum += v;
			}
		}
		}
		}
		const double bmean = ((vsum > 0) ? (bsum / vsum) : 0);
		const double hmean = ((vsum > 0) ? (hsum / vsum) : 0);
		fprintf(fp_log, "  %6d %13.6e %13.6e %13.6e %13.6e %13.6e %5d\n", step, t, cur, lam,
			hmean, bmean, it);
		if (fcsv) fprintf(fcsv, "%.8e,%.8e,%.8e,%.8e,%.8e,%d\n", t, cur, lam, hmean, bmean, it);
		fflush(fp_log);
		memcpy(aold, az, (size_t)n * sizeof(double));
	}

	if (!ierr) {
		fprintf(fp_log, "  energy (per unit length) : input = %.12e, magnetic work = %.12e, "
			"eddy loss = %.12e [J/m]\n", win, wmag, wed);
		if (cur != 0) {
			Mmat[0] = lam * Curr / cur;		// 終了時刻の Λ / I
			HaveM = 1;
		}
		if (FieldOut) {
			field_add_node("Az", az);
			field_add_grad("B", az, 1);
		}
	}
	if (fcsv) fclose(fcsv);
	crs_free(&Ms);
	free(sig); free(comp); free(mv); free(sk); free(dk);
	free(res); free(del); free(r); free(aold); free(da); free(w); free(abest);

	return ierr;
}


static int solve_magnetostatic(FILE *fp_log)
{
	const int n = (int)num_node();
	const int np = NPort;
	int ierr = 0;

	if (np < 1) return 1;

	fprintf(fp_log, "\n=== magnetostatic analysis ===\n");

	// 電流密度の計算に導体の断面積が要る
	for (int p = 0; p <= np; p++) {
		if (CondArea[p] <= 0) {
			fprintf(fp_log, "*** conductor %d has no cell (magnetostatic analysis needs "
				"conductors with a finite cross section)\n", p);
			return 1;
		}
	}
	fprintf(fp_log, "  current = %.6e [A], cross sections [m^2] :", Curr);
	for (int p = 0; p <= np; p++) {
		fprintf(fp_log, " %d:%.4e", p, CondArea[p]);
	}
	fprintf(fp_log, "\n");
	fprintf(fp_log, "  %-10s %8s %13s\n", "port", "iter", "residual");
	fflush(fp_log);

	crs_t A;
	crs_alloc_dist(&A);		// mpisolve = 2 なら自分の行ブロックだけ
	{
		int64_t r0, r1;
		if (crs_dist_info(&A, &r0, &r1)) {
			fprintf(fp_log, "\n=== MPI : the matrix is distributed by rows (mpisolve = 2) : rank %d holds rows %lld..%lld "
				"of %d (%lld nonzeros) ===\n", mpi_rank(), (long long)r0, (long long)(r1 - 1), n, (long long)A.nnz);
		}
	}

	// 鉄損 (bertotti) を評価する材料があるか
	int havebe = 0;
	for (int m = 0; m < NMaterial; m++) {
		if ((Material[m].be_kh > 0) || (Material[m].be_ke > 0) || (Material[m].be_d > 0)) havebe = 1;
	}

	// 要素毎の磁気抵抗率 ν。非線形材料 (B-H) があると反復で更新する
	const int tri2d = (MeshMode && (MeshDim == 2));
	const int64_t ncell = elem_count();
	int nonlinear = 0;
	for (int m = 0; m < NMaterial; m++) {
		if (Material[m].nbh[0] > 0) nonlinear = 1;
	}
	double *nucell = (double *)malloc((size_t)ncell * sizeof(double));
	for (int64_t c = 0; c < ncell; c++) {
		elem_t el;
		elem_get(c, &el);
		const material_t *mt = &Material[el.mat];
		// 線形材料の μr にも温度依存を掛ける (B-H は bh_b が入力時に補正済み)
		const double fm = temp_factor(mt, 1);
		nucell[c] = ((mt->nbh[0] > 0) ? bh_nu(mt, 0, 0) : (1 / (MU0 * mt->mu6[0] * fm)));
	}
	// 2 次元格子では異方性テンソルを面内 2x2 で使えるので nucell は渡さない
	// 線形なら異方性テンソルを使う (nucell は等方性なので非線形反復専用)
	if (tri2d)          assemble_nu_tri(&A, nonlinear ? nucell : NULL);
	else if (nonlinear) assemble_nu(&A, nucell);
	else                assemble(&A, 4);

	// 定数分の不定性を除くために 1 節点だけ固定する
	unsigned char *fix = (unsigned char *)malloc(n * sizeof(unsigned char));
	memset(fix, 0, n * sizeof(unsigned char));
	fix[0] = 1;

	// ポート毎の電流密度ベクトル (要素の形状関数で節点に配分)
	double **b = (double **)malloc((size_t)np * sizeof(double *));
	for (int k = 1; k <= np; k++) {
		double *bk = (double *)malloc(n * sizeof(double));
		memset(bk, 0, n * sizeof(double));
		for (int64_t e = 0; e < ncell; e++) {
			elem_t el;
			elem_get(e, &el);
			double jz = 0;
			if      (el.cond == k) jz = +Curr / CondArea[k];
			else if (el.cond == 0) jz = -Curr / CondArea[0];
			else continue;
			for (int l = 0; l < el.nen; l++) bk[el.nd[l]] += jz * el.w[l];
		}
		b[k - 1] = bk;
	}

	double *rhs = (double *)malloc(n * sizeof(double));
	double *az  = (double *)malloc(n * sizeof(double));
	const double scale = 1 / TlineLength;

	// 過渡解析 (渦電流 + ヒステリシス) : mag_transient が時間発展を解く
	if (TransientOn) {
		ierr = mag_transient(fp_log, &A, b[0], fix, nucell, az);
		free(b[0]); free(b); free(rhs); free(az); free(fix); free(nucell);
		crs_free(&A);
		return ierr;
	}

	// ヒステリシス (Jiles-Atherton) : 電流掃引に沿って履歴を追う
	//
	// ν = H/B はヒステリシスでは第 2 象限で負になるため行列には使えない。
	// 行列は常に正の微分磁気抵抗率 dH/dB で組み (SPD なので Jacobi-PCG が使える)、
	// 残差側だけで正確な H(B) を扱う固定点反復にする。
	// 各掃引点で収束したら J-A の状態を確定させ、次の点はそこから続ける。
	{
		int hysteresis = 0;
		for (int m = 0; m < NMaterial; m++) {
			if (Material[m].ja.on) hysteresis = 1;
		}
		if (hysteresis) {
			fprintf(fp_log, "  Jiles-Atherton hysteresis : %d sweep points, %d sub-steps\n",
				NSweep, JaSub);
			fprintf(fp_log, "  %5s %13s %13s %13s %13s %6s\n",
				"step", "I [A]", "H [A/m]", "B [T]", "L [H/m]", "iter");
			fflush(fp_log);

			double *res  = (double *)malloc(n * sizeof(double));
			double *del  = (double *)malloc(n * sizeof(double));
			double *rhs2 = (double *)malloc(n * sizeof(double));
			memset(az, 0, n * sizeof(double));

			// 残差の正規化に使う基準 (掃引中の最大電流の右辺ノルム)
			double amax = 0;
			for (int q = 0; q < NSweep; q++) {
				if (fabs(Sweep[q]) > amax) amax = fabs(Sweep[q]);
			}
			double bref = 0;
			for (int i = 0; i < n; i++) {
				if (!fix[i]) bref += b[0][i] * b[0][i];
			}
			bref = sqrt(bref) * ((amax > 0) ? (amax / fabs(Curr)) : 1);
			if (bref <= 0) bref = 1;

			// 2 ポートの掃引 (回転磁界) の右辺の正規化は 2 つのポートの最大電流の和で取る
			if (NSweep2 > 0) {
				double a2 = 0, b2 = 0;
				for (int q = 0; q < NSweep2; q++) {
					if (fabs(Sweep2[q]) > a2) a2 = fabs(Sweep2[q]);
				}
				for (int i = 0; i < n; i++) {
					if (!fix[i]) b2 += b[1][i] * b[1][i];
				}
				bref += sqrt(b2) * ((a2 > 0) ? (a2 / fabs(Curr)) : 1);
			}
			// 仕事の収支と観測点
			double win = 0, wmag = 0, iprev[2] = {0, 0}, lprev[2] = {0, 0};
			int64_t probe = -1;
			if (JaProbeOn) {
				double dbest = 1e300;
				const double gq = 1 / sqrt(3.0);
				for (int i = 0; i < Nx; i++) {
				for (int j = 0; j < Ny; j++) {
				for (int k = 0; k < Nz; k++) {
					const int64_t c = ((int64_t)i * Ny + j) * Nz + k;
					if (!Material[CellMaterial[c]].ja.on) continue;
					for (int g = 0; g < 8; g++) {
						const double px = Xn[i] + ((Xn[i + 1] - Xn[i]) * (1 + ((((g >> 2) & 1) ? 1 : -1) * gq)) / 2);
						const double py = Yn[j] + ((Yn[j + 1] - Yn[j]) * (1 + ((((g >> 1) & 1) ? 1 : -1) * gq)) / 2);
						const double pz = Zn[k] + ((Zn[k + 1] - Zn[k]) * (1 + (((g & 1) ? 1 : -1) * gq)) / 2);
						const double dd = ((px - JaProbe[0]) * (px - JaProbe[0])) + ((py - JaProbe[1]) * (py - JaProbe[1]))
						                + ((pz - JaProbe[2]) * (pz - JaProbe[2]));
						if (dd < dbest) { dbest = dd; probe = (c * 8) + g; }
					}
				}
				}
				}
				if (probe < 0) {
					fprintf(fp_log, "*** warning : japrobe : no hysteretic cell, the probe is ignored\n");
				}
				else {
					fprintf(fp_log, "  probe : Gauss point %lld (distance %.3e [m]); columns : step, %s\n",
						(long long)probe, sqrt(dbest), (JaVector ? "H_p H_q [A/m], B_p B_q [T] (in-plane, p, q follow tline)"
						                                           : "H [A/m], B [T] (signed along the first direction)"));
				}
			}
			for (int is = 0; is < NSweep; is++) {
				const double cur = Sweep[is];
				const double fac = cur / Curr;
				// 前のステップの場を捨てる (最後の 1 点だけが後で VTK に残る)
				field_free();
				const double fac2 = ((NSweep2 > 0) ? (Sweep2[is] / Curr) : 0);
				for (int i = 0; i < n; i++) {
					rhs[i] = (fix[i] ? 0 : ((fac * b[0][i]) + ((NSweep2 > 0) ? (fac2 * b[1][i]) : 0)));
				}

				int it;
				double rnorm = 0;
				int converged = 0;
				for (it = 1; it <= NlMaxiter; it++) {
					assemble_ja(az, res, nucell);
					rnorm = 0;
					for (int i = 0; i < n; i++) {
						rhs2[i] = (fix[i] ? 0 : (rhs[i] - res[i]));
						rnorm += rhs2[i] * rhs2[i];
					}
					rnorm = sqrt(rnorm) / bref;
					if (rnorm < NlTol) {
						converged = 1;
						break;
					}
					if (JaNu6 != NULL) assemble_nu6(&A, JaNu6);
					else assemble_nu(&A, nucell);
					if (solver_cg(&A, rhs2, del, fix, Solver.maxiter, Solver.nout,
							Solver.converg, NULL, "ja") < 0) {
						fprintf(fp_log, "*** linear solver did not converge (sweep %d)\n", is + 1);
						ierr = 1;
						break;
					}
					for (int i = 0; i < n; i++) {
						az[i] += NlRelax * del[i];
					}
				}
				if (!converged) {
					fprintf(fp_log, "*** warning : hysteresis iteration did not converge "
						"(step %d, residual = %.3e)\n", is + 1, rnorm);
				}

				// 磁気の仕事 ∫ H・ΔB (台形則、確定前の状態と比べる)。ベクトル J-A は面内 2 成分
				for (int i = 0; i < Nx; i++) {
				for (int j = 0; j < Ny; j++) {
				for (int k = 0; k < Nz; k++) {
					const int64_t c = ((int64_t)i * Ny + j) * Nz + k;
					if (!Material[CellMaterial[c]].ja.on) continue;
					const double w = (Xn[i + 1] - Xn[i]) * (Yn[j + 1] - Yn[j]) * (Zn[k + 1] - Zn[k]) / 8;
					for (int g = 0; g < 8; g++) {
						const int64_t id = (c * 8) + g;
						if (JaVector) {
							for (int q = 0; q < 2; q++) {
								wmag += w * scale * (JaVH[(id * 2) + q] + JaVHn[(id * 2) + q]) / 2
								      * (JaVBn[(id * 2) + q] - JaVB[(id * 2) + q]);
							}
						}
						else {
							wmag += w * scale * (JaH[id] + JaHn[id]) / 2 * (JaBn[id] - JaB[id]);
						}
					}
				}
				}
				}
				// 状態を確定させる (方向は最初に磁化したときだけ決める)
				ja_commit(ncell);
				// 端子の仕事 Σ_k Ī_k ΔΛ_k (台形則)。Λ_k = b_k・A / I_ref (単位長あたり)
				for (int pp = 0; pp < np; pp++) {
					const double ik = ((pp == 0) ? cur : Sweep2[is]);
					double lam = 0;
					for (int i = 0; i < n; i++) lam += b[pp][i] * az[i];
					lam *= scale / Curr;
					win += (ik + iprev[pp]) / 2 * (lam - lprev[pp]);
					iprev[pp] = ik;
					lprev[pp] = lam;
				}
				if (JaProbeOn && (probe >= 0)) {
					if (JaVector) {
						fprintf(fp_log, "  probe %5d %13.6e %13.6e %13.6e %13.6e\n", is + 1,
							JaVH[(probe * 2) + 0], JaVH[(probe * 2) + 1], JaVB[(probe * 2) + 0], JaVB[(probe * 2) + 1]);
					}
					else {
						fprintf(fp_log, "  probe %5d %13.6e %13.6e\n", is + 1, JaH[probe], JaB[probe]);
					}
				}

				// ヒステリシス材料の体積平均 B, H
				double bsum = 0, hsum = 0, vsum = 0;
				for (int i = 0; i < Nx; i++) {
				for (int j = 0; j < Ny; j++) {
				for (int k = 0; k < Nz; k++) {
					const int64_t c = ((int64_t)i * Ny + j) * Nz + k;
					if (!Material[CellMaterial[c]].ja.on) continue;
					const double w = (Xn[i + 1] - Xn[i]) * (Yn[j + 1] - Yn[j])
					               * (Zn[k + 1] - Zn[k]) / 8;
					for (int g = 0; g < 8; g++) {
						const int64_t id = (c * 8) + g;
						bsum += JaB[id] * w;
						hsum += JaH[id] * w;
						vsum += w;
					}
				}
				}
				}
				const double bmean = ((vsum > 0) ? (bsum / vsum) : 0);
				const double hmean = ((vsum > 0) ? (hsum / vsum) : 0);

				// 磁束鎖交インダクタンス (電流 0 の点では定義できない)
				double lval = 0;
				if (cur != 0) {
					double w = 0;
					for (int i = 0; i < n; i++) {
						w += fac * b[0][i] * az[i];
					}
					lval = w / (cur * cur) * scale;
					Mmat[0] = lval;
					HaveM = 1;
				}
				fprintf(fp_log, "  %5d %13.6e %13.6e %13.6e %13.6e %6d\n",
					is + 1, cur, hmean, bmean, lval, it);
				fflush(fp_log);

				// 履歴を系列に追記する。**順序そのものが物理**なので
				// (B は履歴依存)、読み手が並べ替えないよう属性にも書いてある
				ierr |= h5_add_hysteresis(is + 1, cur, hmean, bmean, lval, it);
				if (FieldOut) {
					field_add_node("Az", az);
					field_add_grad("B", az, 1);
					ierr |= h5_add_field(cur);
				}
				if (ierr) break;
			}

			// 仕事の収支 : 収束した各点で Σ I ΔΛ = ∫ H・ΔB (Galerkin の残差 0) なので、台形則でも
			// 反復の打ち切りを除いて一致する。閉じた履歴なら差し引きがヒステリシス損 (単位長あたり)
			// 線形のセル (空気・導体) の分は最後の状態の磁気エネルギー (台形則が厳密に telescoping)
			const double wlin = scale * ja_linear_energy(az);
			fprintf(fp_log, "  hysteresis work : terminals = %.12e, magnetic (int H.dB) = %.12e "
				"(hysteretic %.12e + linear %.12e) [J/m]\n", win, wmag + wlin, wmag, wlin);
			free(res);
			free(del);
			free(rhs2);
			for (int k = 0; k < np; k++) free(b[k]);
			free(b);
			free(rhs);
			free(az);
			free(fix);
			free(nucell);
			crs_free(&A);

			return ierr;
		}
	}

	// 非線形 (B-H) 反復 (Newton-Raphson)
	//
	//   R(A) = K(ν(|∇A|)) A - b = 0
	//   J δ = -R,  A += w δ
	//
	// ν を Gauss 点毎に評価したヤコビアンを使うので 2 次収束する。
	// H(B) が単調なら J は正定値なので Jacobi-PCG で解ける。
	// ν の逐次代入 (successive substitution) は飽和領域で振動するため使わない。
	// 電流駆動なので右辺 b は ν に依らない。重ね合わせが成り立たないので単一ポート。
	if (nonlinear) {
		fprintf(fp_log, "  nonlinear (B-H) Newton iteration : maxiter=%d, tol=%.1e, damping=%.2f\n",
			NlMaxiter, NlTol, NlRelax);

		// 初期値 : 初期透磁率での線形解
		memcpy(rhs, b[0], n * sizeof(double));
		rhs[0] = 0;
		if (solver_cg(&A, rhs, az, fix, Solver.maxiter, Solver.nout,
				Solver.converg, NULL, "bh0") < 0) {
			fprintf(fp_log, "*** linear solver did not converge (B-H initial guess)\n");
			ierr = 1;
		}

		crs_t J;
		crs_alloc(&J);
		crs_distribute(&J);		// mpisolve = 2 : ヤコビアンも自分の行ブロックだけ
		double *res = (double *)malloc(n * sizeof(double));
		double *del = (double *)malloc(n * sizeof(double));

		double bnorm = 0;
		for (int i = 0; i < n; i++) {
			if (!fix[i]) bnorm += b[0][i] * b[0][i];
		}
		bnorm = sqrt(bnorm);

		int nlit;
		double rnorm = 0;
		int converged = 0;
		for (nlit = 1; nlit <= NlMaxiter; nlit++) {
			if (tri2d) assemble_newton_tri(&J, az, res);
			else       assemble_newton(&J, az, res);

			// R = K(ν)A - b (固定節点は 0)
			rnorm = 0;
			for (int i = 0; i < n; i++) {
				res[i] = (fix[i] ? 0 : (b[0][i] - res[i]));	// -R (右辺)
				rnorm += res[i] * res[i];
			}
			rnorm = sqrt(rnorm) / ((bnorm > 0) ? bnorm : 1);
			fprintf(fp_log, "  %-10s %8d %13.5e\n", "B-H", nlit, rnorm);
			fflush(fp_log);
			if (rnorm < NlTol) {
				converged = 1;
				break;
			}

			if (solver_cg(&J, res, del, fix, Solver.maxiter, Solver.nout,
					Solver.converg, NULL, "bh") < 0) {
				fprintf(fp_log, "*** Newton step did not converge (B-H iteration %d)\n", nlit);
				ierr = 1;
				break;
			}
			for (int i = 0; i < n; i++) {
				az[i] += NlRelax * del[i];
			}
		}
		if (!converged) {
			fprintf(fp_log, "*** warning : B-H Newton iteration did not converge "
				"(%d iterations, residual = %.3e)\n", nlit - 1, rnorm);
		}

		// 収束した A から L を求める (ポートループは通さない)
		double w = 0;
		for (int i = 0; i < n; i++) {
			w += b[0][i] * az[i];
		}
		Mmat[0] = w / (Curr * Curr) * scale;
		HaveM = 1;
		if (FieldOut) {
			field_add_node("Az_M", az);
			field_add_grad("B_M", az, 1);
		}

		crs_free(&J);
		free(res);
		free(del);

		free(b[0]);
		free(b);
		free(rhs);
		free(az);
		free(fix);
		free(nucell);
		crs_free(&A);

		return ierr;
	}

	for (int k = 1; k <= np; k++) {
		// 固定節点の右辺は 0 にする (行を恒等行として扱うため)
		memcpy(rhs, b[k - 1], n * sizeof(double));
		rhs[0] = 0;

		char label[BUFSIZ];
		sprintf(label, "port%d", k);
		const int iter = solver_cg(&A, rhs, az, fix,
			Solver.maxiter, Solver.nout, Solver.converg, fp_log, label);
		if (iter < 0) {
			fprintf(fp_log, "*** solver did not converge (port %d)\n", k);
			ierr = 1;
		}

		// L_kj = ∫ J^(j)・A^(k) dV / (I_k I_j) / 線路長
		for (int j = 1; j <= np; j++) {
			double w = 0;
			for (int i = 0; i < n; i++) {
				w += b[j - 1][i] * az[i];
			}
			Mmat[((k - 1) * np) + (j - 1)] = w / (Curr * Curr) * scale;
		}
		if (FieldOut) {
			char nm[BUFSIZ];
			sprintf(nm, "Az_M_port%d", k);
			field_add_node(nm, az);
			sprintf(nm, "B_M_port%d", k);
			field_add_grad(nm, az, 1);
		}

		// 鉄損 (Bertotti の損失分離)。要素毎の |B| を波高値として使う
		if (havebe) {
			Pfemat[((k - 1) * np) + (k - 1)] = iron_loss(az, fp_log, k);
			HavePfe = 1;
		}
	}

	symmetrize(Mmat, np, fp_log, "magnetostatic");
	HaveM = 1;

	for (int k = 0; k < np; k++) {
		free(b[k]);
	}
	free(b);
	free(rhs);
	free(az);
	free(fix);
	free(nucell);
	crs_free(&A);

	return ierr;
}


// 時間調和渦電流解析 (断面 2 次元、複素 Az 定式化)
//
//   ∇・(ν ∇Az) - jωσ Az + σ V' = 0     (導体内)
//   ∇・(ν ∇Az) = 0                      (それ以外)
//
// V' = -dV/dz は導体毎の軸方向電界 (一定)。基準導体を V'=0 に固定し、
// ポート j だけを V'=1 [V/m] で励振して解き、各ポートの電流
//   I_k = ∫_k σ (V'_k - jω Az) dS
// から単位長あたりのアドミタンス行列 Y[k][j] = I_k を作る。
// Z_loop = inv(Y) が基準導体を帰路とするループインピーダンスで、
//   R(f) = Re Z_loop,  L(f) = Im Z_loop / ω
// 導体内部の電流分布を解くので表皮効果・近接効果が入る。
//
// 境界条件は自然境界条件 (磁気壁 : 接線方向 B = 0)。これは領域を貫く
// 正味電流を 0 にする条件、すなわち「帰路電流はモデル内の導体を流れる」
// という伝送線路の前提そのものなので、基準導体が自動的に帰路になる。
// K + jωM は K が特異でも σ>0 の質量項があるため正則になる。
#define AAM (5)		// 非線形渦電流の Anderson 加速の深さ

static int solve_eddy(FILE *fp_log)
{
	const int n = (int)num_node();
	const int np = NPort;
	const int nc = np + 1;			// 基準導体を含む導体数
	int ierr = 0;

	if (np < 1) return 1;

	fprintf(fp_log, "\n=== eddy current (time-harmonic) analysis ===\n");

	if (Freq <= 0) {
		fprintf(fp_log, "*** analysis F requires the frequency key\n");
		return 1;
	}
	for (int p = 0; p < nc; p++) {
		if (CondSigma[p] <= 0) {
			fprintf(fp_log, "*** analysis F requires conductorsigma for conductor %d\n", p);
			return 1;
		}
	}

	const double omega = 2 * PI * Freq;
	const double delta = sqrt(2 / (omega * MU0 * CondSigma[1]));
	fprintf(fp_log, "  frequency = %.6e [Hz], skin depth (conductor 1) = %.4e [m]\n",
		Freq, delta);

	// 表皮深さを格子で刻めているか確認する (刻めていないと R(f) を過小評価する)
	double hmax = 0;
	if (MeshMode && (MeshDim == 2)) {
		// 三角形では辺の最大長を要素寸法とする
		int pa, qa;
		tri_axes(&pa, &qa);
		const double *cp = ((pa == 0) ? Xp : (pa == 1) ? Yp : Zp);
		const double *cq = ((qa == 0) ? Xp : (qa == 1) ? Yp : Zp);
		for (int e = 0; e < NTri; e++) {
			if (TriCond[e] < 0) continue;
			for (int l = 0; l < 3; l++) {
				const int32_t a = Tri[(e * 3) + l], bb = Tri[(e * 3) + ((l + 1) % 3)];
				const double dp = cp[bb] - cp[a], dq = cq[bb] - cq[a];
				const double len = sqrt((dp * dp) + (dq * dq));
				if (len > hmax) hmax = len;
			}
		}
	}
	else {
	for (int i = 0; i < Nx; i++) {
	for (int j = 0; j < Ny; j++) {
	for (int k = 0; k < Nz; k++) {
		if (CellConductor[((int64_t)i * Ny + j) * Nz + k] < 0) continue;
		if ((Tline != 'X') && ((Xn[i + 1] - Xn[i]) > hmax)) hmax = Xn[i + 1] - Xn[i];
		if ((Tline != 'Y') && ((Yn[j + 1] - Yn[j]) > hmax)) hmax = Yn[j + 1] - Yn[j];
		if ((Tline != 'Z') && ((Zn[k + 1] - Zn[k]) > hmax)) hmax = Zn[k + 1] - Zn[k];
	}
	}
	}
	}
	if (hmax > (delta / 2)) {
		fprintf(fp_log, "*** warning : conductor cell size %.4e [m] is coarser than "
			"half the skin depth; R(f) will be underestimated\n", hmax);
		printf("*** warning : conductor mesh does not resolve the skin depth "
			"(cell %.3e > delta/2 = %.3e)\n", hmax, delta / 2);
	}

	fprintf(fp_log, "  %-10s %8s %13s\n", "excite", "iter", "residual");
	fflush(fp_log);

	const int tri2d = (MeshMode && (MeshDim == 2));
	crs_t K, M;
	crs_alloc(&K);
	crs_alloc(&M);
	crs_distribute_pair(fp_log, &K, &M);		// mpisolve = 2 : 自分の行ブロックだけ (同じパターンを同じ分け方で)
	if (tri2d) assemble_mass_tri(&M);		// σ
	else       assemble_mass(&M);

	// 非線形 (B-H) の磁性体があれば**実効磁気抵抗率法**: 要素毎の ν = H(|B|)/|B| を
	// 波高値 |B| = sqrt(|Bx|² + |By|²) (複素振幅の大きさ) で評価して反復で更新する。
	// 単一ポート・与えた電流振幅 (current) での R(f) / L(f)。各回は線形の複素系なので
	// V' = 1 で解いて電流が current になるよう複素数倍し、その場で ν を作り直す
	// (時間調和の場では |B| が解析的でないので Newton は使わず、緩和付きの固定点反復)
	const int64_t nelc = elem_count();
	int nlf = 0;
	for (int m = 0; m < NMaterial; m++) {
		if (Material[m].nbh[0] > 0) nlf = 1;
	}
	double *nucell = NULL, *bvr = NULL, *bvi = NULL;
	double *aay = NULL, *aaf = NULL, *aayp = NULL, *aafp = NULL, *aady = NULL, *aadf = NULL;
	int aan = 0;			// 履歴の数
	if (nlf) {
		nucell = (double *)malloc((size_t)nelc * sizeof(double));
		bvr = (double *)malloc((size_t)nelc * 3 * sizeof(double));
		bvi = (double *)malloc((size_t)nelc * 3 * sizeof(double));
		for (int64_t c = 0; c < nelc; c++) {
			elem_t el;
			elem_get(c, &el);
			const material_t *mt = &Material[el.mat];
			nucell[c] = ((mt->nbh[0] > 0) ? bh_nu(mt, 0, 0)
				: (1 / (MU0 * mt->mu6[0] * temp_factor(mt, 1))));
		}
		// Anderson 加速の履歴 (ln ν の差分を最大 AAM 個)
		aay  = (double *)malloc((size_t)nelc * sizeof(double));
		aaf  = (double *)malloc((size_t)nelc * sizeof(double));
		aayp = (double *)malloc((size_t)nelc * sizeof(double));
		aafp = (double *)malloc((size_t)nelc * sizeof(double));
		aady = (double *)malloc((size_t)nelc * AAM * sizeof(double));
		aadf = (double *)malloc((size_t)nelc * AAM * sizeof(double));
		fprintf(fp_log, "  nonlinear (B-H) : effective reluctivity at current amplitude "
			"%.6e [A], maxiter = %d, tol = %.1e, relax = %.2f\n", Curr, NlMaxiter, NlTol, NlRelax);
		fprintf(fp_log, "  %-6s %16s %16s %12s\n", "nl-it", "R(f) [ohm/m]", "L(f) [H/m]", "max dnu/nu");
	}

	// 自然境界条件 (磁気壁) なので固定節点は無い
	unsigned char *fix = (unsigned char *)malloc(n * sizeof(unsigned char));
	memset(fix, 0, n * sizeof(unsigned char));

	// 導体毎の ∫σ dV
	const int64_t nel = elem_count();
	double *svol = (double *)calloc((size_t)nc, sizeof(double));
	for (int64_t e = 0; e < nel; e++) {
		elem_t el;
		elem_get(e, &el);
		if (el.cond < 0) continue;
		svol[el.cond] += CondSigma[el.cond] * el.vol;
	}

	double *br = (double *)malloc(n * sizeof(double));
	double *bi = (double *)malloc(n * sizeof(double));
	double *ar = (double *)malloc(n * sizeof(double));
	double *ai = (double *)malloc(n * sizeof(double));
	double *yr = (double *)calloc((size_t)np * np, sizeof(double));
	double *yi = (double *)calloc((size_t)np * np, sizeof(double));
	const double tlen = TlineLength;

	// 線形なら 1 回で終わり。非線形は ν が収束するまで回し、最後にもう 1 回
	// (場の出力つき) 解いて Z を作る
	int final = !nlf, nlconv = 0;
	for (int nlit = 1; !ierr; nlit++) {
	if (tri2d) assemble_nu_tri(&K, nlf ? nucell : NULL);	// ν = 1/(μ0 μr)
	else if (nlf) assemble_nu(&K, nucell);
	else       assemble(&K, 4);			// ν を ∇Az の基底で使う

	// 基準導体は V'=0。ポート 1..np を順に V'=1 [V/m] で励振する
	for (int jc = 1; jc <= np; jc++) {
		// b_i = ∫ σ N_i V' dV  (V' = 1 [V/m] を導体 jc にのみ与える)
		memset(br, 0, n * sizeof(double));
		memset(bi, 0, n * sizeof(double));
		for (int64_t e = 0; e < nel; e++) {
			elem_t el;
			elem_get(e, &el);
			if (el.cond != jc) continue;
			for (int l = 0; l < el.nen; l++) br[el.nd[l]] += CondSigma[jc] * el.w[l];
		}
		for (int i = 0; i < n; i++) {
			if (fix[i]) br[i] = 0;
		}

		char label[BUFSIZ];
		sprintf(label, "cond%d", jc);
		const int iter = solver_cocg(&K, &M, omega, br, bi, ar, ai, fix,
			Solver.maxiter, Solver.nout, Solver.converg, (final ? fp_log : NULL), label);
		if (iter < 0) {
			fprintf(fp_log, "*** solver did not converge (conductor %d)\n", jc);
			ierr = 1;
		}

		// I_k = (1/t) [ V'_k ∫_k σ dV - jω ∫_k σ Az dV ]
		double *qr = (double *)calloc((size_t)nc, sizeof(double));
		double *qi = (double *)calloc((size_t)nc, sizeof(double));
		for (int64_t e = 0; e < nel; e++) {
			elem_t el;
			elem_get(e, &el);
			if (el.cond < 0) continue;
			for (int l = 0; l < el.nen; l++) {
				const double w = CondSigma[el.cond] * el.w[l];
				qr[el.cond] += w * ar[el.nd[l]];
				qi[el.cond] += w * ai[el.nd[l]];
			}
		}
		for (int kc = 1; kc <= np; kc++) {
			// I = (V' ∫σdV - jω ∫σAz dV) / t、 -jω(qr + j qi) = ω qi - jω qr
			yr[((kc - 1) * np) + (jc - 1)] = (((kc == jc) ? svol[kc] : 0) + (omega * qi[kc])) / tlen;
			yi[((kc - 1) * np) + (jc - 1)] = (-omega * qr[kc]) / tlen;
		}
		free(qr);
		free(qi);

		// 非線形 : 電流が current になるよう Az を s = I / I_V 倍する (V'_k も s 倍)。
		// 場の出力も与えた電流での物理的な場にする
		double vr = 1, vi = 0;
		if (nlf) {
			const double ivr = yr[0], ivi = yi[0];		// 端子電流 [A] (V' = 1 [V/m])
			const double den = (ivr * ivr) + (ivi * ivi);
			vr = Curr * ivr / den;					// s = Curr / I_V = Curr conj(I_V) / |I_V|²
			vi = -Curr * ivi / den;
			for (int i = 0; i < n; i++) {
				const double xr = ar[i], xi = ai[i];
				ar[i] = (vr * xr) - (vi * xi);
				ai[i] = (vr * xi) + (vi * xr);
			}
		}

		// 場の出力 : 複素 Az と、そこから作る B = ∇×(Az ẑ_t) と J_z。
		// J_z = σ (V'_k - jω Az)、V' は励振した導体のセルだけ 1 [V/m] (非線形では s)
		if (FieldOut && final) {
			char nm[BUFSIZ];
			sprintf(nm, "Az_F_re_port%d", jc);  field_add_node(nm, ar);
			sprintf(nm, "Az_F_im_port%d", jc);  field_add_node(nm, ai);
			sprintf(nm, "B_F_re_port%d", jc);   field_add_grad(nm, ar, 1);
			sprintf(nm, "B_F_im_port%d", jc);   field_add_grad(nm, ai, 1);

			const int64_t ncell = num_cell();
			const int taxis = ((Tline == 'X') ? 0 : (Tline == 'Y') ? 1 : 2);
			double *jvr = (double *)calloc((size_t)ncell * 3, sizeof(double));
			double *jvi = (double *)calloc((size_t)ncell * 3, sizeof(double));
			for (int64_t c = 0; c < ncell; c++) {
				elem_t el;
				elem_get(c, &el);
				if (el.cond < 0) continue;
				// 要素中心の Az。形状関数の重心での値で重みづける
				// (2 次三角形は単純平均ではない)
				double zr8 = 0, zi8 = 0;
				for (int l = 0; l < el.nen; l++) {
					zr8 += el.wc[l] * ar[el.nd[l]];
					zi8 += el.wc[l] * ai[el.nd[l]];
				}
				const double sg = CondSigma[el.cond];
				const double vpr = ((el.cond == jc) ? vr : 0.0);
				const double vpi = ((el.cond == jc) ? vi : 0.0);
				// -jω(zr8 + j zi8) = ω zi8 - j ω zr8
				jvr[(c * 3) + taxis] = sg * (vpr + (omega * zi8));
				jvi[(c * 3) + taxis] = (nlf ? (sg * (vpi - (omega * zr8))) : (sg * (-omega * zr8)));
			}
			sprintf(nm, "J_F_re_port%d", jc);  field_add_cellvec(nm, jvr);
			sprintf(nm, "J_F_im_port%d", jc);  field_add_cellvec(nm, jvi);
			free(jvr);
			free(jvi);
		}
	}
	if (final || ierr) break;

	// ν の更新 : |B|² = |Re B|² + |Im B|²。**H を独立変数にして曲線を辿る**:
	// 要素の |H| = ν |B| (いまの線形解の H) に対して ν_new = |H| / B(|H|)。
	// B を独立変数にした ν_new = H(|B|)/|B| の逐次代入は電流で駆動する問題では
	// 振動して収束しない (実測: 1 次元の板間で 2 値を往復した)。電流駆動では
	// H の方がほぼ決まっているので、こちらは 1 次元では 1 回で厳密に収束する。
	// 緩和 ν += w (ν_new - ν)
	field_cell_grad(ar, 1, bvr);
	field_cell_grad(ai, 1, bvi);
	double dmax = 0;
	for (int64_t c = 0; c < nelc; c++) {
		elem_t el;
		elem_get(c, &el);
		const material_t *mt = &Material[el.mat];
		aay[c] = log(nucell[c]);
		if (mt->nbh[0] <= 0) {
			aaf[c] = 0;
			continue;
		}
		double b2 = 0;
		for (int d = 0; d < 3; d++) {
			b2 += (bvr[(c * 3) + d] * bvr[(c * 3) + d]) + (bvi[(c * 3) + d] * bvi[(c * 3) + d]);
		}
		const double hm = nucell[c] * sqrt(b2);
		const double nun = ((hm > 0) ? (hm / bh_b_of_h(mt, 0, hm)) : bh_nu(mt, 0, 0));
		const double dn = fabs(nun - nucell[c]) / nun;
		if (dn > dmax) dmax = dn;
		aaf[c] = log(nun) - aay[c];		// 固定点写像の残差 (ln ν で)
	}
	// Anderson 加速 (深さ AAM)。ln ν を変数にする (ν は桁で動き、正でなければならない)。
	// 線形収束の固定点反復 (2 次元の格子で 77 回、実測) を数回に縮める
	anderson_step(nelc, aay, aaf, aayp, aafp, aady, aadf, AAM, &aan, NlRelax);
	for (int64_t c = 0; c < nelc; c++) nucell[c] = exp(aay[c]);
	{
		// Z = 1 / Y (単一ポート) を反復の経過としてログに出す
		const double y2 = (yr[0] * yr[0]) + (yi[0] * yi[0]);
		fprintf(fp_log, "  %6d %16.9e %16.9e %12.3e\n", nlit, yr[0] / y2, (-yi[0] / y2) / omega, dmax);
		fflush(fp_log);
	}
	if (dmax < NlTol) {
		nlconv = 1;
		final = 1;
	}
	else if (nlit >= NlMaxiter) {
		fprintf(fp_log, "*** nonlinear (B-H) eddy current : the effective reluctivity did not "
			"converge in %d iterations (max dnu/nu = %.3e)\n", NlMaxiter, dmax);
		ierr = 1;
	}
	}
	if (nlf && nlconv && !ierr) {
		fprintf(fp_log, "  nonlinear (B-H) : converged\n");
	}
	free(nucell); free(bvr); free(bvi);
	free(aay); free(aaf); free(aayp); free(aafp); free(aady); free(aadf);

	// Z_loop = inv(Y)
	double *zr = (double *)malloc((size_t)np * np * sizeof(double));
	double *zi = (double *)malloc((size_t)np * np * sizeof(double));
	if (mat_inverse_c(yr, yi, zr, zi, np)) {
		fprintf(fp_log, "*** admittance matrix is singular; R(f)/L(f) is not available\n");
		ierr = 1;
	}
	else {
		for (int i = 0; i < np * np; i++) {
			Rfmat[i] = zr[i];
			Lfmat[i] = zi[i] / omega;
		}
		symmetrize(Rfmat, np, fp_log, "R(f)");
		symmetrize(Lfmat, np, fp_log, "L(f)");
		HaveF = 1;
	}

	free(zr); free(zi);
	free(yr); free(yi);
	free(br); free(bi); free(ar); free(ai);
	free(svol);
	free(fix);
	crs_free(&K);
	crs_free(&M);

	return ierr;
}


static int solve_one(FILE *fp_log)
{
	const int n = (int)num_node();
	const int np = NPort;
	int ierr = 0;

	// 辺要素の自己検証 (3 次元渦電流の基盤。他の解析とは独立に走る)
	if (Analysis & ANALYSIS_E) {
		// 3 次以上の四面体と 2 次以上の (完全 Lagrange の) 六面体・角柱は eddy3d_pr.c
		ierr |= ((MeshMode && ((TetOrder >= 3) || ((TetOrder >= 2) && ((NHex + NPrism) > 0)))) ? solve_edge_test_pr(fp_log)
		      : (MeshMode && (TetOrder == 2)) ? solve_edge_test_p2(fp_log) : solve_edge_test(fp_log));
		if (Analysis == ANALYSIS_E) return ierr;
	}

	// 節点要素 (P1 / P2) の自己検証。多項式再現の恒等式なので機械精度で合う
	if (Analysis & ANALYSIS_P) {
		ierr |= solve_nodal_test(fp_log);
		if (Analysis == ANALYSIS_P) return ierr;
	}

	crs_t A;
	crs_alloc_dist(&A);		// mpisolve = 2 なら自分の行ブロックだけ (それ以外は crs_alloc と同じ)
	{
		int64_t r0, r1;
		if (crs_dist_info(&A, &r0, &r1)) {
			fprintf(fp_log, "\n=== MPI : the matrix is distributed by rows (mpisolve = 2) : rank %d holds rows %lld..%lld "
				"of %d (%lld nonzeros) ===\n", mpi_rank(), (long long)r0, (long long)(r1 - 1), n, (long long)A.nnz);
		}
	}

	unsigned char *fix = (unsigned char *)malloc(n * sizeof(unsigned char));
	double *phi  = (double *)malloc(n * sizeof(double));
	double *work = (double *)malloc(n * sizeof(double));
	double *diag = (double *)malloc(n * sizeof(double));
	double *q = (double *)malloc((np + 1) * sizeof(double));

	const double scale = ((Tline && (TlineLength > 0)) ? (1 / TlineLength) : 1);

	// mode 0 : 静電界 (容量), 1 : 真空静電界 (インダクタンス), 2 : 定常電流界 (抵抗)
	for (int mode = 0; mode < 3; mode++) {
		if ((mode == 0) && !(Analysis & ANALYSIS_C)) continue;
		if ((mode == 1) && !(Analysis & ANALYSIS_L)) continue;
		if ((mode == 2) && !(Analysis & ANALYSIS_R)) continue;

		// 熱連成 : 定常電流界は σ(T) を要素ごとに更新しながら thermal.c が解く。
		// 熱連成でないときの経路には一切触れない (既存の答えをビット単位で保つ)
		if ((mode == 2) && ThermalOn) {
			double g = 0, x = 0;
			const int e = solve_thermal(fp_log, &g, &x);
			ierr |= e;
			if (!e) {
				Gmat[0] = g;
				HaveR = 1;
			}
			continue;
		}

		const char *modename = ((mode == 0) ? "electrostatic" :
		                        (mode == 1) ? "electrostatic (vacuum)" : "conduction");
		fprintf(fp_log, "\n=== %s analysis ===\n", modename);
		fprintf(fp_log, "  %-10s %8s %13s\n", "port", "iter", "residual");
		fflush(fp_log);

		assemble(&A, mode);

		// Dirichlet 節点 : 導体節点 + (電流界で) 行が空の節点
		crs_diag(&A, diag);
		double dmax = 0;
		for (int i = 0; i < n; i++) {
			if (diag[i] > dmax) dmax = diag[i];
		}
		if (dmax <= 0) {
			fprintf(fp_log, "*** no active element (check material sigma / geometry)\n");
			ierr = 1;
			break;
		}
		int64_t nfix = 0, ninactive = 0;
		for (int i = 0; i < n; i++) {
			const int active = (diag[i] > (1e-12 * dmax));
			if (NodeConductor[i] >= 0) {
				fix[i] = 1;
				nfix++;
			}
			else if (!active) {
				fix[i] = 1;			// 要素の付いていない節点 (σ=0 領域)
				ninactive++;
			}
			else {
				fix[i] = 0;
			}
		}
		fprintf(fp_log, "  conductor nodes = %lld, inactive nodes = %lld, unknowns = %lld\n",
			(long long)nfix, (long long)ninactive, (long long)(n - nfix - ninactive));
		fflush(fp_log);

		double *mat = ((mode == 0) ? Cmat : (mode == 1) ? Lmat : Gmat);

		for (int k = 1; k <= np; k++) {
			char label[BUFSIZ];
			sprintf(label, "port%d", k);
			double energy = 0;
			const int iter = solve_port(&A, fix, k, phi, work, q, &energy, fp_log, label);
			if (iter < 0) {
				fprintf(fp_log, "*** solver did not converge (port %d)\n", k);
				ierr = 1;
			}
			for (int j = 1; j <= np; j++) {
				mat[((k - 1) * np) + (j - 1)] = q[j] / Volt * scale;
			}
			// 場の出力 (fieldout = 1 のときだけ溜める)
			if (FieldOut) {
				char nm[BUFSIZ];
				const char *mn = ((mode == 0) ? "C" : (mode == 1) ? "L" : "R");
				sprintf(nm, "phi_%s_port%d", mn, k);
				field_add_node(nm, phi);
				sprintf(nm, "E_%s_port%d", mn, k);
				field_add_grad(nm, phi, 0);
			}

			// 自己項のエネルギーによる検算 (2W/V^2 = C[k][k])
			if (mode != 2) {
				const double cself = 2 * energy / (Volt * Volt) * scale;
				fprintf(fp_log, "  port%d : diag = %13.6e, 2W/V^2 = %13.6e (check)\n",
					k, mat[((k - 1) * np) + (k - 1)], cself);
			}
			fflush(fp_log);
		}

		symmetrize(mat, np, fp_log, modename);

		if      (mode == 0) HaveC = 1;
		else if (mode == 1) HaveL = 1;
		else                HaveR = 1;
	}

	// L = μ0 ε0 inv(C0)   (C0 は真空の単位長容量行列、Lmat には C0 が入っている)
	if (HaveL) {
		double *c0 = (double *)malloc((size_t)np * np * sizeof(double));
		memcpy(c0, Lmat, (size_t)np * np * sizeof(double));
		if (mat_inverse(c0, Lmat, np)) {
			fprintf(fp_log, "*** vacuum capacitance matrix is singular; L is not available\n");
			HaveL = 0;
		}
		else {
			for (int i = 0; i < np * np; i++) {
				Lmat[i] *= MU0 * EPS0;
			}
		}
		free(c0);
	}

	// R = inv(G)
	if (HaveR) {
		if (mat_inverse(Gmat, Rmat, np)) {
			fprintf(fp_log, "*** conductance matrix is singular; R is not available\n");
			memset(Rmat, 0, (size_t)np * np * sizeof(double));
		}
	}

	free(fix);
	free(phi);
	free(work);
	free(diag);
	free(q);
	crs_free(&A);

	// 静磁場解析 (DC インダクタンス)。行列を作り直すので静電界の後始末をしてから呼ぶ
	if (!ierr && (Analysis & ANALYSIS_M)) {
		ierr |= solve_magnetostatic(fp_log);
	}

	// 時間調和渦電流解析 (表皮効果を含む R(f), L(f))
	if (!ierr && (Analysis & ANALYSIS_F)) {
		if (ThermalOn) {
			// 交流の発熱 (断面 2 次元) : F を σ(T) を要素ごとに更新しながら thermal.c が解く
			// (熱連成でないときの経路には触れない)
			double g = 0, x = 0;
			const int e = solve_thermal(fp_log, &g, &x);
			ierr |= e;
			if (!e) {
				Rfmat[0] = 1 / g;
				Lfmat[0] = x / (2 * PI * Freq);
				HaveF = 1;
			}
		}
		else {
			ierr |= solve_eddy(fp_log);
		}
	}

	// 3 次元渦電流解析 (A-φ、辺要素)
	if (!ierr && (Analysis & ANALYSIS_A)) {
		if (ThermalOn) {
			// 交流の発熱 : A-φ を σ(T) を要素ごとに更新しながら thermal.c が解く
			// (熱連成でないときの経路には触れない)
			double g = 0, x = 0;
			const int e = solve_thermal(fp_log, &g, &x);
			ierr |= e;
			if (!e) {
				Rfmat[0] = 1 / g;
				Lfmat[0] = x / (2 * PI * Freq);
				HaveF = 1;
			}
		}
		else if (TransientOn && MeshMode && (TetOrder == 2) && ((NHex + NPrism) == 0)) {
			ierr |= solve_eddy3d_p2_transient(fp_log);
		}
		else if (TransientOn && MeshMode && ((TetOrder >= 2) || pr_p1_nonlinear())) {
			// 3 次以上の四面体と、2 次以上の (完全 Lagrange の) 六面体・角柱。1 次の六面体・角柱の格子の非線形 (J-A / B-H) も
			// 次数 1 の同じ経路 (積分点ごとの状態)
			ierr |= solve_eddy3d_pr_transient(fp_log);
		}
		else if (TransientOn) {
			ierr |= solve_eddy3d_transient(fp_log);
		}
		else if (MeshMode && ((TetOrder >= 3) || ((TetOrder >= 2) && ((NHex + NPrism) > 0)))) {
			// 3 次以上の四面体と、2 次以上の (完全 Lagrange の) 六面体・角柱
			ierr |= solve_eddy3d_pr(fp_log);
		}
		else if (MeshMode && (TetOrder == 2)) {
			ierr |= solve_eddy3d_p2(fp_log);
		}
		else {
			ierr |= solve_eddy3d(fp_log);
		}
	}

	return ierr;
}


// ---- 周波数掃引 ----
//
// 分散材料の展開 (epsr / tand / eps6) は Freq に依存するので、各点で
// material_freq() を呼び直してから解き直す。結果は 1 行 1 周波数で
// ofe_sweep.csv に出す。ofe.out / rlc.csv には**最後の周波数**の結果が入る
// (掃引が無い従来の挙動と同じ形にするため)。
// ヘッダと行は**同じ列構成**でなければならない。掃引の途中で解析が失敗して
// Have* が変わると CSV が不揃いになるので、最初の点の構成を固定して使い、
// 以降の点で構成が変わったらその列を空にして警告する。
typedef struct { int c, r, f; } sweepcol_t;

static void sweep_header(FILE *fp, sweepcol_t k)
{
	const int np = NPort;

	fprintf(fp, "frequency [Hz]");
	if (k.c) for (int i = 0; i < np; i++) for (int j = 0; j < np; j++) fprintf(fp, ",C(%d;%d)", i + 1, j + 1);
	if (k.r) for (int i = 0; i < np; i++) for (int j = 0; j < np; j++) fprintf(fp, ",G(%d;%d)", i + 1, j + 1);
	if (k.r) for (int i = 0; i < np; i++) for (int j = 0; j < np; j++) fprintf(fp, ",R(%d;%d)", i + 1, j + 1);
	if (k.f) for (int i = 0; i < np; i++) for (int j = 0; j < np; j++) fprintf(fp, ",Rf(%d;%d)", i + 1, j + 1);
	if (k.f) for (int i = 0; i < np; i++) for (int j = 0; j < np; j++) fprintf(fp, ",Lf(%d;%d)", i + 1, j + 1);
	fprintf(fp, "\n");
}


static void sweep_block(FILE *fp, int want, int have, const double *m, int nn)
{
	if (!want) return;
	for (int i = 0; i < nn; i++) {
		if (have) fprintf(fp, ",%.8e", m[i]);
		else      fprintf(fp, ",");
	}
}


static void sweep_row(FILE *fp, sweepcol_t k, FILE *fp_log)
{
	const int nn = NPort * NPort;

	if (((k.c != 0) != (HaveC != 0)) || ((k.r != 0) != (HaveR != 0))
	 || ((k.f != 0) != (HaveF != 0))) {
		fprintf(fp_log, "*** warning : the available matrices changed at %.6e Hz; "
			"those columns are left empty in %s\n", Freq, FN_sweep);
	}
	fprintf(fp, "%.8e", Freq);
	sweep_block(fp, k.c, HaveC, Cmat, nn);
	sweep_block(fp, k.r, HaveR, Gmat, nn);
	sweep_block(fp, k.r, HaveR, Rmat, nn);
	sweep_block(fp, k.f, HaveF, Rfmat, nn);
	sweep_block(fp, k.f, HaveF, Lfmat, nn);
	fprintf(fp, "\n");
}


int solve(FILE *fp_log)
{
	int ierr = 0;
	const int msize = mpi_size();
	const int mrank = mpi_rank();

	/*
	MPI (msize > 1) が並列化するのは**周波数掃引の点の分散だけ**。
	各点の演算は直列と完全に同一なので、mpirun -np N の出力
	(ofe_sweep.csv / rlc.csv / ofe.out / ofe_series.h5) は直列と
	バイト単位で一致する — これが rlc_check.sh の検証の恒等式になる。

	掃引の無い入力と fieldout = 1 は分散する仕事が無い / 場が解いた rank に
	しか無いので、黙って直列で走らず入力エラーにする。
	*/
	if ((msize > 1) && (MpiSolve == 2)
	 && (Direct || (Analysis & ~(ANALYSIS_C | ANALYSIS_L | ANALYSIS_R | ANALYSIS_A | ANALYSIS_M | ANALYSIS_F)))) {
		fprintf(fp_log, "*** mpisolve = 2 (the matrix distributed over the ranks) supports analysis C / L / R / M / F / A "
			"and thermal (not the self tests E / P, no direct = 1); use mpisolve = 1\n");
		return 1;
	}
	if ((msize > 1) && MpiSolve) {
		// 行ブロックの行列ベクトル積 (以下は全 rank が直列と同じ流れを重複して走る)
		mpi_rows_enable(1);
		fprintf(fp_log, "\n=== MPI : matrix-vector products split by rows over %d ranks "
			"(everything else is repeated on every rank in the serial order) ===\n", msize);
	}
	else if (msize > 1) {
		if (NFreqSweep < 2) {
			fprintf(fp_log, "*** MPI parallelism distributes the frequency-sweep "
				"points, but this input has fewer than two; run without mpirun, "
				"or add sweep points\n");
			return 1;
		}
		if (FieldOut) {
			fprintf(fp_log, "*** fieldout = 1 cannot be combined with MPI "
				"(the fields of a sweep point exist only on the rank that "
				"solved it); run with -np 1\n");
			return 1;
		}
	}

	// 系列の出力 (hdf5 = 1)。掃引しない実行でも 1 点の系列として書く
	// (ヒステリシスの履歴は solve_one() の中で追記されるので、ここで開けておく)。
	// MPI では rank 0 だけが書く
	if ((mrank == 0) && h5_open(fp_log)) return 1;

	if (NFreqSweep < 1) {
		ierr = solve_one(fp_log);
		if (!ierr) {
			ierr |= h5_add_sweep(Freq);
			// 過渡解析の時系列 (ofe_transient.csv の列) を /transient に
			if (TransientOn && (mrank == 0)) ierr |= h5_add_transient();
			// 電流掃引 (ヒステリシス) は既にステップ毎に追記済み。
			// ここで足すと最後のステップが二重に入るので、まだ 1 点も
			// 書いていないときだけ「1 点の系列」として書く
			if (h5_nfield() == 0) ierr |= h5_add_field(Freq);
			ierr |= field_write(fp_log);
		}
		h5_close(fp_log);
		field_free();
		return ierr;
	}

	FILE *fp = (((mrank == 0) || MpiSolve) ? mpi_fopen_out(FN_sweep, "w") : NULL);
	if (((mrank == 0) || MpiSolve) && (fp == NULL)) {
		fprintf(fp_log, "*** %s open error\n", FN_sweep);
		h5_close(NULL);
		return 1;
	}
	fprintf(fp_log, "\n=== frequency sweep (%d points%s) ===\n", NFreqSweep,
		(((msize > 1) && !MpiSolve) ? ", MPI" : ""));

	if ((msize > 1) && !MpiSolve) {
		/*
		点 q を rank (q % msize) に round-robin で配る。各 rank は自分の点を
		解いて結果 (Have* と結果行列 8 本) を 1 本の pack にして rank 0 へ送り、
		rank 0 が**点の順に**受けて ofe_sweep.csv / ofe_series.h5 に書く。
		最後の点の結果をグローバルに戻すので、rlc.csv / ofe.out は直列と同じく
		「最後の周波数」の値になる。

		**エラーでも打ち切らないこと。** 解けなかった点も pack (ierr つき) を
		必ず送る。ある rank だけが途中で止まると、rank 0 の受信と食い違って
		デッドロックする (送る側と受ける側の回数は常に一致させる)。
		*/
		const int nn = NPort * NPort;
		const int npk = 8 + (8 * nn);
		double *packs = (double *)malloc((size_t)NFreqSweep * npk * sizeof(double));
		memset(packs, 0, (size_t)NFreqSweep * npk * sizeof(double));

		for (int q = 0; q < NFreqSweep; q++) {
			if ((q % msize) != mrank) continue;
			double *pk = &packs[(size_t)q * npk];
			Freq = FreqSweep[q];
			int perr = 0;
			if (material_freq()) {
				fprintf(fp_log, "*** material expansion failed at %.6e Hz\n", Freq);
				perr = 1;
			}
			else {
				fprintf(fp_log, "\n--- sweep point %d/%d : %.6e [Hz] (rank %d) ---\n",
					q + 1, NFreqSweep, Freq, mrank);
				fflush(fp_log);
				field_free();
				HaveC = HaveL = HaveR = HaveM = HaveF = 0;
				perr = solve_one(fp_log);
			}
			pk[0] = perr;
			pk[1] = (HaveC != 0);
			pk[2] = (HaveL != 0);
			pk[3] = (HaveR != 0);
			pk[4] = (HaveM != 0);
			pk[5] = (HaveF != 0);
			pk[6] = (HavePfe != 0);
			pk[7] = Freq;
			memcpy(&pk[8 + (0 * nn)], Cmat,   (size_t)nn * sizeof(double));
			memcpy(&pk[8 + (1 * nn)], Lmat,   (size_t)nn * sizeof(double));
			memcpy(&pk[8 + (2 * nn)], Gmat,   (size_t)nn * sizeof(double));
			memcpy(&pk[8 + (3 * nn)], Rmat,   (size_t)nn * sizeof(double));
			memcpy(&pk[8 + (4 * nn)], Mmat,   (size_t)nn * sizeof(double));
			memcpy(&pk[8 + (5 * nn)], Rfmat,  (size_t)nn * sizeof(double));
			memcpy(&pk[8 + (6 * nn)], Lfmat,  (size_t)nn * sizeof(double));
			memcpy(&pk[8 + (7 * nn)], Pfemat, (size_t)nn * sizeof(double));
			if (mrank != 0) mpi_send_dbl(pk, npk, 0, q);
			if (perr) ierr = 1;
		}

		if (mrank == 0) {
			sweepcol_t kcol = {0, 0, 0};
			int nhead = 0;
			int nfail = 0;
			for (int q = 0; q < NFreqSweep; q++) {
				double *pk = &packs[(size_t)q * npk];
				const int owner = q % msize;
				if (owner != 0) mpi_recv_dbl(pk, npk, owner, q);
				// **失敗した点の報告は ierr の判定より前に置く。** rank 0 自身の
				// 点が先に失敗して ierr が立っていると、どの点・どの周波数で
				// 落ちたのかが一度も出ない (レビューの実測)
				if (pk[0] != 0) {
					if (!nfail) {
						// rank 0 は ofe.log に直接書いているので、別ファイルを
						// 案内するのは owner != 0 のときだけ
						if (owner != 0) {
							fprintf(fp_log, "*** solve failed at sweep point %d "
								"(%.6e Hz, rank %d); see ofe_rank%d.log for the "
								"details\n", q + 1, pk[7], owner, owner);
						} else {
							fprintf(fp_log, "*** solve failed at sweep point %d "
								"(%.6e Hz, rank 0); see the messages above\n",
								q + 1, pk[7]);
						}
					}
					nfail++;
					ierr = 1;
					continue;
				}
				if (ierr) continue;			// 受信は最後まで続ける (上の注意)
				Freq = pk[7];
				HaveC = (pk[1] != 0);
				HaveL = (pk[2] != 0);
				HaveR = (pk[3] != 0);
				HaveM = (pk[4] != 0);
				HaveF = (pk[5] != 0);
				HavePfe = (pk[6] != 0);
				memcpy(Cmat,   &pk[8 + (0 * nn)], (size_t)nn * sizeof(double));
				memcpy(Lmat,   &pk[8 + (1 * nn)], (size_t)nn * sizeof(double));
				memcpy(Gmat,   &pk[8 + (2 * nn)], (size_t)nn * sizeof(double));
				memcpy(Rmat,   &pk[8 + (3 * nn)], (size_t)nn * sizeof(double));
				memcpy(Mmat,   &pk[8 + (4 * nn)], (size_t)nn * sizeof(double));
				memcpy(Rfmat,  &pk[8 + (5 * nn)], (size_t)nn * sizeof(double));
				memcpy(Lfmat,  &pk[8 + (6 * nn)], (size_t)nn * sizeof(double));
				memcpy(Pfemat, &pk[8 + (7 * nn)], (size_t)nn * sizeof(double));
				if (!nhead) {
					kcol.c = (HaveC != 0);
					kcol.r = (HaveR != 0);
					kcol.f = (HaveF != 0);
					sweep_header(fp, kcol);
					nhead = 1;
				}
				sweep_row(fp, kcol, fp_log);
				fflush(fp);
				ierr |= h5_add_sweep(Freq);
				// h5_add_field は場が無ければ何も書かない
				// (fieldout = 1 との併用は上で弾いてある)
				ierr |= h5_add_field(Freq);
			}
		}
		mpi_wait_sends();		// Isend の完了前に packs を解放しないこと
		free(packs);
		if (fp != NULL) fclose(fp);

		if (!ierr && (mrank == 0)) {
			fprintf(fp_log, "\n  frequency sweep written to %s (%d points, "
				"%d MPI processes)\n", FN_sweep, NFreqSweep, msize);
		}
		if (mrank == 0) h5_close(fp_log);
		field_free();

		return ierr;
	}

	sweepcol_t kcol = {0, 0, 0};
	int nhead = 0;
	for (int q = 0; q < NFreqSweep; q++) {
		Freq = FreqSweep[q];
		if (material_freq()) {
			fprintf(fp_log, "*** material expansion failed at %.6e Hz\n", Freq);
			ierr = 1;
			break;
		}
		fprintf(fp_log, "\n--- sweep point %d/%d : %.6e [Hz] ---\n",
			q + 1, NFreqSweep, Freq);
		fflush(fp_log);

		// 各点で溜め直す (場は最後の周波数のものだけ残す)。
		// HaveS は setup() が 1 度だけ設定する DC 直列抵抗なので消さない
		field_free();
		HaveC = HaveL = HaveR = HaveM = HaveF = 0;

		ierr = solve_one(fp_log);
		if (ierr) break;

		if (!nhead) {
			kcol.c = (HaveC != 0);
			kcol.r = (HaveR != 0);
			kcol.f = (HaveF != 0);
			sweep_header(fp, kcol);
			nhead = 1;
		}
		sweep_row(fp, kcol, fp_log);
		fflush(fp);

		// **場はこの時点でだけ手元にある** (次の点の頭で field_free() される)。
		// 溜めずにその場で追記するので、メモリの使い方は従来と変わらない
		ierr |= h5_add_sweep(Freq);
		ierr |= h5_add_field(Freq);
		if (ierr) break;
	}
	fclose(fp);

	if (!ierr) {
		fprintf(fp_log, "\n  frequency sweep written to %s (%d points)\n",
			FN_sweep, NFreqSweep);
		ierr |= field_write(fp_log);
	}
	h5_close(fp_log);
	field_free();

	return ierr;
}


/*
交流の発熱 (thermal + analysis = F、断面 2 次元の 1 次三角形)。

要素ごとの σ (sigelem[NTri]、熱連成の σ(T)) で F を単一ポートで解き、電流が cur [A]
(波高値) になるよう縮尺した場の時間平均の損失 (単位長あたり) を節点荷重 f にする:
	p(x) = ½ σ |V' - jω Az(x)|²,   f_i = ∫ p N_i dA
V' は導体 1 で s (縮尺)、基準導体 0 で 0。被積分関数は λ の 3 次なので Duffy の 9 点で
厳密。系の行から Σ_i f_i = ½ Re(V' I*) = ½ |I|² R(f) が厳密に出る (aᴴ(K + jωM)a = V' aᴴb
の実部と虚部。損失を組み立てと同じ質量行列の上で積分するのでそのまま成り立つ)。
*zr, *zi = Z = R + jX [Ω/m]。
*/
// 交流の発熱の非線形 (B-H) : 要素ごとの ν と Anderson の作業域。熱連成の回をまたいで持ち越す (最初の呼び出しで作る)
static double *H2Nu, *H2Aay, *H2Aaf, *H2Aayp, *H2Aafp, *H2Aady, *H2Aadf, *H2Bvr, *H2Bvi;
static int H2NlTotal;


int eddy2d_heat_solve(FILE *fp_log, const double *sigelem, double cur, double *f,
	double *zr, double *zi, int it)
{
	const int n = (int)num_node();
	const double omega = 2 * PI * Freq;
	int ierr = 0;
	// 非線形 (B-H) : 3 次元の eddy3d_heat_solve と同じく、熱連成の各回の σ(T) で実効磁気抵抗率の反復 (solve_eddy の F と
	// 同じ H を独立変数にした更新と ln ν の Anderson 加速) を収束するまで解く。ν は前の回の値から始める
	int nlf = 0;
	for (int m = 0; m < NMaterial; m++) {
		if (Material[m].nbh[0] > 0) nlf = 1;
	}
	const int64_t nelc = num_cell();
	if (nlf && (H2Nu == NULL)) {
		H2Nu  = (double *)malloc((size_t)nelc * sizeof(double));
		H2Aay = (double *)malloc((size_t)nelc * sizeof(double));
		H2Aaf = (double *)malloc((size_t)nelc * sizeof(double));
		H2Aayp = (double *)malloc((size_t)nelc * sizeof(double));
		H2Aafp = (double *)malloc((size_t)nelc * sizeof(double));
		H2Aady = (double *)malloc((size_t)nelc * AAM * sizeof(double));
		H2Aadf = (double *)malloc((size_t)nelc * AAM * sizeof(double));
		H2Bvr = (double *)malloc((size_t)nelc * 3 * sizeof(double));
		H2Bvi = (double *)malloc((size_t)nelc * 3 * sizeof(double));
		for (int64_t c = 0; c < nelc; c++) {
			elem_t el;
			elem_get(c, &el);
			const material_t *mt = &Material[el.mat];
			H2Nu[c] = ((mt->nbh[0] > 0) ? bh_nu(mt, 0, 0) : (1 / (MU0 * mt->mu6[0] * temp_factor(mt, 1))));
		}
		H2NlTotal = 0;
		fprintf(fp_log, "  nonlinear (B-H) : effective reluctivity per element, solved to convergence in every thermal "
			"iteration (maxiter = %d, tol = %.1e); the thermal tolerance is at least 10 x that\n", NlMaxiter, NlTol);
	}

	crs_t K, M;
	crs_alloc(&K);
	crs_alloc(&M);
	assemble_mass_tri_sig(&M, sigelem);

	unsigned char *fix = (unsigned char *)calloc((size_t)n, sizeof(unsigned char));
	double *br = (double *)calloc((size_t)n, sizeof(double));
	double *bi = (double *)calloc((size_t)n, sizeof(double));
	double *ar = (double *)calloc((size_t)n, sizeof(double));
	double *ai = (double *)calloc((size_t)n, sizeof(double));

	// b_i = ∫ σ N_i V' dA (V' = 1 [V/m] を導体 1 に)、S_1 = ∫_1 σ dA
	double s1 = 0;
	for (int e = 0; e < NTri; e++) {
		if (TriCond[e] != 1) continue;
		const double sg = sigelem[e];
		for (int l = 0; l < 3; l++) br[Tri[(e * 3) + l]] += sg * TriArea[e] / 3;
		s1 += sg * TriArea[e];
	}
	char label[64];
	sprintf(label, "eddy%d", it);
	int aan = 0;		// Anderson の履歴は σ(T) が変わるたびに捨てる
	double sr = 0, si = 0, ivr = 0, ivi = 0, den = 0;
	for (int nlit = 1; !ierr; nlit++) {
	assemble_nu_tri(&K, (nlf ? H2Nu : NULL));
	for (int i = 0; i < n; i++) ar[i] = ai[i] = 0;
	if (solver_cocg(&K, &M, omega, br, bi, ar, ai, fix,
		Solver.maxiter, Solver.nout, Solver.converg, NULL, label) < 0) {
		fprintf(fp_log, "*** thermal : the eddy-current solve did not converge (iteration %d)\n", it);
		ierr = 1;
	}

	// I_V = S_1 - jω ∫_1 σ Az dA  (V' = 1)
	double qr = 0, qi = 0;
	for (int e = 0; e < NTri; e++) {
		if (TriCond[e] != 1) continue;
		const double w = sigelem[e] * TriArea[e] / 3;
		for (int l = 0; l < 3; l++) {
			qr += w * ar[Tri[(e * 3) + l]];
			qi += w * ai[Tri[(e * 3) + l]];
		}
	}
	ivr = s1 + (omega * qi);
	ivi = -omega * qr;
	den = (ivr * ivr) + (ivi * ivi);
	if (!ierr && (den <= 0)) {
		fprintf(fp_log, "*** thermal : no current flows through the port\n");
		ierr = 1;
	}
	if (ierr) break;
	// s = cur / I_V (V' の縮尺)。Az も s 倍
	sr = cur * ivr / den;
	si = -cur * ivi / den;
	for (int i = 0; i < n; i++) {
		const double xr = ar[i], xi = ai[i];
		ar[i] = (sr * xr) - (si * xi);
		ai[i] = (sr * xi) + (si * xr);
	}
	if (!nlf) break;
	// ν の更新 (solve_eddy の F と同じ式 : |B|² = |Re B|² + |Im B|²、H を独立変数にした ν = |H|/B(|H|))
	field_cell_grad(ar, 1, H2Bvr);
	field_cell_grad(ai, 1, H2Bvi);
	double dmax = 0;
	for (int64_t c = 0; c < nelc; c++) {
		elem_t el;
		elem_get(c, &el);
		const material_t *mt = &Material[el.mat];
		H2Aay[c] = log(H2Nu[c]);
		if (mt->nbh[0] <= 0) {
			H2Aaf[c] = 0;
			continue;
		}
		double b2 = 0;
		for (int d = 0; d < 3; d++) {
			b2 += (H2Bvr[(c * 3) + d] * H2Bvr[(c * 3) + d]) + (H2Bvi[(c * 3) + d] * H2Bvi[(c * 3) + d]);
		}
		const double hm = H2Nu[c] * sqrt(b2);
		const double nun = ((hm > 0) ? (hm / bh_b_of_h(mt, 0, hm)) : bh_nu(mt, 0, 0));
		const double dn = fabs(nun - H2Nu[c]) / nun;
		if (dn > dmax) dmax = dn;
		H2Aaf[c] = log(nun) - H2Aay[c];
	}
	H2NlTotal++;
	if (dmax < NlTol) break;
	if (nlit >= NlMaxiter) {
		fprintf(fp_log, "*** thermal : the effective reluctivity did not converge in %d iterations "
			"(thermal iteration %d, max dnu/nu = %.3e)\n", NlMaxiter, it, dmax);
		ierr = 1;
		break;
	}
	anderson_step(nelc, H2Aay, H2Aaf, H2Aayp, H2Aafp, H2Aady, H2Aadf, AAM, &aan, NlRelax);
	for (int64_t c = 0; c < nelc; c++) H2Nu[c] = exp(H2Aay[c]);
	}
	if (!ierr && nlf) {
		fprintf(fp_log, "  thermal iteration %d : effective reluctivity converged (%d B-H iterations so far)\n",
			it, H2NlTotal);
	}
	if (!ierr) {
		*zr = ivr / den;						// Z = 1 / Y (V' = 1 のとき Y = I_V)
		*zi = -ivi / den;
		// 損失の節点荷重 : Duffy 変換 + 各方向 3 点 (9 点、λ の 5 次まで厳密)
		const double gp[3] = {0.5 - (0.5 * 0.77459666924148337704), 0.5,
		                      0.5 + (0.5 * 0.77459666924148337704)};
		const double gw[3] = {5.0 / 18, 8.0 / 18, 5.0 / 18};
		for (int i = 0; i < n; i++) f[i] = 0;
		for (int e = 0; e < NTri; e++) {
			const double sg = sigelem[e];
			if ((TriCond[e] < 0) || (sg <= 0)) continue;
			const int32_t *nd = &Tri[e * 3];
			const double vr = ((TriCond[e] == 1) ? sr : 0), vi = ((TriCond[e] == 1) ? si : 0);
			for (int a = 0; a < 3; a++) {
			for (int b = 0; b < 3; b++) {
				const double u = gp[a], v = gp[b];
				const double lam[3] = {1 - u - (v * (1 - u)), u, v * (1 - u)};
				const double w = gw[a] * gw[b] * (1 - u) * 2 * TriArea[e];	// ∫ dA = 2A ∫∫ (1-u) du dv
				double zr8 = 0, zi8 = 0;
				for (int l = 0; l < 3; l++) {
					zr8 += lam[l] * ar[nd[l]];
					zi8 += lam[l] * ai[nd[l]];
				}
				// V' - jω Az = (vr + ω zi8) + j (vi - ω zr8)
				const double jr = vr + (omega * zi8), ji = vi - (omega * zr8);
				const double pq = 0.5 * sg * ((jr * jr) + (ji * ji)) * w;
				for (int l = 0; l < 3; l++) f[nd[l]] += pq * lam[l];
			}
			}
		}
	}
	crs_free(&K);
	crs_free(&M);
	free(fix); free(br); free(bi); free(ar); free(ai);
	return ierr;
}
