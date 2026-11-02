/*
assemble.c

8 節点 3 次補間 (trilinear) 6 面体要素による剛性行列の作成。

弱形式 : ∫ c ∇w・∇φ dV = 0
	c = ε0 εr  -> 静電界 (容量)
	c = ε0     -> 真空静電界 (TEM インダクタンス)
	c = σ      -> 定常電流界 (抵抗・コンダクタンス)

要素は軸に平行な直方体なのでヤコビアンは対角一定。2x2x2 Gauss 積分は
この要素に対して厳密。
*/

#include "fem.h"
#include "fem_prototype.h"

// 局所節点番号 l = ((li * 2) + lj) * 2 + lk  (li, lj, lk = 0 or 1)
#define LI(l) (((l) >> 2) & 1)
#define LJ(l) (((l) >> 1) & 1)
#define LK(l) ((l) & 1)


// 直方体要素 (dx x dy x dz) の ∫ (∇Ni)^T C (∇Nj) dV
// C は対称テンソル (成分順 xx, yy, zz, xy, yz, zx)。等方性なら対角のみ
void element_matrix(double dx, double dy, double dz, const double c[6], double ke[8][8])
{
	const double gp = 1 / sqrt(3.0);		// 2 点 Gauss 積分点
	const double detj = (dx * dy * dz) / 8;
	const double rx = 2 / dx;
	const double ry = 2 / dy;
	const double rz = 2 / dz;

	for (int l = 0; l < 8; l++) {
		for (int m = 0; m < 8; m++) {
			ke[l][m] = 0;
		}
	}

	for (int g = 0; g < 8; g++) {
		const double xi  = (LI(g) ? +gp : -gp);
		const double eta = (LJ(g) ? +gp : -gp);
		const double zta = (LK(g) ? +gp : -gp);

		double dn[8][3];
		for (int l = 0; l < 8; l++) {
			const double si = (LI(l) ? +1.0 : -1.0);
			const double sj = (LJ(l) ? +1.0 : -1.0);
			const double sk = (LK(l) ? +1.0 : -1.0);
			dn[l][0] = si * (1 + (sj * eta)) * (1 + (sk * zta)) / 8 * rx;
			dn[l][1] = sj * (1 + (sk * zta)) * (1 + (si * xi )) / 8 * ry;
			dn[l][2] = sk * (1 + (si * xi )) * (1 + (sj * eta)) / 8 * rz;
		}

		for (int l = 0; l < 8; l++) {
			for (int m = 0; m < 8; m++) {
				ke[l][m] += ((c[0] * dn[l][0] * dn[m][0])
				          +  (c[1] * dn[l][1] * dn[m][1])
				          +  (c[2] * dn[l][2] * dn[m][2])
				          +  (c[3] * ((dn[l][0] * dn[m][1]) + (dn[l][1] * dn[m][0])))
				          +  (c[4] * ((dn[l][1] * dn[m][2]) + (dn[l][2] * dn[m][1])))
				          +  (c[5] * ((dn[l][2] * dn[m][0]) + (dn[l][0] * dn[m][2])))) * detj;
			}
		}
	}
}


// 材料 m のセルに掛かる係数
//   mode = 0 : ε0 εr           静電界 (容量)
//        = 1 : ε0              真空静電界 (TEM インダクタンス)
//        = 2 : σ + ω ε0 εr tanδ  定常電流界 (導電損 + 誘電損)
//        = 3 : 1 / (μ0 μr)     静磁場 (DC インダクタンス)
// 温度 T での材料係数の倍率。kind = 0 : εr、1 : μr、2 : B-H、3 : 抵抗率 (σ = σ0 / f)。
// 表 (epstemptable 等) があれば (T, f) の折れ線 (範囲の外は端の値)、
// 無ければ 1 次係数 1 + α (T - T0)
// (tt, ff) の折れ線の x での値 (範囲の外は端の値)。tt は狭義単調増加
double table_interp(int nt, const double *tt, const double *ff, double x)
{
	if (x <= tt[0]) return ff[0];
	for (int k = 0; k < (nt - 1); k++) {
		if (x <= tt[k + 1]) {
			const double w = (x - tt[k]) / (tt[k + 1] - tt[k]);
			return ff[k] + (w * (ff[k + 1] - ff[k]));
		}
	}
	return ff[nt - 1];
}


double temp_factor_at(const material_t *mt, int kind, double t)
{
	const int nt = mt->ntemp[kind];
	if (nt > 0) return table_interp(nt, mt->tempT[kind], mt->tempF[kind], t);
	if (kind == 0) return 1 + (mt->epstempco * (t - mt->epstemp0));
	if (kind == 1) return 1 + (mt->mutempco  * (t - mt->mutemp0));
	if (kind == 2) return 1 + (mt->bhtempco * (t - mt->bhtemp0));
	return 1 + (mt->tempco * (t - mt->temp0));
}


// 動作温度 Temperature での倍率
double temp_factor(const material_t *mt, int kind)
{
	return temp_factor_at(mt, kind, Temperature);
}


void material_coef_pub(int m, int mode, double c[6])
{
	const material_t *mt = &Material[m];

	// εr / μr の温度依存はここ (読み出し時) で掛ける。
	// σ のように入力解釈で一度だけ掛ける方式は使えない : material_freq() が
	// 分散材料の εr を毎回展開し直して上書きするので、周波数掃引で補正が消える。
	// 読み出し時なら何度呼んでも累積しない
	const double fe = temp_factor(mt, 0);
	const double fm = temp_factor(mt, 1);

	for (int d = 0; d < 6; d++) c[d] = 0;

	if      (mode == 0) {
		for (int d = 0; d < 6; d++) c[d] = EPS0 * mt->eps6[d] * fe;
	}
	else if (mode == 1) {
		c[0] = c[1] = c[2] = EPS0;			// 真空 (等方性)
	}
	else if (mode == 2) {
		// 誘電損は等価導電率 σ_d = ω ε0 εr tanδ として扱う (frequency 未指定なら 0)
		const double omega = 2 * PI * Freq;
		for (int d = 0; d < 6; d++) {
			c[d] = omega * EPS0 * mt->eps6[d] * fe * mt->tand;
		}
		c[0] += mt->sigma;					// 導電率は等方性
		c[1] += mt->sigma;
		c[2] += mt->sigma;
	}
	else {
		// 磁気抵抗率テンソル ν = (μ0 μ~)^-1
		double mu[6], nu[6];
		for (int d = 0; d < 6; d++) mu[d] = MU0 * mt->mu6[d] * fm;
		if (tensor6_inverse(mu, nu)) {
			nu[0] = nu[1] = nu[2] = 1 / MU0;
			nu[3] = nu[4] = nu[5] = 0;
		}
		for (int d = 0; d < 6; d++) c[d] = nu[d];

		// mode 4 : 断面 2 次元の Az 定式化で使う「∇Az の基底での ν」
		//
		// エネルギーは ∫ B^T ν B / 2 で、B = ∇×(Az ê_t) は
		//   B_p = +∂Az/∂q,  B_q = -∂Az/∂p
		// なので ∇Az の基底では **面内の 2 成分が入れ替わり、非対角の符号が反転**する:
		//   c_pp = ν_qq,  c_qq = ν_pp,  c_pq = -ν_pq
		// ν をそのまま ∇Az に掛けると μ_p と μ_q を取り違える (実測: μ を
		// (3,2,5) にした平行平板線路で L が +138%)。等方性では一致するので
		// 等方性のケースだけでは絶対に検出できない。
		if (mode == 4) {
			const int t = ((Tline == 'X') ? 0 : (Tline == 'Y') ? 1 : 2);
			const int p = (t + 1) % 3, q = (t + 2) % 3;
			// 6 成分の並び : 0=xx 1=yy 2=zz 3=xy 4=yz 5=zx
			static const int off[3][3] = {{0, 3, 5}, {3, 1, 4}, {5, 4, 2}};
			const double npp = nu[off[p][p]], nqq = nu[off[q][q]], npq = nu[off[p][q]];
			for (int d = 0; d < 6; d++) c[d] = 0;
			c[off[p][p]] = nqq;
			c[off[q][q]] = npp;
			c[off[p][q]] = -npq;
			// 軸方向は Az が一様なので効かないが、行列が特異にならないよう
			// 面内の代表値を入れておく (構造格子は軸方向にも節点を持つ)
			c[off[t][t]] = ((npp + nqq) / 2);
		}
	}
}


// 直方体要素 (dx x dy x dz) の ∫Ni Nj dV
// (Ni Nj は各方向 2 次なので 2 点 Gauss で厳密)
static void element_mass(double dx, double dy, double dz, double me[8][8])
{
	const double gp = 1 / sqrt(3.0);
	const double detj = (dx * dy * dz) / 8;

	for (int l = 0; l < 8; l++) {
		for (int m = 0; m < 8; m++) {
			me[l][m] = 0;
		}
	}

	for (int g = 0; g < 8; g++) {
		const double xi  = (LI(g) ? +gp : -gp);
		const double eta = (LJ(g) ? +gp : -gp);
		const double zta = (LK(g) ? +gp : -gp);

		double nn[8];
		for (int l = 0; l < 8; l++) {
			const double si = (LI(l) ? +1.0 : -1.0);
			const double sj = (LJ(l) ? +1.0 : -1.0);
			const double sk = (LK(l) ? +1.0 : -1.0);
			nn[l] = (1 + (si * xi)) * (1 + (sj * eta)) * (1 + (sk * zta)) / 8;
		}

		for (int l = 0; l < 8; l++) {
			for (int m = 0; m < 8; m++) {
				me[l][m] += nn[l] * nn[m] * detj;
			}
		}
	}
}


// 渦電流項の質量行列 ∫σ N_i N_j dV
// σ は導体セルの導電率 (conductorsigma)。誘電体側の渦電流は無視する
void assemble_mass(crs_t *A)
{
	crs_zero(A);

	for (int i = 0; i < Nx; i++) {
	for (int j = 0; j < Ny; j++) {
	for (int k = 0; k < Nz; k++) {
		const int id = CellConductor[((int64_t)i * Ny + j) * Nz + k];
		if (id < 0) continue;
		const double coef = CondSigma[id];
		if (coef <= 0) continue;

		double me[8][8];
		element_mass(Xn[i + 1] - Xn[i], Yn[j + 1] - Yn[j], Zn[k + 1] - Zn[k], me);

		for (int l = 0; l < 8; l++) {
			const int li = LI(l), lj = LJ(l), lk = LK(l);
			const int64_t row = node_index(i + li, j + lj, k + lk);
			if (!crs_owns(A, row)) continue;		// 分散した行列 (mpisolve = 2) : 他の rank の行
			const int64_t rowstart = A->rowptr[row];
			for (int mm = 0; mm < 8; mm++) {
				const int di = LI(mm) - li;
				const int dj = LJ(mm) - lj;
				const int dk = LK(mm) - lk;
				const int64_t p = crs_offset(rowstart, i + li, j + lj, k + lk, di, dj, dk);
				A->val[p] += coef * me[l][mm];
			}
		}
	}
	}
	}
}


// セルごとの係数の質量行列 ∫c N_i N_j dV (熱連成の過渡解析の熱容量 ρc)
void assemble_mass_cell(crs_t *A, const double *coefcell)
{
	crs_zero(A);

	for (int i = 0; i < Nx; i++) {
	for (int j = 0; j < Ny; j++) {
	for (int k = 0; k < Nz; k++) {
		const double coef = coefcell[((int64_t)i * Ny + j) * Nz + k];
		if (coef <= 0) continue;

		double me[8][8];
		element_mass(Xn[i + 1] - Xn[i], Yn[j + 1] - Yn[j], Zn[k + 1] - Zn[k], me);

		for (int l = 0; l < 8; l++) {
			const int li = LI(l), lj = LJ(l), lk = LK(l);
			const int64_t row = node_index(i + li, j + lj, k + lk);
			if (!crs_owns(A, row)) continue;		// 分散した行列 (mpisolve = 2) : 他の rank の行
			const int64_t rowstart = A->rowptr[row];
			for (int mm = 0; mm < 8; mm++) {
				const int di = LI(mm) - li;
				const int dj = LJ(mm) - lj;
				const int dk = LK(mm) - lk;
				const int64_t p = crs_offset(rowstart, i + li, j + lj, k + lk, di, dj, dk);
				A->val[p] += coef * me[l][mm];
			}
		}
	}
	}
	}
}


// B-H 曲線から磁気抵抗率 ν = H(B)/B を求める
//
// 表は (H, B) の狭義単調増加な点列 (原点は含まない)。
//   B < B[0]      : 初期磁気抵抗率 H[0]/B[0] (原点と第 1 点を結ぶ直線)
//   B[i]..B[i+1]  : H を線形補間
//   B > B[last]   : 最終区間の傾きで外挿 (飽和後も解が発散しないようにする)
// 軸 ax の曲線で H(|b|) と dH/dB を求める (|b| < B[0] は初期透磁率の直線)
static void bh_h_dh(const material_t *mt, int ax, double b, double *h, double *dh)
{
	const int n = mt->nbh[ax];
	const double *hh = mt->bh_h[ax];
	const double *bb = mt->bh_b[ax];
	const double nu0 = hh[0] / bb[0];

	if (b <= bb[0]) {
		*h = nu0 * b;
		*dh = nu0;
		return;
	}

	if (b >= bb[n - 1]) {
		const double slope = (n >= 2)
			? ((hh[n - 1] - hh[n - 2]) / (bb[n - 1] - bb[n - 2])) : nu0;
		*h = hh[n - 1] + (slope * (b - bb[n - 1]));
		*dh = slope;
	}
	else {
		int i = 0;
		while ((i < n - 2) && (b > bb[i + 1])) i++;
		const double slope = (hh[i + 1] - hh[i]) / (bb[i + 1] - bb[i]);
		*h = hh[i] + (slope * (b - bb[i]));
		*dh = slope;
	}
}


// bh_h_dh の逆関数 B(|h|)。H(B) が区分線形 (節点は表の点、最初は原点からの直線、
// 最後は最終区間の傾きで外挿) なので、その逆も同じ節点で区分線形になる
double bh_b_of_h(const material_t *mt, int ax, double h)
{
	const int n = mt->nbh[ax];
	const double *hh = mt->bh_h[ax];
	const double *bb = mt->bh_b[ax];

	if (h <= hh[0]) return h * bb[0] / hh[0];
	if (h >= hh[n - 1]) {
		const double slope = (n >= 2)
			? ((hh[n - 1] - hh[n - 2]) / (bb[n - 1] - bb[n - 2])) : (hh[0] / bb[0]);
		return bb[n - 1] + ((h - hh[n - 1]) / slope);
	}
	int i = 0;
	while ((i < n - 2) && (h > hh[i + 1])) i++;
	return bb[i] + ((h - hh[i]) * (bb[i + 1] - bb[i]) / (hh[i + 1] - hh[i]));
}


double bh_nu(const material_t *mt, int ax, double b)
{
	if (b <= 0) return mt->bh_h[ax][0] / mt->bh_b[ax][0];

	double h, dh;
	bh_h_dh(mt, ax, b, &h, &dh);

	return h / b;
}


// B-H 曲線から ν と dν/d(B^2) を求める (Newton 反復用)
//   ν(B) = H(B)/B,  dν/dB = (H'B - H)/B^2,  dν/d(B^2) = dν/dB / (2B)
// H(B) が単調増加なら ν + 2 B^2 dν/d(B^2) = dH/dB > 0 なのでヤコビアンは正定値になる
static void bh_nu_dnu(const material_t *mt, double b, double *nu, double *dnudb2)
{
	if (mt->nbh[0] <= 0) {
		// 線形材料の μr にも温度依存を掛ける (material_coef_pub と同じ倍率)
		*nu = 1 / (MU0 * mt->mu6[0] * temp_factor(mt, 1));
		*dnudb2 = 0;
		return;
	}
	if (b <= mt->bh_b[0][0]) {
		*nu = mt->bh_h[0][0] / mt->bh_b[0][0];	// 初期透磁率
		*dnudb2 = 0;
		return;
	}

	double h, hp;
	bh_h_dh(mt, 0, b, &h, &hp);

	*nu = h / b;
	*dnudb2 = ((hp * b) - h) / (2 * b * b * b);
}


// 公開版 (3 次元渦電流の過渡解析が四面体の B に使う)
void bh_nu_dnu_pub(const material_t *mt, double b, double *nu, double *dnudb2)
{
	bh_nu_dnu(mt, b, nu, dnudb2);
}

// 公開版 (断面 2 次元の三角形の格子の Newton が軸ごとの曲線に使う)
void bh_h_dh_pub(const material_t *mt, int ax, double b, double *h, double *dh)
{
	bh_h_dh(mt, ax, b, h, dh);
}


// ---- Jiles-Atherton ヒステリシス ----

// Langevin 関数 M_an = Ms(coth(x) - 1/x) と dM_an/dHe
static void ja_langevin(const ja_t *ja, double he, double *man, double *dman)
{
	const double x = he / ja->a;
	const double ax = fabs(x);

	if (ax < 1e-4) {
		// 級数展開 (x -> 0 で coth(x)-1/x -> x/3)
		*man = ja->ms * x / 3;
		*dman = ja->ms / (3 * ja->a);
	}
	else {
		const double sh = sinh(x);
		*man = ja->ms * ((cosh(x) / sh) - (1 / x));
		*dman = (ja->ms / ja->a) * ((1 / (x * x)) - (1 / (sh * sh)));
	}
}


// dM/dH (Jiles-Atherton)。delta = dH の符号
static double ja_dmdh(const ja_t *ja, double h, double m, double delta)
{
	const double he = h + (ja->alpha * m);
	double man, dman;
	ja_langevin(ja, he, &man, &dman);

	const double dm = man - m;
	// δ_M : 非物理な負の帯磁率を防ぐ
	double irr = 0;
	if ((dm * delta) > 0) {
		const double den = (delta * ja->k) - (ja->alpha * dm);
		if (fabs(den) > 0) irr = (1 - ja->c) * dm / den;
	}
	const double num = irr + (ja->c * dman);
	const double den2 = 1 - (ja->alpha * ja->c * dman);

	double chi = ((fabs(den2) > 0) ? (num / den2) : num);
	if (chi < 0) chi = 0;			// 単調性を保つ

	return chi;
}


// 状態 (b0, h0, m0) から磁束密度 b まで積分し、H, M, 微分磁気抵抗率 dH/dB を返す
//
// FEM の未知数は A なので B が駆動変数になる (逆 J-A)。
//   dB/dH = μ0 (1 + dM/dH)
//   -> dH/dB = 1/(μ0(1+χ)),  dM/dB = χ/(μ0(1+χ)),  χ = dM/dH
//
// H は B/μ0 - M で代数的に求めると桁落ちする (高透磁率では B/μ0 と M が
// ほぼ相殺するので、M の相対誤差 1e-3 が H の 100% 誤差になる)。
// そのため H も B に沿って積分する。積分は Heun 法 (2 次)。
void ja_eval(const ja_t *ja, double b0, double h0, double m0, double b,
	double *hout, double *mout, double *dhdb)
{
	const int nsub = ((JaSub > 0) ? JaSub : 20);
	const double db = (b - b0) / nsub;
	const double delta = ((db >= 0) ? +1.0 : -1.0);

	double h = h0, m = m0;
	for (int s = 0; s < nsub; s++) {
		const double chi1 = ja_dmdh(ja, h, m, delta);
		const double f1h = 1 / (MU0 * (1 + chi1));
		const double f1m = chi1 / (MU0 * (1 + chi1));
		const double chi2 = ja_dmdh(ja, h + (f1h * db), m + (f1m * db), delta);
		const double f2h = 1 / (MU0 * (1 + chi2));
		const double f2m = chi2 / (MU0 * (1 + chi2));
		h += (f1h + f2h) / 2 * db;
		m += (f1m + f2m) / 2 * db;
	}

	double chi = ja_dmdh(ja, h, m, delta);
	// H(B) は確定済みの状態で枝が切り替わる (B が b0 の両側で傾きが違う折れ点を持つ)。
	// 過渡解析で B がほとんど動かない点があると反復がその折れ点をまたいで往復し、
	// 小さい方の傾きの行列では固定点が不安定になる (実測)。JaSafeSlope のときは両方の
	// 枝の大きい方の dH/dB (χ の小さい方) を返す。行列だけに使うので解は変わらない
	if (JaSafeSlope) {
		const double chi2 = ja_dmdh(ja, h, m, -delta);
		if (chi2 < chi) chi = chi2;
	}

	*mout = m;
	*hout = h;
	*dhdb = 1 / (MU0 * (1 + chi));
}


// ---- ベクトル Jiles-Atherton (等方性、断面 2 次元の面内 2 成分) ----
//
// スカラーの式をそのまま向きを持つ形にする (場が 1 方向ならスカラーの式に厳密に退化する):
//   M_an = Ms L(|He|/a) He/|He|、He = H + αM
//   χ_an = dM_an/dHe = (Ms/a) L'(x) ê êᵀ + (Ms L(x)/|He|)(I − ê êᵀ)     (対称)
//   f = M_an − M、不可逆成分は f の向きに : (1−c) |f|/(k − α|f|) f̂ f̂ᵀ  (f・dB > 0 のときだけ = δ_M)
//   χ = dM/dH = (I − α c χ_an)⁻¹ [ (1−c) |f|/(k − α|f|) f̂ f̂ᵀ + c χ_an ]
// 1 次元では (1−c)|f|/(k − α|f|) = (1−c) f/(δk − αf) (δ = sign f) なので ja_dmdh と同じ式になる。
// B が駆動変数なので dH = (μ0 (I + χ))⁻¹ dB、dM = χ dH を Heun 法で積分する (ja_eval と同じ)。


// ---- 異方性のベクトル J-A (jaaxis) ----
//
// 成分 i (全体の軸 ax[i]) が自分の (Ms_i, a_i, α_i, k_i, c_i) を持つ。等方性の式を軸ごとに書いたもので、
//   He = H + α M (α は対角)、M_an,i = Ms_i L(|He| / a_i) e_i   (e = He / |He|)
//   χ_an,ij = dM_an,i/dHe_j = M'_i e_i e_j + (M_an,i/|He|) (δ_ij − e_i e_j)
//   f = M_an − M、不可逆 : |f| / (k_f − α_f |f|) f̂ f̂ᵀ、k_f = Σ k_i f̂_i²、α_f = Σ α_i f̂_i² (f の向きへの射影)
//   χ = (I − C χ_an α)⁻¹ [ (I − C) 不可逆 + C χ_an ]   (C = diag c)
// **場が 1 つの軸に沿えばその軸のパラメータのスカラーの J-A に厳密に退化し、全軸が同じなら等方性の式と
// 一致する** (どちらも検証に使う)。等方性の材料は従来の関数をそのまま通るので答えは変わらない。
static void ja_langevin_ma(double ms, double a, double he, double *man, double *dman)
{
	const double x = he / a;
	if (fabs(x) < 1e-4) {
		*man = ms * x / 3;
		*dman = ms / (3 * a);
	}
	else {
		const double sh = sinh(x);
		*man = ms * ((cosh(x) / sh) - (1 / x));
		*dman = (ms / a) * ((1 / (x * x)) - (1 / (sh * sh)));
	}
}


static int inv3j(const double a[3][3], double r[3][3]);

// n = 2 / 3 成分。成分 i のパラメータは ja->axp[ax[i]]
static void ja_chi_an(const ja_t *ja, int n, const int ax[3], const double *h, const double *m, const double *db,
	double chi[3][3])
{
	double ms[3], aa[3], al[3], kk[3], cc[3], he[3], e[3] = {1, 0, 0};
	for (int i = 0; i < n; i++) {
		const double *p = ja->axp[ax[i]];
		ms[i] = p[0]; aa[i] = p[1]; al[i] = p[2]; kk[i] = p[3]; cc[i] = p[4];
		he[i] = h[i] + (al[i] * m[i]);
	}
	double hn2 = 0;
	for (int i = 0; i < n; i++) hn2 += he[i] * he[i];
	const double hn = sqrt(hn2);
	if (hn > 0) {
		for (int i = 0; i < n; i++) e[i] = he[i] / hn;
	}
	double sec[3], dm[3], can[3][3], f[3];
	for (int i = 0; i < n; i++) {
		double man;
		ja_langevin_ma(ms[i], aa[i], hn, &man, &dm[i]);
		sec[i] = ((hn > (1e-4 * aa[i])) ? (man / hn) : dm[i]);
		f[i] = (sec[i] * he[i]) - m[i];
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) can[i][j] = (dm[i] * e[i] * e[j]) + (sec[i] * (((i == j) ? 1 : 0) - (e[i] * e[j])));
	}
	double fn2 = 0, fd = 0;
	for (int i = 0; i < n; i++) {
		fn2 += f[i] * f[i];
		fd += f[i] * db[i];
	}
	const double fn = sqrt(fn2);
	double irr = 0, fh[3] = {0, 0, 0};
	if ((fn > 0) && (fd > 0)) {
		double kf = 0, af = 0, s1 = 0, cf = 0;
		for (int i = 0; i < n; i++) {
			fh[i] = f[i] / fn;
			kf += kk[i] * fh[i] * fh[i];
			af += al[i] * fh[i] * fh[i];
			s1 += (1 - cc[i]) * fh[i] * fh[i];
		}
		const double den = kf - (af * fn);
		if (fabs(den) > 0) irr = fn / den;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) cf += fh[i] * cc[i] * can[i][j] * fh[j];
		}
		// 単調性 : f の向きの帯磁率を負にしない (等方性の式と同じ)
		const double along = (irr * s1) + cf;
		if ((along < 0) && (s1 > 0)) irr -= along / s1;
	}
	double x[3][3], a[3][3], ai[3][3];
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			x[i][j] = a[i][j] = ((i == j) ? 1 : 0);
			chi[i][j] = 0;
		}
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			x[i][j] = ((1 - cc[i]) * irr * fh[i] * fh[j]) + (cc[i] * can[i][j]);
			a[i][j] = ((i == j) ? 1 : 0) - (cc[i] * can[i][j] * al[j]);
		}
	}
	if (inv3j(a, ai)) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) chi[i][j] = x[i][j];
		}
		return;
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			double v = 0;
			for (int l = 0; l < n; l++) v += ai[i][l] * x[l][j];
			chi[i][j] = v;
		}
	}
}

static int inv2(const double a[2][2], double r[2][2])
{
	const double d = (a[0][0] * a[1][1]) - (a[0][1] * a[1][0]);
	if (d == 0) return 1;
	r[0][0] = a[1][1] / d;
	r[0][1] = -a[0][1] / d;
	r[1][0] = -a[1][0] / d;
	r[1][1] = a[0][0] / d;
	return 0;
}


// 帯磁率テンソル χ (2x2)。db は駆動の向き (δ_M の判定に使う)
static void ja_chi_vec(const ja_t *ja, const double h[2], const double m[2], const double db[2],
	double chi[2][2])
{
	if (ja->aniso) {
		// 面内の成分 (B_p, B_q) は伝送線路軸 t の次の 2 軸 (assemble_ja と同じ)
		const int it = ((Tline == 'X') ? 0 : (Tline == 'Y') ? 1 : 2);
		const int ax[3] = {(it + 1) % 3, (it + 2) % 3, it};
		double c3[3][3];
		ja_chi_an(ja, 2, ax, h, m, db, c3);
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) chi[i][j] = c3[i][j];
		}
		return;
	}
	const double he[2] = {h[0] + (ja->alpha * m[0]), h[1] + (ja->alpha * m[1])};
	const double hn = sqrt((he[0] * he[0]) + (he[1] * he[1]));
	double man, dman;
	ja_langevin(ja, hn, &man, &dman);
	double e[2] = {1, 0};
	if (hn > 0) { e[0] = he[0] / hn; e[1] = he[1] / hn; }
	// 横方向の χ_an = Ms L(x)/|He| (x -> 0 で Ms/(3a) = 縦方向に一致)
	const double sec = ((hn > (1e-4 * ja->a)) ? (man / hn) : dman);
	double can[2][2];
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) can[i][j] = (dman * e[i] * e[j]) + (sec * (((i == j) ? 1 : 0) - (e[i] * e[j])));
	}
	const double f[2] = {(sec * he[0]) - m[0], (sec * he[1]) - m[1]};
	const double fn = sqrt((f[0] * f[0]) + (f[1] * f[1]));
	double irr = 0, fh[2] = {0, 0};
	if ((fn > 0) && (((f[0] * db[0]) + (f[1] * db[1])) > 0)) {
		fh[0] = f[0] / fn;
		fh[1] = f[1] / fn;
		const double den = ja->k - (ja->alpha * fn);
		if (fabs(den) > 0) irr = (1 - ja->c) * fn / den;
		// 単調性 (ja_dmdh の χ < 0 -> 0 と同じ) : f の向きの帯磁率を負にしない
		const double along = irr + (ja->c * ((fh[0] * ((can[0][0] * fh[0]) + (can[0][1] * fh[1])))
		                                   + (fh[1] * ((can[1][0] * fh[0]) + (can[1][1] * fh[1])))));
		if (along < 0) irr -= along;
	}
	double x[2][2], a[2][2], ai[2][2];
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			x[i][j] = (irr * fh[i] * fh[j]) + (ja->c * can[i][j]);
			a[i][j] = ((i == j) ? 1 : 0) - (ja->alpha * ja->c * can[i][j]);
		}
	}
	if (inv2(a, ai)) {
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) chi[i][j] = x[i][j];
		}
		return;
	}
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) chi[i][j] = (ai[i][0] * x[0][j]) + (ai[i][1] * x[1][j]);
	}
}


// 対角でない成分の添字 (xx yy zz xy yz zx の並び)
static int ja_off(int a, int b)
{
	const int lo = ((a < b) ? a : b), hi = ((a < b) ? b : a);
	return ((lo == 0) && (hi == 1)) ? 3 : (((lo == 1) && (hi == 2)) ? 4 : 5);
}


// 2x2 の対称部分の最大固有値
static double ja_symmax(const double a[2][2])
{
	const double p = (a[0][0] + a[1][1]) / 2, q = (a[0][0] - a[1][1]) / 2, r = (a[0][1] + a[1][0]) / 2;
	return p + sqrt((q * q) + (r * r));
}


// dH/dB = (μ0 (I + χ))⁻¹
static void ja_dhdb_vec(const double chi[2][2], double nd[2][2])
{
	double a[2][2];
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) a[i][j] = MU0 * (((i == j) ? 1 : 0) + chi[i][j]);
	}
	if (inv2(a, nd)) {
		nd[0][0] = nd[1][1] = 1 / MU0;
		nd[0][1] = nd[1][0] = 0;
	}
}


// 状態 (b0, h0, m0) から b まで積分して H, M と、行列に使う正の微分磁気抵抗率 (B の向きの
// dH/dB の成分) を返す。dref は db = 0 のときの駆動の向き (スカラーの δ = +1 に当たる)
void ja_eval_vec(const ja_t *ja, const double b0[2], const double h0[2], const double m0[2],
	const double b[2], const double dref[2], double hout[2], double mout[2], double ns[2][2])
{
	const int nsub = ((JaSub > 0) ? JaSub : 20);
	double db[2] = {(b[0] - b0[0]) / nsub, (b[1] - b0[1]) / nsub};
	double dd[2] = {db[0], db[1]};
	if ((dd[0] == 0) && (dd[1] == 0)) { dd[0] = dref[0]; dd[1] = dref[1]; }

	double h[2] = {h0[0], h0[1]}, m[2] = {m0[0], m0[1]};
	for (int s = 0; s < nsub; s++) {
		double c1[2][2], n1[2][2], c2[2][2], n2[2][2];
		ja_chi_vec(ja, h, m, dd, c1);
		ja_dhdb_vec(c1, n1);
		const double dh1[2] = {(n1[0][0] * db[0]) + (n1[0][1] * db[1]), (n1[1][0] * db[0]) + (n1[1][1] * db[1])};
		const double dm1[2] = {(c1[0][0] * dh1[0]) + (c1[0][1] * dh1[1]), (c1[1][0] * dh1[0]) + (c1[1][1] * dh1[1])};
		const double hp[2] = {h[0] + dh1[0], h[1] + dh1[1]};
		const double mp[2] = {m[0] + dm1[0], m[1] + dm1[1]};
		ja_chi_vec(ja, hp, mp, dd, c2);
		ja_dhdb_vec(c2, n2);
		const double dh2[2] = {(n2[0][0] * db[0]) + (n2[0][1] * db[1]), (n2[1][0] * db[0]) + (n2[1][1] * db[1])};
		const double dm2[2] = {(c2[0][0] * dh2[0]) + (c2[0][1] * dh2[1]), (c2[1][0] * dh2[0]) + (c2[1][1] * dh2[1])};
		for (int i = 0; i < 2; i++) {
			h[i] += (dh1[i] + dh2[i]) / 2;
			m[i] += (dm1[i] + dm2[i]) / 2;
		}
	}

	// 行列の傾き : dH/dB の対称部分 (2x2 のテンソル)。B の向きの成分 (スカラーの dH/dB) だけを
	// 等方的に使うと、向きが変わる点で横方向の剛性 (|H|/|B| に近い) を過小評価して固定点反復が
	// 振動・発散した。最大固有値を等方的に使うと 1 方向の板で反復が 6 回 -> 145 回になった (実測)。
	// テンソルなら場が 1 方向のとき縦はスカラーの dH/dB そのもの、横は横方向の剛性になる
	double cf[2][2], nf[2][2];
	ja_chi_vec(ja, h, m, dd, cf);
	ja_dhdb_vec(cf, nf);
	if (JaSafeSlope) {
		// 両方の枝のうち剛性の大きい方 (最大固有値で比べる)
		double c2[2][2], n2[2][2];
		const double dr[2] = {-dd[0], -dd[1]};
		ja_chi_vec(ja, h, m, dr, c2);
		ja_dhdb_vec(c2, n2);
		if (ja_symmax(n2) > ja_symmax(nf)) {
			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < 2; j++) nf[i][j] = n2[i][j];
			}
		}
	}
	ns[0][0] = nf[0][0];
	ns[1][1] = nf[1][1];
	ns[0][1] = ns[1][0] = (nf[0][1] + nf[1][0]) / 2;
	// 正定値でなければ (非対称の帯磁率で起こり得る) 最大固有値の等方テンソルにする
	const double tr = ns[0][0] + ns[1][1], dt = (ns[0][0] * ns[1][1]) - (ns[0][1] * ns[0][1]);
	if (!((tr > 0) && (dt > 0))) {
		double t = ja_symmax(nf);
		if (!(t > 0)) t = 1 / MU0;
		ns[0][0] = ns[1][1] = t;
		ns[0][1] = ns[1][0] = 0;
	}
	hout[0] = h[0]; hout[1] = h[1];
	mout[0] = m[0]; mout[1] = m[1];
}


// ---- ベクトル Jiles-Atherton (3 次元、3 成分) ----
//
// 式は面内 2 成分の版 (ja_chi_vec / ja_eval_vec) と同じものを 3 成分で書いたもの
// (3 次元渦電流の過渡解析 analysis = A + transient の四面体が使う)。面内 2 成分の版は
// 既存の答えをビット単位で変えないためにそのまま残す。

static int inv3j(const double a[3][3], double r[3][3])
{
	const double d = (a[0][0] * ((a[1][1] * a[2][2]) - (a[1][2] * a[2][1])))
	               - (a[0][1] * ((a[1][0] * a[2][2]) - (a[1][2] * a[2][0])))
	               + (a[0][2] * ((a[1][0] * a[2][1]) - (a[1][1] * a[2][0])));
	if (d == 0) return 1;
	r[0][0] = ((a[1][1] * a[2][2]) - (a[1][2] * a[2][1])) / d;
	r[0][1] = ((a[0][2] * a[2][1]) - (a[0][1] * a[2][2])) / d;
	r[0][2] = ((a[0][1] * a[1][2]) - (a[0][2] * a[1][1])) / d;
	r[1][0] = ((a[1][2] * a[2][0]) - (a[1][0] * a[2][2])) / d;
	r[1][1] = ((a[0][0] * a[2][2]) - (a[0][2] * a[2][0])) / d;
	r[1][2] = ((a[0][2] * a[1][0]) - (a[0][0] * a[1][2])) / d;
	r[2][0] = ((a[1][0] * a[2][1]) - (a[1][1] * a[2][0])) / d;
	r[2][1] = ((a[0][1] * a[2][0]) - (a[0][0] * a[2][1])) / d;
	r[2][2] = ((a[0][0] * a[1][1]) - (a[0][1] * a[1][0])) / d;
	return 0;
}


// 帯磁率テンソル χ (3x3)。db は駆動の向き (δ_M の判定に使う)
static void ja_chi_vec3(const ja_t *ja, const double h[3], const double m[3], const double db[3],
	double chi[3][3])
{
	if (ja->aniso) {
		const int ax[3] = {0, 1, 2};
		ja_chi_an(ja, 3, ax, h, m, db, chi);
		return;
	}
	double he[3];
	for (int i = 0; i < 3; i++) he[i] = h[i] + (ja->alpha * m[i]);
	const double hn = sqrt((he[0] * he[0]) + (he[1] * he[1]) + (he[2] * he[2]));
	double man, dman;
	ja_langevin(ja, hn, &man, &dman);
	double e[3] = {1, 0, 0};
	if (hn > 0) {
		for (int i = 0; i < 3; i++) e[i] = he[i] / hn;
	}
	const double sec = ((hn > (1e-4 * ja->a)) ? (man / hn) : dman);
	double can[3][3];
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) can[i][j] = (dman * e[i] * e[j]) + (sec * (((i == j) ? 1 : 0) - (e[i] * e[j])));
	}
	double f[3];
	for (int i = 0; i < 3; i++) f[i] = (sec * he[i]) - m[i];
	const double fn = sqrt((f[0] * f[0]) + (f[1] * f[1]) + (f[2] * f[2]));
	double irr = 0, fh[3] = {0, 0, 0};
	if ((fn > 0) && (((f[0] * db[0]) + (f[1] * db[1]) + (f[2] * db[2])) > 0)) {
		for (int i = 0; i < 3; i++) fh[i] = f[i] / fn;
		const double den = ja->k - (ja->alpha * fn);
		if (fabs(den) > 0) irr = (1 - ja->c) * fn / den;
		double cf = 0;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) cf += fh[i] * can[i][j] * fh[j];
		}
		const double along = irr + (ja->c * cf);
		if (along < 0) irr -= along;
	}
	double x[3][3], a[3][3], ai[3][3];
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			x[i][j] = (irr * fh[i] * fh[j]) + (ja->c * can[i][j]);
			a[i][j] = ((i == j) ? 1 : 0) - (ja->alpha * ja->c * can[i][j]);
		}
	}
	if (inv3j(a, ai)) {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) chi[i][j] = x[i][j];
		}
		return;
	}
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) chi[i][j] = (ai[i][0] * x[0][j]) + (ai[i][1] * x[1][j]) + (ai[i][2] * x[2][j]);
	}
}


// dH/dB = (μ0 (I + χ))⁻¹
static void ja_dhdb_vec3(const double chi[3][3], double nd[3][3])
{
	double a[3][3];
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) a[i][j] = MU0 * (((i == j) ? 1 : 0) + chi[i][j]);
	}
	if (inv3j(a, nd)) {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) nd[i][j] = ((i == j) ? (1 / MU0) : 0);
		}
	}
}


// 3x3 の対称部分の最大固有値 (Jacobi 法。対称 3x3 なので数回で収束する)
static double ja_symmax3(const double a[3][3])
{
	double s[3][3];
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) s[i][j] = (a[i][j] + a[j][i]) / 2;
	}
	for (int sweep = 0; sweep < 50; sweep++) {
		const double off = (s[0][1] * s[0][1]) + (s[0][2] * s[0][2]) + (s[1][2] * s[1][2]);
		const double dg = (s[0][0] * s[0][0]) + (s[1][1] * s[1][1]) + (s[2][2] * s[2][2]);
		if (!(off > (1e-30 * dg))) break;
		for (int p = 0; p < 2; p++) {
			for (int q = p + 1; q < 3; q++) {
				if (s[p][q] == 0) continue;
				const double th = (s[q][q] - s[p][p]) / (2 * s[p][q]);
				const double t = ((th >= 0) ? 1.0 : -1.0) / (fabs(th) + sqrt((th * th) + 1));
				const double c = 1 / sqrt((t * t) + 1), sn = t * c;
				for (int k = 0; k < 3; k++) {
					const double skp = s[k][p], skq = s[k][q];
					s[k][p] = (c * skp) - (sn * skq);
					s[k][q] = (sn * skp) + (c * skq);
				}
				for (int k = 0; k < 3; k++) {
					const double spk = s[p][k], sqk = s[q][k];
					s[p][k] = (c * spk) - (sn * sqk);
					s[q][k] = (sn * spk) + (c * sqk);
				}
			}
		}
	}
	double mx = s[0][0];
	if (s[1][1] > mx) mx = s[1][1];
	if (s[2][2] > mx) mx = s[2][2];
	return mx;
}


// 状態 (b0, h0, m0) から b まで積分して H, M と、行列に使う微分磁気抵抗率のテンソル
// ns (xx yy zz xy yz zx の並び、dH/dB の対称部分、正定値) を返す。
// dref は db = 0 のときの駆動の向き (ja_eval_vec と同じ)
void ja_eval_vec3(const ja_t *ja, const double b0[3], const double h0[3], const double m0[3],
	const double b[3], const double dref[3], double hout[3], double mout[3], double ns[6])
{
	const int nsub = ((JaSub > 0) ? JaSub : 20);
	double db[3], dd[3];
	for (int i = 0; i < 3; i++) dd[i] = db[i] = (b[i] - b0[i]) / nsub;
	if ((dd[0] == 0) && (dd[1] == 0) && (dd[2] == 0)) {
		for (int i = 0; i < 3; i++) dd[i] = dref[i];
	}

	double h[3], m[3];
	for (int i = 0; i < 3; i++) { h[i] = h0[i]; m[i] = m0[i]; }
	for (int s = 0; s < nsub; s++) {
		double c1[3][3], n1[3][3], c2[3][3], n2[3][3], dh1[3], dm1[3], hp[3], mp[3], dh2[3], dm2[3];
		ja_chi_vec3(ja, h, m, dd, c1);
		ja_dhdb_vec3(c1, n1);
		for (int i = 0; i < 3; i++) dh1[i] = (n1[i][0] * db[0]) + (n1[i][1] * db[1]) + (n1[i][2] * db[2]);
		for (int i = 0; i < 3; i++) dm1[i] = (c1[i][0] * dh1[0]) + (c1[i][1] * dh1[1]) + (c1[i][2] * dh1[2]);
		for (int i = 0; i < 3; i++) { hp[i] = h[i] + dh1[i]; mp[i] = m[i] + dm1[i]; }
		ja_chi_vec3(ja, hp, mp, dd, c2);
		ja_dhdb_vec3(c2, n2);
		for (int i = 0; i < 3; i++) dh2[i] = (n2[i][0] * db[0]) + (n2[i][1] * db[1]) + (n2[i][2] * db[2]);
		for (int i = 0; i < 3; i++) dm2[i] = (c2[i][0] * dh2[0]) + (c2[i][1] * dh2[1]) + (c2[i][2] * dh2[2]);
		for (int i = 0; i < 3; i++) {
			h[i] += (dh1[i] + dh2[i]) / 2;
			m[i] += (dm1[i] + dm2[i]) / 2;
		}
	}

	double cf[3][3], nf[3][3];
	ja_chi_vec3(ja, h, m, dd, cf);
	ja_dhdb_vec3(cf, nf);
	if (JaSafeSlope) {
		double c2[3][3], n2[3][3];
		const double dr[3] = {-dd[0], -dd[1], -dd[2]};
		ja_chi_vec3(ja, h, m, dr, c2);
		ja_dhdb_vec3(c2, n2);
		if (ja_symmax3(n2) > ja_symmax3(nf)) {
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) nf[i][j] = n2[i][j];
			}
		}
	}
	double sy[3][3];
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) sy[i][j] = (nf[i][j] + nf[j][i]) / 2;
	}
	// 正定値か (Sylvester の主小行列式)。でなければ最大固有値の等方テンソル
	const double d1 = sy[0][0], d2 = (sy[0][0] * sy[1][1]) - (sy[0][1] * sy[1][0]);
	const double d3 = (sy[0][0] * ((sy[1][1] * sy[2][2]) - (sy[1][2] * sy[2][1])))
	                - (sy[0][1] * ((sy[1][0] * sy[2][2]) - (sy[1][2] * sy[2][0])))
	                + (sy[0][2] * ((sy[1][0] * sy[2][1]) - (sy[1][1] * sy[2][0])));
	if ((d1 > 0) && (d2 > 0) && (d3 > 0)) {
		ns[0] = sy[0][0];
		ns[1] = sy[1][1];
		ns[2] = sy[2][2];
		ns[3] = sy[0][1];
		ns[4] = sy[1][2];
		ns[5] = sy[2][0];
	}
	else {
		double t = ja_symmax3(nf);
		if (!(t > 0)) t = 1 / MU0;
		ns[0] = ns[1] = ns[2] = t;
		ns[3] = ns[4] = ns[5] = 0;
	}
	for (int i = 0; i < 3; i++) { hout[i] = h[i]; mout[i] = m[i]; }
}



// ベクトル J-A (3 成分) の自己検証 : 向きを固定した B の経路 (u = (1, 2, 2)/3 に沿って往復) では
// スカラーの ja_eval と同じ Heun の刻みになるので、H が丸めを除いて一致し、横成分は 0 になる。
// 戻り値は max |H_vec − H_s u| / max |H_s| (式の誤りは刻みの次数の誤りまで含めて 1e-6 以上になる)
static double ja_vec3_selftest_dir(const ja_t *ja, const ja_t *js, const double u[3]);

// 3 成分のベクトル J-A を向きを固定した経路で動かし、スカラーの J-A と比べる。異方性 (jaaxis) なら
// 各軸に沿った経路をその軸のパラメータのスカラーの J-A と比べる (成分ごとのパラメータの取り違えの検査)
double ja_vec3_selftest(const ja_t *ja)
{
	if (!ja->aniso) {
		const double u[3] = {1.0 / 3, 2.0 / 3, 2.0 / 3};
		return ja_vec3_selftest_dir(ja, ja, u);
	}
	double emax = 0;
	for (int d = 0; d < 3; d++) {
		ja_t js = *ja;
		js.aniso = 0;
		js.ms = ja->axp[d][0]; js.a = ja->axp[d][1]; js.alpha = ja->axp[d][2]; js.k = ja->axp[d][3]; js.c = ja->axp[d][4];
		const double u[3] = {((d == 0) ? 1.0 : 0.0), ((d == 1) ? 1.0 : 0.0), ((d == 2) ? 1.0 : 0.0)};
		const double e = ja_vec3_selftest_dir(ja, &js, u);
		if (!(e <= emax)) emax = e;
	}
	return emax;
}


static double ja_vec3_selftest_dir(const ja_t *ja, const ja_t *js, const double u[3])
{
	static const double bp[8] = {0.3, 0.9, 1.4, 0.6, -0.4, -1.3, -0.2, 0.8};
	const int safe = JaSafeSlope;
	JaSafeSlope = 0;
	double bs = 0, hs = 0, ms = 0, bv[3] = {0, 0, 0}, hv[3] = {0, 0, 0}, mv[3] = {0, 0, 0};
	double emax = 0, hmax = 0;
	for (int q = 0; q < 8; q++) {
		double h, m, dhdb, bn[3], hn[3], mn[3], ns[6];
		ja_eval(js, bs, hs, ms, bp[q], &h, &m, &dhdb);
		for (int i = 0; i < 3; i++) bn[i] = bp[q] * u[i];
		ja_eval_vec3(ja, bv, hv, mv, bn, u, hn, mn, ns);
		for (int i = 0; i < 3; i++) {
			const double d = fabs(hn[i] - (h * u[i]));
			if (d > emax) emax = d;
		}
		if (fabs(h) > hmax) hmax = fabs(h);
		bs = bp[q]; hs = h; ms = m;
		for (int i = 0; i < 3; i++) { bv[i] = bn[i]; hv[i] = hn[i]; mv[i] = mn[i]; }
	}
	JaSafeSlope = safe;
	return ((hmax > 0) ? (emax / hmax) : emax);
}

// 線形材料 (J-A でないセル) の磁気エネルギー ½ ∫ ν |∇A|² dV (2x2x2 Gauss、assemble_ja と同じ ν)。
// 線形のセルでは台形則の ∫ H・ΔB が厳密に telescoping するので、掃引の仕事の収支には
// 最後の状態のこの値を足せばよい
static double ja_linear_energy_tri(const double *az);
static void assemble_ja_tri(const double *az, double *res, double *nucell);

double ja_linear_energy(const double *az)
{
	if (MeshMode && (MeshDim == 2)) return ja_linear_energy_tri(az);
	const double gp = 1 / sqrt(3.0);
	double wsum = 0;
	for (int i = 0; i < Nx; i++) {
	for (int j = 0; j < Ny; j++) {
	for (int k = 0; k < Nz; k++) {
		const int64_t c = ((int64_t)i * Ny + j) * Nz + k;
		const material_t *mt = &Material[CellMaterial[c]];
		if (mt->ja.on) continue;
		const double dx = Xn[i + 1] - Xn[i], dy = Yn[j + 1] - Yn[j], dz = Zn[k + 1] - Zn[k];
		const double detj = (dx * dy * dz) / 8;
		const double rx = 2 / dx, ry = 2 / dy, rz = 2 / dz;
		const double nu = 1 / (MU0 * mt->mu6[0] * temp_factor(mt, 1));
		double a[8];
		for (int l = 0; l < 8; l++) a[l] = az[node_index(i + LI(l), j + LJ(l), k + LK(l))];
		for (int g = 0; g < 8; g++) {
			const double xi = (LI(g) ? +gp : -gp), eta = (LJ(g) ? +gp : -gp), zta = (LK(g) ? +gp : -gp);
			double gr[3] = {0, 0, 0};
			for (int l = 0; l < 8; l++) {
				const double si = (LI(l) ? +1.0 : -1.0), sj = (LJ(l) ? +1.0 : -1.0), sk = (LK(l) ? +1.0 : -1.0);
				gr[0] += a[l] * si * (1 + (sj * eta)) * (1 + (sk * zta)) / 8 * rx;
				gr[1] += a[l] * sj * (1 + (sk * zta)) * (1 + (si * xi )) / 8 * ry;
				gr[2] += a[l] * sk * (1 + (si * xi )) * (1 + (sj * eta)) / 8 * rz;
			}
			wsum += 0.5 * nu * ((gr[0] * gr[0]) + (gr[1] * gr[1]) + (gr[2] * gr[2])) * detj;
		}
	}
	}
	}
	return wsum;
}


// 反復中の J-A の状態 (Jn) を確定させる (静磁場の電流掃引と過渡解析の各刻みの終わり)。
// 方向の基準は最初に磁化したときだけ決める
void ja_commit(int64_t ncell)
{
	if (JaB == NULL) return;
	memcpy(JaB, JaBn, (size_t)ncell * 8 * sizeof(double));
	memcpy(JaH, JaHn, (size_t)ncell * 8 * sizeof(double));
	memcpy(JaM, JaMn, (size_t)ncell * 8 * sizeof(double));
	if (JaVB != NULL) {
		memcpy(JaVB, JaVBn, (size_t)ncell * 16 * sizeof(double));
		memcpy(JaVH, JaVHn, (size_t)ncell * 16 * sizeof(double));
		memcpy(JaVM, JaVMn, (size_t)ncell * 16 * sizeof(double));
	}
	for (int64_t id = 0; id < ncell * 8; id++) {
		double *d = &JaD[id * 3];
		if ((d[0] == 0) && (d[1] == 0) && (d[2] == 0)) {
			d[0] = JaDn[(id * 3) + 0];
			d[1] = JaDn[(id * 3) + 1];
			d[2] = JaDn[(id * 3) + 2];
		}
	}
}


// ヒステリシス材料の内力ベクトル res = ∫ ν_ja(|B|) ∇A・∇N_i dV と、
// 固定点行列に使う微分磁気抵抗率 nucell (セル毎、正定値)。
// ヒステリシスでは ν = H/B が負になり得る (第 2 象限) ため、行列には
// 常に正の微分磁気抵抗率 dH/dB を使い、残差側だけで正確な H を扱う。
// 反復中の状態は JaBn / JaMn に書き、収束後に JaB / JaM へ確定させる。
void assemble_ja(const double *az, double *res, double *nucell)
{
	if (MeshMode && (MeshDim == 2)) {
		assemble_ja_tri(az, res, nucell);
		return;
	}
	const double gp = 1 / sqrt(3.0);
	const int64_t nnode = num_node();

	memset(res, 0, (size_t)nnode * sizeof(double));

	for (int i = 0; i < Nx; i++) {
	for (int j = 0; j < Ny; j++) {
	for (int k = 0; k < Nz; k++) {
		const int64_t c = ((int64_t)i * Ny + j) * Nz + k;
		const material_t *mt = &Material[CellMaterial[c]];

		const double dx = Xn[i + 1] - Xn[i];
		const double dy = Yn[j + 1] - Yn[j];
		const double dz = Zn[k + 1] - Zn[k];
		const double detj = (dx * dy * dz) / 8;
		const double rx = 2 / dx, ry = 2 / dy, rz = 2 / dz;

		int64_t nd[8];
		double a[8];
		for (int l = 0; l < 8; l++) {
			nd[l] = node_index(i + LI(l), j + LJ(l), k + LK(l));
			a[l] = az[nd[l]];
		}

		double re[8];
		for (int l = 0; l < 8; l++) re[l] = 0;
		double nusum = 0;
		double ten[6] = {0, 0, 0, 0, 0, 0};

		for (int g = 0; g < 8; g++) {
			const double xi  = (LI(g) ? +gp : -gp);
			const double eta = (LJ(g) ? +gp : -gp);
			const double zta = (LK(g) ? +gp : -gp);

			double dn[8][3];
			for (int l = 0; l < 8; l++) {
				const double si = (LI(l) ? +1.0 : -1.0);
				const double sj = (LJ(l) ? +1.0 : -1.0);
				const double sk = (LK(l) ? +1.0 : -1.0);
				dn[l][0] = si * (1 + (sj * eta)) * (1 + (sk * zta)) / 8 * rx;
				dn[l][1] = sj * (1 + (sk * zta)) * (1 + (si * xi )) / 8 * ry;
				dn[l][2] = sk * (1 + (si * xi )) * (1 + (sj * eta)) / 8 * rz;
			}

			double gr[3] = {0, 0, 0};
			for (int l = 0; l < 8; l++) {
				gr[0] += a[l] * dn[l][0];
				gr[1] += a[l] * dn[l][1];
				gr[2] += a[l] * dn[l][2];
			}
			const double bmag = sqrt((gr[0] * gr[0]) + (gr[1] * gr[1]) + (gr[2] * gr[2]));

			double nu, nud;
			if (mt->ja.on) {
				// ヒステリシスは符号を持つ (残留磁束・保磁力) ので、最初に磁化した
				// 向きを基準ベクトルとして符号付きの B を取る。
				// 場が回転しない問題ではこれが厳密な扱いになる。
				const int64_t id = (c * 8) + g;
				double *dir = &JaD[id * 3];
				const int haved = ((dir[0] != 0) || (dir[1] != 0) || (dir[2] != 0));
				double bs = bmag;
				if (haved) {
					bs = (gr[0] * dir[0]) + (gr[1] * dir[1]) + (gr[2] * dir[2]);
				}
				if (bmag > 1e-12) {
					JaDn[(id * 3) + 0] = gr[0] / bmag;
					JaDn[(id * 3) + 1] = gr[1] / bmag;
					JaDn[(id * 3) + 2] = gr[2] / bmag;
				}

				if (JaVector) {
					// ベクトル J-A : 面内の B = (∂A/∂q, −∂A/∂p) (p, q は伝送線路軸 t の次の 2 軸)
					const int it = ((Tline == 'X') ? 0 : (Tline == 'Y') ? 1 : 2);
					const int ip = (it + 1) % 3, iq = (it + 2) % 3;
					const double bv[2] = {gr[iq], -gr[ip]};
					// 基準の向き (勾配の向き dir を B の向きに回したもの)。スカラーの符号付き B と同じ射影
					const double *dg = (haved ? dir : &JaDn[id * 3]);
					double dref[2] = {dg[iq], -dg[ip]};
					double hv[2], mv[2], ns[2][2];
					ja_eval_vec(&mt->ja, &JaVB[id * 2], &JaVH[id * 2], &JaVM[id * 2], bv, dref, hv, mv, ns);
					// 行列のテンソル (∇A の基底) : c_pp = N_qq、c_qq = N_pp、c_pq = −N_pq (B = (∂A/∂q, −∂A/∂p))
					const double nmax = ja_symmax(ns);
					ten[ip] += ns[1][1];
					ten[iq] += ns[0][0];
					ten[it] += nmax;
					ten[ja_off(ip, iq)] -= ns[0][1];
					const double nus = nmax;
					for (int q = 0; q < 2; q++) {
						JaVBn[(id * 2) + q] = bv[q];
						JaVHn[(id * 2) + q] = hv[q];
						JaVMn[(id * 2) + q] = mv[q];
					}
					JaBn[id] = bs;
					JaHn[id] = (hv[0] * dref[0]) + (hv[1] * dref[1]);
					JaMn[id] = (mv[0] * dref[0]) + (mv[1] * dref[1]);
					// 残差 ∫ H・(∇×N ê_t) = ∫ (H_p ∂N/∂q − H_q ∂N/∂p)
					for (int l = 0; l < 8; l++) {
						re[l] += ((hv[0] * dn[l][iq]) - (hv[1] * dn[l][ip])) * detj;
					}
					nusum += nus;
					continue;
				}
				double h, m, dhdb;
				ja_eval(&mt->ja, JaB[id], JaH[id], JaM[id], bs, &h, &m, &dhdb);
				JaBn[id] = bs;
				JaHn[id] = h;
				JaMn[id] = m;
				// ν = H/B は第 2 象限で負になる (残差側だけで使う)
				nu  = ((fabs(bs) > 1e-12) ? (h / bs) : dhdb);
				nud = dhdb;
			}
			else {
				nu = nud = 1 / (MU0 * mt->mu6[0] * temp_factor(mt, 1));
			}
			nusum += nud;

			for (int l = 0; l < 8; l++) {
				const double u = (dn[l][0] * gr[0]) + (dn[l][1] * gr[1]) + (dn[l][2] * gr[2]);
				re[l] += nu * u * detj;
			}
		}

		nucell[c] = nusum / 8;
		if (JaNu6 != NULL) {
			// ベクトル J-A の行列のテンソル。J-A でないセルは等方の ν
			for (int d = 0; d < 6; d++) JaNu6[(c * 6) + d] = ((mt->ja.on) ? (ten[d] / 8) : ((d < 3) ? (nusum / 8) : 0));
		}
		for (int l = 0; l < 8; l++) {
			res[nd[l]] += re[l];
		}
	}
	}
	}
}


/*
断面 2 次元の三角形の格子 (1 次) のヒステリシス : ∇A が要素内で一定なので、J-A の状態は要素に 1 つ (枠の先頭 id = 8e) で
厳密 (構造格子の 2x2x2 Gauss の代わり)。式は構造格子と同じ (∇A の面内 2 成分 (∂A/∂p, ∂A/∂q)、B = (∂A/∂q, −∂A/∂p))。
JaD の向きは面内 2 成分だけを使う (この経路の中だけで閉じる)
*/
static void assemble_ja_tri(const double *az, double *res, double *nucell)
{
	const int it = ((Tline == 'X') ? 0 : (Tline == 'Y') ? 1 : 2);
	const int ip = (it + 1) % 3, iq = (it + 2) % 3;

	memset(res, 0, (size_t)NNode * sizeof(double));
	for (int e = 0; e < NTri; e++) {
		const int32_t *nd = &Tri[e * 3];
		const material_t *mt = &Material[TriMat[e]];
		double g[3][2], area;
		if (tri_grad(nd, g, &area)) continue;		// 退化した三角形は setup で弾いてある
		const double gr[2] = {(az[nd[0]] * g[0][0]) + (az[nd[1]] * g[1][0]) + (az[nd[2]] * g[2][0]),
		                      (az[nd[0]] * g[0][1]) + (az[nd[1]] * g[1][1]) + (az[nd[2]] * g[2][1])};
		const double bmag = sqrt((gr[0] * gr[0]) + (gr[1] * gr[1]));
		double re[3] = {0, 0, 0};
		double nu, nud;
		double ten[6] = {0, 0, 0, 0, 0, 0};
		if (mt->ja.on) {
			const int64_t id = (int64_t)e * 8;
			double *dir = &JaD[id * 3];
			const int haved = ((dir[0] != 0) || (dir[1] != 0));
			const double bs = (haved ? ((gr[0] * dir[0]) + (gr[1] * dir[1])) : bmag);
			if (bmag > 1e-12) {
				JaDn[(id * 3) + 0] = gr[0] / bmag;
				JaDn[(id * 3) + 1] = gr[1] / bmag;
				JaDn[(id * 3) + 2] = 0;
			}
			if (JaVector) {
				const double bv[2] = {gr[1], -gr[0]};
				const double *dg = (haved ? dir : &JaDn[id * 3]);
				double dref[2] = {dg[1], -dg[0]};
				double hv[2], mv[2], ns[2][2];
				ja_eval_vec(&mt->ja, &JaVB[id * 2], &JaVH[id * 2], &JaVM[id * 2], bv, dref, hv, mv, ns);
				const double nmax = ja_symmax(ns);
				ten[ip] = ns[1][1];
				ten[iq] = ns[0][0];
				ten[it] = nmax;
				ten[ja_off(ip, iq)] = -ns[0][1];
				for (int q = 0; q < 2; q++) {
					JaVBn[(id * 2) + q] = bv[q];
					JaVHn[(id * 2) + q] = hv[q];
					JaVMn[(id * 2) + q] = mv[q];
				}
				JaBn[id] = bs;
				JaHn[id] = (hv[0] * dref[0]) + (hv[1] * dref[1]);
				JaMn[id] = (mv[0] * dref[0]) + (mv[1] * dref[1]);
				for (int l = 0; l < 3; l++) re[l] = ((hv[0] * g[l][1]) - (hv[1] * g[l][0])) * area;
				nucell[e] = nmax;
				if (JaNu6 != NULL) for (int d = 0; d < 6; d++) JaNu6[((int64_t)e * 6) + d] = ten[d];
				for (int l = 0; l < 3; l++) res[nd[l]] += re[l];
				continue;
			}
			double h, m, dhdb;
			ja_eval(&mt->ja, JaB[id], JaH[id], JaM[id], bs, &h, &m, &dhdb);
			JaBn[id] = bs;
			JaHn[id] = h;
			JaMn[id] = m;
			nu  = ((fabs(bs) > 1e-12) ? (h / bs) : dhdb);
			nud = dhdb;
		}
		else {
			nu = nud = 1 / (MU0 * mt->mu6[0] * temp_factor(mt, 1));
		}
		for (int l = 0; l < 3; l++) re[l] = nu * ((g[l][0] * gr[0]) + (g[l][1] * gr[1])) * area;
		nucell[e] = nud;
		if (JaNu6 != NULL) {
			for (int d = 0; d < 6; d++) JaNu6[((int64_t)e * 6) + d] = ((d < 3) ? nud : 0);
		}
		for (int l = 0; l < 3; l++) res[nd[l]] += re[l];
	}
}


// 線形の三角形 (J-A でない要素) の磁気エネルギー ½ ∫ ν |∇A|² dS (∇A は要素内一定で厳密)
static double ja_linear_energy_tri(const double *az)
{
	double wsum = 0;
	for (int e = 0; e < NTri; e++) {
		const material_t *mt = &Material[TriMat[e]];
		if (mt->ja.on) continue;
		const int32_t *nd = &Tri[e * 3];
		double g[3][2], area;
		if (tri_grad(nd, g, &area)) continue;
		double gr[2] = {0, 0};
		for (int l = 0; l < 3; l++) { gr[0] += az[nd[l]] * g[l][0]; gr[1] += az[nd[l]] * g[l][1]; }
		const double nu = 1 / (MU0 * mt->mu6[0] * temp_factor(mt, 1));
		wsum += 0.5 * nu * ((gr[0] * gr[0]) + (gr[1] * gr[1])) * area;
	}
	return wsum;
}


// 非線形静磁場の Newton 反復 : ヤコビアン J と内力ベクトル res = K(ν(|∇A|)) A を作る
//
//   R_i = Σ_g w ν(B) (∇N_i・∇A) detJ
//   J_ij = Σ_g w [ ν (∇N_i・∇N_j) + 2 (dν/dB^2) (∇N_i・∇A)(∇N_j・∇A) ] detJ
//
// ν を Gauss 点毎に評価するので R と J が整合し、2 次収束する。
// 線形材料では dν/dB^2 = 0 なので J = K(ν)、R = K(ν)A に退化する。
void assemble_newton(crs_t *J, const double *az, double *res)
{
	const double gp = 1 / sqrt(3.0);
	const int64_t nnode = num_node();

	crs_zero(J);
	memset(res, 0, (size_t)nnode * sizeof(double));

	for (int i = 0; i < Nx; i++) {
	for (int j = 0; j < Ny; j++) {
	for (int k = 0; k < Nz; k++) {
		const int64_t c = ((int64_t)i * Ny + j) * Nz + k;
		const material_t *mt = &Material[CellMaterial[c]];

		const double dx = Xn[i + 1] - Xn[i];
		const double dy = Yn[j + 1] - Yn[j];
		const double dz = Zn[k + 1] - Zn[k];
		const double detj = (dx * dy * dz) / 8;
		const double rx = 2 / dx, ry = 2 / dy, rz = 2 / dz;

		int64_t nd[8];
		double a[8];
		for (int l = 0; l < 8; l++) {
			nd[l] = node_index(i + LI(l), j + LJ(l), k + LK(l));
			a[l] = az[nd[l]];
		}

		double je[8][8], re[8];
		for (int l = 0; l < 8; l++) {
			re[l] = 0;
			for (int m = 0; m < 8; m++) je[l][m] = 0;
		}

		// 直交異方性 (軸毎の B-H) : 伝送線路軸 t に対し B_p = ∂A/∂q, B_q = -∂A/∂p
		// (p, q は t の次の 2 軸を巡回順に取る)
		const int aniso = (mt->bhaniso && (mt->nbh[0] > 0));
		const int it = ((Tline == 'X') ? 0 : (Tline == 'Y') ? 1 : 2);
		const int ip = (it + 1) % 3;
		const int iq = (it + 2) % 3;

		for (int g = 0; g < 8; g++) {
			const double xi  = (LI(g) ? +gp : -gp);
			const double eta = (LJ(g) ? +gp : -gp);
			const double zta = (LK(g) ? +gp : -gp);

			double dn[8][3];
			for (int l = 0; l < 8; l++) {
				const double si = (LI(l) ? +1.0 : -1.0);
				const double sj = (LJ(l) ? +1.0 : -1.0);
				const double sk = (LK(l) ? +1.0 : -1.0);
				dn[l][0] = si * (1 + (sj * eta)) * (1 + (sk * zta)) / 8 * rx;
				dn[l][1] = sj * (1 + (sk * zta)) * (1 + (si * xi )) / 8 * ry;
				dn[l][2] = sk * (1 + (si * xi )) * (1 + (sj * eta)) / 8 * rz;
			}

			double gr[3] = {0, 0, 0};
			for (int l = 0; l < 8; l++) {
				gr[0] += a[l] * dn[l][0];
				gr[1] += a[l] * dn[l][1];
				gr[2] += a[l] * dn[l][2];
			}

			if (aniso) {
				// エネルギーが軸毎に分離するので rank-1 項は出ない
				//   R_i  = ∫ [ H_p(B_p) ∂N_i/∂q - H_q(B_q) ∂N_i/∂p ] dV
				//   J_ij = ∫ [ (dH_p/dB) ∂N_i/∂q ∂N_j/∂q
				//            + (dH_q/dB) ∂N_i/∂p ∂N_j/∂p ] dV
				const double bp = gr[iq];
				const double bq = -gr[ip];
				double hp, dhp, hq, dhq;
				bh_h_dh(mt, ip, fabs(bp), &hp, &dhp);
				bh_h_dh(mt, iq, fabs(bq), &hq, &dhq);
				if (bp < 0) hp = -hp;
				if (bq < 0) hq = -hq;

				for (int l = 0; l < 8; l++) {
					re[l] += ((hp * dn[l][iq]) - (hq * dn[l][ip])) * detj;
				}
				for (int l = 0; l < 8; l++) {
					for (int m = 0; m < 8; m++) {
						je[l][m] += ((dhp * dn[l][iq] * dn[m][iq])
						          +  (dhq * dn[l][ip] * dn[m][ip])) * detj;
					}
				}
				continue;
			}

			const double bmag = sqrt((gr[0] * gr[0]) + (gr[1] * gr[1]) + (gr[2] * gr[2]));

			double nu, dnudb2;
			bh_nu_dnu(mt, bmag, &nu, &dnudb2);

			double u[8];
			for (int l = 0; l < 8; l++) {
				u[l] = (dn[l][0] * gr[0]) + (dn[l][1] * gr[1]) + (dn[l][2] * gr[2]);
				re[l] += nu * u[l] * detj;
			}
			for (int l = 0; l < 8; l++) {
				for (int m = 0; m < 8; m++) {
					const double kk = (dn[l][0] * dn[m][0])
					                + (dn[l][1] * dn[m][1])
					                + (dn[l][2] * dn[m][2]);
					je[l][m] += ((nu * kk) + (2 * dnudb2 * u[l] * u[m])) * detj;
				}
			}
		}

		for (int l = 0; l < 8; l++) {
			const int li = LI(l), lj = LJ(l), lk = LK(l);
			res[nd[l]] += re[l];
			if (!crs_owns(J, nd[l])) continue;		// 分散した行列 (mpisolve = 2) : 他の rank の行 (残差は全行)
			const int64_t rowstart = J->rowptr[nd[l]];
			for (int m = 0; m < 8; m++) {
				const int64_t p = crs_offset(rowstart, i + li, j + lj, k + lk,
					LI(m) - li, LJ(m) - lj, LK(m) - lk);
				J->val[p] += je[l][m];
			}
		}
	}
	}
	}
}


// 全体行列の作成
//   coefcell != NULL のときはセル毎の係数で上書きする (非線形 B-H 反復で使う)
static void assemble_core(crs_t *A, int mode, const double *coefcell)
{
	crs_zero(A);

	for (int i = 0; i < Nx; i++) {
	for (int j = 0; j < Ny; j++) {
	for (int k = 0; k < Nz; k++) {
		const int64_t cid = ((int64_t)i * Ny + j) * Nz + k;
		const int m = CellMaterial[cid];
		double coef[6];
		if (coefcell != NULL) {
			for (int d = 0; d < 6; d++) coef[d] = 0;
			coef[0] = coef[1] = coef[2] = coefcell[cid];
		}
		else {
			material_coef_pub(m, mode, coef);
		}
		if ((coef[0] <= 0) && (coef[1] <= 0) && (coef[2] <= 0)) continue;

		double ke[8][8];
		element_matrix(Xn[i + 1] - Xn[i], Yn[j + 1] - Yn[j], Zn[k + 1] - Zn[k], coef, ke);

		for (int l = 0; l < 8; l++) {
			const int li = LI(l), lj = LJ(l), lk = LK(l);
			const int64_t row = node_index(i + li, j + lj, k + lk);
			if (!crs_owns(A, row)) continue;		// 分散した行列 (mpisolve = 2) : 他の rank の行
			const int64_t rowstart = A->rowptr[row];
			for (int mm = 0; mm < 8; mm++) {
				const int di = LI(mm) - li;
				const int dj = LJ(mm) - lj;
				const int dk = LK(mm) - lk;
				const int64_t p = crs_offset(rowstart, i + li, j + lj, k + lk, di, dj, dk);
				A->val[p] += ke[l][mm];
			}
		}
	}
	}
	}
}


// セル毎の 6 成分テンソル (xx yy zz xy yz zx) で行列を作る (ベクトル J-A の固定点行列)
void assemble_nu6(crs_t *A, const double *nu6)
{
	crs_zero(A);
	for (int i = 0; i < Nx; i++) {
	for (int j = 0; j < Ny; j++) {
	for (int k = 0; k < Nz; k++) {
		const int64_t cid = ((int64_t)i * Ny + j) * Nz + k;
		double ke[8][8];
		element_matrix(Xn[i + 1] - Xn[i], Yn[j + 1] - Yn[j], Zn[k + 1] - Zn[k], &nu6[cid * 6], ke);
		for (int l = 0; l < 8; l++) {
			const int li = LI(l), lj = LJ(l), lk = LK(l);
			const int64_t row = node_index(i + li, j + lj, k + lk);
			if (!crs_owns(A, row)) continue;		// 分散した行列 (mpisolve = 2) : 他の rank の行
			const int64_t rowstart = A->rowptr[row];
			for (int mm = 0; mm < 8; mm++) {
				const int64_t p = crs_offset(rowstart, i + li, j + lj, k + lk, LI(mm) - li, LJ(mm) - lj, LK(mm) - lk);
				A->val[p] += ke[l][mm];
			}
		}
	}
	}
	}
}


void assemble(crs_t *A, int mode)
{
	if (MeshMode) {
		assemble_elem3d(A, mode);		// 四面体・六面体・角柱・ピラミッド (混在可)
		return;
	}
	assemble_core(A, mode, NULL);
}


// セル毎の等方性係数で行列を作る (熱連成 : σ(T) と熱伝導率 k)。
// 係数が 0 以下のセルは要素ごと飛ばす (その節点は行が空になる)
void assemble_cell(crs_t *A, const double *coefcell)
{
	assemble_core(A, 2, coefcell);
}


/*
Joule 発熱の節点荷重 f_i = Σ_gp w σ_c |∇φ|² N_i (熱連成)。

**剛性と同じ 2x2x2 Gauss・同じ重みで積分する。** そうすると Σ_i N_i = 1 から
	Σ_i f_i = Σ_c Σ_gp w σ_c |∇φ|² = φᵀ K φ
が丸めを除いて厳密に成り立ち、さらに自由節点で (Kφ)_i = 0 なので
	φᵀ K φ = V I = I² R
になる (発熱の総和 = 端子で測った電力)。熱の側でも K_T の行和が 0 なので、
ヒートシンクから出ていく熱の総和が Σ_i f_i に一致する。この 2 つの恒等式で
発熱の積分と熱の組み立てを検査できる (rlc_check.sh)。
要素の中心 1 点で |∇φ|² を取るような手抜きをすると、σ が一様でない
(温度が要素内で変わる) ときに 1 つ目が崩れる。
*/
void joule_load(const double *phi, const double *sigcell, double *f)
{
	const int64_t nn = num_node();
	const double gp = 1 / sqrt(3.0);

	for (int64_t i = 0; i < nn; i++) f[i] = 0;

	for (int i = 0; i < Nx; i++) {
	for (int j = 0; j < Ny; j++) {
	for (int k = 0; k < Nz; k++) {
		const int64_t cid = ((int64_t)i * Ny + j) * Nz + k;
		const double sg = sigcell[cid];
		if (sg <= 0) continue;

		const double dx = Xn[i + 1] - Xn[i];
		const double dy = Yn[j + 1] - Yn[j];
		const double dz = Zn[k + 1] - Zn[k];
		const double detj = (dx * dy * dz) / 8;
		const double rx = 2 / dx, ry = 2 / dy, rz = 2 / dz;

		int64_t nd[8];
		double pv[8];
		for (int l = 0; l < 8; l++) {
			nd[l] = node_index(i + LI(l), j + LJ(l), k + LK(l));
			pv[l] = phi[nd[l]];
		}

		for (int g = 0; g < 8; g++) {
			const double xi  = (LI(g) ? +gp : -gp);
			const double eta = (LJ(g) ? +gp : -gp);
			const double zta = (LK(g) ? +gp : -gp);

			double gr[3] = {0, 0, 0}, nv[8];
			for (int l = 0; l < 8; l++) {
				const double si = (LI(l) ? +1.0 : -1.0);
				const double sj = (LJ(l) ? +1.0 : -1.0);
				const double sk = (LK(l) ? +1.0 : -1.0);
				gr[0] += pv[l] * si * (1 + (sj * eta)) * (1 + (sk * zta)) / 8 * rx;
				gr[1] += pv[l] * sj * (1 + (sk * zta)) * (1 + (si * xi )) / 8 * ry;
				gr[2] += pv[l] * sk * (1 + (si * xi )) * (1 + (sj * eta)) / 8 * rz;
				nv[l] = (1 + (si * xi)) * (1 + (sj * eta)) * (1 + (sk * zta)) / 8;
			}
			const double qg = sg * ((gr[0] * gr[0]) + (gr[1] * gr[1]) + (gr[2] * gr[2])) * detj;
			for (int l = 0; l < 8; l++) f[nd[l]] += qg * nv[l];
		}
	}
	}
	}
}


// 静磁場の行列をセル毎の ν で作る (非線形 B-H 反復)
void assemble_nu(crs_t *A, const double *nucell)
{
	// nucell が無いときは異方性テンソルを ∇Az の基底で使う (mode 4)
	assemble_core(A, ((nucell != NULL) ? 3 : 4), nucell);
}
