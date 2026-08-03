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
void material_coef_pub(int m, int mode, double c[6])
{
	const material_t *mt = &Material[m];

	// εr / μr の温度依存はここ (読み出し時) で掛ける。
	// σ のように入力解釈で一度だけ掛ける方式は使えない : material_freq() が
	// 分散材料の εr を毎回展開し直して上書きするので、周波数掃引で補正が消える。
	// 読み出し時なら何度呼んでも累積しない
	const double fe = 1 + (mt->epstempco * (Temperature - mt->epstemp0));
	const double fm = 1 + (mt->mutempco  * (Temperature - mt->mutemp0));

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
		*nu = 1 / (MU0 * mt->mu6[0]);
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

	const double chi = ja_dmdh(ja, h, m, delta);

	*mout = m;
	*hout = h;
	*dhdb = 1 / (MU0 * (1 + chi));
}


// ヒステリシス材料の内力ベクトル res = ∫ ν_ja(|B|) ∇A・∇N_i dV と、
// 固定点行列に使う微分磁気抵抗率 nucell (セル毎、正定値)。
// ヒステリシスでは ν = H/B が負になり得る (第 2 象限) ため、行列には
// 常に正の微分磁気抵抗率 dH/dB を使い、残差側だけで正確な H を扱う。
// 反復中の状態は JaBn / JaMn に書き、収束後に JaB / JaM へ確定させる。
void assemble_ja(const double *az, double *res, double *nucell)
{
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
				nu = nud = 1 / (MU0 * mt->mu6[0]);
			}
			nusum += nud;

			for (int l = 0; l < 8; l++) {
				const double u = (dn[l][0] * gr[0]) + (dn[l][1] * gr[1]) + (dn[l][2] * gr[2]);
				re[l] += nu * u * detj;
			}
		}

		nucell[c] = nusum / 8;
		for (int l = 0; l < 8; l++) {
			res[nd[l]] += re[l];
		}
	}
	}
	}
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


void assemble(crs_t *A, int mode)
{
	if (MeshMode) {
		if (MeshElem == MESHELEM_HEX) assemble_hex(A, mode);
		else                          assemble_tet(A, mode);
		return;
	}
	assemble_core(A, mode, NULL);
}


// 静磁場の行列をセル毎の ν で作る (非線形 B-H 反復)
void assemble_nu(crs_t *A, const double *nucell)
{
	// nucell が無いときは異方性テンソルを ∇Az の基底で使う (mode 4)
	assemble_core(A, ((nucell != NULL) ? 3 : 4), nucell);
}
