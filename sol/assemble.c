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


// 直方体要素 (dx x dy x dz) の ∫ (cx ∂Ni/∂x ∂Nj/∂x + cy ... + cz ...) dV
// 等方性材料では c[0]=c[1]=c[2] とする
void element_matrix(double dx, double dy, double dz, const double c[3], double ke[8][8])
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
				          +  (c[2] * dn[l][2] * dn[m][2])) * detj;
			}
		}
	}
}


// 材料 m のセルに掛かる係数
//   mode = 0 : ε0 εr           静電界 (容量)
//        = 1 : ε0              真空静電界 (TEM インダクタンス)
//        = 2 : σ + ω ε0 εr tanδ  定常電流界 (導電損 + 誘電損)
//        = 3 : 1 / (μ0 μr)     静磁場 (DC インダクタンス)
static void material_coef(int m, int mode, double c[3])
{
	const material_t *mt = &Material[m];

	for (int d = 0; d < 3; d++) {
		if      (mode == 0) {
			c[d] = EPS0 * mt->epsr3[d];
		}
		else if (mode == 1) {
			c[d] = EPS0;
		}
		else if (mode == 2) {
			// 誘電損は等価導電率 σ_d = ω ε0 εr tanδ として扱う (frequency 未指定なら 0)
			const double omega = 2 * PI * Freq;
			c[d] = mt->sigma + (omega * EPS0 * mt->epsr3[d] * mt->tand);
		}
		else {
			c[d] = 1 / (MU0 * mt->mur3[d]);
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
double bh_nu(const material_t *mt, double b)
{
	const int n = mt->nbh;
	const double nu0 = mt->bh_h[0] / mt->bh_b[0];

	if (b <= mt->bh_b[0]) {
		return nu0;					// b -> 0 でも有限 (初期透磁率)
	}

	double h;
	if (b >= mt->bh_b[n - 1]) {
		const double slope = (n >= 2)
			? ((mt->bh_h[n - 1] - mt->bh_h[n - 2]) / (mt->bh_b[n - 1] - mt->bh_b[n - 2]))
			: nu0;
		h = mt->bh_h[n - 1] + (slope * (b - mt->bh_b[n - 1]));
	}
	else {
		int i = 0;
		while ((i < n - 2) && (b > mt->bh_b[i + 1])) i++;
		const double t = (b - mt->bh_b[i]) / (mt->bh_b[i + 1] - mt->bh_b[i]);
		h = mt->bh_h[i] + (t * (mt->bh_h[i + 1] - mt->bh_h[i]));
	}

	return h / b;
}


// B-H 曲線から ν と dν/d(B^2) を求める (Newton 反復用)
//   ν(B) = H(B)/B,  dν/dB = (H'B - H)/B^2,  dν/d(B^2) = dν/dB / (2B)
// H(B) が単調増加なら ν + 2 B^2 dν/d(B^2) = dH/dB > 0 なのでヤコビアンは正定値になる
static void bh_nu_dnu(const material_t *mt, double b, double *nu, double *dnudb2)
{
	if (mt->nbh <= 0) {
		*nu = 1 / (MU0 * mt->mur);
		*dnudb2 = 0;
		return;
	}

	const int n = mt->nbh;
	const double nu0 = mt->bh_h[0] / mt->bh_b[0];

	if (b <= mt->bh_b[0]) {
		*nu = nu0;					// 原点と第 1 点を結ぶ直線 (初期透磁率)
		*dnudb2 = 0;
		return;
	}

	double h, hp;
	if (b >= mt->bh_b[n - 1]) {
		hp = (n >= 2)
			? ((mt->bh_h[n - 1] - mt->bh_h[n - 2]) / (mt->bh_b[n - 1] - mt->bh_b[n - 2]))
			: nu0;
		h = mt->bh_h[n - 1] + (hp * (b - mt->bh_b[n - 1]));
	}
	else {
		int i = 0;
		while ((i < n - 2) && (b > mt->bh_b[i + 1])) i++;
		hp = (mt->bh_h[i + 1] - mt->bh_h[i]) / (mt->bh_b[i + 1] - mt->bh_b[i]);
		h = mt->bh_h[i] + (hp * (b - mt->bh_b[i]));
	}

	*nu = h / b;
	*dnudb2 = ((hp * b) - h) / (2 * b * b * b);
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

			double nu, dnudb2;
			bh_nu_dnu(mt, bmag, &nu, &dnudb2);
			if (mt->nbh <= 0) nu = 1 / (MU0 * mt->mur3[0]);		// 線形セル (異方性は非対応)

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
		double coef[3];
		if (coefcell != NULL) {
			coef[0] = coef[1] = coef[2] = coefcell[cid];
		}
		else {
			material_coef(m, mode, coef);
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
	assemble_core(A, mode, NULL);
}


// 静磁場の行列をセル毎の ν で作る (非線形 B-H 反復)
void assemble_nu(crs_t *A, const double *nucell)
{
	assemble_core(A, 3, nucell);
}
