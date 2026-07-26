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


// 直方体要素 (dx x dy x dz) の ∫∇Ni・∇Nj dV
void element_matrix(double dx, double dy, double dz, double ke[8][8])
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
				ke[l][m] += ((dn[l][0] * dn[m][0])
				          +  (dn[l][1] * dn[m][1])
				          +  (dn[l][2] * dn[m][2])) * detj;
			}
		}
	}
}


// 材料 m のセルに掛かる係数
//   mode = 0 : ε0 εr           静電界 (容量)
//        = 1 : ε0              真空静電界 (TEM インダクタンス)
//        = 2 : σ + ω ε0 εr tanδ  定常電流界 (導電損 + 誘電損)
//        = 3 : 1 / (μ0 μr)     静磁場 (DC インダクタンス)
static double material_coef(int m, int mode)
{
	if      (mode == 0) {
		return EPS0 * Material[m].epsr;
	}
	else if (mode == 1) {
		return EPS0;
	}
	else if (mode == 2) {
		// 誘電損は等価導電率 σ_d = ω ε0 εr tanδ として扱う (frequency 未指定なら 0)
		const double omega = 2 * PI * Freq;
		return Material[m].sigma + (omega * EPS0 * Material[m].epsr * Material[m].tand);
	}
	else {
		return 1 / (MU0 * Material[m].mur);
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


// 全体行列の作成
void assemble(crs_t *A, int mode)
{
	crs_zero(A);

	for (int i = 0; i < Nx; i++) {
	for (int j = 0; j < Ny; j++) {
	for (int k = 0; k < Nz; k++) {
		const int m = CellMaterial[((int64_t)i * Ny + j) * Nz + k];
		const double coef = material_coef(m, mode);
		if (coef <= 0) continue;

		double ke[8][8];
		element_matrix(Xn[i + 1] - Xn[i], Yn[j + 1] - Yn[j], Zn[k + 1] - Zn[k], ke);

		for (int l = 0; l < 8; l++) {
			const int li = LI(l), lj = LJ(l), lk = LK(l);
			const int64_t row = node_index(i + li, j + lj, k + lk);
			const int64_t rowstart = A->rowptr[row];
			for (int mm = 0; mm < 8; mm++) {
				const int di = LI(mm) - li;
				const int dj = LJ(mm) - lj;
				const int dk = LK(mm) - lk;
				const int64_t p = crs_offset(rowstart, i + li, j + lj, k + lk, di, dj, dk);
				A->val[p] += coef * ke[l][mm];
			}
		}
	}
	}
	}
}
