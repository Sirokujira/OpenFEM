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

// 1 ポート励振の求解 -> 各導体の反作用 q[0..NPort]
static int solve_port(const crs_t *A, const unsigned char *fix, int kport,
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
	crs_alloc(&A);

	// セル毎の磁気抵抗率 ν。非線形材料 (B-H) があると反復で更新する
	const int64_t ncell = (int64_t)Nx * Ny * Nz;
	int nonlinear = 0;
	for (int m = 0; m < NMaterial; m++) {
		if (Material[m].nbh > 0) nonlinear = 1;
	}
	double *nucell = (double *)malloc((size_t)ncell * sizeof(double));
	for (int64_t c = 0; c < ncell; c++) {
		const material_t *mt = &Material[CellMaterial[c]];
		nucell[c] = ((mt->nbh > 0) ? bh_nu(mt, 0) : (1 / (MU0 * mt->mur)));
	}
	assemble_nu(&A, nucell);

	// 定数分の不定性を除くために 1 節点だけ固定する
	unsigned char *fix = (unsigned char *)malloc(n * sizeof(unsigned char));
	memset(fix, 0, n * sizeof(unsigned char));
	fix[0] = 1;

	// ポート毎の電流密度ベクトル (要素の形状関数で節点に配分)
	double **b = (double **)malloc((size_t)np * sizeof(double *));
	for (int k = 1; k <= np; k++) {
		double *bk = (double *)malloc(n * sizeof(double));
		memset(bk, 0, n * sizeof(double));
		for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
		for (int m = 0; m < Nz; m++) {
			const int id = CellConductor[((int64_t)i * Ny + j) * Nz + m];
			double jz = 0;
			if      (id == k) jz = +Curr / CondArea[k];
			else if (id == 0) jz = -Curr / CondArea[0];
			else continue;
			const double vol = (Xn[i + 1] - Xn[i]) * (Yn[j + 1] - Yn[j]) * (Zn[m + 1] - Zn[m]);
			const double w = jz * vol / 8;
			for (int l = 0; l < 8; l++) {
				bk[node_index(i + ((l >> 2) & 1), j + ((l >> 1) & 1), m + (l & 1))] += w;
			}
		}
		}
		}
		b[k - 1] = bk;
	}

	double *rhs = (double *)malloc(n * sizeof(double));
	double *az  = (double *)malloc(n * sizeof(double));
	const double scale = 1 / TlineLength;

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
			assemble_newton(&J, az, res);

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
	if (hmax > (delta / 2)) {
		fprintf(fp_log, "*** warning : conductor cell size %.4e [m] is coarser than "
			"half the skin depth; R(f) will be underestimated\n", hmax);
		printf("*** warning : conductor mesh does not resolve the skin depth "
			"(cell %.3e > delta/2 = %.3e)\n", hmax, delta / 2);
	}

	fprintf(fp_log, "  %-10s %8s %13s\n", "excite", "iter", "residual");
	fflush(fp_log);

	crs_t K, M;
	crs_alloc(&K);
	assemble(&K, 3);				// ν = 1/(μ0 μr)
	crs_alloc(&M);
	assemble_mass(&M);				// σ

	// 自然境界条件 (磁気壁) なので固定節点は無い
	unsigned char *fix = (unsigned char *)malloc(n * sizeof(unsigned char));
	memset(fix, 0, n * sizeof(unsigned char));

	// 導体毎の ∫σ dV
	double *svol = (double *)calloc((size_t)nc, sizeof(double));
	for (int i = 0; i < Nx; i++) {
	for (int j = 0; j < Ny; j++) {
	for (int k = 0; k < Nz; k++) {
		const int id = CellConductor[((int64_t)i * Ny + j) * Nz + k];
		if (id < 0) continue;
		svol[id] += CondSigma[id]
			* (Xn[i + 1] - Xn[i]) * (Yn[j + 1] - Yn[j]) * (Zn[k + 1] - Zn[k]);
	}
	}
	}

	double *br = (double *)malloc(n * sizeof(double));
	double *bi = (double *)malloc(n * sizeof(double));
	double *ar = (double *)malloc(n * sizeof(double));
	double *ai = (double *)malloc(n * sizeof(double));
	double *yr = (double *)calloc((size_t)np * np, sizeof(double));
	double *yi = (double *)calloc((size_t)np * np, sizeof(double));
	const double tlen = TlineLength;

	// 基準導体は V'=0。ポート 1..np を順に V'=1 [V/m] で励振する
	for (int jc = 1; jc <= np; jc++) {
		// b_i = ∫ σ N_i V' dV  (V' = 1 [V/m] を導体 jc にのみ与える)
		memset(br, 0, n * sizeof(double));
		memset(bi, 0, n * sizeof(double));
		for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
		for (int k = 0; k < Nz; k++) {
			if (CellConductor[((int64_t)i * Ny + j) * Nz + k] != jc) continue;
			const double vol = (Xn[i + 1] - Xn[i]) * (Yn[j + 1] - Yn[j]) * (Zn[k + 1] - Zn[k]);
			const double w = CondSigma[jc] * vol / 8;
			for (int l = 0; l < 8; l++) {
				br[node_index(i + ((l >> 2) & 1), j + ((l >> 1) & 1), k + (l & 1))] += w;
			}
		}
		}
		}
		for (int i = 0; i < n; i++) {
			if (fix[i]) br[i] = 0;
		}

		char label[BUFSIZ];
		sprintf(label, "cond%d", jc);
		const int iter = solver_cocg(&K, &M, omega, br, bi, ar, ai, fix,
			Solver.maxiter, Solver.nout, Solver.converg, fp_log, label);
		if (iter < 0) {
			fprintf(fp_log, "*** solver did not converge (conductor %d)\n", jc);
			ierr = 1;
		}

		// I_k = (1/t) [ V'_k ∫_k σ dV - jω ∫_k σ Az dV ]
		double *qr = (double *)calloc((size_t)nc, sizeof(double));
		double *qi = (double *)calloc((size_t)nc, sizeof(double));
		for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
		for (int k = 0; k < Nz; k++) {
			const int id = CellConductor[((int64_t)i * Ny + j) * Nz + k];
			if (id < 0) continue;
			const double vol = (Xn[i + 1] - Xn[i]) * (Yn[j + 1] - Yn[j]) * (Zn[k + 1] - Zn[k]);
			const double w = CondSigma[id] * vol / 8;
			for (int l = 0; l < 8; l++) {
				const int64_t nd = node_index(i + ((l >> 2) & 1), j + ((l >> 1) & 1), k + (l & 1));
				qr[id] += w * ar[nd];
				qi[id] += w * ai[nd];
			}
		}
		}
		}
		for (int kc = 1; kc <= np; kc++) {
			// I = (V' ∫σdV - jω ∫σAz dV) / t、 -jω(qr + j qi) = ω qi - jω qr
			yr[((kc - 1) * np) + (jc - 1)] = (((kc == jc) ? svol[kc] : 0) + (omega * qi[kc])) / tlen;
			yi[((kc - 1) * np) + (jc - 1)] = (-omega * qr[kc]) / tlen;
		}
		free(qr);
		free(qi);
	}

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


int solve(FILE *fp_log)
{
	const int n = (int)num_node();
	const int np = NPort;
	int ierr = 0;

	crs_t A;
	crs_alloc(&A);

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
		ierr |= solve_eddy(fp_log);
	}

	return ierr;
}
