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
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i = 0; i < n; i++) {
		const int id = NodeConductor[i];
		phi[i] = ((id == kport) ? Volt : 0);
	}

	// 右辺 b = -(A φ_D) (自由節点のみ、固定節点は 0)
	crs_spmv(A, phi, work, NULL);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i = 0; i < n; i++) {
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
	for (int i = 0; i < n; i++) {
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

	return ierr;
}
