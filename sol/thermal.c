/*
thermal.c

電気-熱連成 (定常、構造格子の analysis = R)。

ポートに電流 I を流したときの Joule 発熱 q = σ|∇φ|² で温度分布を解き、
導電率を要素ごとに σ(T) = σ0 / (1 + α (T - T0)) (tempco) で更新して、
抵抗 R と最高温度が動かなくなるまで繰り返す (Picard 反復):

	∇・(σ(T) ∇φ) = 0           φ = V (ポート), 0 (基準導体)
	∇・(k ∇T)    = -σ(T)|∇φ|²  T = T_sink (heatsink を付けた導体)、他は断熱

電気の系は σ を固定すれば線形なので、V = Volt で 1 回解いて電流が I になるよう
φ を一様に縮尺すれば足りる (R = V / I_V は縮尺に依らない)。

**要素の温度は 8 節点の平均** (直方体の 3 線形補間の体積平均に厳密に一致)。
要素内で σ を一定にするので、温度が場所で変わる問題では h² の離散化誤差が出る
(1 次元の閉形式との比較で誤差比が 4 になることを rlc_check.sh で見る)。

厳密に成り立つ恒等式 (rlc_check.sh で見る):
	Σ_i f_i (場の発熱の積分) = I² R (端子の電力)            [joule_load の注釈]
	ヒートシンクから出ていく熱 = Σ_i f_i                    (K_T の行和が 0)
α = 0 なら σ が温度に依らないので R は熱連成なしの値と一致し、一様断面の棒では
1 次要素の解が節点で厳密になる (定数の発熱源の 1 次元 Poisson)。

熱暴走: ρ が T の 1 次式 (α > 0) の 1 次元の棒では mL > π (m² = J² ρ0 α / k) で
定常解が存在しない。そのとき Picard 反復は温度が増え続けて止まらないので、
上限回数で打ち切って入力エラーにする (黙って最後の値を出さない)。
*/

#include "fem.h"
#include "fem_prototype.h"

#define TI(l) (((l) >> 2) & 1)
#define TJ(l) (((l) >> 1) & 1)
#define TK(l) ((l) & 1)


int solve_thermal(FILE *fp_log, double *gport)
{
	const int n = (int)num_node();
	const int64_t nc = (int64_t)Nx * Ny * Nz;
	int ierr = 0;

	fprintf(fp_log, "\n=== electro-thermal analysis ===\n");
	fprintf(fp_log, "  drive current = %.6e [A], ambient (initial) = %.4f [degC]\n",
		ThermalCurrent, Temperature);
	for (int p = 0; p <= NPort; p++) {
		if (HeatSinkOn[p]) {
			fprintf(fp_log, "  heat sink : conductor %d at %.4f [degC]\n", p, HeatSinkT[p]);
		}
	}
	fflush(fp_log);

	double *tcell = (double *)malloc((size_t)nc * sizeof(double));
	double *scell = (double *)malloc((size_t)nc * sizeof(double));
	double *kcell = (double *)malloc((size_t)nc * sizeof(double));
	for (int64_t c = 0; c < nc; c++) {
		tcell[c] = Temperature;
		kcell[c] = Material[CellMaterial[c]].kth;
	}

	unsigned char *fix  = (unsigned char *)malloc((size_t)n * sizeof(unsigned char));
	unsigned char *fixt = (unsigned char *)malloc((size_t)n * sizeof(unsigned char));
	double *phi  = (double *)malloc((size_t)n * sizeof(double));
	double *work = (double *)malloc((size_t)n * sizeof(double));
	double *diag = (double *)malloc((size_t)n * sizeof(double));
	double *f    = (double *)malloc((size_t)n * sizeof(double));
	double *tn   = (double *)malloc((size_t)n * sizeof(double));
	double *rhs  = (double *)malloc((size_t)n * sizeof(double));
	double *u    = (double *)malloc((size_t)n * sizeof(double));
	double *q    = (double *)malloc(((size_t)NPort + 1) * sizeof(double));

	// 熱の行列は温度に依らないので 1 回だけ作る
	crs_t A, B;
	crs_alloc(&A);
	crs_alloc(&B);
	assemble_cell(&B, kcell);
	crs_diag(&B, diag);
	double dmaxt = 0;
	for (int i = 0; i < n; i++) {
		if (diag[i] > dmaxt) dmaxt = diag[i];
	}
	int64_t nsink = 0;
	for (int i = 0; i < n; i++) {
		const int id = NodeConductor[i];
		if ((id >= 0) && HeatSinkOn[id]) {
			fixt[i] = 1;
			tn[i] = HeatSinkT[id];
			nsink++;
		}
		else if (diag[i] <= (1e-12 * dmaxt)) {
			fixt[i] = 1;					// 熱を通さない領域の節点
			tn[i] = Temperature;
		}
		else {
			fixt[i] = 0;
			tn[i] = 0;
		}
	}
	if (nsink < 1) {
		fprintf(fp_log, "*** thermal : the heat-sink conductors have no node\n");
		ierr = 1;
	}

	double rprev = 0, tmprev = 0, r = 0, g = 0, tmax = 0;
	double ptot = 0, qsink = 0, pterm = 0;
	int it = 0, conv = 0;
	fprintf(fp_log, "  %4s %16s %14s %16s\n", "iter", "R [ohm]", "T_max [degC]", "dR / R");
	for (it = 1; !ierr && (it <= ThermalMaxIter); it++) {
		// σ(T) を要素ごとに。入力解釈で σ には周囲温度での補正
		// σ_s = σ0 / (1 + α(T_amb - T0)) が掛かっているので、それを外して掛け直す。
		// α = 0 なら因子は厳密に 1 で、熱連成なしと同じ σ になる
		for (int64_t c = 0; c < nc; c++) {
			const material_t *mt = &Material[CellMaterial[c]];
			if (mt->sigma <= 0) {
				scell[c] = 0;
				continue;
			}
			const double num = 1 + (mt->tempco * (Temperature - mt->temp0));
			const double den = 1 + (mt->tempco * (tcell[c] - mt->temp0));
			if (den <= 0) {
				// α < 0 で温度が上がりすぎると抵抗率が負になる (モデルの適用範囲外)
				fprintf(fp_log, "*** thermal : the resistivity factor 1 + alpha (T - T0) "
					"became non-positive (T = %.4e degC); the linear sigma(T) model "
					"does not apply here\n", tcell[c]);
				ierr = 1;
				break;
			}
			scell[c] = mt->sigma * num / den;
		}
		if (ierr) break;

		// 電気 : V = Volt で解いて、電流が I になるよう縮尺する
		assemble_cell(&A, scell);
		crs_diag(&A, diag);
		double dmax = 0;
		for (int i = 0; i < n; i++) {
			if (diag[i] > dmax) dmax = diag[i];
		}
		for (int i = 0; i < n; i++) {
			fix[i] = ((NodeConductor[i] >= 0) || (diag[i] <= (1e-12 * dmax))) ? 1 : 0;
		}
		char label[64];
		sprintf(label, "phi%d", it);
		double energy = 0;
		if (solve_port(&A, fix, 1, phi, work, q, &energy, NULL, label) < 0) {
			fprintf(fp_log, "*** thermal : the conduction solve did not converge "
				"(iteration %d)\n", it);
			ierr = 1;
			break;
		}
		const double iv = q[1];
		if (iv <= 0) {
			fprintf(fp_log, "*** thermal : no current flows through the port\n");
			ierr = 1;
			break;
		}
		g = iv / Volt;
		r = 1 / g;
		const double s = ThermalCurrent / iv;
		for (int i = 0; i < n; i++) phi[i] *= s;
		pterm = ThermalCurrent * ThermalCurrent * r;

		// 発熱と熱 : T = T_D + u、自由節点で B u = f - B T_D
		joule_load(phi, scell, f);
		ptot = 0;
		for (int i = 0; i < n; i++) ptot += f[i];
		for (int i = 0; i < n; i++) {
			if (!fixt[i]) tn[i] = 0;
		}
		crs_spmv(&B, tn, work, NULL);
		for (int i = 0; i < n; i++) rhs[i] = (fixt[i] ? 0 : (f[i] - work[i]));
		sprintf(label, "T%d", it);
		if (solver_cg(&B, rhs, u, fixt, Solver.maxiter, Solver.nout, Solver.converg,
				NULL, label) < 0) {
			fprintf(fp_log, "*** thermal : the heat-conduction solve did not converge "
				"(iteration %d)\n", it);
			ierr = 1;
			break;
		}
		for (int i = 0; i < n; i++) {
			if (!fixt[i]) tn[i] += u[i];
		}

		// ヒートシンクへ出ていく熱 (K_T の行和が 0 なので Σ f に一致するはず)
		qsink = 0;
		for (int i = 0; i < n; i++) {
			const int id = NodeConductor[i];
			if ((id >= 0) && HeatSinkOn[id]) qsink += f[i] - crs_row_dot(&B, i, tn);
		}

		// 要素の温度 = 8 節点の平均 (3 線形の体積平均)
		tmax = -1e300;
		for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
		for (int k = 0; k < Nz; k++) {
			const int64_t c = ((int64_t)i * Ny + j) * Nz + k;
			if (kcell[c] <= 0) continue;
			double t = 0;
			for (int l = 0; l < 8; l++) {
				const double tl = tn[node_index(i + TI(l), j + TJ(l), k + TK(l))];
				t += tl / 8;
				if (tl > tmax) tmax = tl;
			}
			tcell[c] = t;
		}
		}
		}

		const double dr = ((rprev > 0) ? fabs(r - rprev) / r : 1);
		const double dt = fabs(tmax - tmprev) / ((fabs(tmax) > 1) ? fabs(tmax) : 1);
		fprintf(fp_log, "  %4d %16.9e %14.6f %16.3e\n", it, r, tmax, dr);
		fflush(fp_log);
		if ((it > 1) && (dr < ThermalTol) && (dt < ThermalTol)) {
			conv = 1;
			break;
		}
		// 熱暴走 : 定常解が無いと温度が増え続ける。物理的に意味の無い温度まで
		// 追いかけずに打ち切る
		if (tmax > 1e5) {
			fprintf(fp_log, "*** thermal runaway : T_max exceeded 1e5 degC at iteration "
				"%d (no steady state for this current)\n", it);
			ierr = 1;
			break;
		}
		rprev = r;
		tmprev = tmax;
	}
	if (!ierr && !conv) {
		fprintf(fp_log, "*** thermal : the coupling did not converge in %d iterations "
			"(thermal runaway: no steady state for this current?)\n", ThermalMaxIter);
		ierr = 1;
	}

	if (!ierr) {
		// 体積加重の平均温度 (導電領域)
		double vsum = 0, tsum = 0;
		for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
		for (int k = 0; k < Nz; k++) {
			const int64_t c = ((int64_t)i * Ny + j) * Nz + k;
			if (scell[c] <= 0) continue;
			const double v = (Xn[i + 1] - Xn[i]) * (Yn[j + 1] - Yn[j]) * (Zn[k + 1] - Zn[k]);
			vsum += v;
			tsum += v * tcell[c];
		}
		}
		}
		fprintf(fp_log, "  converged in %d iterations\n", it);
		fprintf(fp_log, "  R = %.10e [ohm] at I = %.6e [A] (V = %.10e [V])\n",
			r, ThermalCurrent, ThermalCurrent * r);
		fprintf(fp_log, "  T_max = %.10e [degC], T_mean (conductor) = %.10e [degC]\n",
			tmax, ((vsum > 0) ? (tsum / vsum) : Temperature));
		fprintf(fp_log, "  joule power : I^2 R = %.12e, field integral = %.12e, "
			"to heat sinks = %.12e [W]\n", pterm, ptot, qsink);
		*gport = g;

		if (FieldOut) {
			field_add_node("phi_R_port1", phi);
			field_add_grad("E_R_port1", phi, 0);
			field_add_node("T", tn);
		}
	}
	fflush(fp_log);

	crs_free(&A);
	crs_free(&B);
	free(tcell); free(scell); free(kcell);
	free(fix); free(fixt); free(phi); free(work); free(diag);
	free(f); free(tn); free(rhs); free(u); free(q);

	return ierr;
}
