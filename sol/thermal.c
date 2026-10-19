/*
thermal.c

電気-熱連成 (定常、構造格子の analysis = R)。

ポートに電流 I を流したときの Joule 発熱 q = σ|∇φ|² で温度分布を解き、
導電率を要素ごとに σ(T) = σ0 / (1 + α (T - T0)) (tempco) で更新して、
抵抗 R と最高温度が動かなくなるまで繰り返す (Picard 反復):

	∇・(σ(T) ∇φ) = 0           φ = V (ポート), 0 (基準導体)
	∇・(k ∇T)    = -σ(T)|∇φ|²  T = T_sink (heatsink を付けた導体)
	-k ∂T/∂n = h (T - T_amb)    露出面 (convection)、指定が無ければ断熱

電気の系は σ を固定すれば線形なので、V = Volt で 1 回解いて電流が I になるよう
φ を一様に縮尺すれば足りる (R = V / I_V は縮尺に依らない)。

**要素の温度は 8 節点の平均** (直方体の 3 線形補間の体積平均に厳密に一致)。
要素内で σ を一定にするので、温度が場所で変わる問題では h² の離散化誤差が出る
(1 次元の閉形式との比較で誤差比が 4 になることを rlc_check.sh で見る)。

厳密に成り立つ恒等式 (rlc_check.sh で見る):
	Σ_i f_i (場の発熱の積分) = I² R (端子の電力)            [joule_load の注釈]
	ヒートシンク + 対流で出ていく熱 = Σ_i f_i               (K_T の列和が 0)
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


// 対流 -k dT/dn = h (T - T_amb) の面。熱を通すセルの 6 面のうち、隣が格子の外か
// 熱を通さないセルのものが露出面。電極の面 (4 節点がすべて導体の節点) は端子に
// 接しているので除く。
//   B != NULL : 面の質量 h ∫N_i N_j dA を足す (長方形なので (a/6)(1+δ) x (b/6)(1+δ)、厳密)
//   g != NULL : 右辺 h T_amb ∫N_i dA = h T_amb ab/4 を足す
//   t != NULL : 対流で出ていく熱 h ∫(T - T_amb) dA を返す (双 1 次の面の平均 = 4 節点の平均)
// 戻り値は t があればその熱、なければ面の数
static double conv_faces(crs_t *B, double *g, const double *t, const double *kcell)
{
	double qsum = 0;
	int64_t nface = 0;

	for (int i = 0; i < Nx; i++) {
	for (int j = 0; j < Ny; j++) {
	for (int k = 0; k < Nz; k++) {
		const int64_t c = ((int64_t)i * Ny + j) * Nz + k;
		if (kcell[c] <= 0) continue;
		double h[3];
		h[0] = Xn[i + 1] - Xn[i];
		h[1] = Yn[j + 1] - Yn[j];
		h[2] = Zn[k + 1] - Zn[k];
		for (int d = 0; d < 3; d++) {
		for (int s = 0; s < 2; s++) {
			// 隣のセル
			int ni = i, nj = j, nk = k, out = 0;
			if (d == 0) { ni += (s ? 1 : -1); out = (ni < 0) || (ni >= Nx); }
			if (d == 1) { nj += (s ? 1 : -1); out = (nj < 0) || (nj >= Ny); }
			if (d == 2) { nk += (s ? 1 : -1); out = (nk < 0) || (nk >= Nz); }
			if (!out && (kcell[((int64_t)ni * Ny + nj) * Nz + nk] > 0)) continue;

			// 面の 4 節点 (面内の 2 軸 p < q の局所座標 (a, b))
			const int p = ((d == 0) ? 1 : 0);
			const int q = ((d == 2) ? 1 : 2);
			int64_t nd[4];
			int li[4], lj[4], lk[4], ncond = 0;
			for (int m = 0; m < 4; m++) {
				int l3[3];
				l3[d] = s;
				l3[p] = (m >> 1) & 1;
				l3[q] = m & 1;
				li[m] = i + l3[0];
				lj[m] = j + l3[1];
				lk[m] = k + l3[2];
				nd[m] = node_index(li[m], lj[m], lk[m]);
				if (NodeConductor[nd[m]] >= 0) ncond++;
			}
			if (ncond == 4) continue;	// 電極の面
			nface++;

			const double area = h[p] * h[q];
			if (B) {
				for (int m = 0; m < 4; m++) {
					const int64_t rowstart = B->rowptr[nd[m]];
					for (int mm = 0; mm < 4; mm++) {
						const double fa = (((m >> 1) == (mm >> 1)) ? 2 : 1) * h[p] / 6;
						const double fb = (((m & 1) == (mm & 1)) ? 2 : 1) * h[q] / 6;
						const int64_t pos = crs_offset(rowstart, li[m], lj[m], lk[m],
							li[mm] - li[m], lj[mm] - lj[m], lk[mm] - lk[m]);
						B->val[pos] += ConvH * fa * fb;
					}
				}
			}
			if (g) {
				for (int m = 0; m < 4; m++) g[nd[m]] += ConvH * ConvT * area / 4;
			}
			if (t) {
				const double tm = (t[nd[0]] + t[nd[1]] + t[nd[2]] + t[nd[3]]) / 4;
				qsum += ConvH * area * (tm - ConvT);
			}
		}
		}
	}
	}
	}

	return (t ? qsum : (double)nface);
}


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
	if (ConvH > 0) {
		fprintf(fp_log, "  convection : h = %.6e [W/(m^2 K)], T_amb = %.4f [degC]\n",
			ConvH, ConvT);
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
	double *gc   = (double *)malloc((size_t)n * sizeof(double));
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
	// 対流は熱伝導の行列に足す (温度に依らない)。節点の固定は熱伝導だけで決める
	for (int i = 0; i < n; i++) gc[i] = 0;
	int64_t nface = 0;
	if (ConvH > 0) {
		nface = (int64_t)conv_faces(&B, gc, NULL, kcell);
		fprintf(fp_log, "  convection : %lld exposed faces\n", (long long)nface);
	}
	if ((nsink < 1) && (nface < 1)) {
		fprintf(fp_log, "*** thermal : the heat has no way out (the heat-sink conductors "
			"have no node and no face is exposed to convection)\n");
		ierr = 1;
	}

	double rprev = 0, tmprev = 0, r = 0, g = 0, tmax = 0;
	double ptot = 0, qsink = 0, qconv = 0, pterm = 0;
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
		for (int i = 0; i < n; i++) rhs[i] = (fixt[i] ? 0 : (f[i] + gc[i] - work[i]));
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

		// ヒートシンクへ出ていく熱と対流で出ていく熱。K_T の列和が 0 なので
		// 和は Σ f に一致するはず (対流の行列と右辺の和は対流の熱そのもの)
		qsink = 0;
		for (int i = 0; i < n; i++) {
			const int id = NodeConductor[i];
			if ((id >= 0) && HeatSinkOn[id]) qsink += f[i] + gc[i] - crs_row_dot(&B, i, tn);
		}
		qconv = ((ConvH > 0) ? conv_faces(NULL, NULL, tn, kcell) : 0);

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
			"to heat sinks = %.12e, to convection = %.12e [W]\n", pterm, ptot, qsink, qconv);
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
	free(f); free(tn); free(rhs); free(u); free(gc); free(q);

	return ierr;
}
