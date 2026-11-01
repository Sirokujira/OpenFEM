/*
thermal.c

電気-熱連成 (定常、構造格子の analysis = R)。

ポートに電流 I を流したときの Joule 発熱 q = σ|∇φ|² で温度分布を解き、
導電率を要素ごとに σ(T) = σ0 / (1 + α (T - T0)) (tempco、または tempcotable の
折れ線 σ0 / f(T)) で更新して、
抵抗 R と最高温度が動かなくなるまで繰り返す (Picard 反復):

	∇・(σ(T) ∇φ) = 0           φ = V (ポート), 0 (基準導体)
	∇・(k ∇T)    = -σ(T)|∇φ|²  T = T_sink (heatsink を付けた導体)
	-k ∂T/∂n = h (T - T_amb)    露出面 (convection)
	         + ε σ_SB (T_K^4 - T_amb,K^4)  (radiation、絶対温度)、指定が無ければ断熱

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


/*
格子の種類に依る部分 (非構造格子は 1 次四面体だけ)。

熱連成の反復 (th_steady / th_transient) は格子に依らず、要素の数・材料・係数の
組み立て・発熱・質量・要素の温度・露出面だけがここで分かれる。**構造格子の
分岐は従来のコードそのもの** (既存の答えはビット単位で変わらない)。

1 次四面体では ∇N が要素内で一定なので、
	剛性   K_ij = c ∇N_i・∇N_j V                    (厳密)
	発熱   f_i  = σ|∇φ|² V / 4  (Σ_i f_i = φᵀKφ が丸めを除いて厳密)
	質量   M_ij = ρc V (1 + δ_ij) / 20               (厳密)
	要素の温度 = 4 節点の平均 (1 次補間の体積平均に厳密に一致)
露出面は熱を通す四面体の面のうち、熱を通す四面体 1 つにしか属さないもの
(格子の外か k = 0 の要素に接する)。電極の面 (3 節点がすべて導体) は除く。
*/
static int64_t th_nelem(void)
{
	return (MeshMode ? (int64_t)NTet : ((int64_t)Nx * Ny * Nz));
}


static int th_emat(int64_t c)
{
	return (MeshMode ? elem3d_mat((int)c) : (int)CellMaterial[c]);
}


// 露出面の直方体の許容幅 (形状判定 setup.c と同じく解析領域の大きさの EPS 倍)
static double SurfEps;

static void th_extent(void)
{
	if (!MeshMode) {
		SurfEps = EPS * (fabs(Xn[Nx] - Xn[0]) + fabs(Yn[Ny] - Yn[0]) + fabs(Zn[Nz] - Zn[0]));
		return;
	}
	double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
	for (int i = 0; i < NNode; i++) {
		const double p[3] = {Xp[i], Yp[i], Zp[i]};
		for (int d = 0; d < 3; d++) {
			if (p[d] < lo[d]) lo[d] = p[d];
			if (p[d] > hi[d]) hi[d] = p[d];
		}
	}
	SurfEps = EPS * ((hi[0] - lo[0]) + (hi[1] - lo[1]) + (hi[2] - lo[2]));
}


// 四面体の露出面 [3*TFaceN] (熱を通す四面体の境界の三角形、電極の面を除く)
static int32_t *TFace;
static int TFaceN;

static int tface_cmp(const void *a, const void *b)
{
	const int32_t *p = (const int32_t *)a, *q = (const int32_t *)b;
	for (int d = 0; d < 3; d++) {
		if (p[d] != q[d]) return ((p[d] < q[d]) ? -1 : 1);
	}
	return 0;
}

static void tet_faces_build(const double *kcell)
{
	static const int fl[4][3] = {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}};
	int nf = 0;
	for (int e = 0; e < NTet; e++) {
		if (kcell[e] > 0) nf += 4;
	}
	int32_t *all = (int32_t *)malloc(((size_t)nf * 3 + 1) * sizeof(int32_t));
	int m = 0;
	for (int e = 0; e < NTet; e++) {
		if (kcell[e] <= 0) continue;
		for (int f = 0; f < 4; f++) {
			int32_t v[3];
			for (int l = 0; l < 3; l++) v[l] = Tet[(e * 4) + fl[f][l]];
			// 昇順に並べて面の鍵にする
			for (int a = 0; a < 2; a++) {
				for (int b = a + 1; b < 3; b++) {
					if (v[b] < v[a]) { const int32_t w = v[a]; v[a] = v[b]; v[b] = w; }
				}
			}
			for (int l = 0; l < 3; l++) all[(m * 3) + l] = v[l];
			m++;
		}
	}
	qsort(all, (size_t)nf, 3 * sizeof(int32_t), tface_cmp);
	TFace = (int32_t *)malloc(((size_t)nf * 3 + 1) * sizeof(int32_t));
	TFaceN = 0;
	for (int i = 0; i < nf; ) {
		int j = i + 1;
		while ((j < nf) && !tface_cmp(&all[i * 3], &all[j * 3])) j++;
		if ((j - i) == 1) {
			const int32_t *v = &all[i * 3];
			const int ncond = (NodeConductor[v[0]] >= 0) + (NodeConductor[v[1]] >= 0)
			                + (NodeConductor[v[2]] >= 0);
			if (ncond < 3) {
				for (int l = 0; l < 3; l++) TFace[(TFaceN * 3) + l] = v[l];
				TFaceN++;
			}
		}
		i = j;
	}
	free(all);
}


// 面の中心 (xc, yc, zc) が行 bc の直方体に入るか (直方体なしの行はどこでも入る)。
// 許容幅は形状判定 (setup.c) と同じく解析領域の大きさの EPS 倍
static int surf_in(const surfbc_t *bc, double xc, double yc, double zc)
{
	if (!bc->box) return 1;
	const double eps = SurfEps;
	return (xc >= (bc->g[0] - eps)) && (xc <= (bc->g[1] + eps))
	    && (yc >= (bc->g[2] - eps)) && (yc <= (bc->g[3] + eps))
	    && (zc >= (bc->g[4] - eps)) && (zc <= (bc->g[5] + eps));
}


// 面の中心での対流 (rad = 0) / 輻射 (rad = 1) の行。あとに書いた行ほど優先し、
// どの行にも入らなければ NULL (断熱)
static const surfbc_t *surf_bc(int rad, double xc, double yc, double zc, int *idx)
{
	const surfbc_t *list = (rad ? RadBC : ConvBC);
	const int nb = (rad ? NRad : NConv);
	for (int e = nb - 1; e >= 0; e--) {
		if (surf_in(&list[e], xc, yc, zc)) {
			if (idx) *idx = e;
			return &list[e];
		}
	}
	return NULL;
}


// 露出面の境界条件。熱を通すセルの 6 面のうち、隣が格子の外か熱を通さないセルの
// ものが露出面。電極の面 (4 節点がすべて導体の節点) は端子に接しているので除く。
//
// rad = 0 : 対流 -k dT/dn = h (T - T_amb) (線形)
//   B != NULL : 面の質量 h ∫N_i N_j dA を足す (長方形なので (a/6)(1+δ) x (b/6)(1+δ)、厳密)
//   g != NULL : 右辺 h T_amb ∫N_i dA = h T_amb ab/4 を足す
//   t != NULL : 対流で出ていく熱 h ∫(T - T_amb) dA を返す (双 1 次の面の平均 = 4 節点の平均)
//   戻り値は t があればその熱、なければ面の数
// rad = 1 : 輻射 -k dT/dn = q(T) = ε σ_SB (T_K^4 - T_amb,K^4) (非線形、t は必須)
//   t のまわりで Newton 線形化する: q(T) ≈ q(t) + q'(t) (T - t)
//   B != NULL : 接線 ∫q'(t) N_i N_j dA を足す
//   g != NULL : 右辺 ∫(q'(t) t - q(t)) N_i dA を足す
//   戻り値は輻射で出ていく熱 ∫q(t) dA。積分は面の 2x2 Gauss で、接線と右辺と
//   熱で評価点を揃える (揃えないと Newton が 2 次収束しない)。
//   絶対温度が正でない点があれば -1e300 を返す
// rad = 2 : 数えるだけ。対流か輻射の係数が正の面の数を返し、hits があれば
//   各行 (対流 0..NConv-1、輻射 NConv..) の直方体に中心が入る露出面の数を足す
//
// 面の係数と周囲温度は surf_bc() で面の中心から引く (あとに書いた行ほど優先)
static double surf_faces_grid(crs_t *B, double *g, const double *t, const double *kcell,
	int rad, int *hits)
{
	const double gp[2] = {0.5 - 0.5 / sqrt(3.0), 0.5 + 0.5 / sqrt(3.0)};
	double qsum = 0;
	int64_t nface = 0;
	int bad = 0;

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

			// 面の中心で係数を引く
			double xc = 0, yc = 0, zc = 0;
			for (int m = 0; m < 4; m++) {
				xc += Xn[li[m]] / 4;
				yc += Yn[lj[m]] / 4;
				zc += Zn[lk[m]] / 4;
			}
			if (rad == 2) {
				const surfbc_t *cb = surf_bc(0, xc, yc, zc, NULL);
				const surfbc_t *rb = surf_bc(1, xc, yc, zc, NULL);
				if ((cb && (cb->c > 0)) || (rb && (rb->c > 0))) nface++;
				if (hits) {
					for (int e = 0; e < (NConv + NRad); e++) {
						const surfbc_t *bc = ((e < NConv) ? &ConvBC[e] : &RadBC[e - NConv]);
						if (surf_in(bc, xc, yc, zc)) hits[e]++;
					}
				}
				continue;
			}
			const surfbc_t *bc = surf_bc(rad, xc, yc, zc, NULL);
			if ((bc == NULL) || (bc->c <= 0)) continue;
			nface++;
			const double hc = bc->c, ta = bc->t;

			const double area = h[p] * h[q];
			double me[4][4];
			double ge[4] = {0, 0, 0, 0};
			if (!rad) {
				for (int m = 0; m < 4; m++) {
					for (int mm = 0; mm < 4; mm++) {
						const double fa = (((m >> 1) == (mm >> 1)) ? 2 : 1) * h[p] / 6;
						const double fb = (((m & 1) == (mm & 1)) ? 2 : 1) * h[q] / 6;
						me[m][mm] = hc * fa * fb;
					}
					ge[m] = hc * ta * area / 4;
				}
				if (t) {
					const double tm = (t[nd[0]] + t[nd[1]] + t[nd[2]] + t[nd[3]]) / 4;
					qsum += hc * area * (tm - ta);
				}
			}
			else {
				const double ta4 = pow(ta + KELVIN, 4);
				for (int m = 0; m < 4; m++) {
					for (int mm = 0; mm < 4; mm++) me[m][mm] = 0;
				}
				for (int ga = 0; ga < 2; ga++) {
				for (int gb = 0; gb < 2; gb++) {
					double nn[4], tg = 0;
					for (int m = 0; m < 4; m++) {
						nn[m] = (((m >> 1) & 1) ? gp[ga] : (1 - gp[ga]))
						      * ((m & 1) ? gp[gb] : (1 - gp[gb]));
						tg += nn[m] * t[nd[m]];
					}
					const double tk = tg + KELVIN;
					if (tk <= 0) bad = 1;
					const double w = area / 4;
					const double qv = hc * STEFAN * (tk * tk * tk * tk - ta4);
					const double dq = 4 * hc * STEFAN * tk * tk * tk;
					qsum += w * qv;
					for (int m = 0; m < 4; m++) {
						ge[m] += w * (dq * tg - qv) * nn[m];
						for (int mm = 0; mm < 4; mm++) me[m][mm] += w * dq * nn[m] * nn[mm];
					}
				}
				}
			}

			if (B) {
				for (int m = 0; m < 4; m++) {
					const int64_t rowstart = B->rowptr[nd[m]];
					for (int mm = 0; mm < 4; mm++) {
						const int64_t pos = crs_offset(rowstart, li[m], lj[m], lk[m],
							li[mm] - li[m], lj[mm] - lj[m], lk[mm] - lk[m]);
						B->val[pos] += me[m][mm];
					}
				}
			}
			if (g) {
				for (int m = 0; m < 4; m++) g[nd[m]] += ge[m];
			}
		}
		}
	}
	}
	}

	if (bad) return -1e300;
	return ((t || (rad == 1)) ? qsum : (double)nface);
}


// surf_faces_grid の四面体版 (露出面は tet_faces_build で作った三角形)。
// 対流の面の質量は h A (1 + δ_ij) / 12、右辺は h T_amb A / 3、熱は
// h A (3 節点の平均 - T_amb) (いずれも 1 次の三角形で厳密)。
// 輻射は三角形の 3 点則 (重心座標 (2/3, 1/6, 1/6) の巡回、重み A/3、2 次まで厳密) で、
// 接線と右辺と熱で評価点を揃える
static double surf_faces_tet(crs_t *B, double *g, const double *t, int rad, int *hits)
{
	double qsum = 0;
	int64_t nface = 0;
	int bad = 0;

	for (int f = 0; f < TFaceN; f++) {
		const int32_t *nd = &TFace[f * 3];
		const double xc = (Xp[nd[0]] + Xp[nd[1]] + Xp[nd[2]]) / 3;
		const double yc = (Yp[nd[0]] + Yp[nd[1]] + Yp[nd[2]]) / 3;
		const double zc = (Zp[nd[0]] + Zp[nd[1]] + Zp[nd[2]]) / 3;
		if (rad == 2) {
			const surfbc_t *cb = surf_bc(0, xc, yc, zc, NULL);
			const surfbc_t *rb = surf_bc(1, xc, yc, zc, NULL);
			if ((cb && (cb->c > 0)) || (rb && (rb->c > 0))) nface++;
			if (hits) {
				for (int e = 0; e < (NConv + NRad); e++) {
					const surfbc_t *bc = ((e < NConv) ? &ConvBC[e] : &RadBC[e - NConv]);
					if (surf_in(bc, xc, yc, zc)) hits[e]++;
				}
			}
			continue;
		}
		const surfbc_t *bc = surf_bc(rad, xc, yc, zc, NULL);
		if ((bc == NULL) || (bc->c <= 0)) continue;
		nface++;
		const double hc = bc->c, ta = bc->t;

		const double ax = Xp[nd[1]] - Xp[nd[0]], ay = Yp[nd[1]] - Yp[nd[0]], az = Zp[nd[1]] - Zp[nd[0]];
		const double bx = Xp[nd[2]] - Xp[nd[0]], by = Yp[nd[2]] - Yp[nd[0]], bz = Zp[nd[2]] - Zp[nd[0]];
		const double cx = (ay * bz) - (az * by);
		const double cy = (az * bx) - (ax * bz);
		const double cz = (ax * by) - (ay * bx);
		const double area = sqrt((cx * cx) + (cy * cy) + (cz * cz)) / 2;

		double me[3][3];
		double ge[3] = {0, 0, 0};
		if (!rad) {
			for (int m = 0; m < 3; m++) {
				for (int mm = 0; mm < 3; mm++) me[m][mm] = hc * area * ((m == mm) ? 2 : 1) / 12;
				ge[m] = hc * ta * area / 3;
			}
			if (t) {
				const double tm = (t[nd[0]] + t[nd[1]] + t[nd[2]]) / 3;
				qsum += hc * area * (tm - ta);
			}
		}
		else {
			const double ta4 = pow(ta + KELVIN, 4);
			for (int m = 0; m < 3; m++) {
				for (int mm = 0; mm < 3; mm++) me[m][mm] = 0;
			}
			for (int gq = 0; gq < 3; gq++) {
				double nn[3], tg = 0;
				for (int m = 0; m < 3; m++) {
					nn[m] = ((m == gq) ? (2.0 / 3) : (1.0 / 6));
					tg += nn[m] * t[nd[m]];
				}
				const double tk = tg + KELVIN;
				if (tk <= 0) bad = 1;
				const double w = area / 3;
				const double qv = hc * STEFAN * (tk * tk * tk * tk - ta4);
				const double dq = 4 * hc * STEFAN * tk * tk * tk;
				qsum += w * qv;
				for (int m = 0; m < 3; m++) {
					ge[m] += w * (dq * tg - qv) * nn[m];
					for (int mm = 0; mm < 3; mm++) me[m][mm] += w * dq * nn[m] * nn[mm];
				}
			}
		}

		if (B) {
			for (int m = 0; m < 3; m++) {
				for (int mm = 0; mm < 3; mm++) {
					const int64_t pos = crs_find_tri(B, nd[m], nd[mm]);
					if (pos >= 0) B->val[pos] += me[m][mm];
				}
			}
		}
		if (g) {
			for (int m = 0; m < 3; m++) g[nd[m]] += ge[m];
		}
	}

	if (bad) return -1e300;
	return ((t || (rad == 1)) ? qsum : (double)nface);
}


static double surf_faces(crs_t *B, double *g, const double *t, const double *kcell, int rad,
	int *hits)
{
	if (MeshMode) return surf_faces_tet(B, g, t, rad, hits);
	return surf_faces_grid(B, g, t, kcell, rad, hits);
}


// ---- 要素ごとの係数の組み立て・発熱・質量 (四面体) ----

static void tet_assemble(crs_t *A, const double *coef)
{
	crs_zero(A);
	for (int e = 0; e < NTet; e++) {
		const double c = coef[e];
		if (c <= 0) continue;
		const int32_t *nd = &Tet[e * 4];
		double gr[4][3], vol;
		if (tet_grad_pub(nd, gr, &vol)) continue;
		for (int l = 0; l < 4; l++) {
			for (int m = 0; m < 4; m++) {
				// assemble_elem3d の等方性の式と同じ順 (α = 0 で熱連成なしと同じ行列)
				const double v = (c * gr[l][0] * gr[m][0])
				               + (c * gr[l][1] * gr[m][1])
				               + (c * gr[l][2] * gr[m][2]);
				const int64_t p = crs_find_tri(A, nd[l], nd[m]);
				if (p >= 0) A->val[p] += v * vol;
			}
		}
	}
}


static void tet_joule(const double *phi, const double *sig, double *f)
{
	for (int i = 0; i < NNode; i++) f[i] = 0;
	for (int e = 0; e < NTet; e++) {
		const double sg = sig[e];
		if (sg <= 0) continue;
		const int32_t *nd = &Tet[e * 4];
		double gr[4][3], vol;
		if (tet_grad_pub(nd, gr, &vol)) continue;
		double ex = 0, ey = 0, ez = 0;
		for (int l = 0; l < 4; l++) {
			ex += phi[nd[l]] * gr[l][0];
			ey += phi[nd[l]] * gr[l][1];
			ez += phi[nd[l]] * gr[l][2];
		}
		const double qe = sg * ((ex * ex) + (ey * ey) + (ez * ez)) * vol / 4;
		for (int l = 0; l < 4; l++) f[nd[l]] += qe;
	}
}


static void tet_mass(crs_t *A, const double *coef)
{
	crs_zero(A);
	for (int e = 0; e < NTet; e++) {
		const double c = coef[e];
		if (c <= 0) continue;
		const int32_t *nd = &Tet[e * 4];
		double gr[4][3], vol;
		if (tet_grad_pub(nd, gr, &vol)) continue;
		for (int l = 0; l < 4; l++) {
			for (int m = 0; m < 4; m++) {
				const int64_t p = crs_find_tri(A, nd[l], nd[m]);
				if (p >= 0) A->val[p] += c * vol * ((l == m) ? 2 : 1) / 20;
			}
		}
	}
}


static void th_assemble(crs_t *A, const double *coef)
{
	if (MeshMode) tet_assemble(A, coef);
	else          assemble_cell(A, coef);
}


static void th_joule(const double *phi, const double *sig, double *f)
{
	if (MeshMode) tet_joule(phi, sig, f);
	else          joule_load(phi, sig, f);
}


static void th_mass(crs_t *A, const double *coef)
{
	if (MeshMode) tet_mass(A, coef);
	else          assemble_mass_cell(A, coef);
}


// σ(T) を要素ごとに。入力解釈で σ には周囲温度での補正
// σ_s = σ0 / (1 + α(T_amb - T0)) が掛かっているので、それを外して掛け直す。
// α = 0 なら因子は厳密に 1 で、熱連成なしと同じ σ になる
static int th_sigma(FILE *fp_log, int64_t nc, const double *tcell, double *scell)
{
	for (int64_t c = 0; c < nc; c++) {
		const material_t *mt = &Material[th_emat(c)];
		if (mt->sigma <= 0) {
			scell[c] = 0;
			continue;
		}
		const double num = temp_factor_at(mt, 3, Temperature);
		const double den = temp_factor_at(mt, 3, tcell[c]);
		if (den <= 0) {
			// α < 0 で温度が上がりすぎると抵抗率が負になる (モデルの適用範囲外)
			fprintf(fp_log, "*** thermal : the resistivity factor 1 + alpha (T - T0) "
				"became non-positive (T = %.4e degC); the linear sigma(T) model "
				"does not apply here\n", tcell[c]);
			return 1;
		}
		scell[c] = mt->sigma * num / den;
	}
	return 0;
}


// 電気 : V = Volt で解いて、電流が cur になるよう φ を縮尺する。*gout = I_V / V
// (cur = 0 でも R = V / I_V は決まる。発熱は 0)
static int th_conduct(FILE *fp_log, crs_t *A, const double *scell, unsigned char *fix,
	double *phi, double *work, double *diag, double *q, int it, double cur, double *gout)
{
	const int n = (int)A->n;

	th_assemble(A, scell);
	crs_diag(A, diag);
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
	if (solve_port(A, fix, 1, phi, work, q, &energy, NULL, label) < 0) {
		fprintf(fp_log, "*** thermal : the conduction solve did not converge "
			"(iteration %d)\n", it);
		return 1;
	}
	const double iv = q[1];
	if (iv <= 0) {
		fprintf(fp_log, "*** thermal : no current flows through the port\n");
		return 1;
	}
	*gout = iv / Volt;
	const double s = cur / iv;
	for (int i = 0; i < n; i++) phi[i] *= s;
	return 0;
}


// 交流の発熱 (analysis = A) : 電気の段を A-φ で解く (eddy3d.c)。電流は波高値で、
// 発熱は時間平均 ½|I|²R
static double ThX;		// 最後に解いた Z の虚部 (交流のとき)

static int th_ac(void)
{
	return ((Analysis & ANALYSIS_A) != 0);
}


static double th_pfac(void)
{
	return (th_ac() ? 0.5 : 1.0);
}


// 電気の段 : 電流 cur での発熱の節点荷重 f と *gout = 1 / R を作る。
// DC は th_conduct + th_joule (従来の経路そのもの)、交流は eddy3d_heat_solve
static int th_electric(FILE *fp_log, crs_t *A, const double *scell, unsigned char *fix,
	double *phi, double *work, double *diag, double *q, int it, double cur, double *gout,
	double *f)
{
	if (th_ac()) {
		double zr = 0, zi = 0;
		if (eddy3d_heat_solve(fp_log, scell, cur, f, &zr, &zi, it)) return 1;
		if (zr <= 0) {
			fprintf(fp_log, "*** thermal : non-passive impedance R = %.6e (iteration %d)\n",
				zr, it);
			return 1;
		}
		*gout = 1 / zr;
		ThX = zi;
		return 0;
	}
	const int e = th_conduct(fp_log, A, scell, fix, phi, work, diag, q, it, cur, gout);
	if (e) return e;
	th_joule(phi, scell, f);
	return 0;
}


// 要素の温度 = 8 節点の平均 (3 線形の体積平均)。戻り値は熱を通す領域の節点の最高温度
static double th_celltemp(const double *tn, const double *kcell, double *tcell)
{
	double tmax = -1e300;
	if (MeshMode) {
		// 四面体 : 4 節点の平均 (1 次補間の体積平均)
		for (int e = 0; e < NTet; e++) {
			if (kcell[e] <= 0) continue;
			double t = 0;
			for (int l = 0; l < 4; l++) {
				const double tl = tn[Tet[(e * 4) + l]];
				t += tl / 4;
				if (tl > tmax) tmax = tl;
			}
			tcell[e] = t;
		}
		return tmax;
	}
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
	return tmax;
}


// 導電領域の体積加重の平均温度
static double th_tmean(const double *tcell, const double *scell)
{
	double vsum = 0, tsum = 0;
	if (MeshMode) {
		for (int e = 0; e < NTet; e++) {
			if (scell[e] <= 0) continue;
			double gr[4][3], v;
			if (tet_grad_pub(&Tet[e * 4], gr, &v)) continue;
			vsum += v;
			tsum += v * tcell[e];
		}
		return ((vsum > 0) ? (tsum / vsum) : Temperature);
	}
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
	return ((vsum > 0) ? (tsum / vsum) : Temperature);
}


// 輻射 : いまの温度のまわりで Newton 線形化して、伝導 + 対流の行列 b0 に
// 接線を足したものを B に、右辺を gr に入れる
static int th_radlin(FILE *fp_log, crs_t *B, const double *b0, double *gr,
	const double *tn, const double *kcell, int it)
{
	const int64_t nnz = B->rowptr[B->n];
	for (int64_t e = 0; e < nnz; e++) B->val[e] = b0[e];
	for (int64_t i = 0; i < B->n; i++) gr[i] = 0;
	if (surf_faces(B, gr, tn, kcell, 1, NULL) < -1e299) {
		fprintf(fp_log, "*** thermal : the absolute temperature became non-positive "
			"on a radiating face (iteration %d)\n", it);
		return 1;
	}
	return 0;
}


// 時刻 t の駆動電流 = I x (waveform の折れ線の値)。waveform が無ければ一定 (t = 0 で
// 入れる段差)。折れ線の範囲の外は端の値のまま。同じ時刻が 2 回あれば段差で、
// left = 1 ならその時刻の左極限 (前の値)、0 なら右極限 (あとの値) を返す。
// 刻み [t_n, t_{n+1}] の中では f(t_n^+) と f(t_{n+1}^-) を使う。こうすると段差が
// 刻みの境目にあるとき、θ 法の台形側が段差をまたいで平均することがなく、
// Crank-Nicolson が 2 次のまま (段差をまたぐと 1 次に落ちる)
static double th_current(double t, int left)
{
	if (WaveN < 1) return ThermalCurrent;
	int k = -1;
	for (int m = 0; m < WaveN; m++) {
		if (left ? (WaveT[m] < t) : (WaveT[m] <= t)) k = m;
	}
	if (k < 0) return ThermalCurrent * WaveF[0];
	if (k == (WaveN - 1)) return ThermalCurrent * WaveF[k];
	const double w = (t - WaveT[k]) / (WaveT[k + 1] - WaveT[k]);
	return ThermalCurrent * (WaveF[k] + (w * (WaveF[k + 1] - WaveF[k])));
}


// 熱連成の過渡解析 (transient)。θ 法で
//   M (T^{n+1} - T^n)/Δt + θ R(T^{n+1}) + (1-θ) R(T^n) = 0,
//   R(T) = K T + (対流・輻射) - f(T)      (f は Joule 発熱、σ(T) を通して T に依る)
// を各時刻で解く。各時刻の中は定常と同じく σ(T) の Picard 反復 + 輻射の Newton
// 線形化で、R と T_max が動かなくなるまで繰り返す (完全陰的)。
// M は熱容量の質量行列 ∫ρc N_i N_j dV (2x2x2 Gauss で直方体に厳密)。
// 時刻ごとの R・T_max・平均温度を ofe_transient.csv に、熱の収支
//   ∫(発熱) dt = ΔU (蓄熱) + ∫(ヒートシンク + 対流 + 輻射) dt
// をログに出す。M の行和が ∫ρc N_i dV なので、これは丸めと内側の反復の
// 打ち切りを除いて厳密に成り立つ (rlc_check.sh で見る)。
static int th_transient(FILE *fp_log, crs_t *A, crs_t *B, const double *b0,
	const double *kcell, double *tcell, double *scell, unsigned char *fix,
	const unsigned char *fixt, double *phi, double *work, double *diag, double *f,
	double *tn, double *rhs, double *u, const double *gc, double *gr, double *q,
	double *gport)
{
	const int n = (int)B->n;
	const int64_t nc = th_nelem();
	const int64_t nnz = B->rowptr[n];
	const double dts = TransientTime / TransientSteps;
	const double th = TransientTheta;
	int ierr = 0;

	fprintf(fp_log, "  transient : t_end = %.6e [s], %d steps (dt = %.6e [s]), "
		"theta = %.3f\n", TransientTime, TransientSteps, dts, th);

	double *ccell  = (double *)malloc((size_t)nc * sizeof(double));
	double *told   = (double *)malloc((size_t)n * sizeof(double));
	double *t0     = (double *)malloc((size_t)n * sizeof(double));
	double *mtold  = (double *)malloc((size_t)n * sizeof(double));
	double *resold = (double *)malloc((size_t)n * sizeof(double));
	double *f1     = (double *)malloc((size_t)n * sizeof(double));
	for (int64_t c = 0; c < nc; c++) {
		ccell[c] = ((kcell[c] > 0) ? Material[th_emat(c)].cth : 0);
	}
	crs_t M, S;
	crs_alloc(&M);
	crs_alloc(&S);
	th_mass(&M, ccell);

	FILE *fcsv = fopen(FN_transient, "w");
	if (fcsv == NULL) {
		fprintf(fp_log, "*** thermal : cannot open %s\n", FN_transient);
		ierr = 1;
	}
	else {
		fprintf(fcsv, "t [s],R [ohm],V [V],T_max [degC],T_mean [degC],I [A]\n");
	}

	// 時刻 0 : 一様な初期温度 (ヒートシンクの節点だけその温度)。R(T^0) を作るため
	// 電気を 1 回解き、θ 法の (1-θ) 側に要る R(T^0) を持っておく
	double g = 0, r = 0, tmax = 0, ptot = 0, qconv = 0, qrad = 0;
	if (!ierr) {
		tmax = th_celltemp(tn, kcell, tcell);
		ierr = th_sigma(fp_log, nc, tcell, scell);
	}
	// 電気は単位電流で解き、発熱は f = I² f1 (σ を固定すれば発熱は I² に比例するので、
	// 段差の両側の発熱を解き直さずに作れる)
	double cur = th_current(0, 0), p1 = 0;
	if (!ierr) ierr = th_electric(fp_log, A, scell, fix, phi, work, diag, q, 0, 1.0, &g, f1);
	if (!ierr && (RadEps > 0)) ierr = th_radlin(fp_log, B, b0, gr, tn, kcell, 0);
	if (!ierr) {
		r = 1 / g;
		p1 = 0;
		for (int i = 0; i < n; i++) {
			f[i] = cur * cur * f1[i];
			p1 += f1[i];
		}
		ptot = cur * cur * p1;
		crs_spmv(B, tn, work, NULL);
		for (int i = 0; i < n; i++) resold[i] = work[i] - f[i] - gc[i] - gr[i];
		qconv = ((ConvH > 0) ? surf_faces(NULL, NULL, tn, kcell, 0, NULL) : 0);
		qrad = ((RadEps > 0) ? surf_faces(NULL, NULL, tn, kcell, 1, NULL) : 0);
		for (int i = 0; i < n; i++) t0[i] = tn[i];
		fprintf(fcsv, "%.8e,%.8e,%.8e,%.8e,%.8e,%.8e\n", 0.0, r, cur * r, tmax,
			th_tmean(tcell, scell), cur);
		fprintf(fp_log, "  %6s %14s %16s %14s %5s %14s\n", "step", "t [s]", "R [ohm]",
			"T_max [degC]", "iter", "I [A]");
		fprintf(fp_log, "  %6d %14.6e %16.9e %14.6f %5d %14.6e\n", 0, 0.0, r, tmax, 0, cur);
	}

	double ejoule = 0, esink = 0, econv = 0, erad = 0;
	int nittot = 0, step = 0;
	for (step = 1; !ierr && (step <= TransientSteps); step++) {
		const double pold = ptot, qcold = qconv, qrold = qrad;
		// この刻みの中の電流は終わりの時刻の左極限
		cur = th_current(step * dts, 1);
		for (int i = 0; i < n; i++) told[i] = tn[i];
		crs_spmv(&M, told, mtold, NULL);

		double rprev = 0, tmprev = 0, qsink = 0;
		int it = 0, conv = 0;
		for (it = 1; it <= ThermalMaxIter; it++) {
			const int tag = step * 1000 + it;
			if ((ierr = th_sigma(fp_log, nc, tcell, scell)) != 0) break;
			if ((ierr = th_electric(fp_log, A, scell, fix, phi, work, diag, q, tag, 1.0, &g, f1)) != 0) break;
			r = 1 / g;
			p1 = 0;
			for (int i = 0; i < n; i++) {
				f[i] = cur * cur * f1[i];
				p1 += f1[i];
			}
			ptot = cur * cur * p1;
			if ((RadEps > 0) && ((ierr = th_radlin(fp_log, B, b0, gr, tn, kcell, tag)) != 0)) break;

			// 系の行列 S = M/Δt + θ B、右辺 = M T^n/Δt + θ (f + g) - (1-θ) R(T^n)
			for (int64_t e = 0; e < nnz; e++) S.val[e] = M.val[e] / dts + th * B->val[e];
			for (int i = 0; i < n; i++) {
				if (!fixt[i]) tn[i] = 0;
			}
			crs_spmv(&S, tn, work, NULL);
			for (int i = 0; i < n; i++) {
				rhs[i] = (fixt[i] ? 0 : ((mtold[i] / dts) + (th * (f[i] + gc[i] + gr[i]))
					- ((1 - th) * resold[i]) - work[i]));
			}
			char label[64];
			sprintf(label, "T%d", tag);
			if (solver_cg(&S, rhs, u, fixt, Solver.maxiter, Solver.nout, Solver.converg,
					NULL, label) < 0) {
				fprintf(fp_log, "*** thermal : the heat-conduction solve did not converge "
					"(step %d, iteration %d)\n", step, it);
				ierr = 1;
				break;
			}
			for (int i = 0; i < n; i++) {
				if (!fixt[i]) tn[i] += u[i];
			}
			tmax = th_celltemp(tn, kcell, tcell);

			const double dr = ((rprev > 0) ? fabs(r - rprev) / r : 1);
			const double dt = fabs(tmax - tmprev) / ((fabs(tmax) > 1) ? fabs(tmax) : 1);
			if ((it > 1) && (dr < ThermalTol) && (dt < ThermalTol)) {
				conv = 1;
				break;
			}
			if (tmax > 1e5) break;
			rprev = r;
			tmprev = tmax;
		}
		if (ierr) break;
		if (tmax > 1e5) {
			fprintf(fp_log, "*** thermal runaway : T_max exceeded 1e5 degC at t = %.6e [s] "
				"(step %d)\n", step * dts, step);
			ierr = 1;
			break;
		}
		if (!conv) {
			fprintf(fp_log, "*** thermal : the coupling did not converge in %d iterations "
				"at step %d\n", ThermalMaxIter, step);
			ierr = 1;
			break;
		}
		nittot += it;

		// この刻みの熱の収支。ヒートシンクの節点では刻みの式の残差が
		// (符号を変えて) ヒートシンクへ出ていく熱の流れになる
		qsink = 0;
		for (int i = 0; i < n; i++) {
			const int id = NodeConductor[i];
			if ((id >= 0) && HeatSinkOn[id]) {
				const double dm = crs_row_dot(&M, i, tn) - mtold[i];
				qsink -= (dm / dts) + (th * (crs_row_dot(B, i, tn) - f[i] - gc[i] - gr[i]))
					+ ((1 - th) * resold[i]);
			}
		}
		crs_spmv(B, tn, work, NULL);
		for (int i = 0; i < n; i++) resold[i] = work[i] - f[i] - gc[i] - gr[i];
		qconv = ((ConvH > 0) ? surf_faces(NULL, NULL, tn, kcell, 0, NULL) : 0);
		qrad = ((RadEps > 0) ? surf_faces(NULL, NULL, tn, kcell, 1, NULL) : 0);
		ejoule += dts * ((th * ptot) + ((1 - th) * pold));
		econv += dts * ((th * qconv) + ((1 - th) * qcold));
		erad += dts * ((th * qrad) + ((1 - th) * qrold));
		esink += dts * qsink;

		// 段差 : 次の刻みの (1-θ) 側は右極限の電流での R(T^{n+1})。発熱だけが変わる
		const double cr = th_current(step * dts, 0);
		if (cr != cur) {
			for (int i = 0; i < n; i++) {
				resold[i] -= ((cr * cr) - (cur * cur)) * f1[i];
				f[i] = cr * cr * f1[i];
			}
			ptot = cr * cr * p1;
			cur = cr;
		}

		const double tmean = th_tmean(tcell, scell);
		fprintf(fcsv, "%.8e,%.8e,%.8e,%.8e,%.8e,%.8e\n", step * dts, r, cur * r,
			tmax, tmean, cur);
		fprintf(fp_log, "  %6d %14.6e %16.9e %14.6f %5d %14.6e\n", step, step * dts, r, tmax,
			it, cur);
		fflush(fp_log);
	}

	if (!ierr) {
		// 蓄熱 ΔU = Σ_i (M (T - T^0))_i = ∫ρc (T - T^0) dV
		for (int i = 0; i < n; i++) work[i] = tn[i] - t0[i];
		crs_spmv(&M, work, rhs, NULL);
		double du = 0;
		for (int i = 0; i < n; i++) du += rhs[i];
		fprintf(fp_log, "  %d steps, %d iterations in total\n", TransientSteps, nittot);
		fprintf(fp_log, "  R = %.10e [ohm] at I = %.6e [A] (V = %.10e [V]), t = %.6e [s]\n",
			r, cur, cur * r, TransientTime);
		fprintf(fp_log, "  T_max = %.10e [degC], T_mean (conductor) = %.10e [degC]\n",
			tmax, th_tmean(tcell, scell));
		if (th_ac()) {
			fprintf(fp_log, "  X = %.10e [ohm], L = %.10e [H] (I is the peak amplitude; "
				"the joule power below is the time average I^2 R / 2)\n", ThX, ThX / (2 * PI * Freq));
		}
		fprintf(fp_log, "  joule power : I^2 R = %.12e, field integral = %.12e [W] (t = t_end)\n",
			th_pfac() * cur * cur * r, ptot);
		fprintf(fp_log, "  heat balance : joule = %.12e, stored = %.12e, to heat sinks = %.12e, "
			"to convection = %.12e, to radiation = %.12e [J]\n", ejoule, du, esink, econv, erad);
		*gport = g;

		if (FieldOut) {
			if (!th_ac()) {
				for (int i = 0; i < n; i++) phi[i] *= cur;	// 単位電流の φ を終了時刻の電流に
				field_add_node("phi_R_port1", phi);
				field_add_grad("E_R_port1", phi, 0);
			}
			field_add_node("T", tn);
		}
	}
	fflush(fp_log);

	if (fcsv) fclose(fcsv);
	crs_free(&M);
	crs_free(&S);
	free(ccell); free(told); free(t0); free(mtold); free(resold); free(f1);

	return ierr;
}


// 熱連成の定常解析。σ(T) の Picard 反復 (+ 輻射の Newton 線形化) を R と T_max が
// 動かなくなるまで回す
static int th_steady(FILE *fp_log, crs_t *A, crs_t *B, const double *b0,
	const double *kcell, double *tcell, double *scell, unsigned char *fix,
	const unsigned char *fixt, double *phi, double *work, double *diag, double *f,
	double *tn, double *rhs, double *u, const double *gc, double *gr, double *q,
	double *gport)
{
	const int n = (int)B->n;
	const int64_t nc = th_nelem();
	int ierr = 0;

	double rprev = 0, tmprev = 0, r = 0, g = 0, tmax = 0;
	double ptot = 0, qsink = 0, qconv = 0, qrad = 0, pterm = 0;
	int it = 0, conv = 0;
	fprintf(fp_log, "  %4s %16s %14s %16s\n", "iter", "R [ohm]", "T_max [degC]", "dR / R");
	for (it = 1; !ierr && (it <= ThermalMaxIter); it++) {
		if ((ierr = th_sigma(fp_log, nc, tcell, scell)) != 0) break;
		if ((ierr = th_electric(fp_log, A, scell, fix, phi, work, diag, q, it, ThermalCurrent, &g, f)) != 0) break;
		r = 1 / g;
		pterm = th_pfac() * ThermalCurrent * ThermalCurrent * r;

		// 発熱と熱 : T = T_D + u、自由節点で B u = f - B T_D (f は th_electric で作った)
		ptot = 0;
		for (int i = 0; i < n; i++) ptot += f[i];
		// 輻射 : いまの温度のまわりで Newton 線形化して行列と右辺に入れる。
		// σ(T) の Picard 反復と同じ回で回すので、収束すれば両方の不動点になる
		if ((RadEps > 0) && ((ierr = th_radlin(fp_log, B, b0, gr, tn, kcell, it)) != 0)) break;
		for (int i = 0; i < n; i++) {
			if (!fixt[i]) tn[i] = 0;
		}
		crs_spmv(B, tn, work, NULL);
		for (int i = 0; i < n; i++) rhs[i] = (fixt[i] ? 0 : (f[i] + gc[i] + gr[i] - work[i]));
		char label[64];
		sprintf(label, "T%d", it);
		if (solver_cg(B, rhs, u, fixt, Solver.maxiter, Solver.nout, Solver.converg,
				NULL, label) < 0) {
			fprintf(fp_log, "*** thermal : the heat-conduction solve did not converge "
				"(iteration %d)\n", it);
			ierr = 1;
			break;
		}
		for (int i = 0; i < n; i++) {
			if (!fixt[i]) tn[i] += u[i];
		}

		// ヒートシンク・対流・輻射で出ていく熱。K_T の列和が 0 なので
		// 和は Σ f に一致するはず (対流の行列と右辺の和は対流の熱そのもの)。
		// 輻射は線形化した式で解いているので、出ていく熱を新しい温度で
		// 評価し直した値とは O(ΔT²) だけ違う (収束すれば丸めの程度)
		qsink = 0;
		for (int i = 0; i < n; i++) {
			const int id = NodeConductor[i];
			if ((id >= 0) && HeatSinkOn[id]) {
				qsink += f[i] + gc[i] + gr[i] - crs_row_dot(B, i, tn);
			}
		}
		qconv = ((ConvH > 0) ? surf_faces(NULL, NULL, tn, kcell, 0, NULL) : 0);
		qrad = ((RadEps > 0) ? surf_faces(NULL, NULL, tn, kcell, 1, NULL) : 0);

		tmax = th_celltemp(tn, kcell, tcell);

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
		fprintf(fp_log, "  converged in %d iterations\n", it);
		fprintf(fp_log, "  R = %.10e [ohm] at I = %.6e [A] (V = %.10e [V])\n",
			r, ThermalCurrent, ThermalCurrent * r);
		fprintf(fp_log, "  T_max = %.10e [degC], T_mean (conductor) = %.10e [degC]\n",
			tmax, th_tmean(tcell, scell));
		if (th_ac()) {
			fprintf(fp_log, "  X = %.10e [ohm], L = %.10e [H] (I is the peak amplitude; "
				"the joule power below is the time average I^2 R / 2)\n", ThX, ThX / (2 * PI * Freq));
		}
		fprintf(fp_log, "  joule power : I^2 R = %.12e, field integral = %.12e, "
			"to heat sinks = %.12e, to convection = %.12e, to radiation = %.12e [W]\n",
			pterm, ptot, qsink, qconv, qrad);
		*gport = g;

		if (FieldOut) {
			if (!th_ac()) {
				field_add_node("phi_R_port1", phi);
				field_add_grad("E_R_port1", phi, 0);
			}
			field_add_node("T", tn);
		}
	}


	return ierr;
}


int solve_thermal(FILE *fp_log, double *gport, double *xport)
{
	const int n = (int)num_node();
	const int64_t nc = th_nelem();
	int ierr = 0;

	fprintf(fp_log, "\n=== electro-thermal analysis ===\n");
	// 非構造格子は 1 次四面体だけ (発熱・質量・露出面を四面体でしか組んでいない)。
	// 他の種別は黙って誤答にせず入力エラーにする
	if (MeshMode && ((NTet < 1) || (NHex > 0) || (NPrism > 0) || (NPyr > 0) || (TetOrder >= 2))) {
		fprintf(fp_log, "*** thermal on an unstructured mesh supports linear tetrahedra "
			"only (this mesh has %d tetrahedra (order %d), %d hexahedra, %d prisms, "
			"%d pyramids)\n", NTet, TetOrder, NHex, NPrism, NPyr);
		return 1;
	}
	th_extent();
	fprintf(fp_log, "  drive current = %.6e [A]%s, ambient (initial) = %.4f [degC]\n",
		ThermalCurrent, (th_ac() ? " (peak amplitude, AC)" : ""), Temperature);
	ThX = 0;
	for (int p = 0; p <= NPort; p++) {
		if (HeatSinkOn[p]) {
			fprintf(fp_log, "  heat sink : conductor %d at %.4f [degC]\n", p, HeatSinkT[p]);
		}
	}
	for (int e = 0; e < (NConv + NRad); e++) {
		const surfbc_t *bc = ((e < NConv) ? &ConvBC[e] : &RadBC[e - NConv]);
		if (e < NConv) {
			fprintf(fp_log, "  convection : h = %.6e [W/(m^2 K)], T_amb = %.4f [degC]",
				bc->c, bc->t);
		}
		else {
			fprintf(fp_log, "  radiation : emissivity = %.6f, T_amb = %.4f [degC]", bc->c, bc->t);
		}
		if (bc->box) {
			fprintf(fp_log, " on faces in [%g, %g] x [%g, %g] x [%g, %g]",
				bc->g[0], bc->g[1], bc->g[2], bc->g[3], bc->g[4], bc->g[5]);
		}
		fprintf(fp_log, "\n");
	}
	fflush(fp_log);

	double *tcell = (double *)malloc((size_t)nc * sizeof(double));
	double *scell = (double *)malloc((size_t)nc * sizeof(double));
	double *kcell = (double *)malloc((size_t)nc * sizeof(double));
	for (int64_t c = 0; c < nc; c++) {
		tcell[c] = Temperature;
		kcell[c] = Material[th_emat(c)].kth;
	}
	TFace = NULL;
	TFaceN = 0;
	if (MeshMode) tet_faces_build(kcell);

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
	double *gr   = (double *)malloc((size_t)n * sizeof(double));
	double *q    = (double *)malloc(((size_t)NPort + 1) * sizeof(double));

	// 熱の行列は温度に依らないので 1 回だけ作る
	crs_t A, B;
	crs_alloc(&A);
	crs_alloc(&B);
	th_assemble(&B, kcell);
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
			tn[i] = Temperature;			// 輻射の最初の線形化の点・過渡解析の初期値
		}
	}
	// 対流は熱伝導の行列に足す (温度に依らない)。節点の固定は熱伝導だけで決める
	for (int i = 0; i < n; i++) {
		gc[i] = 0;
		gr[i] = 0;
		u[i] = 0;
	}
	int64_t nface = 0;
	if (ConvH > 0) {
		surf_faces(&B, gc, NULL, kcell, 0, NULL);
	}
	if ((NConv + NRad) > 0) {
		// 係数が正の露出面の数。直方体を付けた行は、露出面に 1 つも当たらなければ
		// 書き間違い (座標の桁違いなど) なので黙って無視せず入力エラーにする
		int *hits = (int *)calloc((size_t)(NConv + NRad), sizeof(int));
		nface = (int64_t)surf_faces(NULL, NULL, NULL, kcell, 2, hits);
		for (int e = 0; e < (NConv + NRad); e++) {
			const surfbc_t *bc = ((e < NConv) ? &ConvBC[e] : &RadBC[e - NConv]);
			if (bc->box && (hits[e] < 1)) {
				fprintf(fp_log, "*** thermal : the %s box [%g, %g] x [%g, %g] x [%g, %g] "
					"contains no exposed face\n", ((e < NConv) ? "convection" : "radiation"),
					bc->g[0], bc->g[1], bc->g[2], bc->g[3], bc->g[4], bc->g[5]);
				ierr = 1;
			}
		}
		free(hits);
		fprintf(fp_log, "  %s : %lld exposed faces\n",
			((NConv > 0) ? ((NRad > 0) ? "convection / radiation" : "convection") : "radiation"),
			(long long)nface);
	}
	if (!TransientOn && (nsink < 1) && (nface < 1)) {
		fprintf(fp_log, "*** thermal : the heat has no way out (the heat-sink conductors "
			"have no node and no face is exposed to convection / radiation)\n");
		ierr = 1;
	}
	// 輻射は温度に依るので、反復ごとに伝導 + 対流の行列 (b0) に接線を足し直す
	const int64_t nnz = B.rowptr[n];
	double *b0 = NULL;
	if (RadEps > 0) {
		b0 = (double *)malloc((size_t)nnz * sizeof(double));
		for (int64_t e = 0; e < nnz; e++) b0[e] = B.val[e];
	}

	// 交流 : A-φ の辺の情報と行列は反復の間持っておく
	int ehopen = 0;
	if (!ierr && th_ac()) {
		ierr = eddy3d_heat_open(fp_log);
		ehopen = 1;
	}
	if (!ierr) {
		ierr = (TransientOn ? th_transient : th_steady)(fp_log, &A, &B, b0, kcell, tcell,
			scell, fix, fixt, phi, work, diag, f, tn, rhs, u, gc, gr, q, gport);
	}
	if (ehopen) eddy3d_heat_close();
	*xport = ThX;
	fflush(fp_log);

	crs_free(&A);
	crs_free(&B);
	free(tcell); free(scell); free(kcell);
	free(fix); free(fixt); free(phi); free(work); free(diag);
	free(f); free(tn); free(rhs); free(u); free(gc); free(gr); free(q);
	free(b0);
	free(TFace);
	TFace = NULL;

	return ierr;
}
