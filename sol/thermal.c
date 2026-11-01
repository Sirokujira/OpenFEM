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
	if (MeshMode && (MeshDim == 2)) return (int64_t)NTri;
	return (MeshMode ? (int64_t)elem3d_count() : ((int64_t)Nx * Ny * Nz));
}


static int th_emat(int64_t c)
{
	if (MeshMode && (MeshDim == 2)) return (int)TriMat[c];
	return (MeshMode ? elem3d_mat((int)c) : (int)CellMaterial[c]);
}


// 対流の h の温度依存 (convectiontemptable)。h(T) = h f(T) の f と df/dT (折れ線、
// 範囲の外は端の値で傾き 0)。表が無ければ対流は線形 (従来の経路)
static int th_convnl(void)
{
	return ((ConvH > 0) && (ConvTabN > 0));
}


static int th_nonlin(void)
{
	return ((RadEps > 0) || th_convnl());
}


static double tab_eval(int nt, const double *tt, const double *ff, double t, double *slope)
{
	*slope = 0;
	if (t <= tt[0]) return ff[0];
	for (int k = 0; k < (nt - 1); k++) {
		if (t <= tt[k + 1]) {
			*slope = (ff[k + 1] - ff[k]) / (tt[k + 1] - tt[k]);
			return ff[k] + ((t - tt[k]) * (*slope));
		}
	}
	return ff[nt - 1];
}


static double conv_tab(double t, double *slope)
{
	return tab_eval(ConvTabN, ConvTabT, ConvTabF, t, slope);
}


// 面の熱流 q(T) とその接線 (rad = 1 : 輻射、rad = 3 : h(T) の対流)。
// 絶対温度が正でなければ *bad を立てる (輻射のみ)
static void surf_q(int rad, double hc, double ta, double tg, double *qv, double *dq, int *bad)
{
	if (rad == 1) {
		const double tk = tg + KELVIN, ta4 = pow(ta + KELVIN, 4);
		if (tk <= 0) *bad = 1;
		if (RadTabN > 0) {
			// ε(T) = ε f(T) : q' = ε σ (4 f T³ + f' (T⁴ - T_r⁴))
			double fp;
			const double fv = tab_eval(RadTabN, RadTabT, RadTabF, tg, &fp);
			*qv = hc * fv * STEFAN * (tk * tk * tk * tk - ta4);
			*dq = hc * STEFAN * ((4 * fv * tk * tk * tk) + (fp * (tk * tk * tk * tk - ta4)));
			return;
		}
		*qv = hc * STEFAN * (tk * tk * tk * tk - ta4);
		*dq = 4 * hc * STEFAN * tk * tk * tk;
		return;
	}
	double fp;
	const double fv = conv_tab(tg, &fp);
	*qv = hc * fv * (tg - ta);
	*dq = hc * (fv + (fp * (tg - ta)));
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


// 非構造格子の露出面 [4*TFaceN] (熱を通す要素の境界の三角形・四角形、電極の面を除く)。
// 三角形は昇順の 3 節点 + -1 (四面体だけの格子では従来と同じ並び)、四角形は要素の局所の
// 巡回順 (双 1 次の写像に使う)
static int32_t *TFace;
static int TFaceN;
static int32_t *TFaceEl;		// [TFaceN] 露出面を持つ要素 (外向きの法線を決める)
static const signed char *ThNodeCond;	// ヒートシンクを選ぶ節点の導体番号 (3 次元は NodeConductor、
										// 断面 2 次元は導体の三角形の節点)
static int32_t *TOcc;			// [4*TOccN] 熱を通す要素の境界面すべて (電極の面も含む。相互放射の遮蔽)
static int32_t *TOccEl;			// [TOccN] その面を持つ要素 (高次の四面体の曲がった遮蔽物に使う)
static int TOccN;

// 面の記録 : 昇順の鍵 4 個 (三角形は -1 埋め) + 巡回順の節点 4 個
static int tface_cmp(const void *a, const void *b)
{
	const int32_t *p = (const int32_t *)a, *q = (const int32_t *)b;
	for (int d = 0; d < 4; d++) {
		if (p[d] != q[d]) return ((p[d] < q[d]) ? -1 : 1);
	}
	return 0;
}

// 要素の局所の面 (Gmsh の節点の並び)。nv = 3 / 4
static int elem_faces(int e, int fv[6][4], int nfv[6])
{
	static const int tf[4][3] = {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}};
	static const int hf[6][4] = {{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4},
	                             {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}};
	static const int pf[5][4] = {{0, 2, 1, -1}, {3, 4, 5, -1}, {0, 1, 4, 3},
	                             {1, 2, 5, 4}, {2, 0, 3, 5}};
	static const int yf[5][4] = {{0, 3, 2, 1}, {0, 1, 4, -1}, {1, 2, 4, -1},
	                             {2, 3, 4, -1}, {3, 0, 4, -1}};
	const int kind = elem3d_kind(e);
	int nf = 0;
	if (kind == MESHELEM_TET) {
		nf = 4;
		for (int f = 0; f < 4; f++) {
			nfv[f] = 3;
			for (int l = 0; l < 3; l++) fv[f][l] = tf[f][l];
		}
	}
	else if (kind == MESHELEM_HEX) {
		nf = 6;
		for (int f = 0; f < 6; f++) {
			nfv[f] = 4;
			for (int l = 0; l < 4; l++) fv[f][l] = hf[f][l];
		}
	}
	else {
		nf = 5;
		for (int f = 0; f < 5; f++) {
			const int *src = ((kind == MESHELEM_PRISM) ? pf[f] : yf[f]);
			nfv[f] = ((src[3] < 0) ? 3 : 4);
			for (int l = 0; l < 4; l++) fv[f][l] = src[l];
		}
	}
	return nf;
}

static int edge_cmp(const void *a, const void *b)
{
	const int32_t *p = (const int32_t *)a, *q = (const int32_t *)b;
	for (int d = 0; d < 2; d++) {
		if (p[d] != q[d]) return ((p[d] < q[d]) ? -1 : 1);
	}
	return 0;
}


// 断面 2 次元 (1 次三角形) の露出辺 : 熱を通す三角形 1 つにしか属さない辺 (k = 0 の領域との
// 界面を含む)。TFace の形 {a, b, -1, -1}。単位長あたりの線積分になる
static void tri_edges_build(const double *kcell)
{
	int ne = 0;
	for (int e = 0; e < NTri; e++) ne += ((kcell[e] > 0) ? 3 : 0);
	int32_t *all = (int32_t *)malloc(((size_t)ne * 5 + 1) * sizeof(int32_t));
	int m = 0;
	for (int e = 0; e < NTri; e++) {
		if (kcell[e] <= 0) continue;
		for (int l = 0; l < 3; l++) {
			const int32_t a = Tri[(e * 3) + l], b = Tri[(e * 3) + ((l + 1) % 3)];
			int32_t *r = &all[m * 5];
			r[0] = ((a < b) ? a : b);
			r[1] = ((a < b) ? b : a);
			r[2] = a;
			r[3] = b;
			r[4] = e;
			m++;
		}
	}
	qsort(all, (size_t)ne, 5 * sizeof(int32_t), edge_cmp);
	TFace = (int32_t *)malloc(((size_t)ne * 4 + 1) * sizeof(int32_t));
	TFaceEl = (int32_t *)malloc(((size_t)ne + 1) * sizeof(int32_t));
	TOcc = (int32_t *)malloc(sizeof(int32_t));
	TOccEl = (int32_t *)malloc(sizeof(int32_t));
	TFaceN = 0;
	TOccN = 0;
	for (int i = 0; i < ne; ) {
		int j = i + 1;
		while ((j < ne) && !edge_cmp(&all[i * 5], &all[j * 5])) j++;
		if ((j - i) == 1) {
			TFace[(TFaceN * 4) + 0] = all[(i * 5) + 2];
			TFace[(TFaceN * 4) + 1] = all[(i * 5) + 3];
			TFace[(TFaceN * 4) + 2] = -1;
			TFace[(TFaceN * 4) + 3] = -1;
			TFaceEl[TFaceN] = all[(i * 5) + 4];
			TFaceN++;
		}
		i = j;
	}
	free(all);
}


static void tet_faces_build(const double *kcell)
{
	if (MeshDim == 2) {
		tri_edges_build(kcell);
		return;
	}
	const int ne = elem3d_count();
	int nf = 0;
	for (int e = 0; e < ne; e++) {
		if (kcell[e] <= 0) continue;
		int fv[6][4], nfv[6];
		nf += elem_faces(e, fv, nfv);
	}
	int32_t *all = (int32_t *)malloc(((size_t)nf * 9 + 1) * sizeof(int32_t));
	int m = 0;
	for (int e = 0; e < ne; e++) {
		if (kcell[e] <= 0) continue;
		int32_t nd[ELEM_NMAX];
		elem3d_nodes(e, nd);
		int fv[6][4], nfv[6];
		const int nfe = elem_faces(e, fv, nfv);
		for (int f = 0; f < nfe; f++) {
			int32_t v[4] = {-1, -1, -1, -1};
			for (int l = 0; l < nfv[f]; l++) v[l] = nd[fv[f][l]];
			int32_t *r = &all[m * 9];
			for (int l = 0; l < 4; l++) r[4 + l] = v[l];
			r[8] = e;
			// 昇順に並べて面の鍵にする (三角形は 3 個、4 個目は -1)
			const int k = nfv[f];
			for (int a = 0; a < (k - 1); a++) {
				for (int b = a + 1; b < k; b++) {
					if (v[b] < v[a]) { const int32_t w = v[a]; v[a] = v[b]; v[b] = w; }
				}
			}
			for (int l = 0; l < 4; l++) r[l] = v[l];
			m++;
		}
	}
	qsort(all, (size_t)nf, 9 * sizeof(int32_t), tface_cmp);
	TFace = (int32_t *)malloc(((size_t)nf * 4 + 1) * sizeof(int32_t));
	TFaceEl = (int32_t *)malloc(((size_t)nf + 1) * sizeof(int32_t));
	TOcc = (int32_t *)malloc(((size_t)nf * 4 + 1) * sizeof(int32_t));
	TOccEl = (int32_t *)malloc(((size_t)nf + 1) * sizeof(int32_t));
	TFaceN = 0;
	TOccN = 0;
	for (int i = 0; i < nf; ) {
		int j = i + 1;
		while ((j < nf) && !tface_cmp(&all[i * 9], &all[j * 9])) j++;
		if ((j - i) == 1) {
			const int32_t *v = &all[i * 9];
			for (int l = 0; l < 4; l++) TOcc[(TOccN * 4) + l] = v[4 + l];
			TOccEl[TOccN] = v[8];
			TOccN++;
			const int k = ((v[3] < 0) ? 3 : 4);
			int ncond = 0;
			for (int l = 0; l < k; l++) ncond += (NodeConductor[v[l]] >= 0);
			if (ncond < k) {
				// 三角形は昇順 (従来の並び)、四角形は巡回順
				for (int l = 0; l < 4; l++) {
					TFace[(TFaceN * 4) + l] = ((k == 3) ? v[l] : v[4 + l]);
				}
				TFaceEl[TFaceN] = v[8];
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
// 位置の折れ線 (convectionprofile / radiationprofile) があれば係数に面の中心での倍率を掛けた
// 写しを返す (対流と輻射で別の置き場。呼び出し側は対流と輻射の行を同時に持つことがある)
static const surfbc_t *surf_bc(int rad, double xc, double yc, double zc, int *idx)
{
	static surfbc_t prof[2];
	const int k = ((rad == 1) ? 1 : 0);
	const surfbc_t *list = (k ? RadBC : ConvBC);
	const int nb = (k ? NRad : NConv);
	for (int e = nb - 1; e >= 0; e--) {
		if (surf_in(&list[e], xc, yc, zc)) {
			if (idx) *idx = e;
			if ((ProfN[k] < 1) && !ExprOn[k]) return &list[e];
			const double r[3] = {xc, yc, zc};
			double sl;
			prof[k] = list[e];
			if (!ExprOn[k]) {
				prof[k].c = list[e].c * tab_eval(ProfN[k], ProfX[k], ProfF[k], r[ProfAxis[k]], &sl);
				return &prof[k];
			}
			// 位置の式 (convectionexpr / radiationexpr)。折れ線もあれば両方の積。負・非有限
			// (輻射は ε f > 1 も) は 0 にして数え、solve_thermal が入力エラーにする
			double f = expr_eval(&ExprF[k], xc, yc, zc);
			if (ProfN[k] > 0) f = tab_eval(ProfN[k], ProfX[k], ProfF[k], r[ProfAxis[k]], &sl) * f;
			if (!((f >= 0) && (f <= 1e300)) || (k && ((list[e].c * f) > 1))) {
				ExprBad[k]++;
				f = 0;
			}
			prof[k].c = list[e].c * f;
			return &prof[k];
		}
	}
	return NULL;
}


// 面どうしの相互放射 (radiationexchange、構造格子)。輻射を課した露出面を 1 枚ずつ
// 放射度の要素とし (灰色・拡散)、ε を Q_ii に、周囲 (T_r) との交換を面ごとの一様な
// 項 c_i に置き換える (rx_build / rx_update の注釈)
static struct {
	int on, n;
	int *d, *sg, *row;				// 面の法線の軸、向き (±1)、輻射の行
	double *lo, *hi;				// [3n] 面の範囲 (軸 d では lo = hi)
	double *area, *eps, *tr;
	double *qii, *c, *e;			// Q の対角、一様な項、面の平均の σT⁴
	double *q, *q0, *finf;			// Q [n*n]、q0、F_i∞
	// 非構造格子 : 露出面 (TFace) -> 放射度の要素の番号 (-1 = 輻射なし)、面の頂点
	// (外向きの法線に対して反時計回り、三角形は 3 個)、外向きの単位法線、中心
	// 高次の四面体の曲がった面は、辺を各 2p 分割の折れ線で近似した境界 (周回積分は境界だけで決まる)
	// 高次の四面体では遮蔽物も曲がった面 (各面を (2p)² 枚の三角形に割ったもの、ocv [9 nocc]) で、
	// 線分は曲がった面の中心どうし。自分と相手の面の三角形は除く (okey は面の頂点の昇順の 3 つ組)
	int *fid, *nv;
	double *pv, *nrm, *ctr;			// [3 RX_MAXV n] / [3n] / [3n]
	int nocc;
	double *ocv;
	int32_t *okey, *fkey;			// [4 nocc] / [4n]
	// 高次の要素 : 面の積分点 (3 x 3) の位置と単位法線 [6 RX_NS n]。曲がった面どうしが向き合うかを
	// 中心だけで決めないため (rx_vis_curved)
	double *spx;
} RX;
#define RX_MAXV (48)				// 面の境界の点の最大数 (7 次の三角形で 3 x 14 = 42、4 次の四角形で 4 x 8 = 32)
#define RX_NS (9)


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
			int ridx = 0;
			const surfbc_t *bc = surf_bc(((rad >= 4) ? 1 : rad), xc, yc, zc, &ridx);
			if ((bc == NULL) || (bc->c <= 0)) continue;
			nface++;
			const double hc = bc->c, ta = bc->t;
			const int fid = (int)(nface - 1);	// 輻射の面の通し番号 (rad = 1/4/5 で同じ順)

			const double area = h[p] * h[q];
			if (rad == 4) {
				// 相互放射の面を集める (rx_build)
				if (fid < RX.n) {
					RX.d[fid] = d;
					RX.sg[fid] = (s ? 1 : -1);
					RX.row[fid] = ridx;
					RX.area[fid] = area;
					RX.eps[fid] = hc;
					RX.tr[fid] = ta;
					for (int a = 0; a < 3; a++) {
						const double *xx = ((a == 0) ? Xn : (a == 1) ? Yn : Zn);
						double lo = 1e300, hi = -1e300;
						for (int m = 0; m < 4; m++) {
							const int ix = ((a == 0) ? li[m] : (a == 1) ? lj[m] : lk[m]);
							if (xx[ix] < lo) lo = xx[ix];
							if (xx[ix] > hi) hi = xx[ix];
						}
						RX.lo[(fid * 3) + a] = lo;
						RX.hi[(fid * 3) + a] = hi;
					}
				}
				continue;
			}
			if (rad == 5) {
				// 面の平均の σT_K⁴ (面の 2x2 Gauss、熱流の積分と同じ点)
				double es = 0;
				for (int ga = 0; ga < 2; ga++) {
				for (int gb = 0; gb < 2; gb++) {
					double tg = 0;
					for (int m = 0; m < 4; m++) {
						tg += (((m >> 1) & 1) ? gp[ga] : (1 - gp[ga]))
						    * ((m & 1) ? gp[gb] : (1 - gp[gb])) * t[nd[m]];
					}
					const double tk = tg + KELVIN;
					es += STEFAN * tk * tk * tk * tk / 4;
				}
				}
				RX.e[fid] = es;
				continue;
			}
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
					const double w = area / 4;
					double qv, dq;
					if ((rad == 1) && RX.on) {
						// 相互放射 : q = Q_ii σT⁴ + c_i (c_i は他の面と周囲からの項、反復ごとに更新)
						const double tk = tg + KELVIN;
						if (tk <= 0) bad = 1;
						qv = (RX.qii[fid] * STEFAN * tk * tk * tk * tk) + RX.c[fid];
						dq = 4 * RX.qii[fid] * STEFAN * tk * tk * tk;
					}
					else {
						surf_q(rad, hc, ta, tg, &qv, &dq, &bad);
					}
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
	return ((t || (rad == 1) || (rad == 3)) ? qsum : (double)nface);
}


// 四角形の露出面 1 枚 (節点は巡回順 nd[0..3]) : 双 1 次の写像 x(ξ, η) = Σ N_a x_a
// ((ξ, η) ∈ [0,1]²) の 2x2 Gauss で、対流の面の質量・右辺・熱、輻射・h(T) の Newton
// 線形化を組む (積分点は熱と行列と右辺で同じ)。戻り値は出ていく熱 (t があるとき)
static double surf_quad(crs_t *B, double *g, const double *t, int rad, const int32_t *nd,
	int fid, int *bad, int64_t *nface)
{
	const double xc = (Xp[nd[0]] + Xp[nd[1]] + Xp[nd[2]] + Xp[nd[3]]) / 4;
	const double yc = (Yp[nd[0]] + Yp[nd[1]] + Yp[nd[2]] + Yp[nd[3]]) / 4;
	const double zc = (Zp[nd[0]] + Zp[nd[1]] + Zp[nd[2]] + Zp[nd[3]]) / 4;
	const surfbc_t *bc = surf_bc(((rad == 3) ? 0 : rad), xc, yc, zc, NULL);
	if ((bc == NULL) || (bc->c <= 0)) return 0;
	(*nface)++;
	const double hc = bc->c, ta = bc->t;
	const double gp[2] = {0.5 - 0.5 / sqrt(3.0), 0.5 + 0.5 / sqrt(3.0)};
	double me[4][4], ge[4] = {0, 0, 0, 0}, qsum = 0;
	for (int m = 0; m < 4; m++) {
		for (int mm = 0; mm < 4; mm++) me[m][mm] = 0;
	}
	for (int ga = 0; ga < 2; ga++) {
	for (int gb = 0; gb < 2; gb++) {
		const double x = gp[ga], y = gp[gb];
		const double nn[4] = {(1 - x) * (1 - y), x * (1 - y), x * y, (1 - x) * y};
		const double dx[4] = {-(1 - y), (1 - y), y, -y};
		const double dy[4] = {-(1 - x), -x, x, (1 - x)};
		double a[3] = {0, 0, 0}, b[3] = {0, 0, 0};
		for (int m = 0; m < 4; m++) {
			const double p[3] = {Xp[nd[m]], Yp[nd[m]], Zp[nd[m]]};
			for (int c = 0; c < 3; c++) {
				a[c] += dx[m] * p[c];
				b[c] += dy[m] * p[c];
			}
		}
		const double cx = (a[1] * b[2]) - (a[2] * b[1]);
		const double cy = (a[2] * b[0]) - (a[0] * b[2]);
		const double cz = (a[0] * b[1]) - (a[1] * b[0]);
		const double w = sqrt((cx * cx) + (cy * cy) + (cz * cz)) / 4;	// 重み 1/4 x |J|
		if (!rad) {
			for (int m = 0; m < 4; m++) {
				for (int mm = 0; mm < 4; mm++) me[m][mm] += hc * w * nn[m] * nn[mm];
				ge[m] += hc * ta * w * nn[m];
			}
			if (t) {
				double tg = 0;
				for (int m = 0; m < 4; m++) tg += nn[m] * t[nd[m]];
				qsum += hc * w * (tg - ta);
			}
		}
		else {
			double tg = 0;
			for (int m = 0; m < 4; m++) tg += nn[m] * t[nd[m]];
			double qv, dq;
			if ((rad == 1) && RX.on && (fid >= 0)) {
				const double tk = tg + KELVIN;
				if (tk <= 0) *bad = 1;
				qv = (RX.qii[fid] * STEFAN * tk * tk * tk * tk) + RX.c[fid];
				dq = 4 * RX.qii[fid] * STEFAN * tk * tk * tk;
			}
			else {
				surf_q(rad, hc, ta, tg, &qv, &dq, bad);
			}
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
			for (int mm = 0; mm < 4; mm++) {
				const int64_t pos = crs_find_tri(B, nd[m], nd[mm]);
				if (pos >= 0) B->val[pos] += me[m][mm];
			}
		}
	}
	if (g) {
		for (int m = 0; m < 4; m++) g[nd[m]] += ge[m];
	}
	return qsum;
}


// 高次の四面体の露出面 f (曲がった三角形) : 親の四面体の形状関数を面に制限して
// (反対側の頂点の λ = 0) 面の Duffy 積分 (各方向 p + 2 点) で対流・輻射・h(T) を組む。
// 面積要素は等パラメトリック写像の |∂x/∂μ1 × ∂x/∂μ2|。接線・右辺・熱の評価点は同じ
static double surf_tri_ho(crs_t *B, double *g, const double *t, int rad, int f,
	double hc, double ta, int *bad)
{
	const int e = TFaceEl[f];
	const int32_t *cn = &TFace[f * 4];
	int32_t nd[TET_NMAX];
	const int nen = tet_nodes(e, nd);
	int fa[3];
	for (int k = 0; k < 3; k++) {
		fa[k] = 0;
		for (int l = 0; l < 4; l++) {
			if (nd[l] == cn[k]) fa[k] = l;
		}
	}
	double mu[QTRI_MAX][3], wq[QTRI_MAX];		// 7 次で 9 x 9 = 81 点
	const int nq = tri_quad_pub(TetOrder + 2, mu, wq);
	static double me[TET_NMAX][TET_NMAX];
	double ge[TET_NMAX], qsum = 0;
	for (int a = 0; a < nen; a++) {
		ge[a] = 0;
		for (int b = 0; b < nen; b++) me[a][b] = 0;
	}
	for (int q = 0; q < nq; q++) {
		double nv[TET_NMAX], t1[3], t2[3];
		tet_face_shape(e, fa, mu[q], nv, t1, t2);
		const double cx = (t1[1] * t2[2]) - (t1[2] * t2[1]);
		const double cy = (t1[2] * t2[0]) - (t1[0] * t2[2]);
		const double cz = (t1[0] * t2[1]) - (t1[1] * t2[0]);
		const double w = wq[q] * sqrt((cx * cx) + (cy * cy) + (cz * cz));
		double tg = 0;
		if (t) {
			for (int a = 0; a < nen; a++) tg += nv[a] * t[nd[a]];
		}
		if (!rad) {
			for (int a = 0; a < nen; a++) {
				ge[a] += hc * ta * w * nv[a];
				for (int b = 0; b < nen; b++) me[a][b] += hc * w * nv[a] * nv[b];
			}
			if (t) qsum += hc * w * (tg - ta);
		}
		else {
			double qv, dq;
			if ((rad == 1) && RX.on && RX.fid && (RX.fid[f] >= 0)) {
				// 相互放射 : q = Q_ii σT⁴ + c_i (1 次の面と同じ。c_i は反復ごとに更新)
				const int fi = RX.fid[f];
				const double tk = tg + KELVIN;
				if (tk <= 0) *bad = 1;
				qv = (RX.qii[fi] * STEFAN * tk * tk * tk * tk) + RX.c[fi];
				dq = 4 * RX.qii[fi] * STEFAN * tk * tk * tk;
			}
			else {
				surf_q(rad, hc, ta, tg, &qv, &dq, bad);
			}
			qsum += w * qv;
			for (int a = 0; a < nen; a++) {
				ge[a] += w * ((dq * tg) - qv) * nv[a];
				for (int b = 0; b < nen; b++) me[a][b] += w * dq * nv[a] * nv[b];
			}
		}
	}
	if (B) {
		for (int a = 0; a < nen; a++) {
			for (int b = 0; b < nen; b++) {
				if (me[a][b] == 0) continue;
				const int64_t pos = crs_find_tri(B, nd[a], nd[b]);
				if (pos >= 0) B->val[pos] += me[a][b];
			}
		}
	}
	if (g) {
		for (int a = 0; a < nen; a++) g[nd[a]] += ge[a];
	}
	return qsum;
}


// 2 次の六面体・角柱の露出面 f (曲がった三角形・四角形) : 親の要素の形状関数を面に制限した
// 積分点 (elem3d_face_ip、9 点) で対流・輻射・h(T) を組む。接線・右辺・熱の評価点は同じ
static double surf_face_ho(crs_t *B, double *g, const double *t, int rad, int f,
	double hc, double ta, int *bad)
{
	const int e = TFaceEl[f];
	const int32_t *cn = &TFace[f * 4];
	const int k = ((cn[3] < 0) ? 3 : 4);
	int32_t nd[ELEM_NMAX];
	const int nen = elem3d_nodes(e, nd);
	static double me[ELEM_NMAX][ELEM_NMAX];
	double ge[ELEM_NMAX], qsum = 0;
	for (int a = 0; a < nen; a++) {
		ge[a] = 0;
		for (int b = 0; b < nen; b++) me[a][b] = 0;
	}
	const int nfq = elem3d_face_nip(e);
	for (int q = 0; q < nfq; q++) {
		double nv[ELEM_NMAX], w;
		if (elem3d_face_ip(e, cn, k, q, nv, &w) != nen) { *bad = 1; return 0; }
		double tg = 0;
		if (t) {
			for (int a = 0; a < nen; a++) tg += nv[a] * t[nd[a]];
		}
		if (!rad) {
			for (int a = 0; a < nen; a++) {
				ge[a] += hc * ta * w * nv[a];
				for (int b = 0; b < nen; b++) me[a][b] += hc * w * nv[a] * nv[b];
			}
			if (t) qsum += hc * w * (tg - ta);
		}
		else {
			double qv, dq;
			if ((rad == 1) && RX.on && RX.fid && (RX.fid[f] >= 0)) {
				// 相互放射 : q = Q_ii σT⁴ + c_i (四面体の面と同じ)
				const int fi = RX.fid[f];
				const double tk = tg + KELVIN;
				if (tk <= 0) *bad = 1;
				qv = (RX.qii[fi] * STEFAN * tk * tk * tk * tk) + RX.c[fi];
				dq = 4 * RX.qii[fi] * STEFAN * tk * tk * tk;
			}
			else {
				surf_q(rad, hc, ta, tg, &qv, &dq, bad);
			}
			qsum += w * qv;
			for (int a = 0; a < nen; a++) {
				ge[a] += w * ((dq * tg) - qv) * nv[a];
				for (int b = 0; b < nen; b++) me[a][b] += w * dq * nv[a] * nv[b];
			}
		}
	}
	if (B) {
		for (int a = 0; a < nen; a++) {
			for (int b = 0; b < nen; b++) {
				if (me[a][b] == 0) continue;
				const int64_t pos = crs_find_tri(B, nd[a], nd[b]);
				if (pos >= 0) B->val[pos] += me[a][b];
			}
		}
	}
	if (g) {
		for (int a = 0; a < nen; a++) g[nd[a]] += ge[a];
	}
	return qsum;
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
		const int32_t *nd = &TFace[f * 4];
		if (nd[2] < 0) {
			// 断面 2 次元の露出辺 (単位長あたり) : 対流は 1 次の線分で厳密 (質量 hL(1+δ)/6、
			// 右辺 h T_a L/2)、輻射・h(T) は 2 点 Gauss で接線・右辺・熱の評価点を揃える
			const double xc = (Xp[nd[0]] + Xp[nd[1]]) / 2, yc = (Yp[nd[0]] + Yp[nd[1]]) / 2;
			const double zc = (Zp[nd[0]] + Zp[nd[1]]) / 2;
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
			const surfbc_t *bc = surf_bc(((rad == 3) ? 0 : rad), xc, yc, zc, NULL);
			if ((bc == NULL) || (bc->c <= 0)) continue;
			nface++;
			const double hc = bc->c, ta = bc->t;
			const double dx = Xp[nd[1]] - Xp[nd[0]], dy = Yp[nd[1]] - Yp[nd[0]];
			const double dz = Zp[nd[1]] - Zp[nd[0]];
			const double len = sqrt((dx * dx) + (dy * dy) + (dz * dz));
			double me[2][2] = {{0, 0}, {0, 0}}, ge[2] = {0, 0};
			if (!rad) {
				for (int m = 0; m < 2; m++) {
					for (int mm = 0; mm < 2; mm++) me[m][mm] = hc * len * ((m == mm) ? 2 : 1) / 6;
					ge[m] = hc * ta * len / 2;
				}
				if (t) qsum += hc * len * (((t[nd[0]] + t[nd[1]]) / 2) - ta);
			}
			else {
				const double gq[2] = {0.5 - 0.5 / sqrt(3.0), 0.5 + 0.5 / sqrt(3.0)};
				for (int g2 = 0; g2 < 2; g2++) {
					const double nn[2] = {1 - gq[g2], gq[g2]};
					const double tg = (nn[0] * t[nd[0]]) + (nn[1] * t[nd[1]]);
					const double w = len / 2;
					double qv, dq;
					surf_q(rad, hc, ta, tg, &qv, &dq, &bad);
					qsum += w * qv;
					for (int m = 0; m < 2; m++) {
						ge[m] += w * (dq * tg - qv) * nn[m];
						for (int mm = 0; mm < 2; mm++) me[m][mm] += w * dq * nn[m] * nn[mm];
					}
				}
			}
			if (B) {
				for (int m = 0; m < 2; m++) {
					for (int mm = 0; mm < 2; mm++) {
						const int64_t pos = crs_find_tri(B, nd[m], nd[mm]);
						if (pos >= 0) B->val[pos] += me[m][mm];
					}
				}
			}
			if (g) {
				for (int m = 0; m < 2; m++) g[nd[m]] += ge[m];
			}
			continue;
		}
		if (nd[3] >= 0) {
			// 四角形の面 (六面体・角柱・ピラミッド)
			if (rad == 2) {
				const double xq = (Xp[nd[0]] + Xp[nd[1]] + Xp[nd[2]] + Xp[nd[3]]) / 4;
				const double yq = (Yp[nd[0]] + Yp[nd[1]] + Yp[nd[2]] + Yp[nd[3]]) / 4;
				const double zq = (Zp[nd[0]] + Zp[nd[1]] + Zp[nd[2]] + Zp[nd[3]]) / 4;
				const surfbc_t *cb = surf_bc(0, xq, yq, zq, NULL);
				const surfbc_t *rb = surf_bc(1, xq, yq, zq, NULL);
				if ((cb && (cb->c > 0)) || (rb && (rb->c > 0))) nface++;
				if (hits) {
					for (int e = 0; e < (NConv + NRad); e++) {
						const surfbc_t *bc = ((e < NConv) ? &ConvBC[e] : &RadBC[e - NConv]);
						if (surf_in(bc, xq, yq, zq)) hits[e]++;
					}
				}
				continue;
			}
			if ((HexNen > 8) || (PrismNen > 6)) {
				// 2 次の六面体・角柱 : 曲がった四角形を親の要素の形状関数で積分する
				const double xq = (Xp[nd[0]] + Xp[nd[1]] + Xp[nd[2]] + Xp[nd[3]]) / 4;
				const double yq = (Yp[nd[0]] + Yp[nd[1]] + Yp[nd[2]] + Yp[nd[3]]) / 4;
				const double zq = (Zp[nd[0]] + Zp[nd[1]] + Zp[nd[2]] + Zp[nd[3]]) / 4;
				const surfbc_t *bq = surf_bc(((rad == 3) ? 0 : rad), xq, yq, zq, NULL);
				if ((bq == NULL) || (bq->c <= 0)) continue;
				nface++;
				qsum += surf_face_ho(B, g, t, rad, f, bq->c, bq->t, &bad);
				continue;
			}
			const double qv = surf_quad(B, g, t, rad, nd,
				((RX.on && RX.fid) ? RX.fid[f] : -1), &bad, &nface);
			qsum += qv;
			continue;
		}
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
		if ((PrismNen > 6) && (elem3d_kind(TFaceEl[f]) == MESHELEM_PRISM)) {
			// 高次の角柱の三角形の面 (TetOrder も 2 以上になるので四面体より先に見る。混在格子では
			// 三角形の面を持つ要素が四面体のこともあるので、面を持つ要素の種別で振り分ける)
			qsum += surf_face_ho(B, g, t, rad, f, hc, ta, &bad);
			continue;
		}
		if (TetOrder >= 2) {
			// 高次の四面体 : 曲がった面を親の要素の形状関数で積分する
			qsum += surf_tri_ho(B, g, t, rad, f, hc, ta, &bad);
			continue;
		}

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
			for (int m = 0; m < 3; m++) {
				for (int mm = 0; mm < 3; mm++) me[m][mm] = 0;
			}
			for (int gq = 0; gq < 3; gq++) {
				double nn[3], tg = 0;
				for (int m = 0; m < 3; m++) {
					nn[m] = ((m == gq) ? (2.0 / 3) : (1.0 / 6));
					tg += nn[m] * t[nd[m]];
				}
				const double w = area / 3;
				double qv, dq;
				if ((rad == 1) && RX.on && (RX.fid[f] >= 0)) {
					// 相互放射 : q = Q_ii σT⁴ + c_i (構造格子と同じ)
					const int fi = RX.fid[f];
					const double tk = tg + KELVIN;
					if (tk <= 0) bad = 1;
					qv = (RX.qii[fi] * STEFAN * tk * tk * tk * tk) + RX.c[fi];
					dq = 4 * RX.qii[fi] * STEFAN * tk * tk * tk;
				}
				else {
					surf_q(rad, hc, ta, tg, &qv, &dq, &bad);
				}
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
	return ((t || (rad == 1) || (rad == 3)) ? qsum : (double)nface);
}


static double surf_faces(crs_t *B, double *g, const double *t, const double *kcell, int rad,
	int *hits)
{
	if (MeshMode) return surf_faces_tet(B, g, t, rad, hits);
	return surf_faces_grid(B, g, t, kcell, rad, hits);
}


/*
面どうしの相互放射 (radiationexchange = 1、構造格子)。

輻射を課した露出面を 1 枚ずつ放射度の要素とする (灰色・拡散、面内で一様な放射度)。
面 i の放射度 J、入射 G、正味の熱流 q (外向き) は
	J_i = ε_i E_i + (1 - ε_i) G_i,   G_i = Σ_j F_ij J_j + F_i∞ E∞_i,   q_i = J_i - G_i
(E = σT⁴ は面の平均、E∞ = σT_r⁴、F_i∞ = 1 - Σ_j F_ij は面の行の周囲温度の黒体へ)。
J は E の 1 次式なので一度だけ密に解いておく:
	M = I - diag(1 - ε) F,  X = M⁻¹ diag(ε),  y = M⁻¹ ((1 - ε) F∞ E∞)
	q = Q E + q0,   Q = (I - F) X,   q0 = (I - F) y - F∞ E∞
面の熱流の点の値を q(x) = Q_ii σT(x)⁴ + c_i (c_i = Σ_{j≠i} Q_ij E_j + q0_i) とし、
対角は輻射と同じく Newton 線形化、他の面の分 c は Picard の各回に更新する。
F = 0 (見える面が無い) なら Q_ii = ε_i、c_i = -ε_i E∞_i で通常の輻射に厳密に戻る。

形態係数は**閉形式** (面が軸に平行な長方形なので、周回積分
	A_i F_ij = (1/2π) ∮_i ∮_j ln r (dl_i・dl_j)
で平行な辺の対だけが残り、各対は ln r の 2 重積分の閉形式 F2 になる)。相互性
A_i F_ij = A_j F_ji は式の対称性から厳密。遮蔽は面の中心どうしを結ぶ線分が熱を
通すセルを通るかで判定する (全部見えるか全く見えないかの 2 値。部分的な遮蔽は近似)。
*/
static double rx_f2(double u, double d)
{
	const double r2 = (u * u) + (d * d);
	if (r2 <= 0) return 0;
	const double lr = 0.5 * log(r2);
	double v = (((u * u) - (d * d)) / 2) * lr - (0.75 * u * u);
	if (d > 0) v += d * u * atan(u / d);
	return v;
}


// 平行な線分の対の ∫∫ ln r ds1 ds2 (s1 ∈ [a1, a2], s2 ∈ [b1, b2]、垂直方向の距離 d)
static double rx_seg(double a1, double a2, double b1, double b2, double d)
{
	return rx_f2(a2 - b1, d) - rx_f2(a1 - b1, d) - rx_f2(a2 - b2, d) + rx_f2(a1 - b2, d);
}


// 面 f の 4 辺 (法線に対して反時計回り)。辺 k の軸 ax[k]、向き o[k] (±1)、
// 軸方向の範囲 [s0, s1]、他の 2 軸の座標 pc[k][3] (ax の成分は未使用)
static void rx_edges(int f, int ax[4], int o[4], double s0[4], double s1[4], double pc[4][3])
{
	const int d = RX.d[f];
	const int p = ((d == 0) ? 1 : 0), q = ((d == 2) ? 1 : 2);
	const double *lo = &RX.lo[f * 3], *hi = &RX.hi[f * 3];
	// e_p x e_q = +e_d (d = 0, 2) / -e_d (d = 1)
	const int ccw = (((d == 1) ? -1 : 1) * RX.sg[f]) > 0;
	// (p1,q1) -> (p2,q1) -> (p2,q2) -> (p1,q2) が e_p x e_q の向きに反時計回り
	const int axs[4] = {p, q, p, q};
	const int ors[4] = {+1, +1, -1, -1};
	for (int k = 0; k < 4; k++) {
		ax[k] = axs[k];
		o[k] = (ccw ? ors[k] : -ors[k]);
		s0[k] = lo[ax[k]];
		s1[k] = hi[ax[k]];
		for (int a = 0; a < 3; a++) pc[k][a] = lo[a];
		pc[k][d] = lo[d];
		// 辺の固定座標 : p 向きの辺は q = q1 (k = 0) / q2 (k = 2)、q 向きの辺は p = p2 (k = 1) / p1 (k = 3)
		if (k == 0) pc[k][q] = lo[q];
		if (k == 2) pc[k][q] = hi[q];
		if (k == 1) pc[k][p] = hi[p];
		if (k == 3) pc[k][p] = lo[p];
	}
}


// 面 i, j の A_i F_ij (遮蔽なし)
static double rx_aif(int i, int j)
{
	int ai[4], oi[4], aj[4], oj[4];
	double si0[4], si1[4], pi[4][3], sj0[4], sj1[4], pj[4][3];
	rx_edges(i, ai, oi, si0, si1, pi);
	rx_edges(j, aj, oj, sj0, sj1, pj);
	double sum = 0;
	for (int k = 0; k < 4; k++) {
		for (int l = 0; l < 4; l++) {
			if (ai[k] != aj[l]) continue;
			const int a = ai[k];
			double d2 = 0;
			for (int b = 0; b < 3; b++) {
				if (b == a) continue;
				d2 += (pi[k][b] - pj[l][b]) * (pi[k][b] - pj[l][b]);
			}
			sum += oi[k] * oj[l] * rx_seg(si0[k], si1[k], sj0[l], sj1[l], sqrt(d2));
		}
	}
	return sum / (2 * PI);
}


// 面の中心 (法線方向に少し浮かせた点) どうしを結ぶ線分が熱を通すセルを通るか
static int rx_blocked(int i, int j, const double *kcell)
{
	double a[3], b[3];
	for (int c = 0; c < 3; c++) {
		a[c] = (RX.lo[(i * 3) + c] + RX.hi[(i * 3) + c]) / 2;
		b[c] = (RX.lo[(j * 3) + c] + RX.hi[(j * 3) + c]) / 2;
	}
	// 格子の面上の点から始めると自分のセルに入るので、法線方向に刻みの 1e-6 だけ浮かせる
	const double *xs[3] = {Xn, Yn, Zn};
	const int nn[3] = {Nx, Ny, Nz};
	const double span = fabs(Xn[Nx] - Xn[0]) + fabs(Yn[Ny] - Yn[0]) + fabs(Zn[Nz] - Zn[0]);
	a[RX.d[i]] += RX.sg[i] * 1e-9 * span;
	b[RX.d[j]] += RX.sg[j] * 1e-9 * span;
	// 線分 a + t (b - a) が格子面を横切る t を集め、区間の中点のセルを見る
	double tt[3 * 4096 + 2];
	int nt = 0;
	tt[nt++] = 0;
	tt[nt++] = 1;
	for (int c = 0; c < 3; c++) {
		const double dd = b[c] - a[c];
		if (fabs(dd) <= 0) continue;
		for (int m = 0; (m <= nn[c]) && (nt < (3 * 4096)); m++) {
			const double t = (xs[c][m] - a[c]) / dd;
			if ((t > 0) && (t < 1)) tt[nt++] = t;
		}
	}
	// 挿入ソート (区間の数は高々数百)
	for (int m = 1; m < nt; m++) {
		const double v = tt[m];
		int k = m - 1;
		while ((k >= 0) && (tt[k] > v)) { tt[k + 1] = tt[k]; k--; }
		tt[k + 1] = v;
	}
	for (int m = 0; m + 1 < nt; m++) {
		const double tm = (tt[m] + tt[m + 1]) / 2;
		int ic[3], out = 0;
		for (int c = 0; c < 3; c++) {
			const double x = a[c] + (tm * (b[c] - a[c]));
			// 二分探索でセル
			int lo = 0, hi = nn[c];
			if ((x < xs[c][0]) || (x > xs[c][nn[c]])) { out = 1; break; }
			while ((hi - lo) > 1) {
				const int md = (lo + hi) / 2;
				if (xs[c][md] <= x) lo = md; else hi = md;
			}
			ic[c] = lo;
		}
		if (out) continue;				// 格子の外は空間
		if (kcell[((int64_t)ic[0] * Ny + ic[1]) * Nz + ic[2]] > 0) return 1;
	}
	return 0;
}


// N 元の密な連立 1 次方程式 (部分ピボット付き Gauss)。a は破壊され、b (n x m) に解が入る
static int rx_solve(int n, double *a, double *b, int m)
{
	for (int k = 0; k < n; k++) {
		int pv = k;
		for (int i = k + 1; i < n; i++) {
			if (fabs(a[(i * n) + k]) > fabs(a[(pv * n) + k])) pv = i;
		}
		if (fabs(a[(pv * n) + k]) <= 0) return 1;
		if (pv != k) {
			for (int c = 0; c < n; c++) {
				const double w = a[(k * n) + c]; a[(k * n) + c] = a[(pv * n) + c]; a[(pv * n) + c] = w;
			}
			for (int c = 0; c < m; c++) {
				const double w = b[(k * m) + c]; b[(k * m) + c] = b[(pv * m) + c]; b[(pv * m) + c] = w;
			}
		}
		const double piv = a[(k * n) + k];
		for (int i = k + 1; i < n; i++) {
			const double f = a[(i * n) + k] / piv;
			if (f == 0) continue;
			for (int c = k; c < n; c++) a[(i * n) + c] -= f * a[(k * n) + c];
			for (int c = 0; c < m; c++) b[(i * m) + c] -= f * b[(k * m) + c];
		}
	}
	for (int k = n - 1; k >= 0; k--) {
		for (int c = 0; c < m; c++) {
			double v = b[(k * m) + c];
			for (int j = k + 1; j < n; j++) v -= a[(k * n) + j] * b[(j * m) + c];
			b[(k * m) + c] = v / a[(k * n) + k];
		}
	}
	return 0;
}


#define RX_MAXFACE (2000)
static int rx_finish(FILE *fp_log, double *f, int nvis);

static void rx_free(void)
{
	free(RX.fid); free(RX.nv); free(RX.pv); free(RX.nrm); free(RX.ctr);
	free(RX.ocv); free(RX.okey); free(RX.fkey); free(RX.spx);
	free(RX.d); free(RX.sg); free(RX.row); free(RX.lo); free(RX.hi);
	free(RX.area); free(RX.eps); free(RX.tr);
	free(RX.qii); free(RX.c); free(RX.e); free(RX.q); free(RX.q0); free(RX.finf);
	memset(&RX, 0, sizeof(RX));
}


static int rx_build(FILE *fp_log, const double *kcell)
{
	memset(&RX, 0, sizeof(RX));
	const int n = (int)surf_faces_grid(NULL, NULL, NULL, kcell, 4, NULL);	// 面の数
	if (n > RX_MAXFACE) {
		fprintf(fp_log, "*** radiation exchange : %d radiating faces exceed the limit %d "
			"(the view-factor matrix is dense)\n", n, RX_MAXFACE);
		return 1;
	}
	RX.n = n;
	RX.d = (int *)malloc(((size_t)n + 1) * sizeof(int));
	RX.sg = (int *)malloc(((size_t)n + 1) * sizeof(int));
	RX.row = (int *)malloc(((size_t)n + 1) * sizeof(int));
	RX.lo = (double *)malloc(((size_t)n * 3 + 1) * sizeof(double));
	RX.hi = (double *)malloc(((size_t)n * 3 + 1) * sizeof(double));
	RX.area = (double *)malloc(((size_t)n + 1) * sizeof(double));
	RX.eps = (double *)malloc(((size_t)n + 1) * sizeof(double));
	RX.tr = (double *)malloc(((size_t)n + 1) * sizeof(double));
	RX.qii = (double *)malloc(((size_t)n + 1) * sizeof(double));
	RX.c = (double *)malloc(((size_t)n + 1) * sizeof(double));
	RX.e = (double *)malloc(((size_t)n + 1) * sizeof(double));
	RX.q0 = (double *)malloc(((size_t)n + 1) * sizeof(double));
	RX.finf = (double *)malloc(((size_t)n + 1) * sizeof(double));
	RX.q = (double *)malloc(((size_t)n * n + 1) * sizeof(double));
	surf_faces_grid(NULL, NULL, NULL, kcell, 4, NULL);		// 面の幾何を集める

	// 形態係数 (面どうしが向き合い、同一平面でなく、遮られないときだけ)
	double *f = (double *)calloc((size_t)n * n + 1, sizeof(double));
	int nvis = 0;
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			double ci[3], cj[3];
			for (int c = 0; c < 3; c++) {
				ci[c] = (RX.lo[(i * 3) + c] + RX.hi[(i * 3) + c]) / 2;
				cj[c] = (RX.lo[(j * 3) + c] + RX.hi[(j * 3) + c]) / 2;
			}
			const double fi = RX.sg[i] * (cj[RX.d[i]] - ci[RX.d[i]]);
			const double fj = RX.sg[j] * (ci[RX.d[j]] - cj[RX.d[j]]);
			if ((fi <= 0) || (fj <= 0)) continue;		// 向き合っていない (同一平面を含む)
			if (rx_blocked(i, j, kcell)) continue;
			const double aif = rx_aif(i, j);
			if (aif <= 0) continue;
			f[(i * n) + j] = aif / RX.area[i];
			f[(j * n) + i] = aif / RX.area[j];			// 相互性 A_i F_ij = A_j F_ji
			nvis++;
		}
	}
	return rx_finish(fp_log, f, nvis);
}


// 形態係数 f (n x n、F_ij) から F∞ と放射度の行列 Q, q0 を作る (構造格子・非構造格子で共通)。
// f は解放する
static int rx_finish(FILE *fp_log, double *f, int nvis)
{
	const int n = RX.n;
	double amax = 0, asum = 0, afinf = 0;
	for (int i = 0; i < n; i++) {
		double rs = 0;
		for (int j = 0; j < n; j++) rs += f[(i * n) + j];
		RX.finf[i] = 1 - rs;
		if (rs > amax) amax = rs;
		asum += RX.area[i];
		afinf += RX.area[i] * RX.finf[i];
	}
	fprintf(fp_log, "  radiation exchange : %d faces, %d visible pairs, max sum_j F_ij = %.12f, "
		"sum A F_inf / sum A = %.12e\n", n, nvis, amax, afinf / ((asum > 0) ? asum : 1));
	for (int r = 0; r < NRad; r++) {
		double ar = 0, af = 0;
		for (int i = 0; i < n; i++) {
			if (RX.row[i] != r) continue;
			ar += RX.area[i];
			af += RX.area[i] * RX.finf[i];
		}
		if (ar > 0) {
			fprintf(fp_log, "  radiation exchange : row %d, area = %.6e [m^2], mean F_inf = %.12f\n",
				r + 1, ar, af / ar);
		}
	}
	// 高次の四面体の曲がった面は境界を折れ線で近似するので、行の和が 1 を僅かに越え得る (粗い球の空隙で
	// 1.6e-3、実測)。1% までは警告にして続ける (F∞ が僅かに負 = 周囲からの僅かな入射になるだけで、熱の収支の
	// 恒等式は保たれる)。平らな面 (1 次) は従来どおり厳密に弾く
	const double rtol = ((MeshMode && (TetOrder >= 2)) ? 1e-2 : 1e-9);
	if (amax > 1 + rtol) {
		fprintf(fp_log, "*** radiation exchange : a view-factor row sum exceeds 1 (%.12f)\n", amax);
		free(f);
		return 1;
	}
	if (amax > 1 + 1e-9) {
		fprintf(fp_log, "*** warning : radiation exchange : a view-factor row sum exceeds 1 by %.3e "
			"(curved faces; refine the mesh or raise the element order)\n", amax - 1);
	}

	// M = I - diag(1-ε) F、右辺 [diag(ε) | (1-ε) F∞ E∞]
	double *m = (double *)malloc(((size_t)n * n + 1) * sizeof(double));
	double *rhs = (double *)malloc(((size_t)n * (n + 1) + 1) * sizeof(double));
	for (int i = 0; i < n; i++) {
		const double einf = STEFAN * pow(RX.tr[i] + KELVIN, 4);
		for (int j = 0; j < n; j++) {
			m[(i * n) + j] = ((i == j) ? 1 : 0) - ((1 - RX.eps[i]) * f[(i * n) + j]);
			rhs[(i * (n + 1)) + j] = ((i == j) ? RX.eps[i] : 0);
		}
		rhs[(i * (n + 1)) + n] = (1 - RX.eps[i]) * RX.finf[i] * einf;
	}
	if (rx_solve(n, m, rhs, n + 1)) {
		fprintf(fp_log, "*** radiation exchange : the radiosity system is singular\n");
		free(f); free(m); free(rhs);
		return 1;
	}
	// Q = (I - F) X、q0 = (I - F) y - F∞ E∞
	for (int i = 0; i < n; i++) {
		const double einf = STEFAN * pow(RX.tr[i] + KELVIN, 4);
		for (int j = 0; j <= n; j++) {
			double v = rhs[(i * (n + 1)) + j];
			for (int k = 0; k < n; k++) v -= f[(i * n) + k] * rhs[(k * (n + 1)) + j];
			if (j < n) RX.q[(i * n) + j] = v;
			else       RX.q0[i] = v - (RX.finf[i] * einf);
		}
		RX.qii[i] = RX.q[(i * n) + i];
		RX.c[i] = RX.q0[i];
	}
	free(f); free(m); free(rhs);
	RX.on = 1;
	return 0;
}


/*
非構造格子の相互放射 (四面体・六面体・角柱・ピラミッドの露出面)。

面は一般の平面多角形 (三角形・四角形) なので、形態係数は同じ周回積分
	A_i F_ij = (1/2π) Σ_k Σ_l (t_k・t_l) ∫∫ ln r ds dt     (t は辺の単位接線)
を辺の対ごとに計算する。辺は外向きの法線に対して反時計回り (構造格子と同じ規約)。
	・平行な辺 (共線を含む) : 構造格子と同じ閉形式 rx_seg (相互性が式の対称性から厳密)
	・直交する辺 : t_k・t_l = 0 で寄与なし (軸に平行な六面体の面はこの 2 種類だけなので、
	  構造格子と同じ閉形式の和になり答えが一致する)
	・それ以外 : 内側の線積分は閉形式 ∫ ln|x - q(t)| dt = G(L - τ) - G(-τ)
	  (G(v) = v ln√(v² + d²) - v + d atan(v/d))、外側は 8 点 Gauss。頂点を共有する辺
	  (隣り合う面) では被積分関数がその端で v ln v の形になるので、端に向けて区間を
	  幾何級数で細かくする (24 段)
遮蔽は構造格子と同じく面の中心どうしを結ぶ線分 (法線方向に少し浮かせる) で 2 値に
判定し、熱を通す要素の境界面 (電極の面も含む) との交差を見る。
*/
static double rx_g(double v, double d)
{
	if (d > 0) return (v * 0.5 * log((v * v) + (d * d))) - v + (d * atan(v / d));
	if (v == 0) return 0;
	return (v * log(fabs(v))) - v;
}


// 点 x と線分 q0 + t e (t ∈ [0, L]) の ∫ ln|x - q(t)| dt
static double rx_pseg(const double x[3], const double q0[3], const double e[3], double l)
{
	double w[3], tau = 0, r2 = 0;
	for (int c = 0; c < 3; c++) {
		w[c] = x[c] - q0[c];
		tau += w[c] * e[c];
		r2 += w[c] * w[c];
	}
	double d2 = r2 - (tau * tau);
	if (d2 < 0) d2 = 0;
	const double d = sqrt(d2);
	return rx_g(l - tau, d) - rx_g(-tau, d);
}


// 点 x と線分 [a, b] の距離
static double rx_pdist(const double x[3], const double a[3], const double b[3])
{
	double ab[3], ax[3], l2 = 0, t = 0;
	for (int c = 0; c < 3; c++) {
		ab[c] = b[c] - a[c];
		ax[c] = x[c] - a[c];
		l2 += ab[c] * ab[c];
		t += ab[c] * ax[c];
	}
	t = ((l2 > 0) ? (t / l2) : 0);
	if (t < 0) t = 0;
	if (t > 1) t = 1;
	double r2 = 0;
	for (int c = 0; c < 3; c++) {
		const double v = ax[c] - (t * ab[c]);
		r2 += v * v;
	}
	return sqrt(r2);
}


// 辺の対 (p0 -> p1, q0 -> q1) の (t_p・t_q) ∫∫ ln r ds dt
static double rx_segpair(const double p0[3], const double p1[3], const double q0[3], const double q1[3])
{
	double ep[3], eq[3], lp = 0, lq = 0;
	for (int c = 0; c < 3; c++) {
		ep[c] = p1[c] - p0[c];
		eq[c] = q1[c] - q0[c];
		lp += ep[c] * ep[c];
		lq += eq[c] * eq[c];
	}
	lp = sqrt(lp);
	lq = sqrt(lq);
	if ((lp <= 0) || (lq <= 0)) return 0;
	double cs = 0;
	for (int c = 0; c < 3; c++) {
		ep[c] /= lp;
		eq[c] /= lq;
		cs += ep[c] * eq[c];
	}
	const double cr[3] = {(ep[1] * eq[2]) - (ep[2] * eq[1]), (ep[2] * eq[0]) - (ep[0] * eq[2]),
	                      (ep[0] * eq[1]) - (ep[1] * eq[0])};
	const double sn = sqrt((cr[0] * cr[0]) + (cr[1] * cr[1]) + (cr[2] * cr[2]));
	if (fabs(cs) <= 1e-14) return 0;				// 直交 : 寄与なし
	if (sn <= 1e-12) {
		// 平行 (共線を含む) : ep に沿った座標で閉形式
		double b0 = 0, b1 = 0;
		for (int c = 0; c < 3; c++) {
			b0 += (q0[c] - p0[c]) * ep[c];
			b1 += (q1[c] - p0[c]) * ep[c];
		}
		double d2 = 0;
		for (int c = 0; c < 3; c++) {
			const double v = (q0[c] - p0[c]) - (b0 * ep[c]);
			d2 += v * v;
		}
		const double lo = ((b0 < b1) ? b0 : b1), hi = ((b0 < b1) ? b1 : b0);
		return ((cs > 0) ? 1 : -1) * rx_seg(0, lp, lo, hi, sqrt(d2));
	}
	// 一般の向き : 外側を 8 点 Gauss。端が相手の辺に近ければ、その端に向けて幾何級数で刻む
	static const double gx[4] = {0.18343464249564980494, 0.52553240991632898582,
	                             0.79666647741362673959, 0.96028985649753623168};
	static const double gw[4] = {0.36268378337836198297, 0.31370664587788728734,
	                             0.22238103445337447054, 0.10122853629037625915};
	const int g0 = (rx_pdist(p0, q0, q1) < (0.1 * lp));
	const int g1 = (rx_pdist(p1, q0, q1) < (0.1 * lp));
	double cut[64];
	int nc = 0;
	cut[nc++] = 0;
	if (g0) {
		for (int k = 24; k >= 2; k--) cut[nc++] = lp * ldexp(1.0, -k);
	}
	if (g0 || g1) cut[nc++] = 0.5 * lp;
	if (g1) {
		for (int k = 2; k <= 24; k++) cut[nc++] = lp * (1 - ldexp(1.0, -k));
	}
	cut[nc++] = lp;
	double sum = 0;
	for (int m = 0; m + 1 < nc; m++) {
		const double a = cut[m], b = cut[m + 1];
		const double h = 0.5 * (b - a), mid = 0.5 * (a + b);
		for (int k = 0; k < 4; k++) {
			for (int sg = -1; sg <= 1; sg += 2) {
				const double sv = mid + (sg * h * gx[k]);
				double x[3];
				for (int c = 0; c < 3; c++) x[c] = p0[c] + (sv * ep[c]);
				sum += h * gw[k] * rx_pseg(x, q0, eq, lq);
			}
		}
	}
	return cs * sum;
}


// 面 i, j (非構造格子) の A_i F_ij (遮蔽なし)
static double rx_aif_poly(int i, int j)
{
	const int ni = RX.nv[i], nj = RX.nv[j];
	const double *pi = &RX.pv[i * 3 * RX_MAXV], *pj = &RX.pv[j * 3 * RX_MAXV];
	double sum = 0;
	for (int k = 0; k < ni; k++) {
		const double *a0 = &pi[k * 3], *a1 = &pi[((k + 1) % ni) * 3];
		for (int l = 0; l < nj; l++) {
			sum += rx_segpair(a0, a1, &pj[l * 3], &pj[((l + 1) % nj) * 3]);
		}
	}
	return sum / (2 * PI);
}


// 線分 a -> b が三角形 (v0, v1, v2) と交わるか (端点を除く。Möller-Trumbore)
static int rx_hit_tri(const double a[3], const double b[3], const double *v0, const double *v1,
	const double *v2)
{
	double d[3], e1[3], e2[3];
	for (int c = 0; c < 3; c++) {
		d[c] = b[c] - a[c];
		e1[c] = v1[c] - v0[c];
		e2[c] = v2[c] - v0[c];
	}
	const double pv[3] = {(d[1] * e2[2]) - (d[2] * e2[1]), (d[2] * e2[0]) - (d[0] * e2[2]),
	                      (d[0] * e2[1]) - (d[1] * e2[0])};
	const double det = (e1[0] * pv[0]) + (e1[1] * pv[1]) + (e1[2] * pv[2]);
	if (fabs(det) <= 0) return 0;					// 面と平行
	const double tv[3] = {a[0] - v0[0], a[1] - v0[1], a[2] - v0[2]};
	// **辺の上は「当たり」にする** (許容 1e-9)。隣り合う三角形の共有辺をちょうど通る
	// 線分は、丸めで両方の三角形を外れることがある (水密でない) ので、境界を厳密に
	// 判定すると壁越しの対が OS によって見えたり見えなかったりする (実測 : CI の
	// Linux / macOS でだけ 2 つの空洞の F∞ が 6e-4 ずれた)。端点は面から浮かせてあるので、
	// 線分の両端の近く (自分の面とその隣) は t の許容で除く
	const double tol = 1e-9;
	const double u = ((tv[0] * pv[0]) + (tv[1] * pv[1]) + (tv[2] * pv[2])) / det;
	if ((u < -tol) || (u > 1 + tol)) return 0;
	const double qv[3] = {(tv[1] * e1[2]) - (tv[2] * e1[1]), (tv[2] * e1[0]) - (tv[0] * e1[2]),
	                      (tv[0] * e1[1]) - (tv[1] * e1[0])};
	const double v = ((d[0] * qv[0]) + (d[1] * qv[1]) + (d[2] * qv[2])) / det;
	if ((v < -tol) || ((u + v) > 1 + tol)) return 0;
	const double t = ((e2[0] * qv[0]) + (e2[1] * qv[1]) + (e2[2] * qv[2])) / det;
	return ((t > 1e-6) && (t < 1 - 1e-6));
}


// 高次の曲がった面 : 中心で向き合わなくても、面の積分点のすべての対が互いを前に見る (両方の cos ≥ 0) なら
// 見える (周回積分はその対で厳密)。実測 : 円筒の穴の z に積んだ 2 枚の面は中心の法線が中心を結ぶ線と
// 直交して「同一平面」と判定され、その分の F が抜けて F∞ が +4.4e-3 (曲がった円環、8 x 2 枚)。
// 一部の対だけが見える (部分的な可視) は従来どおり中心で決める。平らな面は全部の点で ±丸めになり見えない。
// 負の側の許容は cos ≥ −5e-2 (約 3°) : gmsh の高次の格子は面の内部の節点を角度方向にずらして置く (円筒の
// 上には乗るが母線が直線でない、実測 3 次 14.64〜15°・4 次 10.80〜11.25°) ので、同じ角度の点の対の cos が
// 0 でなく 3 次 −2e-3・4 次 −1.5e-2 になる。本当に一部だけ見える対は O(1) の負の cos を持つ。わずかに負の
// 対の周回積分への寄与は cos の積なので無視できる
static int rx_vis_curved(int i, int j)
{
	if (RX.spx == NULL) return 0;
	int pos = 0;
	for (int q = 0; q < RX_NS; q++) {
		const double *pq = &RX.spx[(((size_t)i * RX_NS) + q) * 6];
		for (int r = 0; r < RX_NS; r++) {
			const double *pr = &RX.spx[(((size_t)j * RX_NS) + r) * 6];
			const double d[3] = {pr[0] - pq[0], pr[1] - pq[1], pr[2] - pq[2]};
			const double dl = sqrt((d[0] * d[0]) + (d[1] * d[1]) + (d[2] * d[2]));
			const double c1 = (pq[3] * d[0]) + (pq[4] * d[1]) + (pq[5] * d[2]);
			const double c2 = -((pr[3] * d[0]) + (pr[4] * d[1]) + (pr[5] * d[2]));
			if ((c1 < -5e-2 * dl) || (c2 < -5e-2 * dl)) return 0;
			if ((c1 > 1e-9 * dl) && (c2 > 1e-9 * dl)) pos = 1;
		}
	}
	return pos;
}


static int rx_blocked_tet(int i, int j)
{
	double a[3], b[3];
	const double span = SurfEps / EPS;
	for (int c = 0; c < 3; c++) {
		a[c] = RX.ctr[(i * 3) + c] + (1e-9 * span * RX.nrm[(i * 3) + c]);
		b[c] = RX.ctr[(j * 3) + c] + (1e-9 * span * RX.nrm[(j * 3) + c]);
	}
	if (RX.ocv != NULL) {
		const int32_t *ki = &RX.fkey[i * 4], *kj = &RX.fkey[j * 4];
		for (int o = 0; o < RX.nocc; o++) {
			const int32_t *ko = &RX.okey[o * 4];
			if ((ko[0] == ki[0]) && (ko[1] == ki[1]) && (ko[2] == ki[2]) && (ko[3] == ki[3])) continue;
			if ((ko[0] == kj[0]) && (ko[1] == kj[1]) && (ko[2] == kj[2]) && (ko[3] == kj[3])) continue;
			const double *v = &RX.ocv[o * 9];
			if (rx_hit_tri(a, b, &v[0], &v[3], &v[6])) return 1;
		}
		return 0;
	}
	for (int f = 0; f < TOccN; f++) {
		const int32_t *nd = &TOcc[f * 4];
		double v[4][3];
		const int k = ((nd[3] < 0) ? 3 : 4);
		for (int l = 0; l < k; l++) {
			v[l][0] = Xp[nd[l]]; v[l][1] = Yp[nd[l]]; v[l][2] = Zp[nd[l]];
		}
		if (rx_hit_tri(a, b, v[0], v[1], v[2])) return 1;
		if ((k == 4) && rx_hit_tri(a, b, v[0], v[2], v[3])) return 1;
	}
	return 0;
}


// 高次の露出面の角の媒介変数 (三角形 (0,0) (1,0) (0,1)、四角形 (0,0) (1,0) (1,1) (0,1))
static const double RX_CST[4][2] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
static const double RX_CST3[3][2] = {{0, 0}, {1, 0}, {0, 1}};

// 面の角 cn[k] を昇順にした鍵 (三角形の 4 番目は -1)
static void rx_key(const int32_t *cn, int k, int32_t key[4])
{
	for (int l = 0; l < 4; l++) key[l] = ((l < k) ? cn[l] : -1);
	for (int x = 0; x < k - 1; x++) {
		for (int y = x + 1; y < k; y++) {
			if (key[y] < key[x]) { const int32_t w = key[x]; key[x] = key[y]; key[y] = w; }
		}
	}
}


// 非構造格子の rx_build : 輻射を課した露出面を集め、形態係数を作る
static int rx_build_tet(FILE *fp_log)
{
	memset(&RX, 0, sizeof(RX));
	RX.fid = (int *)malloc(((size_t)TFaceN + 1) * sizeof(int));
	int n = 0;
	for (int f = 0; f < TFaceN; f++) {
		const int32_t *nd = &TFace[f * 4];
		const int k = ((nd[3] < 0) ? 3 : 4);
		double c[3] = {0, 0, 0};
		for (int l = 0; l < k; l++) {
			c[0] += Xp[nd[l]] / k; c[1] += Yp[nd[l]] / k; c[2] += Zp[nd[l]] / k;
		}
		const surfbc_t *bc = surf_bc(1, c[0], c[1], c[2], NULL);
		RX.fid[f] = (((bc != NULL) && (bc->c > 0)) ? n++ : -1);
	}
	if (n > RX_MAXFACE) {
		fprintf(fp_log, "*** radiation exchange : %d radiating faces exceed the limit %d "
			"(the view-factor matrix is dense)\n", n, RX_MAXFACE);
		return 1;
	}
	RX.n = n;
	RX.nv = (int *)malloc(((size_t)n + 1) * sizeof(int));
	RX.pv = (double *)malloc(((size_t)n * 3 * RX_MAXV + 1) * sizeof(double));
	RX.nrm = (double *)malloc(((size_t)n * 3 + 1) * sizeof(double));
	RX.ctr = (double *)malloc(((size_t)n * 3 + 1) * sizeof(double));
	RX.fkey = (int32_t *)malloc(((size_t)n * 4 + 1) * sizeof(int32_t));
	RX.row = (int *)malloc(((size_t)n + 1) * sizeof(int));
	RX.area = (double *)malloc(((size_t)n + 1) * sizeof(double));
	RX.eps = (double *)malloc(((size_t)n + 1) * sizeof(double));
	RX.tr = (double *)malloc(((size_t)n + 1) * sizeof(double));
	RX.qii = (double *)malloc(((size_t)n + 1) * sizeof(double));
	RX.c = (double *)malloc(((size_t)n + 1) * sizeof(double));
	RX.e = (double *)malloc(((size_t)n + 1) * sizeof(double));
	RX.q0 = (double *)malloc(((size_t)n + 1) * sizeof(double));
	RX.finf = (double *)malloc(((size_t)n + 1) * sizeof(double));
	RX.q = (double *)malloc(((size_t)n * n + 1) * sizeof(double));
	if (TetOrder >= 2) RX.spx = (double *)calloc((size_t)n * RX_NS * 6 + 1, sizeof(double));
	for (int f = 0; f < TFaceN; f++) {
		const int i = RX.fid[f];
		if (i < 0) continue;
		const int32_t *nd = &TFace[f * 4];
		int k = ((nd[3] < 0) ? 3 : 4);
		double v[RX_MAXV][3], c[3] = {0, 0, 0}, an[3] = {0, 0, 0}, ar = 0;
		if (TetOrder >= 2) {
			// 高次の要素 (四面体・六面体・角柱) の曲がった面 : 境界は辺を 2p 分割した折れ線 (点は曲がった辺の上)、
			// 面積は曲がった面の積分、中心と法線は面の媒介変数の中央の点
			const int e = TFaceEl[f], msub = 2 * TetOrder;
			const double (*cst)[2] = ((k == 3) ? RX_CST3 : RX_CST);
			double nv[ELEM_NMAX], t1[3], t2[3];
			for (int s3 = 0; s3 < k; s3++) {
				for (int jj = 0; jj < msub; jj++) {
					const double w = (double)jj / msub;
					const double ss = ((1 - w) * cst[s3][0]) + (w * cst[(s3 + 1) % k][0]);
					const double tt = ((1 - w) * cst[s3][1]) + (w * cst[(s3 + 1) % k][1]);
					elem3d_face_point(e, nd, k, ss, tt, nv, v[(s3 * msub) + jj], t1, t2);
				}
			}
			elem3d_face_point(e, nd, k, ((k == 3) ? (1.0 / 3) : 0.5), ((k == 3) ? (1.0 / 3) : 0.5), nv, c, t1, t2);
			an[0] = (t1[1] * t2[2]) - (t1[2] * t2[1]);
			an[1] = (t1[2] * t2[0]) - (t1[0] * t2[2]);
			an[2] = (t1[0] * t2[1]) - (t1[1] * t2[0]);
			double st[QTRI_MAX][2], wq[QTRI_MAX];
			const int nq = face_quad_pub(k, TetOrder + 2, st, wq);
			for (int q = 0; q < nq; q++) {
				double xq[3];
				elem3d_face_point(e, nd, k, st[q][0], st[q][1], nv, xq, t1, t2);
				const double cx = (t1[1] * t2[2]) - (t1[2] * t2[1]);
				const double cy = (t1[2] * t2[0]) - (t1[0] * t2[2]);
				const double cz = (t1[0] * t2[1]) - (t1[1] * t2[0]);
				ar += wq[q] * sqrt((cx * cx) + (cy * cy) + (cz * cz));
			}
			k = k * msub;
			// 折れ線の向きを法線 an に合わせる (面積ベクトルと an が同じ向き)
			double pa[3] = {0, 0, 0};
			for (int l = 0; l < k; l++) {
				const double *pp = v[l], *qq = v[(l + 1) % k];
				pa[0] += 0.5 * ((pp[1] * qq[2]) - (pp[2] * qq[1]));
				pa[1] += 0.5 * ((pp[2] * qq[0]) - (pp[0] * qq[2]));
				pa[2] += 0.5 * ((pp[0] * qq[1]) - (pp[1] * qq[0]));
			}
			if (((pa[0] * an[0]) + (pa[1] * an[1]) + (pa[2] * an[2])) < 0) {
				for (int l = 1; l < (k + 1) / 2; l++) {
					for (int a = 0; a < 3; a++) {
						const double w = v[l][a]; v[l][a] = v[k - l][a]; v[k - l][a] = w;
					}
				}
			}
		}
		else {
			for (int l = 0; l < k; l++) {
				v[l][0] = Xp[nd[l]]; v[l][1] = Yp[nd[l]]; v[l][2] = Zp[nd[l]];
				for (int a = 0; a < 3; a++) c[a] += v[l][a] / k;
			}
			// 面積ベクトル (多角形の外積の和の半分)
			for (int l = 0; l < k; l++) {
				const double *p = v[l], *q = v[(l + 1) % k];
				an[0] += 0.5 * ((p[1] * q[2]) - (p[2] * q[1]));
				an[1] += 0.5 * ((p[2] * q[0]) - (p[0] * q[2]));
				an[2] += 0.5 * ((p[0] * q[1]) - (p[1] * q[0]));
			}
		}
		// 外向き : 要素の重心から面の中心へ向かう側。高次の四面体は面に無い頂点から面の中心へ向かう側
		// (粗く曲がった要素では全節点の重心が面の外側に出ることがある。実測 : 粗い球の空隙で 97 対の
		// 法線が裏返って「向き合わない」と判定された)
		int32_t en[ELEM_NMAX];
		const int nen = elem3d_nodes(TFaceEl[f], en);
		double ec[3] = {0, 0, 0};
		if (TetOrder >= 2) {
			const int kk = ((nd[3] < 0) ? 3 : 4), kind = elem3d_kind(TFaceEl[f]);
			const int ncor = ((kind == MESHELEM_HEX) ? 8 : (kind == MESHELEM_PRISM) ? 6 : 4);
			for (int l = ncor - 1; l >= 0; l--) {
				int on = 0;
				for (int m = 0; m < kk; m++) on |= (en[l] == nd[m]);
				if (!on) { ec[0] = Xp[en[l]]; ec[1] = Yp[en[l]]; ec[2] = Zp[en[l]]; }
			}
		}
		else {
			for (int l = 0; l < nen; l++) {
				ec[0] += Xp[en[l]] / nen; ec[1] += Yp[en[l]] / nen; ec[2] += Zp[en[l]] / nen;
			}
		}
		const double dot = (an[0] * (c[0] - ec[0])) + (an[1] * (c[1] - ec[1])) + (an[2] * (c[2] - ec[2]));
		const double anr = sqrt((an[0] * an[0]) + (an[1] * an[1]) + (an[2] * an[2]));
		if (TetOrder < 2) ar = anr;
		RX.nv[i] = k;
		for (int l = 0; l < k; l++) {
			const int ls = ((dot >= 0) ? l : ((k - l) % k));	// 内向きなら並びを逆に
			for (int a = 0; a < 3; a++) RX.pv[(i * 3 * RX_MAXV) + (l * 3) + a] = v[ls][a];
		}
		for (int a = 0; a < 3; a++) {
			RX.nrm[(i * 3) + a] = ((dot >= 0) ? 1 : -1) * an[a] / anr;
			RX.ctr[(i * 3) + a] = c[a];
		}
		rx_key(nd, ((nd[3] < 0) ? 3 : 4), &RX.fkey[i * 4]);
		RX.area[i] = ar;
		if (RX.spx != NULL) {
			const int ks = ((nd[3] < 0) ? 3 : 4);
			double st[QTRI_MAX][2], wq[QTRI_MAX];
			const int ns = face_quad_pub(ks, 3, st, wq);
			for (int q = 0; (q < ns) && (q < RX_NS); q++) {
				double nv[ELEM_NMAX], t1[3], t2[3], *sp = &RX.spx[(((size_t)i * RX_NS) + q) * 6];
				elem3d_face_point(TFaceEl[f], nd, ks, st[q][0], st[q][1], nv, sp, t1, t2);
				const double cr[3] = {(t1[1] * t2[2]) - (t1[2] * t2[1]), (t1[2] * t2[0]) - (t1[0] * t2[2]),
				                      (t1[0] * t2[1]) - (t1[1] * t2[0])};
				const double cl = sqrt((cr[0] * cr[0]) + (cr[1] * cr[1]) + (cr[2] * cr[2]));
				const double sg = ((((cr[0] * RX.nrm[i * 3]) + (cr[1] * RX.nrm[(i * 3) + 1]) + (cr[2] * RX.nrm[(i * 3) + 2])) >= 0) ? 1 : -1);
				for (int a = 0; a < 3; a++) sp[3 + a] = ((cl > 0) ? (sg * cr[a] / cl) : 0);
			}
		}
		// 係数は面を選んだときと同じ点 (頂点の平均) で引く。曲がった面の中央の点で引くと、直方体の境目の
		// 近くの面が選ばれたのに係数が無い (NULL) ことがある (実測 : 円環の外側の面で異常終了)
		const int kc = ((nd[3] < 0) ? 3 : 4);
		double c0[3] = {0, 0, 0};
		for (int l = 0; l < kc; l++) {
			c0[0] += Xp[nd[l]] / kc; c0[1] += Yp[nd[l]] / kc; c0[2] += Zp[nd[l]] / kc;
		}
		int row = 0;
		const surfbc_t *bc = surf_bc(1, c0[0], c0[1], c0[2], &row);
		RX.row[i] = row;
		RX.eps[i] = bc->c;
		RX.tr[i] = bc->t;
	}

	// 高次の要素 : 遮蔽物を曲がった面にする (各面を媒介変数の格子で平らな三角形に割る。点は曲がった面の上。
	// 三角形の面は (2p)² 枚、四角形の面は 2 (2p)² 枚)
	if (TetOrder >= 2) {
		const int m = 2 * TetOrder;
		RX.nocc = TOccN * 2 * m * m;
		RX.ocv = (double *)malloc(((size_t)RX.nocc * 9 + 1) * sizeof(double));
		RX.okey = (int32_t *)malloc(((size_t)RX.nocc * 4 + 1) * sizeof(int32_t));
		int o = 0;
		for (int fo = 0; fo < TOccN; fo++) {
			const int32_t *cn = &TOcc[fo * 4];
			const int k = ((cn[3] < 0) ? 3 : 4);
			const int e = TOccEl[fo];
			int32_t key[4];
			rx_key(cn, k, key);
			// 格子点 (a, b) の物理座標 : 媒介変数 (a/m, b/m)
			double lp[(2 * 7 + 1) * (2 * 7 + 1)][3];
			for (int aa = 0; aa <= m; aa++) {
				for (int bb = 0; bb <= m; bb++) {
					if ((k == 3) && (aa + bb > m)) continue;
					double nv[ELEM_NMAX], t1[3], t2[3];
					elem3d_face_point(e, cn, k, (double)aa / m, (double)bb / m, nv, lp[(aa * (m + 1)) + bb], t1, t2);
				}
			}
			for (int aa = 0; aa < m; aa++) {
				for (int bb = 0; bb < m; bb++) {
					if ((k == 3) && (aa + bb >= m)) continue;
					const int tri[2][3][2] = {{{aa, bb}, {aa + 1, bb}, {aa, bb + 1}},
					                          {{aa + 1, bb}, {aa + 1, bb + 1}, {aa, bb + 1}}};
					const int nt = (((k == 4) || (aa + bb + 1 < m)) ? 2 : 1);
					for (int t = 0; t < nt; t++) {
						for (int v = 0; v < 3; v++) {
							const double *x = lp[(tri[t][v][0] * (m + 1)) + tri[t][v][1]];
							for (int c = 0; c < 3; c++) RX.ocv[(o * 9) + (v * 3) + c] = x[c];
						}
						for (int c = 0; c < 4; c++) RX.okey[(o * 4) + c] = key[c];
						o++;
					}
				}
			}
		}
		RX.nocc = o;
	}
	double *f = (double *)calloc((size_t)n * n + 1, sizeof(double));
	int nvis = 0;
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			double fi = 0, fj = 0, dd = 0;
			for (int a = 0; a < 3; a++) {
				const double dc = RX.ctr[(j * 3) + a] - RX.ctr[(i * 3) + a];
				fi += RX.nrm[(i * 3) + a] * dc;
				fj -= RX.nrm[(j * 3) + a] * dc;
				dd += dc * dc;
			}
			// 向き合っていない (同一平面を含む)。**判定は相対の許容で** : 同一平面の三角形は
			// 厳密には fi = 0 だが、座標の丸めで ±1e-20 になり、その対の (厳密には 0 の)
			// 形態係数の数値積分の誤差 (相対 1e-9) が F に入る (実測 : 座標を 1 ulp 揺らすと
			// 見える対が 19300 -> 19658 に増えて F∞ が 3e-9 ずれ、CI の OS で結果が割れた)
			const double ftol = 1e-9 * sqrt(dd);
			if (((fi <= ftol) || (fj <= ftol)) && !rx_vis_curved(i, j)) continue;
			if (rx_blocked_tet(i, j)) continue;
			const double aif = rx_aif_poly(i, j);
			if (aif <= 0) continue;
			f[(i * n) + j] = aif / RX.area[i];
			f[(j * n) + i] = aif / RX.area[j];
			nvis++;
		}
	}
	// 曲がった面は凹なら自分を見る (F_ii > 0。球の内側では F_ii = A_i / A_球、平らな面は 0)。周回積分は
	// 自分自身には使えないので、面の積分点の対で直接積分する : A_i F_ii = ∫∫ cosθ1 cosθ2 / (π r²) dA dA
	// (両方の cos が正の対だけ)。被積分関数は有界なので、対角 (q = q') を除いた和を
	// A² / (A² − Σ w²) 倍して補う (被積分関数が一定の球では厳密)
	if (TetOrder >= 2) {
		for (int fc = 0; fc < TFaceN; fc++) {
			const int i = RX.fid[fc];
			if (i < 0) continue;
			const int e = TFaceEl[fc];
			const int32_t *cn = &TFace[fc * 4];
			const int k = ((cn[3] < 0) ? 3 : 4);
			double st[QTRI_MAX][2], wq[QTRI_MAX];
			const int nq = face_quad_pub(k, 7, st, wq);
			double xq[QTRI_MAX][3], nq3[QTRI_MAX][3], wa[QTRI_MAX], aa = 0, w2 = 0;
			for (int q = 0; q < nq; q++) {
				double nv[ELEM_NMAX], t1[3], t2[3];
				elem3d_face_point(e, cn, k, st[q][0], st[q][1], nv, xq[q], t1, t2);
				double cr[3] = {(t1[1] * t2[2]) - (t1[2] * t2[1]), (t1[2] * t2[0]) - (t1[0] * t2[2]),
				                (t1[0] * t2[1]) - (t1[1] * t2[0])};
				const double cl = sqrt((cr[0] * cr[0]) + (cr[1] * cr[1]) + (cr[2] * cr[2]));
				const double sg = ((((cr[0] * RX.nrm[i * 3]) + (cr[1] * RX.nrm[(i * 3) + 1]) + (cr[2] * RX.nrm[(i * 3) + 2])) >= 0) ? 1 : -1);
				for (int a = 0; a < 3; a++) nq3[q][a] = ((cl > 0) ? (sg * cr[a] / cl) : 0);
				wa[q] = wq[q] * cl;
				aa += wa[q];
				w2 += wa[q] * wa[q];
			}
			double sum = 0;
			for (int q = 0; q < nq; q++) {
				for (int r = 0; r < nq; r++) {
					if (r == q) continue;
					const double d[3] = {xq[r][0] - xq[q][0], xq[r][1] - xq[q][1], xq[r][2] - xq[q][2]};
					const double r2 = (d[0] * d[0]) + (d[1] * d[1]) + (d[2] * d[2]);
					if (r2 <= 0) continue;
					const double c1 = (nq3[q][0] * d[0]) + (nq3[q][1] * d[1]) + (nq3[q][2] * d[2]);
					const double c2 = -((nq3[r][0] * d[0]) + (nq3[r][1] * d[1]) + (nq3[r][2] * d[2]));
					if ((c1 <= 0) || (c2 <= 0)) continue;
					sum += wa[q] * wa[r] * c1 * c2 / (PI * r2 * r2);
				}
			}
			const double den = (aa * aa) - w2;
			if (den > 0) f[(i * n) + i] = (sum * (aa * aa) / den) / RX.area[i];
		}
	}
	return rx_finish(fp_log, f, nvis);
}


// 非構造格子の面の平均の σT⁴ (輻射と同じ積分点 : 三角形の 3 点、四角形の 2x2 Gauss)
static void rx_emean_tet(const double *tn)
{
	const double gp[2] = {0.5 - 0.5 / sqrt(3.0), 0.5 + 0.5 / sqrt(3.0)};
	for (int f = 0; f < TFaceN; f++) {
		const int i = RX.fid[f];
		if (i < 0) continue;
		const int32_t *nd = &TFace[f * 4];
		double es = 0, ws = 0;
		if (TetOrder >= 2) {
			// 高次の要素の曲がった面 : 面の媒介変数の積分点 (p + 2 点) で
			const int e = TFaceEl[f];
			const int k = ((nd[3] < 0) ? 3 : 4);
			int32_t en[ELEM_NMAX];
			elem3d_nodes(e, en);
			double st[QTRI_MAX][2], wq[QTRI_MAX];
			const int nq = face_quad_pub(k, TetOrder + 2, st, wq);
			for (int q = 0; q < nq; q++) {
				double nv[ELEM_NMAX], x[3], t1[3], t2[3], tg = 0;
				const int ten = elem3d_face_point(e, nd, k, st[q][0], st[q][1], nv, x, t1, t2);
				const double cx = (t1[1] * t2[2]) - (t1[2] * t2[1]);
				const double cy = (t1[2] * t2[0]) - (t1[0] * t2[2]);
				const double cz = (t1[0] * t2[1]) - (t1[1] * t2[0]);
				const double w = wq[q] * sqrt((cx * cx) + (cy * cy) + (cz * cz));
				for (int a = 0; a < ten; a++) tg += nv[a] * tn[en[a]];
				const double tk = tg + KELVIN;
				es += w * STEFAN * tk * tk * tk * tk;
				ws += w;
			}
		}
		else if (nd[3] < 0) {
			for (int gq = 0; gq < 3; gq++) {
				double tg = 0;
				for (int m = 0; m < 3; m++) tg += ((m == gq) ? (2.0 / 3) : (1.0 / 6)) * tn[nd[m]];
				const double tk = tg + KELVIN;
				es += STEFAN * tk * tk * tk * tk;
				ws += 1;
			}
		}
		else {
			// 四角形は双 1 次の写像の |J| を重みに (surf_quad と同じ点・同じ重み)
			for (int ga = 0; ga < 2; ga++) {
			for (int gb = 0; gb < 2; gb++) {
				const double x = gp[ga], y = gp[gb];
				const double nn[4] = {(1 - x) * (1 - y), x * (1 - y), x * y, (1 - x) * y};
				const double dx[4] = {-(1 - y), (1 - y), y, -y};
				const double dy[4] = {-(1 - x), -x, x, (1 - x)};
				double a[3] = {0, 0, 0}, b[3] = {0, 0, 0}, tg = 0;
				for (int m = 0; m < 4; m++) {
					const double p[3] = {Xp[nd[m]], Yp[nd[m]], Zp[nd[m]]};
					for (int c = 0; c < 3; c++) {
						a[c] += dx[m] * p[c];
						b[c] += dy[m] * p[c];
					}
					tg += nn[m] * tn[nd[m]];
				}
				const double cx = (a[1] * b[2]) - (a[2] * b[1]);
				const double cy = (a[2] * b[0]) - (a[0] * b[2]);
				const double cz = (a[0] * b[1]) - (a[1] * b[0]);
				const double w = sqrt((cx * cx) + (cy * cy) + (cz * cz));
				const double tk = tg + KELVIN;
				es += w * STEFAN * tk * tk * tk * tk;
				ws += w;
			}
			}
		}
		RX.e[i] = ((ws > 0) ? (es / ws) : 0);
	}
}


// いまの温度で面の平均の E を作り、c_i = Σ_{j≠i} Q_ij E_j + q0_i を更新する
static void rx_update(const double *tn, const double *kcell)
{
	const int n = RX.n;
	if (MeshMode) rx_emean_tet(tn);
	else          surf_faces_grid(NULL, NULL, tn, kcell, 5, NULL);
	for (int i = 0; i < n; i++) {
		double v = RX.q0[i];
		for (int j = 0; j < n; j++) {
			if (j != i) v += RX.q[(i * n) + j] * RX.e[j];
		}
		RX.c[i] = v;
	}
}


// ---- 要素ごとの係数の組み立て・発熱・質量 (四面体) ----

// 四面体以外の 1 次要素 (六面体・角柱・ピラミッド) : 剛性と同じ積分点 (elem3d_ip) で
// K_ij = Σ w c ∇N_i・∇N_j、発熱 f_i = Σ w σ|∇φ|² N_i、質量 M_ij = Σ w c N_i N_j、
// 要素の温度 = Σ w N_i T_i / Σ w (1 次補間の体積平均)
// 積分点の経路 (gen_*) で組む最初の要素。1 次の四面体は閉形式の経路 (tet_*) のまま
// (従来の答えとバイト一致)、高次の四面体は積分点の経路に入れる
static int gen_first(void)
{
	return ((TetOrder >= 2) ? 0 : NTet);
}


static void gen_assemble(crs_t *A, const double *coef, int mass)
{
	const int ne = elem3d_count();
	for (int e = gen_first(); e < ne; e++) {
		const double c = coef[e];
		if (c <= 0) continue;
		int32_t nd[ELEM_NMAX];
		const int nen = elem3d_nodes(e, nd);
		const int nq = elem3d_nip(e);
		double ke[ELEM_NMAX][ELEM_NMAX];
		for (int l = 0; l < nen; l++) {
			for (int m = 0; m < nen; m++) ke[l][m] = 0;
		}
		int bad = 0;
		for (int q = 0; q < nq; q++) {
			double nv[ELEM_NMAX], g[ELEM_NMAX][3], w;
			if (elem3d_ip(e, q, nv, g, &w)) { bad = 1; break; }
			for (int l = 0; l < nen; l++) {
				for (int m = 0; m < nen; m++) {
					ke[l][m] += (mass ? (w * c * nv[l] * nv[m])
					           : (w * c * ((g[l][0] * g[m][0]) + (g[l][1] * g[m][1]) + (g[l][2] * g[m][2]))));
				}
			}
		}
		if (bad) continue;
		for (int l = 0; l < nen; l++) {
			for (int m = 0; m < nen; m++) {
				const int64_t p = crs_find_tri(A, nd[l], nd[m]);
				if (p >= 0) A->val[p] += ke[l][m];
			}
		}
	}
}


static void gen_joule(const double *phi, const double *sig, double *f)
{
	const int ne = elem3d_count();
	for (int e = gen_first(); e < ne; e++) {
		const double sg = sig[e];
		if (sg <= 0) continue;
		int32_t nd[ELEM_NMAX];
		const int nen = elem3d_nodes(e, nd);
		const int nq = elem3d_nip(e);
		for (int q = 0; q < nq; q++) {
			double nv[ELEM_NMAX], g[ELEM_NMAX][3], w;
			if (elem3d_ip(e, q, nv, g, &w)) break;
			double ex = 0, ey = 0, ez = 0;
			for (int l = 0; l < nen; l++) {
				ex += phi[nd[l]] * g[l][0];
				ey += phi[nd[l]] * g[l][1];
				ez += phi[nd[l]] * g[l][2];
			}
			const double qg = sg * ((ex * ex) + (ey * ey) + (ez * ez)) * w;
			for (int l = 0; l < nen; l++) f[nd[l]] += qg * nv[l];
		}
	}
}


// 要素 e (四面体以外) の体積と、温度場 tn の要素平均
static double gen_mean(int e, const double *tn, double *vol)
{
	int32_t nd[ELEM_NMAX];
	const int nen = elem3d_nodes(e, nd);
	const int nq = elem3d_nip(e);
	double v = 0, ts = 0;
	for (int q = 0; q < nq; q++) {
		double nv[ELEM_NMAX], g[ELEM_NMAX][3], w;
		if (elem3d_ip(e, q, nv, g, &w)) break;
		v += w;
		if (tn) {
			for (int l = 0; l < nen; l++) ts += w * nv[l] * tn[nd[l]];
		}
	}
	*vol = v;
	return ((v > 0) ? (ts / v) : 0);
}


// 断面 2 次元 (1 次三角形) : 剛性 K_ij = c S ∇λ_i・∇λ_j、質量 M_ij = c S (1 + δ_ij)/12
static void tri2d_assemble(crs_t *A, const double *coef, int mass)
{
	crs_zero(A);
	for (int e = 0; e < NTri; e++) {
		const double c = coef[e];
		if (c <= 0) continue;
		const int32_t *nd = &Tri[e * 3];
		double g[3][2], area;
		if (tri_grad(nd, g, &area)) continue;
		for (int l = 0; l < 3; l++) {
			for (int m = 0; m < 3; m++) {
				const double v = (mass ? (c * area * ((l == m) ? 2 : 1) / 12)
				               : (c * ((g[l][0] * g[m][0]) + (g[l][1] * g[m][1])) * area));
				const int64_t p = crs_find_tri(A, nd[l], nd[m]);
				if (p >= 0) A->val[p] += v;
			}
		}
	}
}


static void tet_assemble(crs_t *A, const double *coef)
{
	if (MeshDim == 2) {
		tri2d_assemble(A, coef, 0);
		return;
	}
	crs_zero(A);
	gen_assemble(A, coef, 0);
	for (int e = 0; (TetOrder < 2) && (e < NTet); e++) {
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
	gen_joule(phi, sig, f);
	for (int e = 0; (TetOrder < 2) && (e < NTet); e++) {
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
	if (MeshDim == 2) {
		tri2d_assemble(A, coef, 1);
		return;
	}
	crs_zero(A);
	gen_assemble(A, coef, 1);
	for (int e = 0; (TetOrder < 2) && (e < NTet); e++) {
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
static double cond_rho_factor(int p, double t)
{
	return ((CondNTemp[p] > 0) ? table_interp(CondNTemp[p], CondTempT[p], CondTempF[p], t)
	        : (1 + (CondTempco[p] * (t - CondTemp0[p]))));
}

static int th_sigma(FILE *fp_log, int64_t nc, const double *tcell, double *scell)
{
	if (MeshMode && (MeshDim == 2)) {
		// 断面 2 次元 (F) の σ は導体ごとの CondSigma (conductorsigma / conductortempco)。
		// 入力解釈で周囲温度の補正が掛かっているので外して掛け直す
		for (int64_t c = 0; c < nc; c++) {
			const int id = TriCond[c];
			if ((id < 0) || (CondSigma[id] <= 0)) {
				scell[c] = 0;
				continue;
			}
			const double den = cond_rho_factor(id, tcell[c]);
			if (den <= 0) {
				fprintf(fp_log, "*** thermal : the resistivity factor of conductor %d became "
					"non-positive (T = %.4e degC)\n", id, tcell[c]);
				return 1;
			}
			scell[c] = CondSigma[id] * cond_rho_factor(id, Temperature) / den;
		}
		return 0;
	}
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
	return ((Analysis & (ANALYSIS_A | ANALYSIS_F)) != 0);
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
		if (MeshMode && (MeshDim == 2)) {
			if (eddy2d_heat_solve(fp_log, scell, cur, f, &zr, &zi, it)) return 1;
		}
		else if (MeshMode && (TetOrder == 2) && ((NHex + NPrism + NPyr) == 0)) {
			if (eddy3d_p2_heat_solve(fp_log, scell, cur, f, &zr, &zi, it)) return 1;
		}
		else if (eddy3d_heat_solve(fp_log, scell, cur, f, &zr, &zi, it)) return 1;
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
	if (MeshMode && (MeshDim == 2)) {
		// 断面 2 次元 : 3 節点の平均 (1 次補間の面積平均)
		for (int e = 0; e < NTri; e++) {
			if (kcell[e] <= 0) continue;
			double t = 0;
			for (int l = 0; l < 3; l++) {
				const double tl = tn[Tri[(e * 3) + l]];
				t += tl / 3;
				if (tl > tmax) tmax = tl;
			}
			tcell[e] = t;
		}
		return tmax;
	}
	if (MeshMode) {
		// 四面体 : 4 節点の平均 (1 次補間の体積平均)
		for (int e = 0; (TetOrder < 2) && (e < NTet); e++) {
			if (kcell[e] <= 0) continue;
			double t = 0;
			for (int l = 0; l < 4; l++) {
				const double tl = tn[Tet[(e * 4) + l]];
				t += tl / 4;
				if (tl > tmax) tmax = tl;
			}
			tcell[e] = t;
		}
		const int ne = elem3d_count();
		for (int e = gen_first(); e < ne; e++) {
			if (kcell[e] <= 0) continue;
			double v;
			tcell[e] = gen_mean(e, tn, &v);
			int32_t nd[ELEM_NMAX];
			const int nen = elem3d_nodes(e, nd);
			for (int l = 0; l < nen; l++) {
				if (tn[nd[l]] > tmax) tmax = tn[nd[l]];
			}
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
	if (MeshMode && (MeshDim == 2)) {
		for (int e = 0; e < NTri; e++) {
			if (scell[e] <= 0) continue;
			vsum += TriArea[e];
			tsum += TriArea[e] * tcell[e];
		}
		return ((vsum > 0) ? (tsum / vsum) : Temperature);
	}
	if (MeshMode) {
		for (int e = 0; (TetOrder < 2) && (e < NTet); e++) {
			if (scell[e] <= 0) continue;
			double gr[4][3], v;
			if (tet_grad_pub(&Tet[e * 4], gr, &v)) continue;
			vsum += v;
			tsum += v * tcell[e];
		}
		const int ne = elem3d_count();
		for (int e = gen_first(); e < ne; e++) {
			if (scell[e] <= 0) continue;
			double v;
			gen_mean(e, NULL, &v);
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
	if (RX.on) rx_update(tn, kcell);	// 相互放射 : 他の面からの項を今の温度で
	if ((RadEps > 0) && (surf_faces(B, gr, tn, kcell, 1, NULL) < -1e299)) {
		fprintf(fp_log, "*** thermal : the absolute temperature became non-positive "
			"on a radiating face (iteration %d)\n", it);
		return 1;
	}
	// h(T) の対流も同じく線形化する (表が無ければ対流は b0 と gc に入っている)
	if (th_convnl()) surf_faces(B, gr, tn, kcell, 3, NULL);
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

	FILE *fcsv = mpi_fopen_out(FN_transient, "w");
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
	if (!ierr && th_nonlin()) ierr = th_radlin(fp_log, B, b0, gr, tn, kcell, 0);
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
		qconv = ((ConvH > 0) ? surf_faces(NULL, NULL, tn, kcell, (th_convnl() ? 3 : 0), NULL) : 0);
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
			if (th_nonlin() && ((ierr = th_radlin(fp_log, B, b0, gr, tn, kcell, tag)) != 0)) break;

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
			const int id = ThNodeCond[i];
			if ((id >= 0) && HeatSinkOn[id]) {
				const double dm = crs_row_dot(&M, i, tn) - mtold[i];
				qsink -= (dm / dts) + (th * (crs_row_dot(B, i, tn) - f[i] - gc[i] - gr[i]))
					+ ((1 - th) * resold[i]);
			}
		}
		crs_spmv(B, tn, work, NULL);
		for (int i = 0; i < n; i++) resold[i] = work[i] - f[i] - gc[i] - gr[i];
		qconv = ((ConvH > 0) ? surf_faces(NULL, NULL, tn, kcell, (th_convnl() ? 3 : 0), NULL) : 0);
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
		if (th_nonlin() && ((ierr = th_radlin(fp_log, B, b0, gr, tn, kcell, it)) != 0)) break;
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
			const int id = ThNodeCond[i];
			if ((id >= 0) && HeatSinkOn[id]) {
				qsink += f[i] + gc[i] + gr[i] - crs_row_dot(B, i, tn);
			}
		}
		qconv = ((ConvH > 0) ? surf_faces(NULL, NULL, tn, kcell, (th_convnl() ? 3 : 0), NULL) : 0);
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
	// 10 節点四面体だけの格子の交流の発熱は 2 次の Nedelec (eddy3d_p2.c) で解く
	const int p2tet = (MeshMode && (MeshDim == 3) && (TetOrder == 2) && ((NHex + NPrism + NPyr) == 0)
		&& (Analysis & ANALYSIS_A));

	fprintf(fp_log, "\n=== electro-thermal analysis ===\n");
	// 非構造格子は 1 次の 3 次元要素 (四面体・六面体・角柱・ピラミッド、混在可)、高次の四面体、
	// 2 次の六面体・角柱 (いずれも DC、相互放射なし)。それ以外は黙って誤答にせず入力エラー
	// 断面 2 次元は交流の発熱 (F) の 1 次三角形だけ (F は断面 2 次元でしか解けない)
	if (MeshMode && (MeshDim == 2) && (Analysis & ANALYSIS_F) && (TetOrder < 2)) {
		if (RadExchange) {
			fprintf(fp_log, "*** radiationexchange is not implemented on a 2-D cross-section mesh\n");
			return 1;
		}
	}
	else if ((Analysis & ANALYSIS_F) || (MeshMode && ((MeshDim != 3)
	         || ((TetOrder >= 2) && ((Analysis & ANALYSIS_A) && !p2tet))))) {
		if (Analysis & ANALYSIS_F) {
			fprintf(fp_log, "*** AC heating with analysis F (thermal + F) needs a 2-D cross-section "
				"mesh of 3-node triangles\n");
			return 1;
		}
		fprintf(fp_log, "*** thermal on an unstructured mesh supports first-order 3-D elements "
			"and higher-order tetrahedra / hexahedra / prisms (DC) only (this mesh is %d-D, element "
			"order %d, %d-node hexahedra, %d-node prisms)\n",
			MeshDim, TetOrder, HexNen, PrismNen);
		return 1;
	}
	// 交流の発熱 (A) は辺要素の種別 (四面体・六面体・角柱、混在可) ならどれでも。損失は
	// edge_elem_matrices の T_e の上の ½ eᴴ T_e e なので種別に依らない。ピラミッドは A 自体が弾く
	if (MeshMode && th_ac() && (MeshDim == 3) && (NPyr > 0)) {
		fprintf(fp_log, "*** AC heating (thermal + analysis A) does not support pyramids "
			"(no polynomial Nedelec basis; this mesh has %d pyramids)\n", NPyr);
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
	for (int k = 0; k < 2; k++) {
		if (ProfN[k] < 1) continue;
		fprintf(fp_log, "  %s(r) = %s f(%c), f :", (k ? "emissivity e" : "convection h"), (k ? "e" : "h"),
			'x' + ProfAxis[k]);
		for (int i = 0; i < ProfN[k]; i++) fprintf(fp_log, " (%g, %g)", ProfX[k][i], ProfF[k][i]);
		fprintf(fp_log, " (at the face centre)\n");
	}
	for (int k = 0; k < 2; k++) {
		if (!ExprOn[k]) continue;
		fprintf(fp_log, "  %s(r) = %s x f(x, y, z), f = %s (at the face centre)\n",
			(k ? "emissivity e" : "convection h"), (k ? "e" : "h"), ExprStr[k]);
	}
	if ((RadEps > 0) && (RadTabN > 0)) {
		fprintf(fp_log, "  emissivity e(T) = e f(T), f :");
		for (int k = 0; k < RadTabN; k++) fprintf(fp_log, " (%g, %g)", RadTabT[k], RadTabF[k]);
		fprintf(fp_log, "\n");
	}
	if (th_convnl()) {
		fprintf(fp_log, "  convection h(T) = h f(T), f :");
		for (int k = 0; k < ConvTabN; k++) fprintf(fp_log, " (%g, %g)", ConvTabT[k], ConvTabF[k]);
		fprintf(fp_log, " (Newton-linearized each iteration)\n");
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
	TFaceEl = NULL;
	TOcc = NULL;
	TOccEl = NULL;
	TFaceN = 0;
	TOccN = 0;
	if (MeshMode) tet_faces_build(kcell);
	signed char *nc2d = NULL;
	ThNodeCond = NodeConductor;
	if (MeshMode && (MeshDim == 2)) {
		// 断面 2 次元 : ヒートシンクの導体は、その導体の三角形の節点を固定する
		nc2d = (signed char *)malloc((size_t)n * sizeof(signed char));
		memset(nc2d, -1, (size_t)n * sizeof(signed char));
		for (int e = 0; e < NTri; e++) {
			if (TriCond[e] < 0) continue;
			for (int l = 0; l < 3; l++) nc2d[Tri[(e * 3) + l]] = TriCond[e];
		}
		ThNodeCond = nc2d;
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
		const int id = ThNodeCond[i];
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
	if ((ConvH > 0) && !th_convnl()) {
		surf_faces(&B, gc, NULL, kcell, 0, NULL);
	}
	if ((NConv + NRad) > 0) {
		// 係数が正の露出面の数。直方体を付けた行は、露出面に 1 つも当たらなければ
		// 書き間違い (座標の桁違いなど) なので黙って無視せず入力エラーにする
		int *hits = (int *)calloc((size_t)(NConv + NRad), sizeof(int));
		ExprBad[0] = ExprBad[1] = 0;		// 数えるのはこの 1 回の走査の面 (組み立ての走査と二重に数えない)
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
		for (int k = 0; k < 2; k++) {
			if (ExprOn[k] && (ExprBad[k] > 0)) {
				fprintf(fp_log, "*** thermal : %s = %s is negative or not finite%s on %d exposed faces\n",
					(k ? "radiationexpr" : "convectionexpr"), ExprStr[k], (k ? " (or e f > 1)" : ""), ExprBad[k]);
				ierr = 1;
			}
		}
		fprintf(fp_log, "  %s : %lld exposed faces\n",
			((NConv > 0) ? ((NRad > 0) ? "convection / radiation" : "convection") : "radiation"),
			(long long)nface);
	}
	if (!TransientOn && (nsink < 1) && (nface < 1)) {
		fprintf(fp_log, "*** thermal : the heat has no way out (the heat-sink conductors "
			"have no node and no face is exposed to convection / radiation)\n");
		ierr = 1;
	}
	// 面どうしの相互放射 : 形態係数と放射度の行列を一度だけ作る
	memset(&RX, 0, sizeof(RX));
	if (!ierr && RadExchange && (RadEps > 0)) ierr = (MeshMode ? rx_build_tet(fp_log) : rx_build(fp_log, kcell));
	// 輻射は温度に依るので、反復ごとに伝導 + 対流の行列 (b0) に接線を足し直す
	const int64_t nnz = B.rowptr[n];
	double *b0 = NULL;
	if (th_nonlin()) {
		b0 = (double *)malloc((size_t)nnz * sizeof(double));
		for (int64_t e = 0; e < nnz; e++) b0[e] = B.val[e];
	}

	// 交流 : A-φ の辺の情報と行列は反復の間持っておく
	int ehopen = 0;
	if (!ierr && th_ac() && !(MeshMode && (MeshDim == 2))) {
		ierr = (p2tet ? eddy3d_p2_heat_open(fp_log) : eddy3d_heat_open(fp_log));
		ehopen = 1;
	}
	if (!ierr) {
		ierr = (TransientOn ? th_transient : th_steady)(fp_log, &A, &B, b0, kcell, tcell,
			scell, fix, fixt, phi, work, diag, f, tn, rhs, u, gc, gr, q, gport);
	}
	if (ehopen) {
		if (p2tet) eddy3d_p2_heat_close();
		else eddy3d_heat_close();
	}
	*xport = ThX;
	fflush(fp_log);

	crs_free(&A);
	crs_free(&B);
	free(tcell); free(scell); free(kcell);
	free(fix); free(fixt); free(phi); free(work); free(diag);
	free(f); free(tn); free(rhs); free(u); free(gc); free(gr); free(q);
	free(b0);
	free(TFace); free(TFaceEl); free(TOcc); free(TOccEl); free(nc2d);
	TFace = NULL;
	rx_free();

	return ierr;
}
