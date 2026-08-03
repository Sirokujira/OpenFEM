/*
setup.c

格子・材料分布・導体 (電極) 節点の設定。

・セル中心が形状の内部にあれば、そのセルに材料番号を与える (後の geometry が優先)。
・節点が conductor 形状の内部にあれば、その節点を導体番号で固定する
  (Dirichlet 境界)。厚さ 0 の指定は電極面として扱える。
*/

#include "fem.h"
#include "fem_prototype.h"

int64_t node_index(int i, int j, int k)
{
	return (((int64_t)i * (Ny + 1)) + j) * (Nz + 1) + k;
}


int64_t num_node(void)
{
	return (MeshMode ? (int64_t)NNode : ((int64_t)(Nx + 1) * (Ny + 1) * (Nz + 1)));
}


// 結果行列の確保 (格子の種類に依らない)
static void alloc_matrices(void)
{
	const size_t msize = (size_t)NPort * NPort * sizeof(double);

	Cmat = (double *)malloc(msize);
	Lmat = (double *)malloc(msize);
	Gmat = (double *)malloc(msize);
	Rmat = (double *)malloc(msize);
	Mmat = (double *)malloc(msize);
	Smat = (double *)malloc(msize);
	Rfmat = (double *)malloc(msize);
	Lfmat = (double *)malloc(msize);
	Pfemat = (double *)malloc(msize);
	memset(Cmat, 0, msize);
	memset(Lmat, 0, msize);
	memset(Gmat, 0, msize);
	memset(Rmat, 0, msize);
	memset(Mmat, 0, msize);
	memset(Smat, 0, msize);
	memset(Rfmat, 0, msize);
	memset(Lfmat, 0, msize);
	memset(Pfemat, 0, msize);
	HaveC = HaveL = HaveR = HaveM = HaveS = HaveF = HavePfe = 0;
}


// 断面 2 次元の非構造格子 (三角形) のセットアップ
//
// M / F は伝送線路軸 t に垂直な断面での 2 次元問題なので、四面体ではなく
// 三角形で切った格子に載る。三角形が体積要素になり、`region` が材料、
// `electrode` が**導体断面**の物理タグを指す (3 次元格子の「電極面」に
// あたるものは 2 次元では存在しない)。
//
// 単位長あたりの量として扱うので TlineLength = 1 とし、以降の
// 「体積」はすべて面積になる。
static int setup_tri2d(void)
{
	// 断面 2 次元の格子で解けるのは M / F だけ
	if (Analysis & ~(ANALYSIS_M | ANALYSIS_F)) {
		printf("%s\n", "*** a 2-D (cross-section) mesh supports only analysis M and F "
			"(C / L / R / E / A / P need a 3-D tetrahedral mesh)");
		return 1;
	}
	// 面は伝送線路軸に垂直でなければならない (その軸の座標が一定)
	if (!Tline) {
		printf("%s\n", "*** a 2-D (cross-section) mesh needs the tline key");
		return 1;
	}
	// 非線形 (B-H) とヒステリシス (J-A) は Gauss 点毎の状態を構造格子の
	// セル配列で持っているので、2 次元格子では使えない
	for (int m = 0; m < NMaterial; m++) {
		if ((Material[m].nbh[0] > 0) || Material[m].ja.on) {
			printf("%s\n", "*** nonlinear (bh) and hysteresis (ja) materials are "
				"available only on a structured mesh");
			return 1;
		}
	}
	{
		const double *c = ((Tline == 'X') ? Xp : (Tline == 'Y') ? Yp : Zp);
		double lo = c[0], hi = c[0], span = 0;
		for (int i = 1; i < NNode; i++) {
			if (c[i] < lo) lo = c[i];
			if (c[i] > hi) hi = c[i];
		}
		int p, q;
		tri_axes(&p, &q);
		const double *cp = ((p == 0) ? Xp : (p == 1) ? Yp : Zp);
		const double *cq = ((q == 0) ? Xp : (q == 1) ? Yp : Zp);
		double plo = cp[0], phi = cp[0], qlo = cq[0], qhi = cq[0];
		for (int i = 1; i < NNode; i++) {
			if (cp[i] < plo) plo = cp[i];
			if (cp[i] > phi) phi = cp[i];
			if (cq[i] < qlo) qlo = cq[i];
			if (cq[i] > qhi) qhi = cq[i];
		}
		span = (phi - plo) + (qhi - qlo);
		if ((hi - lo) > (EPS * span)) {
			printf("*** a 2-D mesh must lie in a plane normal to tline = %c "
				"(the %c coordinate spans %.4e [m])\n", Tline, Tline, hi - lo);
			return 1;
		}
	}

	// 物理タグ -> 材料番号 / 導体番号 (どちらも三角形のタグ)
	TriMat = (unsigned char *)malloc((size_t)NTri * sizeof(unsigned char));
	TriCond = (signed char *)malloc((size_t)NTri * sizeof(signed char));
	TriArea = (double *)malloc((size_t)NTri * sizeof(double));
	for (int e = 0; e < NTri; e++) {
		int m = 0, id = -1;
		for (int q = 0; q < NRegion; q++) {
			if (TriTag[e] == RegionTag[q]) m = RegionMat[q];
		}
		for (int q = 0; q < NElectrode; q++) {
			if (TriTag[e] == ElecTag[q]) id = ElecCond[q];
		}
		TriMat[e] = (unsigned char)m;
		TriCond[e] = (signed char)id;

		double g[3][2], area;
		if (tri_grad(&Tri[e * 3], g, &area)) {
			printf("*** triangle %d is degenerate (zero area)\n", e + 1);
			return 1;
		}
		// 2 次要素では等パラメトリック写像で積分した面積を使う
		// (頂点 3 個から作る面積は曲がった要素では内接多角形の値になる)
		if (TetOrder >= 2) {
			area = tri6_area(e);
			if (area <= 0) {
				printf("*** triangle %d is degenerate or inverted\n", e + 1);
				return 1;
			}
		}
		TriArea[e] = area;
	}

	// 導体の断面積。2 次元では線路長 1 m あたりで扱う
	for (int p = 0; p < MAXPORT; p++) CondArea[p] = 0;
	for (int e = 0; e < NTri; e++) {
		const int id = TriCond[e];
		if (id >= 0) CondArea[id] += TriArea[e];
	}
	for (int p = 0; p <= NPort; p++) {
		if (CondArea[p] <= 0) {
			printf("*** conductor %d has no triangle (check the physical tags)\n", p);
			return 1;
		}
	}

	// 節点の Dirichlet は使わない (M / F は自然境界条件)
	NodeConductor = (signed char *)malloc((size_t)NNode * sizeof(signed char));
	memset(NodeConductor, -1, (size_t)NNode * sizeof(signed char));

	TlineLength = 1;
	if (LineLength <= 0) LineLength = 1;

	alloc_matrices();

	// 導体の DC 直列抵抗 [ohm/m] (3 次元格子と同じ式。面積は単位長あたり)
	{
		int have = 0;
		for (int p = 0; p <= NPort; p++) {
			if (CondSigma[p] > 0) have = 1;
		}
		if (have) {
			const double r0 = ((CondSigma[0] > 0) ? (1 / (CondSigma[0] * CondArea[0])) : 0);
			for (int k = 1; k <= NPort; k++) {
				const double rk = ((CondSigma[k] > 0) ? (1 / (CondSigma[k] * CondArea[k])) : 0);
				for (int j = 1; j <= NPort; j++) {
					Smat[((k - 1) * NPort) + (j - 1)] = r0 + ((k == j) ? rk : 0);
				}
			}
			HaveS = 1;
		}
	}

	return 0;
}


// 非構造格子のセットアップ
static int setup_unstruct(void)
{
	if (mesh_read(MeshFile)) return 1;

	// 辺要素 (E / A) は 1 次四面体の Whitney 形状関数に基づくので、2 次格子を
	// 渡されると辺上の中間節点がどの要素にも現れず、節点ブロックが特異になる。
	// 黙って誤答を出すより弾く
	// 辺要素 (Whitney) も節点要素の自己検証も四面体の形状関数に基づく。
	// 六面体格子を渡されたら黙って誤答を出さずに弾く
	if ((MeshElem != MESHELEM_TET) && (Analysis & (ANALYSIS_E | ANALYSIS_A))) {
		printf("%s\n", "*** analysis E / A (edge elements) need a tetrahedral mesh "
			"(this mesh has hexahedra or prisms)");
		return 1;
	}
	if ((TetOrder >= 2) && (Analysis & (ANALYSIS_E | ANALYSIS_A))) {
		printf("%s\n", "*** analysis E / A (edge elements) need a first-order mesh "
			"(this mesh has 10-node tetrahedra)");
		return 1;
	}

	if (MeshDim == 2) return setup_tri2d();

	// 物理タグ -> 材料番号
	if (NHex > 0) {
		HexMat = (unsigned char *)malloc((size_t)NHex * sizeof(unsigned char));
		for (int e = 0; e < NHex; e++) {
			int m = 0;
			for (int q = 0; q < NRegion; q++) {
				if (HexTag[e] == RegionTag[q]) m = RegionMat[q];
			}
			HexMat[e] = (unsigned char)m;
		}
	}
	if (NPrism > 0) {
		PrismMat = (unsigned char *)malloc((size_t)NPrism * sizeof(unsigned char));
		for (int e = 0; e < NPrism; e++) {
			int m = 0;
			for (int q = 0; q < NRegion; q++) {
				if (PrismTag[e] == RegionTag[q]) m = RegionMat[q];
			}
			PrismMat[e] = (unsigned char)m;
		}
	}
	if (NTet > 0) {
		TetMat = (unsigned char *)malloc((size_t)NTet * sizeof(unsigned char));
		for (int e = 0; e < NTet; e++) {
			int m = 0;
			for (int q = 0; q < NRegion; q++) {
				if (TetTag[e] == RegionTag[q]) m = RegionMat[q];
			}
			TetMat[e] = (unsigned char)m;
		}
	}

	// 物理タグ -> 電極 (三角形の節点を Dirichlet にする)
	// 2 次格子では辺上の中間節点も固定する。1 次では Tri2 に頂点が入っているので
	// 同じ値を二度塗るだけになり無害
	NodeConductor = (signed char *)malloc((size_t)NNode * sizeof(signed char));
	memset(NodeConductor, -1, (size_t)NNode * sizeof(signed char));
	// 六面体格子の境界面は四角形になる
	for (int t = 0; t < NQuad; t++) {
		for (int q = 0; q < NElectrode; q++) {
			if (QuadTag[t] != ElecTag[q]) continue;
			for (int l = 0; l < 4; l++) {
				NodeConductor[Quad[(t * 4) + l]] = (signed char)ElecCond[q];
			}
		}
	}
	for (int t = 0; t < NTri; t++) {
		for (int q = 0; q < NElectrode; q++) {
			if (TriTag[t] != ElecTag[q]) continue;
			for (int l = 0; l < 3; l++) {
				NodeConductor[Tri[(t * 3) + l]] = (signed char)ElecCond[q];
				NodeConductor[Tri2[(t * 3) + l]] = (signed char)ElecCond[q];
			}
		}
	}

	int64_t count[MAXPORT];
	for (int p = 0; p < MAXPORT; p++) count[p] = 0;
	for (int i = 0; i < NNode; i++) {
		const int id = NodeConductor[i];
		if (id >= 0) count[id]++;
	}
	for (int p = 0; p <= NPort; p++) {
		if (count[p] == 0) {
			printf("*** electrode %d has no node (check the physical tags)\n", p);
			return 1;
		}
	}

	// 伝送線路長 (節点の外接直方体から取る)
	TlineLength = 0;
	if (Tline) {
		const double *c = ((Tline == 'X') ? Xp : (Tline == 'Y') ? Yp : Zp);
		double lo = c[0], hi = c[0];
		for (int i = 1; i < NNode; i++) {
			if (c[i] < lo) lo = c[i];
			if (c[i] > hi) hi = c[i];
		}
		TlineLength = hi - lo;
	}
	if (LineLength <= 0) LineLength = TlineLength;

	for (int p = 0; p < MAXPORT; p++) CondArea[p] = 0;
	alloc_matrices();

	return 0;
}


int setup(void)
{
	if (MeshMode) return setup_unstruct();

	// 節点番号は int32 で扱う (CRS の列番号・OpenMP のループ変数)
	if (num_node() >= INT_MAX) {
		printf("*** too many nodes (%lld); reduce the mesh division\n", (long long)num_node());
		return 1;
	}

	// セル中心

	Xc = (double *)malloc(Nx * sizeof(double));
	Yc = (double *)malloc(Ny * sizeof(double));
	Zc = (double *)malloc(Nz * sizeof(double));
	for (int i = 0; i < Nx; i++) Xc[i] = (Xn[i] + Xn[i + 1]) / 2;
	for (int j = 0; j < Ny; j++) Yc[j] = (Yn[j] + Yn[j + 1]) / 2;
	for (int k = 0; k < Nz; k++) Zc[k] = (Zn[k] + Zn[k + 1]) / 2;

	// 形状判定の許容幅 (解析領域の大きさに対する相対値)
	const double size = fabs(Xn[Nx] - Xn[0]) + fabs(Yn[Ny] - Yn[0]) + fabs(Zn[Nz] - Zn[0]);
	const double eps = EPS * size;

	// セル材料

	const int64_t ncell = (int64_t)Nx * Ny * Nz;
	CellMaterial = (unsigned char *)malloc(ncell * sizeof(unsigned char));
	memset(CellMaterial, 0, ncell * sizeof(unsigned char));

	for (int n = 0; n < NGeometry; n++) {
		const int    m     = Geometry[n].m;
		const int    shape = Geometry[n].shape;
		const double *g    = Geometry[n].g;
		int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
		for (int k = 0; k < Nz; k++) {
			if (ingeometry(Xc[i], Yc[j], Zc[k], shape, g, eps)) {
				CellMaterial[((int64_t)i * Ny + j) * Nz + k] = (unsigned char)m;
			}
		}
		}
		}
	}

	// 導体セル (断面積・電流密度の計算に使う。節点と同じ順序で重ね塗りする)

	CellConductor = (signed char *)malloc(ncell * sizeof(signed char));
	memset(CellConductor, -1, ncell * sizeof(signed char));

	for (int n = 0; n < NConductor; n++) {
		const int    id    = Conductor[n].id;
		const int    shape = Conductor[n].shape;
		const double *g    = Conductor[n].g;
		int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
		for (int k = 0; k < Nz; k++) {
			if (ingeometry(Xc[i], Yc[j], Zc[k], shape, g, eps)) {
				CellConductor[((int64_t)i * Ny + j) * Nz + k] = (signed char)id;
			}
		}
		}
		}
	}

	// 導体節点

	const int64_t nnode = num_node();
	NodeConductor = (signed char *)malloc(nnode * sizeof(signed char));
	memset(NodeConductor, -1, nnode * sizeof(signed char));

	for (int n = 0; n < NConductor; n++) {
		const int    id    = Conductor[n].id;
		const int    shape = Conductor[n].shape;
		const double *g    = Conductor[n].g;
		int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (i = 0; i <= Nx; i++) {
		for (int j = 0; j <= Ny; j++) {
		for (int k = 0; k <= Nz; k++) {
			if (ingeometry(Xn[i], Yn[j], Zn[k], shape, g, eps)) {
				NodeConductor[node_index(i, j, k)] = (signed char)id;
			}
		}
		}
		}
	}

	// 導体毎の節点数を数え、空の導体があればエラーにする

	int64_t count[MAXPORT];
	for (int p = 0; p < MAXPORT; p++) count[p] = 0;
	for (int64_t n = 0; n < nnode; n++) {
		const int id = NodeConductor[n];
		if (id >= 0) count[id]++;
	}
	for (int p = 0; p <= NPort; p++) {
		if (count[p] == 0) {
			printf("*** conductor %d has no node (check its position and the mesh)\n", p);
			return 1;
		}
	}

	// 伝送線路長 (単位長あたりの出力に使う)

	TlineLength = 0;
	if      (Tline == 'X') TlineLength = Xn[Nx] - Xn[0];
	else if (Tline == 'Y') TlineLength = Yn[Ny] - Yn[0];
	else if (Tline == 'Z') TlineLength = Zn[Nz] - Zn[0];

	// 等価回路の線路長 (linelength 省略時は解析領域長)
	if (LineLength <= 0) LineLength = TlineLength;

	// 導体の断面積 (体積 / 線路長)。静磁場の電流密度と DC 直列抵抗に使う
	for (int p = 0; p < MAXPORT; p++) {
		CondArea[p] = 0;
	}
	if (TlineLength > 0) {
		for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
		for (int k = 0; k < Nz; k++) {
			const int id = CellConductor[((int64_t)i * Ny + j) * Nz + k];
			if (id < 0) continue;
			CondArea[id] += (Xn[i + 1] - Xn[i]) * (Yn[j + 1] - Yn[j]) * (Zn[k + 1] - Zn[k]);
		}
		}
		}
		for (int p = 0; p <= NPort; p++) {
			CondArea[p] /= TlineLength;
		}
	}

	// Jiles-Atherton の状態 (Gauss 点毎)

	int hysteresis = 0;
	for (int m = 0; m < NMaterial; m++) {
		if (Material[m].ja.on) hysteresis = 1;
	}
	JaB = JaH = JaM = JaBn = JaHn = JaMn = JaD = JaDn = NULL;
	if (hysteresis) {
		const size_t njasz = (size_t)ncell * 8 * sizeof(double);
		JaB  = (double *)malloc(njasz);
		JaH  = (double *)malloc(njasz);
		JaM  = (double *)malloc(njasz);
		JaBn = (double *)malloc(njasz);
		JaHn = (double *)malloc(njasz);
		JaMn = (double *)malloc(njasz);
		memset(JaB,  0, njasz);		// 初期状態 : 処女 (B = H = M = 0)
		memset(JaH,  0, njasz);
		memset(JaM,  0, njasz);
		memset(JaBn, 0, njasz);
		memset(JaHn, 0, njasz);
		memset(JaMn, 0, njasz);
		JaD  = (double *)malloc(njasz * 3);
		JaDn = (double *)malloc(njasz * 3);
		memset(JaD,  0, njasz * 3);		// 0 : 方向未確定 (未磁化)
		memset(JaDn, 0, njasz * 3);
	}

	// 結果行列

	alloc_matrices();

	// 導体の DC 直列抵抗 [ohm/m] (conductorsigma 指定時)
	// 帰路を基準導体が共有するので Smat[k][j] = R0 + (k==j ? Rk : 0) とする
	if (TlineLength > 0) {
		int have = 0;
		for (int p = 0; p <= NPort; p++) {
			if ((CondSigma[p] > 0) && (CondArea[p] > 0)) have = 1;
		}
		if (have) {
			const double r0 = ((CondSigma[0] > 0) && (CondArea[0] > 0))
				? (1 / (CondSigma[0] * CondArea[0])) : 0;
			for (int k = 1; k <= NPort; k++) {
				const double rk = ((CondSigma[k] > 0) && (CondArea[k] > 0))
					? (1 / (CondSigma[k] * CondArea[k])) : 0;
				for (int j = 1; j <= NPort; j++) {
					Smat[((k - 1) * NPort) + (j - 1)] = r0 + ((k == j) ? rk : 0);
				}
			}
			HaveS = 1;
		}
	}

	return 0;
}


void memfree(void)
{
	free(Xn); free(Yn); free(Zn);
	free(Xc); free(Yc); free(Zc);
	free(Material);
	free(Geometry);
	free(Conductor);
	free(Xp); free(Yp); free(Zp);
	free(Tet); free(TetTag); free(TetMat); free(Tet2);
	free(Tri); free(TriTag); free(Tri2); free(TriMat); free(TriCond); free(TriArea);
	free(CellMaterial);
	free(CellConductor);
	free(NodeConductor);
	free(Cmat); free(Lmat); free(Gmat); free(Rmat); free(Mmat); free(Smat);
	free(Rfmat); free(Lfmat); free(Pfemat);
	free(JaB); free(JaH); free(JaM); free(JaBn); free(JaHn); free(JaMn);
	free(JaD); free(JaDn);
}
