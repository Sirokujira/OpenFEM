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
	return (int64_t)(Nx + 1) * (Ny + 1) * (Nz + 1);
}


int setup(void)
{
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
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
		for (int k = 0; k < Nz; k++) {
			if (ingeometry(Xc[i], Yc[j], Zc[k], shape, g, eps)) {
				CellMaterial[((int64_t)i * Ny + j) * Nz + k] = (unsigned char)m;
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
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (int i = 0; i <= Nx; i++) {
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

	// 結果行列

	const size_t msize = (size_t)NPort * NPort * sizeof(double);
	Cmat = (double *)malloc(msize);
	Lmat = (double *)malloc(msize);
	Gmat = (double *)malloc(msize);
	Rmat = (double *)malloc(msize);
	memset(Cmat, 0, msize);
	memset(Lmat, 0, msize);
	memset(Gmat, 0, msize);
	memset(Rmat, 0, msize);
	HaveC = HaveL = HaveR = 0;

	return 0;
}


void memfree(void)
{
	free(Xn); free(Yn); free(Zn);
	free(Xc); free(Yc); free(Zc);
	free(Material);
	free(Geometry);
	free(Conductor);
	free(CellMaterial);
	free(NodeConductor);
	free(Cmat); free(Lmat); free(Gmat); free(Rmat);
}
