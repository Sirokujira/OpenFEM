/*
crs.c

CRS 形式の疎行列。

格子が構造格子 (6 面体要素) なので、節点 (i,j,k) の非ゼロ列は
3x3x3 ステンシルを領域内に切り詰めたものに一致する。列は昇順に並ぶので
探索なしに格納位置を計算できる (crs_offset)。
*/

#include "fem.h"
#include "fem_prototype.h"

// 節点 (i,j,k) の隣接範囲
static void stencil_range(int i, int j, int k,
	int *dilo, int *dihi, int *djlo, int *djhi, int *dklo, int *dkhi)
{
	*dilo = (i > 0)  ? -1 : 0;
	*dihi = (i < Nx) ? +1 : 0;
	*djlo = (j > 0)  ? -1 : 0;
	*djhi = (j < Ny) ? +1 : 0;
	*dklo = (k > 0)  ? -1 : 0;
	*dkhi = (k < Nz) ? +1 : 0;
}


// 行 row (= 節点 (i,j,k)) の中で隣接節点 (i+di, j+dj, k+dk) が入る位置
int64_t crs_offset(int64_t rowstart, int i, int j, int k, int di, int dj, int dk)
{
	int dilo, dihi, djlo, djhi, dklo, dkhi;
	stencil_range(i, j, k, &dilo, &dihi, &djlo, &djhi, &dklo, &dkhi);

	const int ndj = djhi - djlo + 1;
	const int ndk = dkhi - dklo + 1;

	return rowstart + ((((int64_t)(di - dilo) * ndj) + (dj - djlo)) * ndk) + (dk - dklo);
}


void crs_alloc(crs_t *A)
{
	if (MeshMode) {
		if (MeshDim == 2) {
			crs_alloc_tri(A);
			return;
		}
		if      (MeshElem == MESHELEM_HEX)   crs_alloc_hex(A);
		else if (MeshElem == MESHELEM_PRISM) crs_alloc_prism(A);
		else                                 crs_alloc_tet(A);
		return;
	}

	const int64_t n = num_node();

	A->n = n;
	A->rowptr = (int64_t *)malloc((n + 1) * sizeof(int64_t));

	// 行毎の非ゼロ数 -> rowptr
	A->rowptr[0] = 0;
	for (int i = 0; i <= Nx; i++) {
	for (int j = 0; j <= Ny; j++) {
	for (int k = 0; k <= Nz; k++) {
		int dilo, dihi, djlo, djhi, dklo, dkhi;
		stencil_range(i, j, k, &dilo, &dihi, &djlo, &djhi, &dklo, &dkhi);
		const int64_t nz = (int64_t)(dihi - dilo + 1) * (djhi - djlo + 1) * (dkhi - dklo + 1);
		const int64_t row = node_index(i, j, k);
		A->rowptr[row + 1] = nz;
	}
	}
	}
	for (int64_t row = 0; row < n; row++) {
		A->rowptr[row + 1] += A->rowptr[row];
	}
	A->nnz = A->rowptr[n];

	A->col = (int32_t *)malloc(A->nnz * sizeof(int32_t));
	A->val = (double *)malloc(A->nnz * sizeof(double));

	// 列番号 (昇順)
	int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (i = 0; i <= Nx; i++) {
	for (int j = 0; j <= Ny; j++) {
	for (int k = 0; k <= Nz; k++) {
		int dilo, dihi, djlo, djhi, dklo, dkhi;
		stencil_range(i, j, k, &dilo, &dihi, &djlo, &djhi, &dklo, &dkhi);
		int64_t p = A->rowptr[node_index(i, j, k)];
		for (int di = dilo; di <= dihi; di++) {
		for (int dj = djlo; dj <= djhi; dj++) {
		for (int dk = dklo; dk <= dkhi; dk++) {
			A->col[p++] = (int32_t)node_index(i + di, j + dj, k + dk);
		}
		}
		}
	}
	}
	}

	crs_zero(A);
}


void crs_zero(crs_t *A)
{
	// IEEE754 の 0.0 は全ビット 0 なので memset で足りる
	memset(A->val, 0, (size_t)A->nnz * sizeof(double));
}


void crs_free(crs_t *A)
{
	free(A->rowptr);
	free(A->col);
	free(A->val);
	A->rowptr = NULL;
	A->col = NULL;
	A->val = NULL;
	A->n = A->nnz = 0;
}


// y = A x (fix[row] != 0 の行は y = x : Dirichlet 節点の恒等行)
void crs_spmv(const crs_t *A, const double *x, double *y, const unsigned char *fix)
{
	// MSVC の OpenMP 2.0 は 64bit のループ変数を受け付けないので int を使う
	// (節点数が INT32_MAX 未満であることは setup() で確認済み)
	const int n = (int)A->n;

	int row;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (row = 0; row < n; row++) {
		if ((fix != NULL) && fix[row]) {
			y[row] = x[row];
			continue;
		}
		double s = 0;
		const int64_t p0 = A->rowptr[row];
		const int64_t p1 = A->rowptr[row + 1];
		for (int64_t p = p0; p < p1; p++) {
			s += A->val[p] * x[A->col[p]];
		}
		y[row] = s;
	}
}


/*
複素対称行列の積 y = (K + jω M) x  (fix[row] != 0 の行は恒等行 y = x)。

COCG (反復) と直接解法の**両方**がこれを使う。同じ関数を通しておかないと、
直接解法が最後に計算する残差が「自分の思う行列」に対する残差になってしまい、
検査として意味を失う。w1..w4 は長さ n の作業配列。
*/
void crs_spmv_c(const crs_t *K, const crs_t *M, double omega,
	const double *xr, const double *xi, double *yr, double *yi,
	const unsigned char *fix, double *w1, double *w2, double *w3, double *w4)
{
	const int n = (int)K->n;

	crs_spmv(K, xr, w1, fix);		// 固定行は w1 = xr
	crs_spmv(K, xi, w2, fix);
	crs_spmv(M, xr, w3, NULL);
	crs_spmv(M, xi, w4, NULL);

	int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (i = 0; i < n; i++) {
		if ((fix != NULL) && fix[i]) {
			yr[i] = w1[i];
			yi[i] = w2[i];
		}
		else {
			yr[i] = w1[i] - (omega * w4[i]);
			yi[i] = w2[i] + (omega * w3[i]);
		}
	}
}


// 対角成分
void crs_diag(const crs_t *A, double *d)
{
	const int n = (int)A->n;

	int row;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (row = 0; row < n; row++) {
		double s = 0;
		const int64_t p0 = A->rowptr[row];
		const int64_t p1 = A->rowptr[row + 1];
		for (int64_t p = p0; p < p1; p++) {
			if (A->col[p] == (int32_t)row) {
				s = A->val[p];
				break;
			}
		}
		d[row] = s;
	}
}


// 1 行分の内積 (A x)_row : 反作用 (電荷・電流) の計算に使う
double crs_row_dot(const crs_t *A, int64_t row, const double *x)
{
	double s = 0;
	const int64_t p0 = A->rowptr[row];
	const int64_t p1 = A->rowptr[row + 1];
	for (int64_t p = p0; p < p1; p++) {
		s += A->val[p] * x[A->col[p]];
	}

	return s;
}
