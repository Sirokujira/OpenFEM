/*
matutil.c

小さい密行列の逆行列 (Gauss-Jordan 法、部分ピボット選択)。
ポート数は高々 MAXPORT なので素朴な実装で足りる。

戻り値 : 0 = 正常、1 = 特異
*/

#include "fem.h"
#include "fem_prototype.h"

int mat_inverse(const double *a, double *inv, int n)
{
	if (n < 1) return 1;

	double *w = (double *)malloc((size_t)n * n * sizeof(double));
	memcpy(w, a, (size_t)n * n * sizeof(double));

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			inv[(i * n) + j] = ((i == j) ? 1.0 : 0.0);
		}
	}

	for (int k = 0; k < n; k++) {
		// ピボット選択
		int ip = k;
		double amax = fabs(w[(k * n) + k]);
		for (int i = k + 1; i < n; i++) {
			const double aa = fabs(w[(i * n) + k]);
			if (aa > amax) {
				amax = aa;
				ip = i;
			}
		}
		if (amax <= 0) {
			free(w);
			return 1;
		}
		if (ip != k) {
			for (int j = 0; j < n; j++) {
				double t = w[(k * n) + j];
				w[(k * n) + j] = w[(ip * n) + j];
				w[(ip * n) + j] = t;
				t = inv[(k * n) + j];
				inv[(k * n) + j] = inv[(ip * n) + j];
				inv[(ip * n) + j] = t;
			}
		}

		const double piv = w[(k * n) + k];
		for (int j = 0; j < n; j++) {
			w[(k * n) + j] /= piv;
			inv[(k * n) + j] /= piv;
		}

		for (int i = 0; i < n; i++) {
			if (i == k) continue;
			const double f = w[(i * n) + k];
			if (f == 0) continue;
			for (int j = 0; j < n; j++) {
				w[(i * n) + j]   -= f * w[(k * n) + j];
				inv[(i * n) + j] -= f * inv[(k * n) + j];
			}
		}
	}

	free(w);

	return 0;
}


// 複素小行列の逆行列 (Gauss-Jordan 法、部分ピボット選択)
// 入力 a = ar + j ai、出力 b = br + j bi
// 戻り値 : 0 = 正常、1 = 特異
int mat_inverse_c(const double *ar, const double *ai, double *br, double *bi, int n)
{
	if (n < 1) return 1;

	const size_t nn = (size_t)n * n;
	double *wr = (double *)malloc(nn * sizeof(double));
	double *wi = (double *)malloc(nn * sizeof(double));
	memcpy(wr, ar, nn * sizeof(double));
	memcpy(wi, ai, nn * sizeof(double));

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			br[(i * n) + j] = ((i == j) ? 1.0 : 0.0);
			bi[(i * n) + j] = 0.0;
		}
	}

	for (int k = 0; k < n; k++) {
		int ip = k;
		double amax = (wr[(k * n) + k] * wr[(k * n) + k]) + (wi[(k * n) + k] * wi[(k * n) + k]);
		for (int i = k + 1; i < n; i++) {
			const double aa = (wr[(i * n) + k] * wr[(i * n) + k]) + (wi[(i * n) + k] * wi[(i * n) + k]);
			if (aa > amax) {
				amax = aa;
				ip = i;
			}
		}
		if (amax <= 0) {
			free(wr);
			free(wi);
			return 1;
		}
		if (ip != k) {
			for (int j = 0; j < n; j++) {
				double t;
				t = wr[(k * n) + j]; wr[(k * n) + j] = wr[(ip * n) + j]; wr[(ip * n) + j] = t;
				t = wi[(k * n) + j]; wi[(k * n) + j] = wi[(ip * n) + j]; wi[(ip * n) + j] = t;
				t = br[(k * n) + j]; br[(k * n) + j] = br[(ip * n) + j]; br[(ip * n) + j] = t;
				t = bi[(k * n) + j]; bi[(k * n) + j] = bi[(ip * n) + j]; bi[(ip * n) + j] = t;
			}
		}

		const double pr = wr[(k * n) + k];
		const double pi = wi[(k * n) + k];
		const double pd = (pr * pr) + (pi * pi);
		for (int j = 0; j < n; j++) {
			double xr, xi;
			xr = ((wr[(k * n) + j] * pr) + (wi[(k * n) + j] * pi)) / pd;
			xi = ((wi[(k * n) + j] * pr) - (wr[(k * n) + j] * pi)) / pd;
			wr[(k * n) + j] = xr;
			wi[(k * n) + j] = xi;
			xr = ((br[(k * n) + j] * pr) + (bi[(k * n) + j] * pi)) / pd;
			xi = ((bi[(k * n) + j] * pr) - (br[(k * n) + j] * pi)) / pd;
			br[(k * n) + j] = xr;
			bi[(k * n) + j] = xi;
		}

		for (int i = 0; i < n; i++) {
			if (i == k) continue;
			const double fr = wr[(i * n) + k];
			const double fi = wi[(i * n) + k];
			if ((fr == 0) && (fi == 0)) continue;
			for (int j = 0; j < n; j++) {
				wr[(i * n) + j] -= (fr * wr[(k * n) + j]) - (fi * wi[(k * n) + j]);
				wi[(i * n) + j] -= (fr * wi[(k * n) + j]) + (fi * wr[(k * n) + j]);
				br[(i * n) + j] -= (fr * br[(k * n) + j]) - (fi * bi[(k * n) + j]);
				bi[(i * n) + j] -= (fr * bi[(k * n) + j]) + (fi * br[(k * n) + j]);
			}
		}
	}

	free(wr);
	free(wi);

	return 0;
}


// 対称 3x3 テンソル (成分順 xx, yy, zz, xy, yz, zx) の逆行列
// 透磁率テンソル μ から磁気抵抗率テンソル ν = μ^-1 を作るのに使う
// 戻り値 : 0 = 正常、1 = 特異
int tensor6_inverse(const double a[6], double b[6])
{
	const double xx = a[0], yy = a[1], zz = a[2];
	const double xy = a[3], yz = a[4], zx = a[5];

	const double c00 = (yy * zz) - (yz * yz);
	const double c11 = (xx * zz) - (zx * zx);
	const double c22 = (xx * yy) - (xy * xy);
	const double c01 = (yz * zx) - (xy * zz);
	const double c12 = (xy * zx) - (xx * yz);
	const double c02 = (xy * yz) - (yy * zx);

	const double det = (xx * c00) + (xy * c01) + (zx * c02);
	if (fabs(det) <= 0) return 1;

	b[0] = c00 / det;
	b[1] = c11 / det;
	b[2] = c22 / det;
	b[3] = c01 / det;
	b[4] = c12 / det;
	b[5] = c02 / det;

	return 0;
}
