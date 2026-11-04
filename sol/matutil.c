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


/*
Anderson 加速 (深さ m) の 1 段。非線形の固定点反復 y <- y + f (f は固定点写像の残差)
を加速する。y, f は長さ n、yp / fp は前回の y / f、dy / df は差分の履歴 (n x m、
循環)、*cnt はこれまでの段の数。y を上書きして次の点を返す:
	y_new = y + w f - Σ γ_j (ΔY_j + w ΔF_j),   γ = argmin |f - ΔF γ|
最小二乗は正規方程式 (m x m、正則化つき) で解く。破綻したら γ = 0 (緩和付きの
素の固定点) に落ちる。
*/
void anderson_step(int64_t n, double *y, const double *f, double *yp, double *fp,
	double *dy, double *df, int m, int *cnt, double relax)
{
	if (*cnt > 0) {
		const int slot = (*cnt - 1) % m;
		for (int64_t c = 0; c < n; c++) {
			dy[(slot * n) + c] = y[c] - yp[c];
			df[(slot * n) + c] = f[c] - fp[c];
		}
	}
	for (int64_t c = 0; c < n; c++) {
		yp[c] = y[c];
		fp[c] = f[c];
	}
	const int mh = ((*cnt < m) ? *cnt : m);
	double gam[ANDERSON_MAX];
	for (int j = 0; j < ANDERSON_MAX; j++) gam[j] = 0;
	if (mh > 0) {
		double am[ANDERSON_MAX][ANDERSON_MAX + 1];
		double tr = 0;
		for (int i = 0; i < mh; i++) {
			for (int j = 0; j < mh; j++) {
				double sum = 0;
				for (int64_t c = 0; c < n; c++) sum += df[(i * n) + c] * df[(j * n) + c];
				am[i][j] = sum;
			}
			double sum = 0;
			for (int64_t c = 0; c < n; c++) sum += df[(i * n) + c] * f[c];
			am[i][mh] = sum;
			tr += am[i][i];
		}
		// 正則化 (履歴がほぼ線形従属のとき)。正規方程式は正定値なのでピボット不要
		for (int i = 0; i < mh; i++) am[i][i] += (1e-12 * tr) + 1e-300;
		int sing = 0;
		for (int k = 0; k < mh; k++) {
			if (am[k][k] <= 0) { sing = 1; break; }
			for (int i = k + 1; i < mh; i++) {
				const double r = am[i][k] / am[k][k];
				for (int j = k; j <= mh; j++) am[i][j] -= r * am[k][j];
			}
		}
		for (int k = mh - 1; !sing && (k >= 0); k--) {
			double sum = am[k][mh];
			for (int j = k + 1; j < mh; j++) sum -= am[k][j] * gam[j];
			gam[k] = sum / am[k][k];
		}
		if (sing) {
			for (int j = 0; j < ANDERSON_MAX; j++) gam[j] = 0;
		}
	}
	for (int64_t c = 0; c < n; c++) {
		double v = y[c] + (relax * f[c]);
		for (int j = 0; j < mh; j++) v -= gam[j] * (dy[(j * n) + c] + (relax * df[(j * n) + c]));
		y[c] = v;
	}
	(*cnt)++;
}
