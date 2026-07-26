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
