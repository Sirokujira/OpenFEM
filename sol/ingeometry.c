/*
ingeometry.c

点 (x, y, z) が形状の内部にあるか (1 : 内部または境界、0 : 外部)。
形状番号と引数 g[] の意味は OpenFDTD の geometry と共通:

	 1 : 直方体   g = xmin xmax ymin ymax zmin zmax
	 2 : 楕円体   g = 外接直方体
	11 : X 円柱   g = x1 x2 (y,z の外接長方形)
	12 : Y 円柱   g = (x の外接) y1 y2 (z の外接)
	13 : Z 円柱   g = (x,y の外接長方形) z1 z2

厚さ 0 の指定 (例 zmin == zmax) は eps の許容幅を持つ平板として扱う。
これにより電極を「面」として指定できる。
*/

#include <math.h>

int ingeometry(double x, double y, double z, int shape, const double *g, double eps)
{
	const double zero = 1e-6;
	const double eps2 = eps * eps;

	// 直方体
	if (shape == 1) {
		if (((x - g[0]) * (x - g[1]) <= eps2) &&
		    ((y - g[2]) * (y - g[3]) <= eps2) &&
		    ((z - g[4]) * (z - g[5]) <= eps2)) {
			return 1;
		}
	}
	// 楕円体
	else if (shape == 2) {
		const double x0 = (g[0] + g[1]) / 2;
		const double y0 = (g[2] + g[3]) / 2;
		const double z0 = (g[4] + g[5]) / 2;
		const double xr = fabs(g[0] - g[1]) / 2;
		const double yr = fabs(g[2] - g[3]) / 2;
		const double zr = fabs(g[4] - g[5]) / 2;
		if ((xr > 0) && (yr > 0) && (zr > 0)) {
			if (((x - x0) * (x - x0) / (xr * xr)
			   + (y - y0) * (y - y0) / (yr * yr)
			   + (z - z0) * (z - z0) / (zr * zr)) < 1 + zero) {
				return 1;
			}
		}
	}
	// X 円柱
	else if (shape == 11) {
		const double y0 = (g[2] + g[3]) / 2;
		const double z0 = (g[4] + g[5]) / 2;
		const double yr = fabs(g[2] - g[3]) / 2;
		const double zr = fabs(g[4] - g[5]) / 2;
		if ((yr > 0) && (zr > 0)) {
			if (((x - g[0]) * (x - g[1]) <= eps2) &&
			    (((y - y0) * (y - y0) / (yr * yr)
			    + (z - z0) * (z - z0) / (zr * zr)) < 1 + zero)) {
				return 1;
			}
		}
	}
	// Y 円柱
	else if (shape == 12) {
		const double z0 = (g[4] + g[5]) / 2;
		const double x0 = (g[0] + g[1]) / 2;
		const double zr = fabs(g[4] - g[5]) / 2;
		const double xr = fabs(g[0] - g[1]) / 2;
		if ((zr > 0) && (xr > 0)) {
			if (((y - g[2]) * (y - g[3]) <= eps2) &&
			    (((z - z0) * (z - z0) / (zr * zr)
			    + (x - x0) * (x - x0) / (xr * xr)) < 1 + zero)) {
				return 1;
			}
		}
	}
	// Z 円柱
	else if (shape == 13) {
		const double x0 = (g[0] + g[1]) / 2;
		const double y0 = (g[2] + g[3]) / 2;
		const double xr = fabs(g[0] - g[1]) / 2;
		const double yr = fabs(g[2] - g[3]) / 2;
		if ((xr > 0) && (yr > 0)) {
			if (((z - g[4]) * (z - g[5]) <= eps2) &&
			    (((x - x0) * (x - x0) / (xr * xr)
			    + (y - y0) * (y - y0) / (yr * yr)) < 1 + zero)) {
				return 1;
			}
		}
	}

	return 0;
}
