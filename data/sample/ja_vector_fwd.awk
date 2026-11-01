# H 駆動の順方向のベクトル J-A (FEM の B 駆動の逆とは独立に書いたもの)。probe の行 (H_p H_q B_p B_q) を読み、
# H を折れ線でつないで dM = χ dH を積分し、B = μ0 (H + M) を FEM の B と比べる
function lang(x) { if (x < 1e-4) return x / 3; return (exp(x) + exp(-x)) / (exp(x) - exp(-x)) - 1 / x }
function dlang(x,   s) { if (x < 1e-4) return 1 / 3; s = (exp(x) - exp(-x)) / 2; return 1 / (x * x) - 1 / (s * s) }
function chi(hx, hy, mx, my, dx, dy,   ex, ey, en, x, man, dman, sec, c11, c12, c22, fx, fy, fn, irr, ux, uy, a11, a12, a21, a22, d, x11, x12, x21, x22) {
	ex = hx + al * mx; ey = hy + al * my; en = sqrt(ex * ex + ey * ey)
	x = en / aa; man = ms * lang(x); dman = ms / aa * dlang(x)
	sec = (en > 1e-4 * aa) ? man / en : dman
	if (en > 0) { ux = ex / en; uy = ey / en } else { ux = 1; uy = 0 }
	c11 = dman * ux * ux + sec * (1 - ux * ux); c12 = (dman - sec) * ux * uy; c22 = dman * uy * uy + sec * (1 - uy * uy)
	fx = sec * ex - mx; fy = sec * ey - my; fn = sqrt(fx * fx + fy * fy); irr = 0
	if ((fn > 0) && (fx * dx + fy * dy > 0)) { fx /= fn; fy /= fn; irr = (1 - cc) * fn / (kk - al * fn) } else { fx = 0; fy = 0 }
	x11 = irr * fx * fx + cc * c11; x12 = irr * fx * fy + cc * c12; x21 = x12; x22 = irr * fy * fy + cc * c22
	a11 = 1 - al * cc * c11; a12 = -al * cc * c12; a21 = a12; a22 = 1 - al * cc * c22; d = a11 * a22 - a12 * a21
	X11 = (a22 * x11 - a12 * x21) / d; X12 = (a22 * x12 - a12 * x22) / d
	X21 = (-a21 * x11 + a11 * x21) / d; X22 = (-a21 * x12 + a11 * x22) / d
}
BEGIN { mu = 4e-7 * atan2(0, -1); hx = 0; hy = 0; mx = 0; my = 0; emax = 0; bmax = 0; ns = 400 }
$1 == "probe" && NF == 6 {
	tx = $3; ty = $4; dx = (tx - hx) / ns; dy = (ty - hy) / ns
	for (s = 0; s < ns; s++) {
		chi(hx, hy, mx, my, dx, dy); k1x = X11 * dx + X12 * dy; k1y = X21 * dx + X22 * dy
		chi(hx + dx, hy + dy, mx + k1x, my + k1y, dx, dy); k2x = X11 * dx + X12 * dy; k2y = X21 * dx + X22 * dy
		mx += (k1x + k2x) / 2; my += (k1y + k2y) / 2; hx += dx; hy += dy
	}
	hx = tx; hy = ty
	bx = mu * (hx + mx); by = mu * (hy + my)
	e = sqrt((bx - $5) ^ 2 + (by - $6) ^ 2); if (e > emax) emax = e
	b = sqrt($5 * $5 + $6 * $6); if (b > bmax) bmax = b
}
END { printf "%.6e %.6e\n", emax / bmax, bmax }
