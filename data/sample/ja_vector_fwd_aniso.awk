# 異方性のベクトル J-A (jaaxis) の H 駆動の順方向 (FEM の B 駆動の逆とは独立に書いたもの)。
# ja_vector_fwd.awk と同じく probe の行 (H_p H_q B_p B_q) を読み、H を折れ線でつないで dM = χ dH を積分し、
# B = μ0 (H + M) を FEM の B と比べる。成分 p は (msp aap alp kkp ccp)、成分 q は (msq aaq alq kkq ccq):
#   He = H + α M (対角)、M_an,i = Ms_i L(|He| / a_i) e_i
#   不可逆 |f| / (k_f − α_f |f|) f̂ f̂ᵀ、k_f = Σ k_i f̂_i²、α_f = Σ α_i f̂_i²
#   χ = (I − C χ_an α)⁻¹ [ (I − C) 不可逆 + C χ_an ]
function lang(x) { if (x < 1e-4) return x / 3; return (exp(x) + exp(-x)) / (exp(x) - exp(-x)) - 1 / x }
function dlang(x,   s) { if (x < 1e-4) return 1 / 3; s = (exp(x) - exp(-x)) / 2; return 1 / (x * x) - 1 / (s * s) }
function chi(hx, hy, mx, my, dx, dy,   ex, ey, en, ux, uy, mpx, dpx, spx, mqy, dqy, sqy, c11, c12, c21, c22, fx, fy, fn, irr, kf, af, s1, cf, al, a11, a12, a21, a22, d, x11, x12, x21, x22) {
	ex = hx + alp * mx; ey = hy + alq * my; en = sqrt(ex * ex + ey * ey)
	if (en > 0) { ux = ex / en; uy = ey / en } else { ux = 1; uy = 0 }
	mpx = msp * lang(en / aap); dpx = msp / aap * dlang(en / aap); spx = (en > 1e-4 * aap) ? mpx / en : dpx
	mqy = msq * lang(en / aaq); dqy = msq / aaq * dlang(en / aaq); sqy = (en > 1e-4 * aaq) ? mqy / en : dqy
	c11 = dpx * ux * ux + spx * (1 - ux * ux); c12 = (dpx - spx) * ux * uy
	c21 = (dqy - sqy) * ux * uy; c22 = dqy * uy * uy + sqy * (1 - uy * uy)
	fx = spx * ex - mx; fy = sqy * ey - my; fn = sqrt(fx * fx + fy * fy); irr = 0
	if ((fn > 0) && (fx * dx + fy * dy > 0)) {
		fx /= fn; fy /= fn
		kf = kkp * fx * fx + kkq * fy * fy; af = alp * fx * fx + alq * fy * fy
		s1 = (1 - ccp) * fx * fx + (1 - ccq) * fy * fy
		irr = fn / (kf - af * fn)
		cf = fx * ccp * (c11 * fx + c12 * fy) + fy * ccq * (c21 * fx + c22 * fy)
		al = irr * s1 + cf; if ((al < 0) && (s1 > 0)) irr -= al / s1
	} else { fx = 0; fy = 0 }
	x11 = (1 - ccp) * irr * fx * fx + ccp * c11; x12 = (1 - ccp) * irr * fx * fy + ccp * c12
	x21 = (1 - ccq) * irr * fy * fx + ccq * c21; x22 = (1 - ccq) * irr * fy * fy + ccq * c22
	a11 = 1 - ccp * c11 * alp; a12 = -ccp * c12 * alq; a21 = -ccq * c21 * alp; a22 = 1 - ccq * c22 * alq
	d = a11 * a22 - a12 * a21
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
