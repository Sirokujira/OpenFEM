# 平行平板線路の導体 (厚さ t、電流密度一様) を J-A の材料にしたときの、導体の平均の B の参照値 (FEM の B 駆動の逆 J-A とは
# 独立に書いた H 駆動の順方向のスカラー J-A)。導体の中の H は Ampère で H(s) = s I/W (s = 0 .. 1 は外側の面から内側の面へ) に
# 決まり、材料に依らない (1 次元)。各 s の点で H を電流の列 cur (単調でなくてよい) に沿って nsub 等分の Heun で積分し、
# 各電流での B = μ0 (H + M) を s について中点則 (nz 点) で平均して出す。
#   awk -v ms= -v aa= -v al= -v kk= -v cc= -v w= -v cur="0.5 1 2" -f ja_plate_cond.awk
# 出力 : 1 行に 1 点 "番号 電流 B_mean"
function lang(x,   ax) { ax = (x < 0) ? -x : x; if (ax < 1e-4) return x / 3; return (exp(x) + exp(-x)) / (exp(x) - exp(-x)) - 1 / x }
function dlang(x,   ax, s) { ax = (x < 0) ? -x : x; if (ax < 1e-4) return 1 / 3; s = (exp(x) - exp(-x)) / 2; return 1 / (x * x) - 1 / (s * s) }
function dmdh(h, m, dl,   he, man, dman, dm, irr, den, num, den2, c) {
	he = h + al * m; man = ms * lang(he / aa); dman = ms / aa * dlang(he / aa); dm = man - m; irr = 0
	if (dm * dl > 0) { den = dl * kk - al * dm; if (den != 0) irr = (1 - cc) * dm / den }
	num = irr + cc * dman; den2 = 1 - al * cc * dman; c = (den2 != 0) ? num / den2 : num
	return (c < 0) ? 0 : c
}
BEGIN {
	mu = 4e-7 * atan2(0, -1); if (nz == "") nz = 400; if (nsub == "") nsub = 400
	np = split(cur, ci, " ")
	for (k = 1; k <= np; k++) bs[k] = 0
	for (iz = 0; iz < nz; iz++) {
		s = (iz + 0.5) / nz; h = 0; m = 0
		for (k = 1; k <= np; k++) {
			hn = s * ci[k] / w; dh = (hn - h) / nsub; dl = (dh >= 0) ? 1 : -1
			for (q = 0; q < nsub; q++) { c1 = dmdh(h, m, dl); c2 = dmdh(h + dh, m + c1 * dh, dl); m += (c1 + c2) / 2 * dh; h += dh }
			h = hn
			bs[k] += mu * (h + m) / nz
		}
	}
	for (k = 1; k <= np; k++) printf "%d %.6e %.9e\n", k, ci[k], bs[k]
}
