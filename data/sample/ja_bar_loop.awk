# 準静的な導体棒のヒステリシス損 (FEM の B 駆動の逆 J-A とは独立に書いた H 駆動の順方向のスカラー J-A)。
# 導体の厚さ方向 z に電流が一様なら H_y(z) = J(t)(z - t/2) (Ampere)。各点の H を J0|z - t/2| sin(2π s)
# (s = 0 .. ncyc 周期) で動かし、最後の周期の ∮ H dB を体積で積分した値 [J] を出す。
#   awk -v ms= -v aa= -v al= -v kk= -v cc= -v j0= -v tt= -v area= -v ncyc= -f ja_bar_loop.awk
#   (-v bh="H1 B1 H2 B2 ..." -v nstep= なら J-A の代わりに履歴の無い B-H 曲線)
# nstep > 0 なら H を 1 周期 nstep 点で取り、仕事を FEM と同じ右端の和にする (時間の打ち切りまで揃える)。
# area = 長さ x 幅 (z に垂直な断面ではなく z 方向に積分する面積)。z は [0, t/2] を nz 点の中点則で積分して 2 倍
function lang(x,   ax) { ax = (x < 0) ? -x : x; if (ax < 1e-4) return x / 3; return (exp(x) + exp(-x)) / (exp(x) - exp(-x)) - 1 / x }
function dlang(x,   ax, s) { ax = (x < 0) ? -x : x; if (ax < 1e-4) return 1 / 3; s = (exp(x) - exp(-x)) / 2; return 1 / (x * x) - 1 / (s * s) }
function dmdh(h, m, dl,   he, man, dman, dm, irr, den, num, den2, c) {
	he = h + al * m; man = ms * lang(he / aa); dman = ms / aa * dlang(he / aa); dm = man - m; irr = 0
	if (dm * dl > 0) { den = dl * kk - al * dm; if (den != 0) irr = (1 - cc) * dm / den }
	num = irr + cc * dman; den2 = 1 - al * cc * dman; c = (den2 != 0) ? num / den2 : num
	return (c < 0) ? 0 : c
}
# B-H 曲線 (履歴なし) の B(H)。bh = "H1 B1 H2 B2 ..." (sol の bh と同じ区分線形 H(B) の逆。最初は原点からの直線、
# 最後は最終区間の傾きで外挿)
function bofh(h,   ah, s, i, sl) {
	ah = (h < 0) ? -h : h; s = (h < 0) ? -1 : 1
	if (ah <= bhh[1]) return s * ah * bhb[1] / bhh[1]
	if (ah >= bhh[nbh]) { sl = (nbh >= 2) ? (bhh[nbh] - bhh[nbh - 1]) / (bhb[nbh] - bhb[nbh - 1]) : bhh[1] / bhb[1]; return s * (bhb[nbh] + (ah - bhh[nbh]) / sl) }
	i = 1; while ((i < nbh - 1) && (ah > bhh[i + 1])) i++
	return s * (bhb[i] + (ah - bhh[i]) * (bhb[i + 1] - bhb[i]) / (bhh[i + 1] - bhh[i]))
}
BEGIN {
	pi = atan2(0, -1); mu = 4e-7 * pi; nz = 64; nper = 2000; tot = 0
	nbh = 0; if (bh != "") { nf = split(bh, bv, " "); for (q = 1; q < nf; q += 2) { nbh++; bhh[nbh] = bv[q] + 0; bhb[nbh] = bv[q + 1] + 0 } }
	for (iz = 0; iz < nz; iz++) {
		z = (iz + 0.5) * (tt / 2) / nz; hamp = j0 * (tt / 2 - z)
		h = 0; m = 0; w = 0
		if (nbh > 0) {
			# B-H 曲線 : FEM と同じ刻みで H を取り、右端の和 Σ H^{n+1} (B^{n+1} - B^n)
			b = 0
			for (k = 1; k <= ncyc * nstep; k++) {
				hn = hamp * sin(2 * pi * k / nstep); bn = bofh(hn)
				if (k > (ncyc - 1) * nstep) w += hn * (bn - b)
				b = bn
			}
		}
		else if (nstep > 0) {
			# FEM と同じ刻み (1 周期 nstep 点) で H を取り、刻みの中は細かく積分して、
			# 仕事は FEM と同じ右端の和 Σ H^{n+1} (B^{n+1} - B^n)
			nsub = int(nper / nstep); if (nsub < 1) nsub = 1
			for (k = 1; k <= ncyc * nstep; k++) {
				hn = hamp * sin(2 * pi * k / nstep); dh = (hn - h) / nsub; dl = (dh >= 0) ? 1 : -1; m0 = m; h0 = h
				for (q = 0; q < nsub; q++) {
					c1 = dmdh(h, m, dl); c2 = dmdh(h + dh, m + c1 * dh, dl); m += (c1 + c2) / 2 * dh; h += dh
				}
				h = hn
				if (k > (ncyc - 1) * nstep) w += hn * mu * ((hn - h0) + (m - m0))
			}
		}
		else {
			for (k = 1; k <= ncyc * nper; k++) {
				hn = hamp * sin(2 * pi * k / nper); dh = hn - h; dl = (dh >= 0) ? 1 : -1
				c1 = dmdh(h, m, dl); c2 = dmdh(hn, m + c1 * dh, dl); mn = m + (c1 + c2) / 2 * dh
				if (k > (ncyc - 1) * nper) w += (h + hn) / 2 * mu * ((hn - h) + (mn - m))
				h = hn; m = mn
			}
		}
		tot += w * (tt / 2) / nz
	}
	printf "%.10e\n", 2 * tot * area
}
