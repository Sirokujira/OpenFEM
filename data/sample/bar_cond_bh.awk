# 非線形の導体 (B-H 曲線) の棒の 1 次元の参照解 (rlc_check.sh の [A order 2] が使う)。
#   awk -v f=<Hz> -v I=<A 波高値> [-v n=<導体の分割>] -f bar_cond_bh.awk  -> R L iters
#   -v lin=<μr> で導体を線形にする (bar_layer.awk の空気の層 mr = 1 と同じ形になる。参照解自身の検査)
# bar_air の形 (導体 0 < z < t、その上に空気の層 t < z < t + g、上面と下面で A = 0)
# で、導体の μ が B-H 曲線 (100,0.5) (1000,1.0) (10000,1.2) (外挿は最後の傾き) に従う。
# 1 次元の式 −(ν A')' + jωσ A = σ E0 を 1 次要素 (要素ごとの ν) で解き、実効磁気抵抗率法と同じ
# 更新 ν = |H| / B(|H|)、|H| = ν |A'| (波高値) を収束するまで繰り返す。端子電流が I になるよう
# E0 を倍率で決める (ν を固定すれば線形)。ソルバー (3 次元・2 次の Nedelec・積分点ごとの ν) とは
# 独立な実装で、分割を細かくすれば連続の解に h² で近づく (n = 4000 で 1e-8 程度)。
function bofh(h) {
	if (h <= 100) return h * 0.5 / 100
	if (h <= 1000) return 0.5 + (h - 100) * 0.5 / 900
	if (h <= 10000) return 1.0 + (h - 1000) * 0.2 / 9000
	return 1.2 + (h - 10000) * 0.2 / 9000
}
BEGIN {
	pi = atan2(0, -1); mu0 = 4 * pi * 1e-7
	Lx = 2e-3; W = 0.25e-3; t = 1e-3; g = 0.5e-3; sg = 5.8e7
	om = 2 * pi * f
	if (n == "") n = 4000
	na = int(n * g / t + 0.5)
	N = n + na					# 要素数、節点 0..N (節点 N が上面で A = 0)
	for (e = 0; e < N; e++) {
		hz[e] = ((e < n) ? t / n : g / na)
		se[e] = ((e < n) ? sg : 0)
		nu[e] = ((e < n) ? 1 / (0.5 / 100) : 1 / mu0)
		if ((lin != "") && (e < n)) nu[e] = 1 / (mu0 * lin)
	}
	for (it = 1; it <= 1000; it++) {
		# 3 重対角 (複素) : 対角 dr + j di、副対角 or + j oi (要素 e は節点 e, e+1)
		for (i = 0; i <= N; i++) { dr[i] = di[i] = 0; br[i] = bi[i] = 0; ofr[i] = ofi[i] = 0 }
		for (e = 0; e < N; e++) {
			k = nu[e] / hz[e]; m = om * se[e] * hz[e] / 6
			dr[e] += k; dr[e + 1] += k; ofr[e] += -k
			di[e] += 2 * m; di[e + 1] += 2 * m; ofi[e] += m
			br[e] += se[e] * hz[e] / 2; br[e + 1] += se[e] * hz[e] / 2	# E0 = 1
		}
		# A(0) = A(N) = 0 : 節点 1..N-1 を Thomas 法で (対称な複素系、ピボットなし)
		M = N - 1
		for (i = 1; i <= M; i++) {
			ar_ = dr[i]; ai_ = di[i]
			if (i > 1) {
				# 対角 -= off(i-1)^2 / dd(i-1)
				xr = ofr[i - 1]; xi = ofi[i - 1]
				sr = xr * xr - xi * xi; si = 2 * xr * xi
				den = ddr[i - 1] * ddr[i - 1] + ddi[i - 1] * ddi[i - 1]
				qr = (sr * ddr[i - 1] + si * ddi[i - 1]) / den; qi = (si * ddr[i - 1] - sr * ddi[i - 1]) / den
				ar_ -= qr; ai_ -= qi
				# 右辺 -= off(i-1) * y(i-1) / dd(i-1)
				den2 = den
				tr = (yr[i - 1] * ddr[i - 1] + yi[i - 1] * ddi[i - 1]) / den2
				ti = (yi[i - 1] * ddr[i - 1] - yr[i - 1] * ddi[i - 1]) / den2
				yr[i] = br[i] - (xr * tr - xi * ti); yi[i] = bi[i] - (xr * ti + xi * tr)
			}
			else { yr[i] = br[i]; yi[i] = bi[i] }
			ddr[i] = ar_; ddi[i] = ai_
		}
		for (i = M; i >= 1; i--) {
			sr = yr[i]; si = yi[i]
			if (i < M) { sr -= ofr[i] * Ar[i + 1] - ofi[i] * Ai[i + 1]; si -= ofr[i] * Ai[i + 1] + ofi[i] * Ar[i + 1] }
			den = ddr[i] * ddr[i] + ddi[i] * ddi[i]
			Ar[i] = (sr * ddr[i] + si * ddi[i]) / den; Ai[i] = (si * ddr[i] - sr * ddi[i]) / den
		}
		Ar[0] = Ai[0] = Ar[N] = Ai[N] = 0
		# 電流 (E0 = 1) : I1 = W Σ σ ∫(1 − jωA) dz
		Kr = 0; Ki = 0
		for (e = 0; e < n; e++) {
			am_r = (Ar[e] + Ar[e + 1]) / 2; am_i = (Ai[e] + Ai[e + 1]) / 2
			Kr += se[e] * hz[e] * (1 + om * am_i); Ki += se[e] * hz[e] * (-om * am_r)
		}
		I1r = W * Kr; I1i = W * Ki
		# Z = E0 Lx / I1
		den = I1r * I1r + I1i * I1i
		Zr = Lx * I1r / den; Zi = -Lx * I1i / den
		# 端子電流を I にする倍率 |s| = I / |I1| (|B| にだけ効く)
		s = I / sqrt(den)
		dmax = 0
		for (e = 0; e < n; e++) {
			bpr = (Ar[e + 1] - Ar[e]) / hz[e]; bpi = (Ai[e + 1] - Ai[e]) / hz[e]
			b = s * sqrt(bpr * bpr + bpi * bpi)
			hm = nu[e] * b
			nn = ((hm > 0) ? hm / bofh(hm) : 1 / (0.5 / 100))
			d = (nn - nu[e]) / nn; if (d < 0) d = -d
			if (d > dmax) dmax = d
			if (lin == "") nu[e] = nn
		}
		if (lin != "") break
		if (dmax < 1e-11) break
	}
	printf "%.12e %.12e %d %.3e\n", Zr, Zi / om, it, dmax
}
