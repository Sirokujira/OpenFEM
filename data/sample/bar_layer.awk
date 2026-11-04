# bar_air.ofe / bar_air_bh.ofe の 1 次元閉形式 (rlc_check.sh が使う)。
#   awk -v f=<Hz> -v mode=lin -v mr=<層の μr> -v I=<A> -f bar_layer.awk  -> R L H
#   awk -v f=<Hz> -v mode=nl -v I=<A> -f bar_layer.awk               -> R L H μr (非線形の層)
# bar_air の 1 次元閉形式を層の μr (ml) で一般化。導体 μc = 1。
# 界面条件 A(t) + κ A'(t) = 0, κ = g ml / μc。Z = γ L / (σ W D)、D = sinh γt - X (cosh γt - 1)
# X = ((cosh γt - 1) + κγ sinh γt) / (sinh γt + κγ cosh γt)
# 層の H = |A(t)| / (g μ0 ml)、A(t) = (E0/(jω)) ((1 - cosh γt) + X sinh γt)、E0 = I γ / (W σ D)
function cmul(ar, ai, br, bi) { RR = ar * br - ai * bi; RI = ar * bi + ai * br }
function cdiv(ar, ai, br, bi,   d) { d = br * br + bi * bi; RR = (ar * br + ai * bi) / d; RI = (ai * br - ar * bi) / d }
function solve(ml,   gr, gi, a, b, chr, chi, shr, shi, kr, ki, nr, ni, dr, di, xr, xi, Dr, Di, zr, zi, e0r, e0i, atr, ati) {
	# γ = sqrt(j ω μ0 σ) = (1 + j) sqrt(ω μ0 σ / 2)
	gr = sqrt(om * mu0 * sg / 2); gi = gr
	a = gr * t; b = gi * t
	chr = cosh_(a) * cos(b); chi = sinh_(a) * sin(b)
	shr = sinh_(a) * cos(b); shi = cosh_(a) * sin(b)
	ka = g * ml
	kr = ka * gr; ki = ka * gi                          # κγ
	cmul(kr, ki, shr, shi); nr = (chr - 1) + RR; ni = chi + RI
	cmul(kr, ki, chr, chi); dr = shr + RR; di = shi + RI
	cdiv(nr, ni, dr, di); xr = RR; xi = RI
	cmul(xr, xi, chr - 1, chi); Dr = shr - RR; Di = shi - RI
	cmul(gr, gi, Lx, 0); cdiv(RR, RI, sg * W * Dr, sg * W * Di); zr = RR; zi = RI
	# E0 = I γ / (W σ D)、A(t) = (E0 / (jω)) ((1 - cosh) + X sinh)
	cdiv(I * gr, I * gi, W * sg * Dr, W * sg * Di); e0r = RR; e0i = RI
	cmul(xr, xi, shr, shi); atr = (1 - chr) + RR; ati = -chi + RI
	cmul(e0r, e0i, atr, ati); cdiv(RR, RI, 0, om)
	H = sqrt(RR * RR + RI * RI) / (g * mu0 * ml)
	ZR = zr; ZL = zi / om
}
function cosh_(x) { return (exp(x) + exp(-x)) / 2 }
function sinh_(x) { return (exp(x) - exp(-x)) / 2 }
function bofh(h,   i) {   # B-H の逆 (曲線 (100,0.5) (1000,1.0) (10000,1.2))
	if (h <= 100) return h * 0.5 / 100
	if (h <= 1000) return 0.5 + (h - 100) * 0.5 / 900
	if (h <= 10000) return 1.0 + (h - 1000) * 0.2 / 9000
	return 1.2 + (h - 10000) * 0.2 / 9000
}
BEGIN {
	mu0 = 4 * atan2(0, -1) * 1e-7; Lx = 2e-3; W = 0.25e-3; t = 1e-3; g = 0.5e-3; sg = 5.8e7; om = 2 * atan2(0, -1) * f
	if (mode == "lin") { solve(mr); printf "%.10e %.10e %.10e\n", ZR, ZL, H; exit }
	# 非線形 : g(ml) = B(H(ml))/(μ0 H(ml)) - ml = 0 を ln ml で二分法 (符号の変化は 1 つ)
	lo = log(1); hi = log(4000)
	for (it = 0; it < 200; it++) {
		ml = exp((lo + hi) / 2); solve(ml)
		if (bofh(H) / (mu0 * H) - ml > 0) lo = log(ml); else hi = log(ml)
	}
	ml = exp((lo + hi) / 2); solve(ml)
	printf "%.10e %.10e %.10e %.10e\n", ZR, ZL, H, ml
}
