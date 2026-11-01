# bar_eddy の 1 次元の厳密解 Z = γL / (2σW tanh(γt/2))、γ = sqrt(jωμ0σ) の R と L を出す (awk -v f=<Hz>)
function cmul(ar, ai, br, bi) { RR = ar * br - ai * bi; II = ar * bi + ai * br }
function cdiv(ar, ai, br, bi,   d) { d = br * br + bi * bi; RR = (ar * br + ai * bi) / d; II = (ai * br - ar * bi) / d }
BEGIN {
	pi = atan2(0, -1); mu = 4e-7 * pi; L = 2e-3; W = 0.25e-3; t = 1e-3; s = 5.8e7; w = 2 * pi * f
	k = sqrt(w * mu * s / 2); gr = k; gi = k
	xr = gr * t / 2; xi = gi * t / 2
	# tanh(x) = sinh(2xr) / (cosh 2xr + cos 2xi) + j sin(2xi) / (cosh 2xr + cos 2xi)
	ch = (exp(2 * xr) + exp(-2 * xr)) / 2; sh = (exp(2 * xr) - exp(-2 * xr)) / 2; den = ch + cos(2 * xi)
	tr = sh / den; ti = sin(2 * xi) / den
	cdiv(gr * L, gi * L, 2 * s * W * tr, 2 * s * W * ti)
	printf "%.9e %.9e\n", RR, II / w
}
