# ac_heat.awk : 交流の発熱 (thermal_ac.ofe の棒、α = 0) の閉形式
#
# 棒 (x ∈ [0, L]、幅 W、z ∈ [0, t]) の z = 0 / t を awall にした 1 次元の表皮効果では
# J(z) ∝ cosh(γ(z - t/2))、γ = (1 + j)/δ なので、時間平均の発熱は
#   q(u) ∝ |cosh(γu)|² = (cosh(2u/δ) + cos(2u/δ)) / 2,   u = z - t/2
# で、平均は q̄ = ½ I² R / (L W t)。両端 (x = 0, L) を T_s、他を断熱にすると
# -k ∇²T = q(z) は z の cos(mπz/t) で分離でき、モード m の解は
#   θ_m(x) = (q_m / (k κ²)) (1 - cosh(κ(x - L/2)) / cosh(κL/2)),   κ = mπ/t
# (m = 0 は放物線)。中央 x = L/2 の表面 (z = 0) と中心 (z = t/2) の温度差は
#   ΔT = Σ_{m = 2p} (q_m / (k κ²)) (1 - sech(κL/2)) (1 - (-1)^p)
# で、z の分布だけで決まる (q̄ の放物線は消える)。体積平均の温度上昇は z の分布に
# 依らず q̄ L² / (12 k)。R は FE の値を渡す (分布の形だけを比べるため)。
#
#   awk -v f=<Hz> -v sig=<S/m> -v k=<W/(m K)> -v I=<A, 波高値> -v R=<ohm> -f ac_heat.awk
#   -> "ΔT [K] 体積平均の温度上昇 [K]"
#   awk -v f=<Hz> -v sig=<S/m> -v k=<W/(m K)> -v I=<A peak> -v R=<ohm> -f series.awk
BEGIN {
	PI = atan2(0, -1); mu0 = 4e-7 * PI
	L = 2e-3; W = 0.25e-3; t = 1e-3
	d = sqrt(2 / (2 * PI * f * mu0 * sig))
	P = 0.5 * I * I * R; V = L * W * t; qb = P / V
	al = 2 / d
	gb = (d / t) * ((exp(t / d) - exp(-t / d)) / 2 + sin(t / d))
	ds = 0
	for (p = 1; p < 4000; p++) {
		m = 2 * p; be = m * PI / t
		ch = (exp(al * t / 2) + exp(-al * t / 2)) / 2; sh = (exp(al * t / 2) - exp(-al * t / 2)) / 2
		i1 = 2 * (al * sh * cos(be * t / 2) + be * ch * sin(be * t / 2)) / (al * al + be * be)
		i2 = sin((al - be) * t / 2) / (al - be) + sin((al + be) * t / 2) / (al + be)
		sgn = ((p % 2) ? -1 : 1)
		qm = (2 / t) * (qb / gb) * (i1 + i2) * sgn
		ka = m * PI / t; x = ka * L / 2
		sech = (x > 700) ? 0 : 2 / (exp(x) + exp(-x))
		ds += (qm / (k * ka * ka)) * (1 - sech) * (1 - ((p % 2) ? -1 : 1))
	}
	printf "%.12e %.12e\n", ds, qb * L * L / (12 * k)
}
