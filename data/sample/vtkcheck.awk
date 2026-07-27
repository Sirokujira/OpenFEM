# vtkcheck.awk — ofe_field.vtk のベクトル場を体積積分して検証に使う
#
# 場の出力そのものを検証するために、書き出した場から**集中定数を作り直して**
# 元の抽出値と比べる (厳密な恒等式なので誤差の言い訳が効かない)。
# python に依存させないため awk で書く。
#
# 使い方 : awk -v arr=<配列名> -f vtkcheck.awk ofe_field.vtk
#
# 出力 (1 行 1 項目):
#   ncell  <セル数>
#   vol    <Σ セル体積>
#   int2   <∫ |v|^2 dV>          … エネルギー / 損失の恒等式に使う
#   vmin   <min |v|>,  vmax <max |v|>
#   amax0  <max |v_x|>, amax1 <max |v_y|>, amax2 <max |v_z|>
#   mean0/1/2  <体積加重平均の各成分>   … **符号**が効くので向きの誤りを捕まえる
#   lo / hi    <xcut で切った両側の v_<axis> の体積加重平均>
#              (-v axis=0|1|2 -v xcut=<座標> を与えたときだけ)
#              セルと座標の対応が崩れると値が入れ替わるので、
#              VTK の並べ替え (i 最内) の誤りを検出できる
#
# 構造格子 (RECTILINEAR_GRID) と非構造格子 (UNSTRUCTURED_GRID、VTK_TETRA) の
# どちらも扱う。VTK のセルの並びは i が最内なので、体積もその順に作る。

BEGIN { mode = ""; nx = ny = nz = 0; np = 0; nc = 0; want = 0; got = 0
        if (xcut == "") xcut = ""; if (axis == "") axis = 0 }

/^DATASET/            { mode = $2; next }
/^X_COORDINATES/      { rd = "x"; n = $2; k = 0; next }
/^Y_COORDINATES/      { rd = "y"; n = $2; k = 0; next }
/^Z_COORDINATES/      { rd = "z"; n = $2; k = 0; next }
/^POINTS/             { rd = "p"; n = $2; k = 0; next }
/^CELLS/              { rd = "c"; n = $2; k = 0; next }
/^CELL_TYPES/         { rd = ""; next }
/^(POINT_DATA|CELL_DATA|LOOKUP_TABLE)/ { rd = ""; next }

/^SCALARS/ { rd = ""; want = 0; next }
/^VECTORS/ { rd = ""; want = ($2 == arr); k = 0; next }

{
	if (rd == "x") { xc[k++] = $1; if (k == n) { nx = n; rd = "" }; next }
	if (rd == "y") { yc[k++] = $1; if (k == n) { ny = n; rd = "" }; next }
	if (rd == "z") { zc[k++] = $1; if (k == n) { nz = n; rd = "" }; next }
	if (rd == "p") { px[k] = $1; py[k] = $2; pz[k] = $3; k++;
	                 if (k == n) { np = n; rd = "" }; next }
	if (rd == "c") { c0[k] = $2; c1[k] = $3; c2[k] = $4; c3[k] = $5; k++;
	                 if (k == n) { nc = n; rd = "" }; next }
	if (want && (NF == 3)) { vx[got] = $1; vy[got] = $2; vz[got] = $3; got++; next }
	if (want && (NF != 3) && (got > 0)) { want = 0 }
}

END {
	if (got == 0) { print "ERROR array " arr " not found"; exit 1 }

	# セル体積
	if (mode == "RECTILINEAR_GRID") {
		m = 0
		for (k = 0; k < nz - 1; k++)
		for (j = 0; j < ny - 1; j++)
		for (i = 0; i < nx - 1; i++) {
			cv[m++] = (xc[i+1] - xc[i]) * (yc[j+1] - yc[j]) * (zc[k+1] - zc[k])
		}
		ncell = m
	}
	else {
		for (e = 0; e < nc; e++) {
			ax = px[c1[e]] - px[c0[e]]; ay = py[c1[e]] - py[c0[e]]; az = pz[c1[e]] - pz[c0[e]]
			bx = px[c2[e]] - px[c0[e]]; by = py[c2[e]] - py[c0[e]]; bz = pz[c2[e]] - pz[c0[e]]
			dx = px[c3[e]] - px[c0[e]]; dy = py[c3[e]] - py[c0[e]]; dz = pz[c3[e]] - pz[c0[e]]
			d = (ay*bz - az*by) * dx + (az*bx - ax*bz) * dy + (ax*by - ay*bx) * dz
			cv[e] = ((d < 0) ? -d : d) / 6
		}
		ncell = nc
	}
	if (ncell != got) { print "ERROR cell count " ncell " != vector count " got; exit 1 }

	# セル重心 (xcut を使うときだけ要る)
	if (xcut != "") {
		if (mode == "RECTILINEAR_GRID") {
			m = 0
			for (k = 0; k < nz - 1; k++)
			for (j = 0; j < ny - 1; j++)
			for (i = 0; i < nx - 1; i++) {
				gx[m] = (xc[i] + xc[i+1]) / 2
				gy[m] = (yc[j] + yc[j+1]) / 2
				gz[m] = (zc[k] + zc[k+1]) / 2
				m++
			}
		}
		else {
			for (e = 0; e < nc; e++) {
				gx[e] = (px[c0[e]] + px[c1[e]] + px[c2[e]] + px[c3[e]]) / 4
				gy[e] = (py[c0[e]] + py[c1[e]] + py[c2[e]] + py[c3[e]]) / 4
				gz[e] = (pz[c0[e]] + pz[c1[e]] + pz[c2[e]] + pz[c3[e]]) / 4
			}
		}
	}

	vol = 0; s2 = 0; vmin = -1; vmax = 0; a0 = a1 = a2 = 0
	m0 = m1 = m2 = 0
	vlo = vhi = 0; wlo = whi = 0
	for (e = 0; e < ncell; e++) {
		q = (vx[e]*vx[e]) + (vy[e]*vy[e]) + (vz[e]*vz[e])
		vol += cv[e]
		s2  += q * cv[e]
		m = sqrt(q)
		if ((vmin < 0) || (m < vmin)) vmin = m
		if (m > vmax) vmax = m
		t = (vx[e] < 0) ? -vx[e] : vx[e]; if (t > a0) a0 = t
		t = (vy[e] < 0) ? -vy[e] : vy[e]; if (t > a1) a1 = t
		t = (vz[e] < 0) ? -vz[e] : vz[e]; if (t > a2) a2 = t
		m0 += vx[e] * cv[e]; m1 += vy[e] * cv[e]; m2 += vz[e] * cv[e]
		if (xcut != "") {
			g = (axis == 0) ? gx[e] : ((axis == 1) ? gy[e] : gz[e])
			c = (axis == 0) ? vx[e] : ((axis == 1) ? vy[e] : vz[e])
			if (g < xcut) { vlo += c * cv[e]; wlo += cv[e] }
			else           { vhi += c * cv[e]; whi += cv[e] }
		}
	}
	printf "ncell %d\n", ncell
	printf "vol %.10e\n", vol
	printf "int2 %.10e\n", s2
	printf "vmin %.10e\n", vmin
	printf "vmax %.10e\n", vmax
	printf "amax0 %.10e\n", a0
	printf "amax1 %.10e\n", a1
	printf "amax2 %.10e\n", a2
	printf "mean0 %.10e\n", m0 / vol
	printf "mean1 %.10e\n", m1 / vol
	printf "mean2 %.10e\n", m2 / vol
	if (xcut != "") {
		printf "lo %.10e\n", ((wlo > 0) ? (vlo / wlo) : 0)
		printf "hi %.10e\n", ((whi > 0) ? (vhi / whi) : 0)
	}
}
