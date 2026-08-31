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
#   vlo / vhi  <両側の体積>、ilo / ihi <両側の ∫ v_<axis> dV>
#              (-v axis=0|1|2 -v xcut=<座標> を与えたときだけ)
#              積分する成分は -v comp=0|1|2 で別に指定できる (既定は axis)。
#              分割軸と成分が違う場合 (例: y で切って z 成分を積分) に要る
#              セルと座標の対応が崩れると値が入れ替わるので、
#              VTK の並べ替え (i 最内) の誤りを検出できる
#
# 構造格子 (RECTILINEAR_GRID) と非構造格子 (UNSTRUCTURED_GRID) の どちらも扱う。
# VTK のセルの並びは i が最内なので、体積もその順に作る。
#
# 非構造格子のセルは **CELL_TYPES で見分ける**:
#   10 / 24 -> 四面体 (VTK_TETRA / VTK_QUADRATIC_TETRA) -> 体積 (頂点 4 個)
#    5 / 22 -> 三角形 (VTK_TRIANGLE / VTK_QUADRATIC_TRIANGLE) -> **面積** (頂点 3 個)
#   12      -> 六面体 (VTK_HEXAHEDRON) -> 2x2x2 Gauss で積分した体積
#   13      -> 角柱   (VTK_WEDGE)      -> 三角形 3 点 x ζ 2 点 Gauss
#   14      -> ピラミッド (VTK_PYRAMID) -> 立方体座標の 2x2x2 Gauss (ζ は [0,1])
#
# **節点数だけでは見分けられない**: 角柱と 2 次三角形はどちらも 6 節点なので、
# 型を読まないと角柱の体積を「三角形の面積」として足してしまう
# (実測: 六面体 576 + 角柱 1152 の格子で 5.0e-10 のところ 6.0e-06 を返した)。
#
# 六面体・角柱・ピラミッドの体積は**ソルバーと同じ求積**で出す。等パラメトリック
# 写像の体積は ∫det J であって、四面体に割った和ではない (面が平面でない要素では
# 一致しない)。同じ求積を使えば、体積の恒等式が「格子の形」ではなく実装を見る。
# 四面体・三角形は先頭の頂点だけで測る (直線要素では厳密。曲がった 2 次要素では
# 内接多角形になるので、その場合は体積そのものを恒等式に使わないこと)
# 断面 2 次元は単位長あたりで扱うので、面積がそのまま「体積」になり
# ∫|J|²/(2σ) dA = ½Re(Y) のような恒等式がそのまま使える。


# --- 求積で使う道具 ---

# 節点 a, b, c, d の四面体の体積
function tetvol(a, b, c, d,   ax, ay, az, bx, by, bz, dx, dy, dz, t) {
	ax = px[b] - px[a]; ay = py[b] - py[a]; az = pz[b] - pz[a]
	bx = px[c] - px[a]; by = py[c] - py[a]; bz = pz[c] - pz[a]
	dx = px[d] - px[a]; dy = py[d] - py[a]; dz = pz[d] - pz[a]
	t = ((ay*bz - az*by) * dx) + ((az*bx - ax*bz) * dy) + ((ax*by - ay*bx) * dz)

	return ((t < 0) ? -t : t) / 6
}

# 形状関数の微分 dn[a,j] (呼び出し側で埋める) と節点から det J を作る
function detj(e, nen,   i, j, a, w, jm, d, d0, d1, d2) {
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) jm[i, j] = 0
	}
	for (a = 0; a < nen; a++) {
		w = cnod[e, a]
		for (j = 0; j < 3; j++) {
			jm[0, j] += px[w] * dn[a, j]
			jm[1, j] += py[w] * dn[a, j]
			jm[2, j] += pz[w] * dn[a, j]
		}
	}
	# **式を複数行に分けないこと。** awk の文は改行で終わるので、次の行の
	# "- (...)" は別の文になり、第 1 項だけの値になる (実測: ピラミッドの
	# 底面が x = 一定の面だと jm[0,0] = 0 で体積が 0 になり、box_pyr の総体積が
	# 解析値の 2/3 になった)。小行列式を変数に落として 1 行で書く
	d0 = (jm[1,1] * jm[2,2]) - (jm[1,2] * jm[2,1])
	d1 = (jm[1,0] * jm[2,2]) - (jm[1,2] * jm[2,0])
	d2 = (jm[1,0] * jm[2,1]) - (jm[1,1] * jm[2,0])
	d = (jm[0,0] * d0) - (jm[0,1] * d1) + (jm[0,2] * d2)

	return ((d < 0) ? -d : d)
}

# 六面体 (3 線形、節点の並びは VTK = Gmsh)。2x2x2 Gauss、重み 1
function hexvol(e,   gp, i, j, k, a, sx, sy, sz, xi, et, ze, v) {
	gp = 1 / sqrt(3)
	v = 0
	for (i = 0; i < 2; i++) {
	for (j = 0; j < 2; j++) {
	for (k = 0; k < 2; k++) {
		xi = (i ? gp : -gp); et = (j ? gp : -gp); ze = (k ? gp : -gp)
		for (a = 0; a < 8; a++) {
			sx = HSX[a]; sy = HSY[a]; sz = HSZ[a]
			dn[a, 0] = sx * (1 + (sy * et)) * (1 + (sz * ze)) / 8
			dn[a, 1] = sy * (1 + (sx * xi)) * (1 + (sz * ze)) / 8
			dn[a, 2] = sz * (1 + (sx * xi)) * (1 + (sy * et)) / 8
		}
		v += detj(e, 8)
	}
	}
	}

	return v
}

# 角柱 (下面 0-2 / 上面 3-5)。三角形 3 点 (重み 1/6) x ζ 2 点 Gauss
function prismvol(e,   gp, q, k, a, t, u, v2, ze, s, h, lam, v) {
	gp = 1 / sqrt(3)
	v = 0
	for (q = 0; q < 3; q++) {
		u = TU[q]; v2 = TV[q]
		lam[0] = 1 - u - v2; lam[1] = u; lam[2] = v2
		for (k = 0; k < 2; k++) {
			ze = (k ? gp : -gp)
			for (a = 0; a < 6; a++) {
				t = a % 3
				s = ((a < 3) ? -1 : 1)
				h = (1 + (s * ze)) / 2
				dn[a, 0] = DLU[t] * h
				dn[a, 1] = DLV[t] * h
				dn[a, 2] = lam[t] * s / 2
			}
			v += detj(e, 6) / 6
		}
	}

	return v
}

# ピラミッド (底面 0-3、頂点 4)。立方体座標 (u, v ∈ [-1,1]、ζ ∈ [0,1]) の
# 2x2x2 Gauss。ζ 方向の重みが 1/2 (区間が [0,1] なので)
function pyrvol(e,   gp, i, j, k, a, su, sv, u, v2, ze, v) {
	gp = 1 / sqrt(3)
	v = 0
	for (i = 0; i < 2; i++) {
	for (j = 0; j < 2; j++) {
	for (k = 0; k < 2; k++) {
		u = (i ? gp : -gp); v2 = (j ? gp : -gp)
		ze = (k ? (0.5 + (0.5 * gp)) : (0.5 - (0.5 * gp)))
		for (a = 0; a < 4; a++) {
			su = PSU[a]; sv = PSV[a]
			dn[a, 0] = su * (1 + (sv * v2)) * (1 - ze) / 4
			dn[a, 1] = sv * (1 + (su * u)) * (1 - ze) / 4
			dn[a, 2] = -(1 + (su * u)) * (1 + (sv * v2)) / 4
		}
		dn[4, 0] = 0; dn[4, 1] = 0; dn[4, 2] = 1
		v += detj(e, 5) / 2
	}
	}
	}

	return v
}

BEGIN { mode = ""; nx = ny = nz = 0; np = 0; nc = 0; want = 0; got = 0
        if (xcut == "") xcut = ""; if (axis == "") axis = 0
        if (comp == "") comp = axis
        # 六面体の頂点の符号 (VTK = Gmsh の並び。sol/unstruct.c の HEX_SGN と同じ)
        HSX[0] = -1; HSY[0] = -1; HSZ[0] = -1
        HSX[1] = +1; HSY[1] = -1; HSZ[1] = -1
        HSX[2] = +1; HSY[2] = +1; HSZ[2] = -1
        HSX[3] = -1; HSY[3] = +1; HSZ[3] = -1
        HSX[4] = -1; HSY[4] = -1; HSZ[4] = +1
        HSX[5] = +1; HSY[5] = -1; HSZ[5] = +1
        HSX[6] = +1; HSY[6] = +1; HSZ[6] = +1
        HSX[7] = -1; HSY[7] = +1; HSZ[7] = +1
        # 角柱の三角形 3 点則 (次数 2) と面積座標の微分
        TU[0] = 1/6; TV[0] = 1/6
        TU[1] = 2/3; TV[1] = 1/6
        TU[2] = 1/6; TV[2] = 2/3
        DLU[0] = -1; DLU[1] = 1; DLU[2] = 0
        DLV[0] = -1; DLV[1] = 0; DLV[2] = 1
        # ピラミッドの底面の符号
        PSU[0] = -1; PSU[1] = 1; PSU[2] = 1; PSU[3] = -1
        PSV[0] = -1; PSV[1] = -1; PSV[2] = 1; PSV[3] = 1 }

/^DATASET/            { mode = $2; next }
/^X_COORDINATES/      { rd = "x"; n = $2; k = 0; next }
/^Y_COORDINATES/      { rd = "y"; n = $2; k = 0; next }
/^Z_COORDINATES/      { rd = "z"; n = $2; k = 0; next }
/^POINTS/             { rd = "p"; n = $2; k = 0; next }
/^CELLS/              { rd = "c"; n = $2; k = 0; next }
/^CELL_TYPES/         { rd = "t"; n = $2; k = 0; next }
/^(POINT_DATA|CELL_DATA|LOOKUP_TABLE)/ { rd = ""; next }

/^SCALARS/ { rd = ""; want = 0; next }
/^VECTORS/ { rd = ""; want = ($2 == arr); k = 0; next }

{
	if (rd == "x") { xc[k++] = $1; if (k == n) { nx = n; rd = "" }; next }
	if (rd == "y") { yc[k++] = $1; if (k == n) { ny = n; rd = "" }; next }
	if (rd == "z") { zc[k++] = $1; if (k == n) { nz = n; rd = "" }; next }
	if (rd == "p") { px[k] = $1; py[k] = $2; pz[k] = $3; k++;
	                 if (k == n) { np = n; rd = "" }; next }
	if (rd == "c") { cn[k] = $1
	                 for (i = 0; i < $1; i++) cnod[k, i] = $(i + 2)
	                 c0[k] = $2; c1[k] = $3; c2[k] = $4; c3[k] = $5; k++
	                 if (k == n) { nc = n; rd = "" }; next }
	if (rd == "t") { ct[k++] = $1; if (k == n) rd = ""; next }
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
			ty = ct[e]
			if ((ty == 5) || (ty == 22)) {
				# 三角形 : 面積 = |a x b| / 2 (単位長あたりなのでこれが「体積」)
				ax = px[c1[e]] - px[c0[e]]; ay = py[c1[e]] - py[c0[e]]; az = pz[c1[e]] - pz[c0[e]]
				bx = px[c2[e]] - px[c0[e]]; by = py[c2[e]] - py[c0[e]]; bz = pz[c2[e]] - pz[c0[e]]
				nx1 = (ay*bz - az*by); ny1 = (az*bx - ax*bz); nz1 = (ax*by - ay*bx)
				cv[e] = sqrt((nx1*nx1) + (ny1*ny1) + (nz1*nz1)) / 2
			}
			else if ((ty == 10) || (ty == 24)) {
				cv[e] = tetvol(c0[e], c1[e], c2[e], c3[e])
			}
			else if (ty == 12) { cv[e] = hexvol(e) }
			else if (ty == 13) { cv[e] = prismvol(e) }
			else if (ty == 14) { cv[e] = pyrvol(e) }
			else {
				printf "ERROR unknown VTK cell type %d\n", ty
				exit 1
			}
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
				# 頂点だけで重心を出す (中間節点は含めない)。切断面のどちら側かを
				# 決めるだけなので、要素中心の厳密な位置は要らない
				ty = ct[e]
				nn = 4
				if ((ty == 5) || (ty == 22)) nn = 3
				if (ty == 12) nn = 8
				if (ty == 13) nn = 6
				if (ty == 14) nn = 5
				sx2 = sy2 = sz2 = 0
				for (a = 0; a < nn; a++) {
					sx2 += px[cnod[e, a]]; sy2 += py[cnod[e, a]]; sz2 += pz[cnod[e, a]]
				}
				gx[e] = sx2 / nn; gy[e] = sy2 / nn; gz[e] = sz2 / nn
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
			c = (comp == 0) ? vx[e] : ((comp == 1) ? vy[e] : vz[e])
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
		printf "vlo %.10e\n", wlo		# 切断面より下側の体積
		printf "vhi %.10e\n", whi		# 上側の体積 (積分 = lo*vlo, hi*vhi)
		printf "ilo %.10e\n", vlo		# ∫ v_<axis> dV (下側)
		printf "ihi %.10e\n", vhi		# ∫ v_<axis> dV (上側)
	}
}
