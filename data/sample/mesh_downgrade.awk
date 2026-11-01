# mesh_downgrade.awk : gmsh の 2.2 ASCII の格子で、座標 (列 col、既定 2 = x) の重心が xc より小さい要素 (型 tf) と面
# (三角形の型 ff、四角形の型 qf) を低い次数に落とす (次数の違う要素が面を共有する格子の検証用。元は gmsh が作った
# 一様な高次の格子)。kind = tet (既定) / hex / prism は体積要素の種別。
#   awk -v col=<2|3|4> -v xc=<しきい値> -v tf=<型> -v ff=<型> [-v qf=<型>] [-v kind=hex|prism] -v to=<1|2> -f mesh_downgrade.awk in.msh > out.msh
#   to = 1 : 頂点だけ (四面体 4 / 六面体 5 / 角柱 6、三角形 2 / 四角形 3)。to = 2 (四面体だけ) : 4 次 (型 30 / 23) の
#   辺の 3 個の真ん中を中間節点にした 10 / 6 節点 (型 11 / 9)
#   使われなくなった節点はそのまま残す (行列の行が空になるだけ)
BEGIN { if (qf == "") qf = -1; nvv = (kind == "hex") ? 8 : (kind == "prism") ? 6 : 4; tv = (kind == "hex") ? 5 : (kind == "prism") ? 6 : 4 }
/^\$Nodes/ { nd = 1; print; getline; print; next }
/^\$EndNodes/ { nd = 0 }
nd { X[$1] = $((col > 0) ? col : 2) }
/^\$Elements/ { el = 1; print; getline; print; next }
/^\$EndElements/ { el = 0 }
el && (($2 == tf) || ($2 == ff) || ($2 == qf)) {
	nt = $3; vol = ($2 == tf); quad = ($2 == qf); nv = vol ? nvv : quad ? 4 : 3; b = 3 + nt; ne = vol ? 6 : 3
	c = 0; for (i = 1; i <= nv; i++) c += X[$(b + i)]; c /= nv
	if (c < xc) {
		s = $1 " " ((to == 1) ? (vol ? tv : quad ? 3 : 2) : (vol ? 11 : 9)) " " nt
		for (i = 1; i <= nt; i++) s = s " " $(3 + i)
		for (i = 1; i <= nv; i++) s = s " " $(b + i)
		if (to == 2) { for (k = 0; k < ne; k++) s = s " " $(b + nv + 3 * k + 2) }
		print s; next
	}
}
{ print }
