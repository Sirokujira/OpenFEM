# mesh_downgrade.awk : gmsh の 2.2 ASCII の格子で、座標 (列 col、既定 2 = x) の重心が xc より小さい四面体 (型 tf) と三角形
# (型 ff) を低い次数に落とす (次数の違う要素が面を共有する格子の検証用。元は gmsh が作った直線 / 曲がった高次の格子)。
#   awk -v col=<2|3|4> -v xc=<しきい値> -v tf=<型> -v ff=<型> -v to=<1|2> -f mesh_downgrade.awk in.msh > out.msh
#   to = 1 : 頂点だけ (型 4 / 2)。to = 2 : 4 次 (型 30 / 23) の辺の 3 個の真ん中を中間節点にした 10 / 6 節点 (型 11 / 9)
#   使われなくなった節点はそのまま残す (行列の行が空になるだけ)
/^\$Nodes/ { nd = 1; print; getline; print; next }
/^\$EndNodes/ { nd = 0 }
nd { X[$1] = $((col > 0) ? col : 2) }
/^\$Elements/ { el = 1; print; getline; print; next }
/^\$EndElements/ { el = 0 }
el && (($2 == tf) || ($2 == ff)) {
	nt = $3; tet = ($2 == tf); nv = tet ? 4 : 3; b = 3 + nt; ne = tet ? 6 : 3
	c = 0; for (i = 1; i <= nv; i++) c += X[$(b + i)]; c /= nv
	if (c < xc) {
		s = $1 " " ((to == 1) ? (tet ? 4 : 2) : (tet ? 11 : 9)) " " nt
		for (i = 1; i <= nt; i++) s = s " " $(3 + i)
		for (i = 1; i <= nv; i++) s = s " " $(b + i)
		if (to == 2) { for (k = 0; k < ne; k++) s = s " " $(b + nv + 3 * k + 2) }
		print s; next
	}
}
{ print }
