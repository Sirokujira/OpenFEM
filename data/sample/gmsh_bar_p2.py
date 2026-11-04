# bar_tet.msh (1 次の四面体) から 10 節点の四面体の格子 bar_tet2.msh を gmsh に作らせる (中間節点は直線の辺の中点)。
#   python gmsh_bar_p2.py bar_tet.msh bar_tet2.msh     (gmsh 4.15.2 の Python API)
# 2 次の A (eddy3d_p2.c) の検証用。物理タグはそのまま引き継がれる。bar_air_p2.msh も同じ (元は
# python3 mkmesh.py bar_air bar_air_c.msh -nx 8 -ny 1 -nz 8 -nza 2)
import gmsh, sys
src, dst = sys.argv[1], sys.argv[2]
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.open(src)
gmsh.model.mesh.setOrder(2)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.option.setNumber("Mesh.SaveAll", 0)
gmsh.write(dst)
gmsh.finalize()
