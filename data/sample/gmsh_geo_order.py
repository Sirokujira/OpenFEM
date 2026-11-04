# .geo から高次の四面体 / 三角形の格子を gmsh (4.15.2、Python API) に作らせる (gmsh のコマンドが無い環境用)。
#   python gmsh_geo_order.py <geo> <次元> <次数> <出力.msh> [<版 2.2|4.1> <バイナリ 0|1>]
#   box_p5.msh       = box_p3.geo 3 5 (box_p5_41.msh = ... 4.1 0、box_p5_bin_41.msh = ... 4.1 1)
#   coax_p5.msh      = coax3.geo の先頭に "h = 0.4e-3;" を足したもの 3 5 (同じ手順の 4 次が coax_p4.msh と全バイト一致)
#   plate2d_p5.msh   = plate2d_p3.geo 2 5
import gmsh, sys
geo, dim, order, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
ver = float(sys.argv[5]) if len(sys.argv) > 5 else 2.2
binary = int(sys.argv[6]) if len(sys.argv) > 6 else 0
gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
gmsh.open(geo); gmsh.model.mesh.generate(dim); gmsh.model.mesh.setOrder(order)
gmsh.option.setNumber("Mesh.MshFileVersion", ver); gmsh.option.setNumber("Mesh.Binary", binary)
gmsh.option.setNumber("Mesh.SaveAll", 0)
gmsh.write(out); gmsh.finalize()
