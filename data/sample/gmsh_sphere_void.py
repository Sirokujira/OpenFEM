# 導体の直方体 (x 方向 10 mm、断面 3 x 3 mm) の中の球の空隙 (半径 1 mm) : gmsh (OCC) が曲がった高次の四面体にする
#   sphere_void_p2.msh / _p3.msh = 2 / 3 0.6e-3 をバイナリ 4.1 に変換したもの (gmsh の open + write、Mesh.Binary = 1)
#   python gmsh_sphere_void.py <次数> <要素の大きさ [m]> <出力.msh>
import gmsh, sys
order, h, out = int(sys.argv[1]), float(sys.argv[2]), sys.argv[3]
gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0); gmsh.model.add("v")
occ = gmsh.model.occ
b = occ.addBox(0, 0, 0, 10e-3, 3e-3, 3e-3)
sp = occ.addSphere(5e-3, 1.5e-3, 1.5e-3, 1e-3)
vol, _ = occ.cut([(3, b)], [(3, sp)])
occ.synchronize()
v = vol[0][1]
gmsh.model.addPhysicalGroup(3, [v], 1)
x0, x1 = [], []
for (d, t) in gmsh.model.getBoundary([(3, v)], oriented=False):
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(d, t)
    xm = 0.5 * (xmin + xmax); flat = (xmax - xmin) < 1e-6
    if flat and abs(xm) < 1e-6: x0.append(t)
    if flat and abs(xm - 10e-3) < 1e-6: x1.append(t)
gmsh.model.addPhysicalGroup(2, x0, 10); gmsh.model.addPhysicalGroup(2, x1, 11)
gmsh.option.setNumber("Mesh.MeshSizeMax", h); gmsh.option.setNumber("Mesh.MeshSizeMin", h)
gmsh.model.mesh.generate(3); gmsh.model.mesh.setOrder(order)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2); gmsh.option.setNumber("Mesh.SaveAll", 0)
gmsh.write(out); gmsh.finalize()
