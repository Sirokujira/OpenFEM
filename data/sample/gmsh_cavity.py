# 開いた立方体の空洞 (thermal_cavity.ofe と同じ形) を gmsh (OCC) で四面体にする (高次の面の相互放射の検査)
#   python gmsh_cavity.py <次数> <要素の大きさ [m]> <出力.msh>   thermal_cavity_tet2.msh = 2 1.0e-3
import gmsh, sys
order, h, out = int(sys.argv[1]), float(sys.argv[2]), sys.argv[3]
gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0); gmsh.model.add("c")
occ = gmsh.model.occ
b = occ.addBox(0, 0, 0, 10e-3, 3e-3, 3e-3)
cv = occ.addBox(4e-3, 0.5e-3, 1e-3, 2e-3, 2e-3, 2e-3)
vol, _ = occ.cut([(3, b)], [(3, cv)])
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
gmsh.option.setNumber("Mesh.MeshSizeMax", h)
gmsh.model.mesh.generate(3); gmsh.model.mesh.setOrder(order)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2); gmsh.option.setNumber("Mesh.SaveAll", 0)
gmsh.write(out); gmsh.finalize()
