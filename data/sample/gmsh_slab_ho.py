# gmsh (4.15.2、Python API) で 2mm x 1mm x 0.2mm の板 (電極は x = 0 / x = 2mm、z は 2 層) を作る
#   python gmsh_slab_ho.py <hex|prism> <次数> slab_hex3.msh
import gmsh, sys
# slab.py <hex|prism> <order> <out> : 2mm x 1mm x 0.2mm、電極は x = 0 (10) / x = 2mm (11)、z は 2 層
kind, order, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
lx, ly, lz, nx, ny, nz = 2e-3, 1e-3, 0.2e-3, 4, 2, 2
gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0); gmsh.model.add("s"); g = gmsh.model.geo
p = [g.addPoint(0,0,0), g.addPoint(lx,0,0), g.addPoint(lx,ly,0), g.addPoint(0,ly,0)]
l = [g.addLine(p[i], p[(i+1)%4]) for i in range(4)]
s = g.addPlaneSurface([g.addCurveLoop(l)])
for i, c in enumerate(l): g.mesh.setTransfiniteCurve(c, (nx if i % 2 == 0 else ny) + 1)
g.mesh.setTransfiniteSurface(s)
if kind == "hex": g.mesh.setRecombine(2, s)
ex = g.extrude([(2, s)], 0, 0, lz, [nz], recombine=True); g.synchronize()
gmsh.model.addPhysicalGroup(3, [e[1] for e in ex if e[0] == 3], 1)
x0, x1 = [], []
for d, t in gmsh.model.getEntities(2):
    bb = gmsh.model.getBoundingBox(d, t)
    if abs(bb[0]) < 1e-12 and abs(bb[3]) < 1e-12: x0.append(t)
    if abs(bb[0] - lx) < 1e-12 and abs(bb[3] - lx) < 1e-12: x1.append(t)
gmsh.model.addPhysicalGroup(2, x0, 10); gmsh.model.addPhysicalGroup(2, x1, 11)
gmsh.model.mesh.generate(3); gmsh.model.mesh.setOrder(order)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2); gmsh.write(out); gmsh.finalize()
