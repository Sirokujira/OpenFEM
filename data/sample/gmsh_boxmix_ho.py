import gmsh, sys
# boxmix.py <order> <out> : 1mm x 1mm x 0.2mm。z = 0..0.1mm を角柱 (三角形を押し出す)、z = 0.1..0.2mm を四面体。
# 角柱の上面の三角形を四面体が共有する (適合)。電極は z = 0 (10) と z = 0.2mm (11)
order, out = int(sys.argv[1]), sys.argv[2]
lx, ly, lz = 1e-3, 1e-3, 0.2e-3
gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0); gmsh.model.add("b"); g = gmsh.model.geo
p = [g.addPoint(0, 0, 0, 0.4e-3), g.addPoint(lx, 0, 0, 0.4e-3), g.addPoint(lx, ly, 0, 0.4e-3), g.addPoint(0, ly, 0, 0.4e-3)]
l = [g.addLine(p[i], p[(i + 1) % 4]) for i in range(4)]
s = g.addPlaneSurface([g.addCurveLoop(l)])
e1 = g.extrude([(2, s)], 0, 0, lz / 2, [1], recombine=True)
mid = e1[0][1]
e2 = g.extrude([(2, mid)], 0, 0, lz / 2)
g.synchronize()
vols = [e[1] for e in e1 + e2 if e[0] == 3]
gmsh.model.addPhysicalGroup(3, vols, 1)
gmsh.model.addPhysicalGroup(2, [s], 10)
gmsh.model.addPhysicalGroup(2, [e2[0][1]], 11)
gmsh.model.mesh.generate(3); gmsh.model.mesh.setOrder(order)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2); gmsh.write(out); gmsh.finalize()
