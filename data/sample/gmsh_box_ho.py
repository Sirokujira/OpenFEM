# gmsh (4.15.2、Python API) で 1mm x 1mm x 0.2mm の箱を高次の六面体・角柱で切る (電極は z = 0 / z = 0.2mm)
#   python gmsh_box_ho.py <hex|prism> <次数> box_hex3.msh   (同じセッションで _41 / _bin_22 / _bin_41 も書く)
import gmsh, sys
def box(kind, order, out, nx=2, ny=2, nz=2, lx=1e-3, ly=1e-3, lz=0.2e-3, fmt=2.2):
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("b")
    p = [gmsh.model.geo.addPoint(0,0,0), gmsh.model.geo.addPoint(lx,0,0), gmsh.model.geo.addPoint(lx,ly,0), gmsh.model.geo.addPoint(0,ly,0)]
    l = [gmsh.model.geo.addLine(p[i], p[(i+1)%4]) for i in range(4)]
    cl = gmsh.model.geo.addCurveLoop(l); s = gmsh.model.geo.addPlaneSurface([cl])
    for i, c in enumerate(l): gmsh.model.geo.mesh.setTransfiniteCurve(c, (nx if i % 2 == 0 else ny) + 1)
    gmsh.model.geo.mesh.setTransfiniteSurface(s)
    if kind == "hex": gmsh.model.geo.mesh.setRecombine(2, s)
    ex = gmsh.model.geo.extrude([(2, s)], 0, 0, lz, [nz], recombine=True)
    gmsh.model.geo.synchronize()
    vol = [e[1] for e in ex if e[0] == 3]
    top = ex[0][1]
    gmsh.model.addPhysicalGroup(3, vol, 1)
    gmsh.model.addPhysicalGroup(2, [s], 10)
    gmsh.model.addPhysicalGroup(2, [top], 11)
    gmsh.model.mesh.generate(3); gmsh.model.mesh.setOrder(order)
    # 同じセッションで 4 形式を書く (節点の番号が揃うので結果がバイト一致するはず)
    for v, b, suf in ((2.2, 0, ""), (4.1, 0, "_41"), (2.2, 1, "_bin_22"), (4.1, 1, "_bin_41")):
        gmsh.option.setNumber("Mesh.MshFileVersion", v); gmsh.option.setNumber("Mesh.Binary", b)
        gmsh.write(out.replace(".msh", suf + ".msh"))
    gmsh.finalize()
kind, order, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
box(kind, order, out)
