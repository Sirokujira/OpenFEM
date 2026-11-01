# gmsh (4.15.2、Python API) で同軸 / 円環の高次の六面体・角柱の格子を作る (rlc_check.sh の [order-3 hex/prism])
#   python gmsh_coax_ho.py <hex|prism>[_axial] <次数> <1/4 周の分割> <半径方向の分割> <出力.msh>
#   coax_hex3a.msh = hex 3 2 1、coax_hex3b.msh = hex 3 4 2、ann_hex3a/b.msh = hex_axial 3 (2 1 / 4 2)、
#   ann_hex1b.msh = hex_axial 1 4 2 (角柱は prism)。_axial は電極を z = 0 / z = lz の面にして z を 2 層にする
import gmsh, sys, math
# coax.py <hex|prism> <order> <nq (1/4 周の分割)> <nr (半径方向)> <out> [<format> <binary>]
kind, order, nq, nr, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
axial = kind.endswith("_axial"); kind = kind.replace("_axial", "")
fmts = sys.argv[6:]
ra, rb, lz = 0.5e-3, 1.5e-3, 0.1e-3
gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0); gmsh.model.add("c")
g = gmsh.model.geo
c0 = g.addPoint(0, 0, 0)
pa = [g.addPoint(ra * math.cos(k * math.pi / 2), ra * math.sin(k * math.pi / 2), 0) for k in range(4)]
pb = [g.addPoint(rb * math.cos(k * math.pi / 2), rb * math.sin(k * math.pi / 2), 0) for k in range(4)]
aa = [g.addCircleArc(pa[k], c0, pa[(k + 1) % 4]) for k in range(4)]
ab = [g.addCircleArc(pb[k], c0, pb[(k + 1) % 4]) for k in range(4)]
rl = [g.addLine(pa[k], pb[k]) for k in range(4)]
surf = []
for k in range(4):
    cl = g.addCurveLoop([rl[k], ab[k], -rl[(k + 1) % 4], -aa[k]]); s = g.addPlaneSurface([cl]); surf.append(s)
for c in aa + ab: g.mesh.setTransfiniteCurve(c, nq + 1)
for c in rl: g.mesh.setTransfiniteCurve(c, nr + 1)
for s in surf:
    g.mesh.setTransfiniteSurface(s)
    if kind == "hex": g.mesh.setRecombine(2, s)
ex = g.extrude([(2, s) for s in surf], 0, 0, lz, [2 if axial else 1], recombine=True)
g.synchronize()
vols = [e[1] for e in ex if e[0] == 3]
gmsh.model.addPhysicalGroup(3, vols, 1)
# 側面 : 内側 (r = a) と外側 (r = b) を中心からの距離で分ける
inner, outer = [], []
for d, t in gmsh.model.getEntities(2):
    bb = gmsh.model.getBoundingBox(d, t)
    if abs(bb[2] - bb[5]) < 1e-12: continue
    rmax = max(math.hypot(bb[0], bb[1]), math.hypot(bb[3], bb[4]), math.hypot(bb[0], bb[4]), math.hypot(bb[3], bb[1]))
    if rmax < 0.5 * (ra + rb): inner.append(t)
    else:
        # 半径方向の面 (x = 0 や y = 0 の平面) を除く : 曲がった面だけ
        if (abs(bb[0] - bb[3]) < 1e-12) or (abs(bb[1] - bb[4]) < 1e-12): continue
        outer.append(t)
if axial:
    # 軸方向の電流 : 電極は z = 0 (10) と z = lz (11)。曲がった側面は露出面
    top = [e[1] for e in ex if e[0] == 2 and abs(gmsh.model.getBoundingBox(2, e[1])[2] - lz) < 1e-12 and abs(gmsh.model.getBoundingBox(2, e[1])[5] - lz) < 1e-12]
    gmsh.model.addPhysicalGroup(2, surf, 10)
    gmsh.model.addPhysicalGroup(2, top, 11)
else:
    gmsh.model.addPhysicalGroup(2, outer, 10)
    gmsh.model.addPhysicalGroup(2, inner, 11)
gmsh.model.mesh.generate(3); gmsh.model.mesh.setOrder(order)
if not fmts: fmts = ["2.2", "0"]
for i in range(0, len(fmts), 2):
    gmsh.option.setNumber("Mesh.MshFileVersion", float(fmts[i])); gmsh.option.setNumber("Mesh.Binary", int(fmts[i + 1]))
    o = out if i == 0 else out.replace(".msh", "_%s%s.msh" % (fmts[i].replace(".", ""), "_bin" if fmts[i + 1] == "1" else ""))
    gmsh.write(o)
gmsh.finalize()
