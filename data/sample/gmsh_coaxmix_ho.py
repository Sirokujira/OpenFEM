import gmsh, sys, math
# coaxmix.py <order> <nq> <nr (each ring)> <out> [tet] : 同軸の内側のリング (r = a..m) を三角形 -> 角柱、
# 外側のリング (r = m..b) を四角形 -> 六面体 (tet を付けると外側を四面体に)。界面 r = m で面を共有
order, nq, nr, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
outer_tet = (len(sys.argv) > 5 and sys.argv[5] == "tet")
ra, rm, rb, lz = 0.5e-3, 0.866e-3, 1.5e-3, 0.1e-3
gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0); gmsh.model.add("c")
g = gmsh.model.geo
c0 = g.addPoint(0, 0, 0)
def ring(r): return [g.addPoint(r * math.cos(k * math.pi / 2), r * math.sin(k * math.pi / 2), 0) for k in range(4)]
pa, pm, pb = ring(ra), ring(rm), ring(rb)
aa = [g.addCircleArc(pa[k], c0, pa[(k + 1) % 4]) for k in range(4)]
am = [g.addCircleArc(pm[k], c0, pm[(k + 1) % 4]) for k in range(4)]
ab = [g.addCircleArc(pb[k], c0, pb[(k + 1) % 4]) for k in range(4)]
r1 = [g.addLine(pa[k], pm[k]) for k in range(4)]
r2 = [g.addLine(pm[k], pb[k]) for k in range(4)]
s_in, s_out = [], []
for k in range(4):
    s_in.append(g.addPlaneSurface([g.addCurveLoop([r1[k], am[k], -r1[(k + 1) % 4], -aa[k]])]))
    s_out.append(g.addPlaneSurface([g.addCurveLoop([r2[k], ab[k], -r2[(k + 1) % 4], -am[k]])]))
for c in aa + am + ab: g.mesh.setTransfiniteCurve(c, nq + 1)
for c in r1 + r2: g.mesh.setTransfiniteCurve(c, nr + 1)
for s in s_out:
    g.mesh.setTransfiniteSurface(s)
    if not outer_tet: g.mesh.setRecombine(2, s)
for s in s_in: g.mesh.setTransfiniteSurface(s)	# 三角形 (再結合しない)
ex_in = g.extrude([(2, s) for s in s_in], 0, 0, lz, [1], recombine=True)
if outer_tet:
    ex_out = g.extrude([(2, s) for s in s_out], 0, 0, lz, [1], recombine=False)
else:
    ex_out = g.extrude([(2, s) for s in s_out], 0, 0, lz, [1], recombine=True)
g.synchronize()
vols = [e[1] for e in ex_in + ex_out if e[0] == 3]
gmsh.model.addPhysicalGroup(3, vols, 1)
inner, outer = [], []
for d, t in gmsh.model.getEntities(2):
    bb = gmsh.model.getBoundingBox(d, t)
    if abs(bb[2] - bb[5]) < 1e-12: continue
    if (abs(bb[0] - bb[3]) < 1e-12) or (abs(bb[1] - bb[4]) < 1e-12): continue
    rr = max(abs(bb[0]), abs(bb[3]), abs(bb[1]), abs(bb[4]))	# 1/4 周の円弧の面の半径
    if abs(rr - ra) < 1e-9: inner.append(t)
    elif abs(rr - rb) < 1e-9: outer.append(t)
gmsh.model.addPhysicalGroup(2, outer, 10)
gmsh.model.addPhysicalGroup(2, inner, 11)
gmsh.model.mesh.generate(3); gmsh.model.mesh.setOrder(order)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2); gmsh.write(out); gmsh.finalize()
