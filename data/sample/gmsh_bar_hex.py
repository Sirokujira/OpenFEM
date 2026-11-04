# bar_eddy と同じ導体棒 (L(x) = 2 mm, W(y) = 0.25 mm, t(z) = 1 mm) を gmsh の幾何から次数 p の六面体 (完全 Lagrange) で切る
# (高次の六面体の A、eddy3d_pr.c の部分辺の Nedelec の検証用)。物理タグは gmsh_bar_ho.py と同じ :
#   体積 1、x = 0 の面 10 (電極 0)、x = L の面 11 (電極 1)、z = 0 / z = t / x = 0 / x = L の面 20 (awall)
#   python gmsh_bar_hex.py <出力> <次数> <nx> <ny> <nz> [air [yawall]]     (gmsh 4.15.2 の Python API)
# air で gmsh_bar_ho.py と同じ空気層 (z = t .. t + 0.5 mm、体積 2、z に nz/2 区間)、yawall で電極の面を awall にせず
# y 面を awall にする (ゲージ固定の検査用)。twist (6 番目の引数) で要素ごとに局所の座標系を参照立方体の回転 (90° を
# 要素ごとに順に x / y / z 軸のまわり、恒等も含む) で振り直す : transfinite の格子は全要素の局所の向きが揃っているので、
# 隣と向きが食い違う格子でないと辺の向きの符号の誤りが検出できない (gmsh の API で節点の並びを付け替える)
import gmsh, sys
dst, order = sys.argv[1], int(sys.argv[2])
nx, ny, nz = int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
air = (len(sys.argv) > 6) and (sys.argv[6] == "air")
twist = (len(sys.argv) > 6) and (sys.argv[6] == "twist")
yawall = air and (len(sys.argv) > 7) and (sys.argv[7] == "yawall")
L, W, T = 2e-3, 0.25e-3, 1e-3
G = 0.5e-3 if air else 0
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.add("bar")
v = gmsh.model.occ.addBox(0, 0, 0, L, W, T)
if air:
    a = gmsh.model.occ.addBox(0, 0, T, L, W, G)
    gmsh.model.occ.fragment([(3, v)], [(3, a)])
gmsh.model.occ.synchronize()
def faces(c, x, ztop=None, zbot=0):
    e = 2e-7
    lo = [-e, -e, zbot - e]; hi = [L + e, W + e, (T + G if ztop is None else ztop) + e]
    lo[c] = x - e; hi[c] = x + e
    return [s for (d, s) in gmsh.model.getEntitiesInBoundingBox(lo[0], lo[1], lo[2], hi[0], hi[1], hi[2], 2)]
def vols(z0, z1):
    e = 2e-7
    return [s for (d, s) in gmsh.model.getEntitiesInBoundingBox(-e, -e, z0 - e, L + e, W + e, z1 + e, 3)]
# 構造格子 (曲線を区間数で切り、面・体積を transfinite + recombine で六面体に)
for (d, c) in gmsh.model.getEntities(1):
    x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(1, c)
    dx, dy, dz = x1 - x0, y1 - y0, z1 - z0
    if dx >= max(dy, dz): n = nx
    elif dy >= max(dx, dz): n = ny
    else: n = nz if z1 <= T + 1e-6 else max(1, nz // 2)
    gmsh.model.mesh.setTransfiniteCurve(c, n + 1)
for (d, f) in gmsh.model.getEntities(2):
    gmsh.model.mesh.setTransfiniteSurface(f)
    gmsh.model.mesh.setRecombine(2, f)
for (d, b) in gmsh.model.getEntities(3):
    gmsh.model.mesh.setTransfiniteVolume(b)
    gmsh.model.mesh.setRecombine(3, b)
gmsh.model.addPhysicalGroup(3, vols(0, T), 1)
if air:
    gmsh.model.addPhysicalGroup(3, vols(T, T + G), 2)
gmsh.model.addPhysicalGroup(2, faces(0, 0, T), 10)
gmsh.model.addPhysicalGroup(2, faces(0, L, T), 11)
if yawall:
    gmsh.model.addPhysicalGroup(2, faces(2, 0) + faces(2, T + G) + faces(0, 0, None, T) + faces(0, L, None, T)
                                + faces(1, 0) + faces(1, W), 20)
else:
    gmsh.model.addPhysicalGroup(2, faces(2, 0) + faces(2, T + G) + faces(0, 0) + faces(0, L), 20)
gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 0)
gmsh.model.mesh.generate(3)
gmsh.model.mesh.setOrder(order)
if twist:
    rots = [lambda x, y, z: (x, y, z), lambda x, y, z: (x, -z, y), lambda x, y, z: (z, y, -x), lambda x, y, z: (-y, x, z)]
    for (d, b) in gmsh.model.getEntities(3):
        types, etags, ntags = gmsh.model.mesh.getElements(3, b)
        for (t, et, nt) in zip(types, etags, ntags):
            name, dim, eo, nn, lc, npn = gmsh.model.mesh.getElementProperties(t)
            pos = {(round(lc[3 * a] * 1e6), round(lc[3 * a + 1] * 1e6), round(lc[3 * a + 2] * 1e6)): a for a in range(nn)}
            new = []
            for k in range(len(et)):
                r = rots[k % len(rots)]
                old = nt[k * nn:(k + 1) * nn]
                for a in range(nn):
                    c = r(lc[3 * a], lc[3 * a + 1], lc[3 * a + 2])
                    new.append(old[pos[(round(c[0] * 1e6), round(c[1] * 1e6), round(c[2] * 1e6))]])
            gmsh.model.mesh.removeElements(3, b, list(et))
            gmsh.model.mesh.addElements(3, b, [t], [list(et)], [new])
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.option.setNumber("Mesh.SaveAll", 0)
gmsh.write(dst)
gmsh.finalize()
