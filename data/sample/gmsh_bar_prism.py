# bar_eddy と同じ導体棒 (L(x) = 2 mm, W(y) = 0.25 mm, t(z) = 1 mm) を gmsh が次数 p の角柱 (完全 Lagrange) で切る
# (高次の角柱の A、eddy3d_pr.c の検証用)。物理タグは gmsh_bar_ho.py と同じ : 体積 1 (空気は 2)、x = 0 の面 10 (電極 0)、
# x = L の面 11 (電極 1)、z = 0 / z = t (空気つきは t + 0.5 mm) / x = 0 / x = L の面 20 (awall)
#   python gmsh_bar_prism.py <出力> <次数> <h [m]> <n> <z|y> [air] [twist]     (gmsh 4.15.2 の Python API)
#   z : x-y の長方形を三角形 (大きさ h) で切って z に n 層で押し出す (場は ζ 方向にだけ変わる : 面内の高次の関数は効かない)
#   y : x-z の長方形 (空気つきは導体と空気の 2 枚) を三角形で切って y に n 層で押し出す (場が三角形の面内で変わる)
#   x : y-z の長方形を三角形で切って x (電流の向き) に n 層で押し出す (A_x と ∇φ が ζ 方向の部分辺に載る。電極は三角形の面)
#   twist で要素ごとに局所の座標系を参照角柱の対称 (三角形の巡回 / ζ の反転 + λ1 と λ2 の入れ替え) で振り直す
#   (1 次の 6 節点の角柱で付け替えてから gmsh が次数を上げる。gmsh 4.15.2 は 3 次以上の角柱の局所座標を返せない)
import gmsh, sys, random
dst, order, h, n, ax = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
air = "air" in sys.argv[6:]
twist = "twist" in sys.argv[6:]
L, W, T = 2e-3, 0.25e-3, 1e-3
G = 0.5e-3
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.add("bar")
geo = gmsh.model.geo
def rect(x0, y0, z0, dx, dy, dz):
    # 長方形 (dy = 0 なら x-z、dz = 0 なら x-y) の 4 点と 4 辺 (辺は (x0 側, ...) の順)
    if dz == 0:
        c = [(x0, y0, z0), (x0 + dx, y0, z0), (x0 + dx, y0 + dy, z0), (x0, y0 + dy, z0)]
    else:
        c = [(x0, y0, z0), (x0 + dx, y0, z0), (x0 + dx, y0, z0 + dz), (x0, y0, z0 + dz)]
    return c
pts = {}
def pt(c):
    if c not in pts: pts[c] = geo.addPoint(c[0], c[1], c[2], h)
    return pts[c]
lines = {}
def ln(a, b):
    if (a, b) in lines: return lines[(a, b)]
    if (b, a) in lines: return -lines[(b, a)]
    lines[(a, b)] = geo.addLine(pt(a), pt(b))
    return lines[(a, b)]
def surf(c):
    return geo.addPlaneSurface([geo.addCurveLoop([ln(c[i], c[(i + 1) % 4]) for i in range(4)])])
vols = {1: [], 2: []}
faces = {10: [], 11: [], 20: []}
def classify(sd, x0, x1, z0, z1, mat):
    pass
if ax == "z":
    s = surf(rect(0, 0, 0, L, W, 0))
    ex = geo.extrude([(2, s)], 0, 0, T, [n], recombine=True)
    vols[1].append(ex[1][1])
    if air:
        ea = geo.extrude([ex[0]], 0, 0, G, [max(1, n // 2)], recombine=True)
        vols[2].append(ea[1][1])
elif ax == "x":
    def rect_yz(y0, z0, dy, dz):
        return [(0, y0, z0), (0, y0 + dy, z0), (0, y0 + dy, z0 + dz), (0, y0, z0 + dz)]
    ss = [(surf(rect_yz(0, 0, W, T)), 1)]
    if air:
        ss.append((surf(rect_yz(0, T, W, G)), 2))
    for (sf, m) in ss:
        ex = geo.extrude([(2, sf)], L, 0, 0, [n], recombine=True)
        vols[m].append(ex[1][1])
else:
    s1 = surf(rect(0, 0, 0, L, 0, T))
    ss = [(s1, 1)]
    if air:
        ss.append((surf(rect(0, 0, T, L, 0, G)), 2))
    for (sf, m) in ss:
        ex = geo.extrude([(2, sf)], 0, W, 0, [n], recombine=True)
        vols[m].append(ex[1][1])
geo.synchronize()
# 境界の面を座標で選ぶ (OCC と違い geo の境界箱は正確)
e = 1e-9
ztop = T + (G if air else 0)
def box(x0, y0, z0, x1, y1, z1):
    return [t for (d, t) in gmsh.model.getEntitiesInBoundingBox(x0 - e, y0 - e, z0 - e, x1 + e, y1 + e, z1 + e, 2)]
faces[10] = box(0, 0, 0, 0, W, T)
faces[11] = box(L, 0, 0, L, W, T)
faces[20] = box(0, 0, 0, L, W, 0) + box(0, 0, ztop, L, W, ztop) + box(0, 0, 0, 0, W, ztop) + box(L, 0, 0, L, W, ztop)
for m in (1, 2):
    if vols[m]: gmsh.model.addPhysicalGroup(3, vols[m], m)
for t in (10, 11, 20):
    gmsh.model.addPhysicalGroup(2, sorted(set(faces[t])), t)
gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 0)
gmsh.model.mesh.generate(3)
if twist:
    rots = [lambda u, v, z: (u, v, z), lambda u, v, z: (v, 1 - u - v, z), lambda u, v, z: (v, u, -z), lambda u, v, z: (1 - u - v, u, z)]
    for (d, b) in gmsh.model.getEntities(3):
        types, etags, ntags = gmsh.model.mesh.getElements(3, b)
        for (t, et, nt) in zip(types, etags, ntags):
            name, dim, eo, nn, lc, npn = gmsh.model.mesh.getElementProperties(t)
            pos = {(round(lc[3 * a] * 1e6), round(lc[3 * a + 1] * 1e6), round(lc[3 * a + 2] * 1e6)): a for a in range(nn)}
            new = []
            rng = random.Random(12345)
            for k in range(len(et)):
                # 擬似乱数で選ぶ (要素の順に巡回させると、押し出しの列の同じ層の隣どうしが同じ回転になって向きが食い違わない)
                r = rots[rng.randrange(len(rots))]
                old = nt[k * nn:(k + 1) * nn]
                for a in range(nn):
                    c = r(lc[3 * a], lc[3 * a + 1], lc[3 * a + 2])
                    new.append(old[pos[(round(c[0] * 1e6), round(c[1] * 1e6), round(c[2] * 1e6))]])
            gmsh.model.mesh.removeElements(3, b, list(et))
            gmsh.model.mesh.addElements(3, b, [t], [list(et)], [new])
gmsh.model.mesh.setOrder(order)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.option.setNumber("Mesh.SaveAll", 0)
gmsh.write(dst)
gmsh.finalize()
