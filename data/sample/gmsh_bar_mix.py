# bar_eddy と同じ導体棒 (L(x) = 2 mm, W(y) = 0.25 mm, t(z) = 1 mm) を種別の混在した次数 p の格子で切る (高次の辺要素の種別混在の
# 検証用、eddy3d_pr.c の「種別の混在」)。物理タグは gmsh_bar_ho.py と同じ (体積 1、x = 0 の面 10、x = L の面 11、z = 0 / z = t /
# x = 0 / x = L の面 20)。
#   python gmsh_bar_mix.py <出力> <次数> <h [m]> hp [twist]   z < t/2 を六面体 (x-z の四角形を y に押し出す)、z > t/2 を角柱
#                                                         (x-z の三角形を y に押し出す)。界面は z = t/2 の四角形の面
#   python gmsh_bar_mix.py <出力> <次数> <h [m]> pt [twist]   z < t/2 を角柱 (x-y の三角形を z に 2 層で押し出す)、z > t/2 を
#                                                         四面体。界面は z = t/2 の三角形の面
#   python gmsh_bar_mix.py <出力> <次数> <h [m]> tri [twist]  四面体・六面体・角柱の 3 種別 (六面体と四面体は角柱を挟んで接しない)
#   air (5 番目以降の引数) で gmsh_bar_ho.py と同じ空気層 (z = t .. t + 0.5 mm、体積 2) を載せる (hp は角柱、pt は四面体。
#   導体と空気で種別の組が違うので、種別ごとの材料の添字の誤りが効く)
#   yawall で電極の面を awall にせず y の面を awall にする (ゲージの「電極の節点を 1 つの根に」の検査用)
#   twist で要素ごとに局所の座標系を擬似乱数の回転で振り直す (1 次の要素で付け替えてから gmsh が次数を上げる)
import gmsh, sys, random
dst, order, h, mode = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
twist = "twist" in sys.argv[5:]
air = "air" in sys.argv[5:]
yawall = "yawall" in sys.argv[5:]
L, W, T = 2e-3, 0.25e-3, 1e-3
G = 0.5e-3 if air else 0
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.add("bar")
geo = gmsh.model.geo
vols = []
if mode == "hp":
    # z < t/2 を六面体、z > t/2 を角柱 (どちらも x-z を y に 1 層で押し出す)。界面 z = t/2 は四角形の面で、電流 (x 向き) が
    # 界面の接線方向に流れる (界面の接線のトレースが効く。x = L/2 で割ると A が界面の法線方向で検査にならない)
    nx, nz = int(round(L / h)), max(1, int(round(T / 2 / h)))
    P = [geo.addPoint(x, 0, z, h) for (x, z) in [(0, 0), (L, 0), (L, T / 2), (0, T / 2), (L, T), (0, T)]]
    l = [geo.addLine(P[0], P[1]), geo.addLine(P[1], P[2]), geo.addLine(P[2], P[3]), geo.addLine(P[3], P[0]),
         geo.addLine(P[2], P[4]), geo.addLine(P[4], P[5]), geo.addLine(P[5], P[3])]
    sl = geo.addPlaneSurface([geo.addCurveLoop([l[0], l[1], l[2], l[3]])])
    su = geo.addPlaneSurface([geo.addCurveLoop([-l[2], l[4], l[5], l[6]])])
    for c in (l[0], l[2]): geo.mesh.setTransfiniteCurve(c, nx + 1)
    for c in (l[1], l[3]): geo.mesh.setTransfiniteCurve(c, nz + 1)
    geo.mesh.setTransfiniteSurface(sl)
    geo.mesh.setRecombine(2, sl)
    sfs = [sl, su]
    if air:
        Pa = [geo.addPoint(L, 0, T + G, h), geo.addPoint(0, 0, T + G, h)]
        sfs.append(geo.addPlaneSurface([geo.addCurveLoop([-l[5], geo.addLine(P[4], Pa[0]), geo.addLine(Pa[0], Pa[1]), geo.addLine(Pa[1], P[5])])]))
    for sf in sfs:
        ex = geo.extrude([(2, sf)], 0, W, 0, [1], recombine=True)
        vols.append(ex[1][1])
elif mode == "tri":
    # 3 種別 : x > L/2, z < t/2 を角柱 (x-y の三角形を z に押し出す)、その上 (z > t/2) を四面体、x < L/2 を x に押し出した六面体
    # (z < t/2、角柱の x = L/2 の側面の四角形から) と角柱 (z > t/2、四面体の x = L/2 の三角形の面から)。六面体と四面体は
    # 接しない (角柱の側面と上下の面を挟む)
    nx2, ny, nz2 = max(1, int(round(L / 2 / h))), max(1, int(round(W / h))), max(1, int(round(T / 2 / h)))
    p0 = [geo.addPoint(L / 2, 0, 0, h), geo.addPoint(L, 0, 0, h), geo.addPoint(L, W, 0, h), geo.addPoint(L / 2, W, 0, h)]
    s0 = geo.addPlaneSurface([geo.addCurveLoop([geo.addLine(p0[i], p0[(i + 1) % 4]) for i in range(4)])])
    eq = geo.extrude([(2, s0)], 0, 0, T / 2, [nz2], recombine=True)
    et = geo.extrude([eq[0]], 0, 0, T / 2)
    geo.synchronize()
    se = [s for (d, s) in et[2:] if abs(gmsh.model.getBoundingBox(2, s)[0] - L / 2) < 1e-9 and abs(gmsh.model.getBoundingBox(2, s)[3] - L / 2) < 1e-9][0]
    geo.extrude([(2, eq[5][1]), (2, se)], -L / 2, 0, 0, [nx2], recombine=True)
    vols = None
else:
    p0 = [geo.addPoint(0, 0, 0, h), geo.addPoint(L, 0, 0, h), geo.addPoint(L, W, 0, h), geo.addPoint(0, W, 0, h)]
    s0 = geo.addPlaneSurface([geo.addCurveLoop([geo.addLine(p0[i], p0[(i + 1) % 4]) for i in range(4)])])
    ex = geo.extrude([(2, s0)], 0, 0, T / 2, [2], recombine=True)
    vols.append(ex[1][1])
    et = geo.extrude([ex[0]], 0, 0, T / 2)
    vols.append(et[1][1])
    if air:
        ea = geo.extrude([et[0]], 0, 0, G)
        vols.append(ea[1][1])
geo.synchronize()
if vols is None:
    vols = [t for (d, t) in gmsh.model.getEntities(3)]
    if air: raise SystemExit("air is not available with tri")
e = 1e-9
def box(x0, y0, z0, x1, y1, z1):
    return [t for (d, t) in gmsh.model.getEntitiesInBoundingBox(x0 - e, y0 - e, z0 - e, x1 + e, y1 + e, z1 + e, 2)]
gmsh.model.addPhysicalGroup(3, vols[:2] if air else vols, 1)
if air:
    gmsh.model.addPhysicalGroup(3, vols[2:], 2)
gmsh.model.addPhysicalGroup(2, box(0, 0, 0, 0, W, T), 10)
gmsh.model.addPhysicalGroup(2, box(L, 0, 0, L, W, T), 11)
if yawall:
    gmsh.model.addPhysicalGroup(2, sorted(set(box(0, 0, 0, L, W, 0) + box(0, 0, T + G, L, W, T + G) + box(0, 0, 0, L, 0, T + G) + box(0, W, 0, L, W, T + G)
                                              + (box(0, 0, T, 0, W, T + G) + box(L, 0, T, L, W, T + G) if air else []))), 20)
else:
    gmsh.model.addPhysicalGroup(2, sorted(set(box(0, 0, 0, L, W, 0) + box(0, 0, T + G, L, W, T + G) + box(0, 0, 0, 0, W, T + G) + box(L, 0, 0, L, W, T + G))), 20)
gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 0)
gmsh.model.mesh.generate(3)
if twist:
    rng = random.Random(12345)
    rot = {5: [lambda u, v, z: (u, v, z), lambda u, v, z: (u, -z, v), lambda u, v, z: (z, v, -u), lambda u, v, z: (-v, u, z)],
           6: [lambda u, v, z: (u, v, z), lambda u, v, z: (v, 1 - u - v, z), lambda u, v, z: (v, u, -z), lambda u, v, z: (1 - u - v, u, z)],
           4: [lambda u, v, z: (u, v, z), lambda u, v, z: (v, z, u), lambda u, v, z: (z, u, v)]}
    for (d, b) in gmsh.model.getEntities(3):
        types, etags, ntags = gmsh.model.mesh.getElements(3, b)
        for (t, et_, nt) in zip(types, etags, ntags):
            name, dim, eo, nn, lc, npn = gmsh.model.mesh.getElementProperties(t)
            rots = rot[t]
            pos = {(round(lc[3 * a] * 1e6), round(lc[3 * a + 1] * 1e6), round(lc[3 * a + 2] * 1e6)): a for a in range(nn)}
            new = []
            for k in range(len(et_)):
                r = rots[rng.randrange(len(rots))]
                old = nt[k * nn:(k + 1) * nn]
                for a in range(nn):
                    c = r(lc[3 * a], lc[3 * a + 1], lc[3 * a + 2])
                    new.append(old[pos[(round(c[0] * 1e6), round(c[1] * 1e6), round(c[2] * 1e6))]])
            gmsh.model.mesh.removeElements(3, b, list(et_))
            gmsh.model.mesh.addElements(3, b, [t], [list(et_)], [new])
gmsh.model.mesh.setOrder(order)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.option.setNumber("Mesh.SaveAll", 0)
gmsh.write(dst)
gmsh.finalize()
