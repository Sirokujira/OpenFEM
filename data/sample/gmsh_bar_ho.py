# bar_eddy と同じ導体棒 (L(x) = 2 mm, W(y) = 0.25 mm, t(z) = 1 mm) を gmsh の幾何から次数 p の四面体で切る
# (3 次以上の A、eddy3d_pr.c の検証用。gmsh は読み込んだ離散の格子を 3 次以上に上げられないので、
# gmsh_bar_p2.py と違い幾何から作る)。物理タグは mkmesh.py の bar と同じ :
#   体積 1、x = 0 の面 10 (電極 0)、x = L の面 11 (電極 1)、z = 0 / z = t / x = 0 / x = L の面 20 (awall)
#   python gmsh_bar_ho.py <出力> <次数> <要素の大きさ [m]>     (gmsh 4.15.2 の Python API)
#   python gmsh_bar_ho.py bar_ho2.msh 2 0.25e-3    (同じ頂点の 2 次 : 3 次との比較用)
#   python gmsh_bar_ho.py bar_ho3.msh 3 0.25e-3
# 4 番目の引数 air で bar_air.ofe と同じ空気層 (z = t .. t + 0.5 mm、体積 2) を載せる (電極は導体の断面だけ、
# awall は z = 0 / z = t + g と導体・空気の両方の x 面) :
#   python gmsh_bar_ho.py bar_air_ho3.msh 3 0.25e-3 air
#   python gmsh_bar_ho.py bar_air_ho3c.msh 3 0.5e-3 air   (粗い格子 : 3 次の B-H の検査、同じ頂点の 2 次は bar_air_ho2c.msh)
# 5 番目の引数 yawall でゲージ固定の検査用に電極の面を awall にせず、代わりに y 面 (導体・空気) を awall にする
# (mkmesh.py bar_air -yawall 1 と同じ。空気の x 面は awall のまま。電極の面は内部の頂点が無い構造格子) :
#   python gmsh_bar_ho.py bar_air_yw_ho3c.msh 3 0.5e-3 air yawall
import gmsh, sys
dst, order, h = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
air = (len(sys.argv) > 4) and (sys.argv[4] == "air")
yawall = air and (len(sys.argv) > 5) and (sys.argv[5] == "yawall")
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
    # 座標 c が x の平面に乗る面 (OCC の境界箱は 1e-7 だけ膨らむので、それより広い箱で拾う)。ztop / zbot で z の範囲
    e = 2e-7
    lo = [-e, -e, zbot - e]; hi = [L + e, W + e, (T + G if ztop is None else ztop) + e]
    lo[c] = x - e; hi[c] = x + e
    return [s for (d, s) in gmsh.model.getEntitiesInBoundingBox(lo[0], lo[1], lo[2], hi[0], hi[1], hi[2], 2)]
def vols(z0, z1):
    e = 2e-7
    return [s for (d, s) in gmsh.model.getEntitiesInBoundingBox(-e, -e, z0 - e, L + e, W + e, z1 + e, 3)]
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
if yawall:
    # 電極の面に内部の頂点を作らない (どの頂点も y 面の awall の辺に乗る。内部の頂点は awall の辺で
    # つながらないのでゲージ固定が弾く) : 電極の面を y に 1 区間・z に h の区間の構造格子にする
    for f in faces(0, 0, T) + faces(0, L, T):
        for (d, c) in gmsh.model.getBoundary([(2, f)], oriented=False):
            x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(1, c)
            gmsh.model.mesh.setTransfiniteCurve(c, 2 if (y1 - y0) > (z1 - z0) else int(round(T / h)) + 1)
        gmsh.model.mesh.setTransfiniteSurface(f)
gmsh.option.setNumber("Mesh.MeshSizeMax", h)
gmsh.option.setNumber("Mesh.MeshSizeMin", h)
gmsh.model.mesh.generate(3)
gmsh.model.mesh.setOrder(order)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.option.setNumber("Mesh.SaveAll", 0)
gmsh.write(dst)
gmsh.finalize()
