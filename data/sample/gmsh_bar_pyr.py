# bar_eddy と同じ導体棒 (L(x) = 2 mm, W(y) = 0.25 mm, t(z) = 1 mm) を、x < L/2 を六面体 (transfinite + recombine)、x > L/2 を
# 四面体で gmsh に切らせる。四面体の側の境界の四角形の面 (x = L/2 の界面) に gmsh が**ピラミッド**を自動で挟む
# (ピラミッドの辺要素 = edge.c の有理の Nedelec 基底の検証用)。物理タグは gmsh_bar_hex.py と同じ :
#   体積 1、x = 0 の面 10 (電極 0)、x = L の面 11 (電極 1)、z = 0 / z = t / x = 0 / x = L の面 20 (awall)
#   python gmsh_bar_pyr.py <出力> <h [m]>      (gmsh 4.15.2 の Python API、1 次)
import gmsh, sys
dst, h = sys.argv[1], float(sys.argv[2])
L, W, T = 2e-3, 0.25e-3, 1e-3
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.add("bar")
a = gmsh.model.occ.addBox(0, 0, 0, L / 2, W, T)
b = gmsh.model.occ.addBox(L / 2, 0, 0, L / 2, W, T)
gmsh.model.occ.fragment([(3, a)], [(3, b)])
gmsh.model.occ.synchronize()
e = 2e-7
def ents(dim, x0, x1, y0=-e, y1=W + e, z0=-e, z1=T + e):
    return [s for (d, s) in gmsh.model.getEntitiesInBoundingBox(x0, y0, z0, x1, y1, z1, dim)]
def faces(c, x):
    lo = [-e, -e, -e]; hi = [L + e, W + e, T + e]
    lo[c] = x - e; hi[c] = x + e
    return [s for (d, s) in gmsh.model.getEntitiesInBoundingBox(lo[0], lo[1], lo[2], hi[0], hi[1], hi[2], 2)]
hexv = ents(3, -e, L / 2 + e)
# 六面体の側 : 曲線を h で区切り、面・体積を transfinite + recombine
for c in ents(1, -e, L / 2 + e):
    x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(1, c)
    n = max(1, int(round(max(x1 - x0, y1 - y0, z1 - z0) / h)))
    gmsh.model.mesh.setTransfiniteCurve(c, n + 1)
for f in ents(2, -e, L / 2 + e):
    gmsh.model.mesh.setTransfiniteSurface(f)
    gmsh.model.mesh.setRecombine(2, f)
for v in hexv:
    gmsh.model.mesh.setTransfiniteVolume(v)
    gmsh.model.mesh.setRecombine(3, v)
gmsh.option.setNumber("Mesh.MeshSizeMax", h)
gmsh.option.setNumber("Mesh.MeshSizeMin", h)
gmsh.model.addPhysicalGroup(3, [s for (d, s) in gmsh.model.getEntities(3)], 1)
gmsh.model.addPhysicalGroup(2, faces(0, 0), 10)
gmsh.model.addPhysicalGroup(2, faces(0, L), 11)
gmsh.model.addPhysicalGroup(2, faces(2, 0) + faces(2, T) + faces(0, 0) + faces(0, L), 20)
gmsh.model.mesh.generate(3)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.option.setNumber("Mesh.SaveAll", 0)
gmsh.write(dst)
gmsh.finalize()
