# gmsh の高次の四面体 / 三角形の局所節点の並びを「再帰の規則」で作り、gmsh (4.15.2、Python API) の出力を重心座標で
# 分類したものと全項で比べる (unstruct.c の tet_rule / tri_rule と同じ規則。2〜7 次で一致することを確かめてある)。
#   python gmsh_ho_order.py 2 3 4 5 6 7
# 規則 : 頂点、辺 (0,1) (1,2) (2,0) (3,0) (3,2) (3,1) に p−1 個ずつ (始点の側から)、面 (0,2,1) (0,1,3) (0,3,2) (3,1,2) に
# (1,1,1) + 次数 p−3 の三角形 (同じ規則、面の頂点の順に写す)、内部に (1,1,1,1) + 次数 p−4 の四面体。
# 三角形は頂点、辺 (0,1) (1,2) (2,0)、内部 (1,1,1) + 次数 p−3
import gmsh, sys
def tri_order(p):
    if p == 0: return [(0, 0, 0)]
    out = [(p, 0, 0), (0, p, 0), (0, 0, p)]
    for (i, j) in ((0, 1), (1, 2), (2, 0)):
        for k in range(1, p):
            a = [0, 0, 0]; a[i] = p - k; a[j] = k; out.append(tuple(a))
    if p >= 3:
        for s in tri_order(p - 3): out.append(tuple(x + 1 for x in s))
    return out
def tet_order(p):
    if p == 0: return [(0, 0, 0, 0)]
    out = [tuple(p if v == w else 0 for w in range(4)) for v in range(4)]
    for (i, j) in ((0, 1), (1, 2), (2, 0), (3, 0), (3, 2), (3, 1)):
        for k in range(1, p):
            a = [0] * 4; a[i] = p - k; a[j] = k; out.append(tuple(a))
    if p >= 3:
        for f in ((0, 2, 1), (0, 1, 3), (0, 3, 2), (3, 1, 2)):
            for s in tri_order(p - 3):
                a = [0] * 4
                for m in range(3): a[f[m]] = s[m] + 1
                out.append(tuple(a))
    if p >= 4:
        for s in tet_order(p - 4): out.append(tuple(x + 1 for x in s))
    return out
def solve3(A, b):
    # 3 x 3 の最小二乗 (列が 2 本なら正規方程式) : 重心座標
    import itertools
    n = len(A[0])
    M = [[sum(A[r][i] * A[r][j] for r in range(3)) for j in range(n)] for i in range(n)]
    v = [sum(A[r][i] * b[r] for r in range(3)) for i in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            f = M[j][i] / M[i][i]
            for k in range(n): M[j][k] -= f * M[i][k]
            v[j] -= f * v[i]
    x = [0] * n
    for i in reversed(range(n)):
        x[i] = (v[i] - sum(M[i][k] * x[k] for k in range(i + 1, n))) / M[i][i]
    return x
def classify(etype, dim, p):
    tags, nodes = gmsh.model.mesh.getElementsByType(etype)
    nn = len(nodes) // len(tags)
    first = None; bad = 0
    for e in range(len(tags)):
        X = [gmsh.model.mesh.getNode(int(n))[0] for n in nodes[e * nn:(e + 1) * nn]]
        V = X[:dim + 1]
        A = [[V[k][c] - V[0][c] for k in range(1, dim + 1)] for c in range(3)]
        res = []
        for x in X:
            l = solve3(A, [x[c] - V[0][c] for c in range(3)])
            res.append(tuple(int(round(p * v)) for v in [1 - sum(l)] + list(l)))
        if first is None: first = res
        elif res != first: bad += 1
    return first, bad, len(tags)
gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
ok = True
for p in [int(a) for a in sys.argv[1:]]:
    gmsh.clear(); gmsh.model.add("b"); gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1); gmsh.model.occ.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.6); gmsh.model.mesh.generate(3); gmsh.model.mesh.setOrder(p)
    for dim, name, gen in ((3, "Tetrahedron", tet_order), (2, "Triangle", tri_order)):
        et = gmsh.model.mesh.getElementType(name, p)
        got, bad, ne = classify(et, dim, p)
        m = (got == gen(p)) and (bad == 0)
        ok = ok and m
        print("order %d %s : gmsh type %d, %d nodes, %d elements, rule %s" % (p, name, et, len(got), ne, "MATCHES" if m else "DIFFERS"))
gmsh.finalize()
sys.exit(0 if ok else 1)
