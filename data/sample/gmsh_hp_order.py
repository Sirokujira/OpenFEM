# 高次の六面体・角柱の局所節点の並びを gmsh (4.15.2、Python API) の出力から分類し、unstruct.c の実測の表 (HEX3_I / HEX4_I /
# HEX5_I / PRISM3_A / PRISM4_A / PRISM5_A) と全項で比べる。六面体は角の局所座標を p 倍して丸めた格子位置 (0..p)、
# 角柱は三角形の重心座標の多重指数 (a, b, c) と ζ の格子位置。5 次の表はこの分類で作った (3 次・4 次の表を全項再現する
# ことを先に確かめた)。   python gmsh_hp_order.py   (unstruct.c は ../../sol/unstruct.c を読む)
import gmsh, sys, re
src = open(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..', '..', 'sol', 'unstruct.c')).read()
def ctab(name):
    m = re.search(r'static const signed char %s\[\d+\]\[\d\] = \{(.*?)\n\};' % name, src, re.S)
    return [tuple(int(x) for x in t.split(',')) for t in re.findall(r'\{([-\d, ]+)\}', m.group(1))]
def mesh(kind, p):
    gmsh.clear(); gmsh.model.add("m")
    g = gmsh.model.geo
    pts = [g.addPoint(*c) for c in ((0,0,0),(1,0,0),(1,1,0),(0,1,0))]
    ls = [g.addLine(pts[i], pts[(i+1)%4]) for i in range(4)]
    s = g.addPlaneSurface([g.addCurveLoop(ls)])
    for l in ls: g.mesh.setTransfiniteCurve(l, 3)
    g.mesh.setTransfiniteSurface(s)
    if kind == "hex": g.mesh.setRecombine(2, s)
    g.extrude([(2, s)], 0, 0, 1, [2], recombine=True)
    g.synchronize(); gmsh.model.mesh.generate(3); gmsh.model.mesh.setOrder(p)
def classify(kind, p):
    et = gmsh.model.mesh.getElementType("Hexahedron" if kind == "hex" else "Prism", p)
    tags, nodes = gmsh.model.mesh.getElementsByType(et)
    nn = len(nodes) // len(tags); first = None; bad = 0
    for e in range(len(tags)):
        X = [gmsh.model.mesh.getNode(int(n))[0] for n in nodes[e*nn:(e+1)*nn]]
        res = []
        if kind == "hex":
            o = X[0]; ax = [[X[1][c]-o[c] for c in range(3)], [X[3][c]-o[c] for c in range(3)], [X[4][c]-o[c] for c in range(3)]]
            for x in X:
                d = [x[c]-o[c] for c in range(3)]
                res.append(tuple(int(round(p * sum(d[c]*a[c] for c in range(3)) / sum(a[c]*a[c] for c in range(3)))) for a in ax))
        else:
            o = X[0]; e1 = [X[1][c]-o[c] for c in range(3)]; e2 = [X[2][c]-o[c] for c in range(3)]; e3 = [X[3][c]-o[c] for c in range(3)]
            # 角柱は平行移動の三角形 : 面内 (u, v) は e1, e2 の 2x2 の正規方程式、ζ は e3
            def dot(a, b): return sum(a[c]*b[c] for c in range(3))
            for x in X:
                d = [x[c]-o[c] for c in range(3)]
                t = dot(d, e3) / dot(e3, e3); dd = [d[c] - t*e3[c] for c in range(3)]
                a11, a12, a22 = dot(e1,e1), dot(e1,e2), dot(e2,e2); b1, b2 = dot(dd,e1), dot(dd,e2)
                det = a11*a22 - a12*a12; u = (b1*a22 - b2*a12)/det; v = (a11*b2 - a12*b1)/det
                res.append((int(round(p*(1-u-v))), int(round(p*u)), int(round(p*v)), int(round(p*t))))
        if first is None: first = res
        elif res != first: bad += 1
    return first, bad, len(tags), et
gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
out = {}; ok = True
for kind in ("hex", "prism"):
    for p in (3, 4, 5):
        mesh(kind, p); got, bad, ne, et = classify(kind, p); out[(kind, p)] = got
        msg = ""
        if True:
            ref = ctab(("HEX%d_I" if kind == "hex" else "PRISM%d_A") % p)
            msg = "C table " + ("MATCHES" if ref == got else "DIFFERS"); ok = ok and (ref == got) and (bad == 0)
        print(kind, p, "type", et, len(got), "nodes", ne, "elements", bad, "disagree", msg)
gmsh.finalize()
sys.exit(0 if ok else 1)
