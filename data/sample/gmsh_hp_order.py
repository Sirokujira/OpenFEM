# 高次の六面体・角柱の局所節点の並びを gmsh (4.15.2、Python API) の出力から分類し、unstruct.c の実測の表 (HEX3_I / HEX4_I /
# HEX5_I / PRISM3_A / PRISM4_A / PRISM5_A) と全項で比べる。六面体は角の局所座標を p 倍して丸めた格子位置 (0..p)、
# 角柱は三角形の重心座標の多重指数 (a, b, c) と ζ の格子位置。5 次の表はこの分類で作った (3 次・4 次の表を全項再現する
# ことを先に確かめた)。6 次以上の表を作る再帰の規則 (hex_rule / prism_rule) も 2〜8 次の全項で比べる。
#   python gmsh_hp_order.py   (unstruct.c は ../../sol/unstruct.c を読む)
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
# 6 次以上の表を作る再帰の規則 (unstruct.c の hex_rule / prism_rule / quad_rule と同じ形)。2〜8 次の gmsh の出力と全項で比べる
def line_pts(q):  # nodes 0..q in gmsh line order
    if q==0: return [0]
    return [0,q]+list(range(1,q))
def quad(q):  # list of (i,j) in 0..q, corners c0=(0,0),c1=(q,0),c2=(q,q),c3=(0,q) in gmsh quad order
    if q==0: return [(0,0)]
    pts=[(0,0),(q,0),(q,q),(0,q)]
    for t in range(1,q): pts.append((t,0))
    for t in range(1,q): pts.append((q,t))
    for t in range(1,q): pts.append((q-t,q))
    for t in range(1,q): pts.append((0,q-t))
    if q>=2: pts += [(i+1,j+1) for (i,j) in quad(q-2)]
    return pts
def tri(q):  # barycentric multi-index (a,b,c), a+b+c=q, vertices v0=(q,0,0),v1,v2
    if q<0: return []
    if q==0: return [(0,0,0)]
    pts=[(q,0,0),(0,q,0),(0,0,q)]
    for t in range(1,q): pts.append((q-t,t,0))
    for t in range(1,q): pts.append((0,q-t,t))
    for t in range(1,q): pts.append((t,0,q-t))
    pts += [(a+1,b+1,c+1) for (a,b,c) in tri(q-3)]
    return pts
HV=[(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)]
HE=[(0,1),(0,3),(0,4),(1,2),(1,5),(2,3),(2,6),(3,7),(4,5),(4,7),(5,6),(6,7)]
HF=[(0,3,2,1),(0,1,5,4),(0,4,7,3),(1,2,6,5),(2,3,7,6),(4,5,6,7)]
def hexr(p):
    if p==0: return [(0,0,0)]
    V=[tuple(p*c for c in v) for v in HV]
    out=list(V)
    for a,b in HE:
        for t in range(1,p): out.append(tuple(V[a][c]+(V[b][c]-V[a][c])*t//p for c in range(3)))
    for f in (HF if p>=2 else []):
        c0,c1,c3=V[f[0]],V[f[1]],V[f[3]]
        for (i,j) in quad(p-2):
            out.append(tuple(c0[c]+(c1[c]-c0[c])*(i+1)//p+(c3[c]-c0[c])*(j+1)//p for c in range(3)))
    if p>=2: out += [(x+1,y+1,z+1) for (x,y,z) in hexr(p-2)]
    return out
PV=[(1,0,0,0),(0,1,0,0),(0,0,1,0),(1,0,0,1),(0,1,0,1),(0,0,1,1)]
PE=[(0,1),(0,2),(0,3),(1,2),(1,4),(2,5),(3,4),(3,5),(4,5)]
PT=[(0,2,1),(3,4,5)]
PQ=[(0,1,4,3),(0,3,5,2),(1,2,5,4)]
def prismr(p):
    V=[tuple(p*c for c in v) for v in PV]
    out=list(V)
    for a,b in PE:
        for t in range(1,p): out.append(tuple(V[a][c]+(V[b][c]-V[a][c])*t//p for c in range(4)))
    for f in PT:
        for (a,b,c) in tri(p-3):
            w=(a+1,b+1,c+1)
            out.append(tuple(sum(w[k]*PV[f[k]][c2] for k in range(3)) if c2<3 else PV[f[0]][3]*p for c2 in range(4)))
    for f in PQ:
        c0,c1,c3=V[f[0]],V[f[1]],V[f[3]]
        for (i,j) in quad(p-2):
            out.append(tuple(c0[c]+(c1[c]-c0[c])*(i+1)//p+(c3[c]-c0[c])*(j+1)//p for c in range(4)))
    for (a,b,c) in tri(p-3):
        for k in line_pts(p-2):
            out.append((a+1,b+1,c+1,k+1))
    return out
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
for kind in ("hex", "prism"):
    for p in range(2, 9):
        mesh(kind, p); got, bad, ne, et = classify(kind, p)
        r = [tuple(x) for x in (hexr(p) if kind == "hex" else prismr(p))]
        same = ([tuple(x) for x in got] == r) and (bad == 0); ok = ok and same
        print(kind, p, "type", et, len(got), "nodes : recursive rule", "MATCHES" if same else "DIFFERS")
gmsh.finalize()
sys.exit(0 if ok else 1)
