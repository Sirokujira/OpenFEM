#!/usr/bin/env python3
"""mkmesh.py — 検証用の四面体メッシュ (Gmsh ASCII 2.2) を生成する

外部のメッシュ生成器に依存せずに検証できるよう、解析解のある 2 つの形状を
自前で切る。どちらも 6 面体/角柱を四面体に分割した構造的なメッシュだが、
OpenFEM 側から見れば一般の非構造格子 (節点の並びも隣接関係も任意) になる。

  box  : 直方体を六面体 -> 6 四面体に分割 (平行平板コンデンサ)
         物理タグ 1 = 体積、10 = z 下面 (電極 0)、11 = z 上面 (電極 1)

  coax : 円環を極座標で切って三角形柱にし、1 層押し出して 3 四面体に分割
         (同軸線路)。円形境界に**適合**するので階段近似の誤差が出ない。
         物理タグ 1 = 体積、10 = 外側 r=b (電極 0)、11 = 内側 r=a (電極 1)

  bar  : 導体棒 (3 次元渦電流 A-φ)。表皮効果の 1 次元厳密解と比較する。
         物理タグ 1 = 体積、10 = x=0 面 (電極 0)、11 = x=lx 面 (電極 1)、
         20 = A_t = 0 の面 (z=0, z=lz, x=0, x=lx)

  bar_air : bar に非導電層 (空気) を載せたもの。解は 1 次元のままで
         閉形式が残るので、空気を含む系の検証に使える。
         物理タグ 1 = 導体、2 = 空気、10/11 = 導体断面の電極、20 = A_t = 0

使い方:
  python3 mkmesh.py box  box_tet.msh
  python3 mkmesh.py coax coax_tet.msh
  python3 mkmesh.py bar  bar_tet.msh
  python3 mkmesh.py bar_air bar_air.msh

  -order 2 を付けると 2 次要素 (tet10 / tri6) にする。coax では中間節点を
  円筒面に載せるので、境界が折れ線でなく円になる (等パラメトリック要素):
  python3 mkmesh.py coax coax_p2.msh -order 2 -nr 4 -nt 12

  分割数は -<キーワード引数> <値> で上書きできる (例 -nt 12 -nr 4)。
"""

import math
import sys


def write_msh(path, nodes, tets, tris):
    """要素の節点数から次数を決めて書く (4/3 なら 1 次、10/6 なら 2 次)"""
    ttype = {4: 4, 10: 11}
    stype = {3: 2, 6: 9}
    with open(path, "w") as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        f.write("$Nodes\n%d\n" % len(nodes))
        for i, (x, y, z) in enumerate(nodes):
            f.write("%d %.16g %.16g %.16g\n" % (i + 1, x, y, z))
        f.write("$EndNodes\n")
        f.write("$Elements\n%d\n" % (len(tets) + len(tris)))
        e = 0
        for lst, tmap in ((tris, stype), (tets, ttype)):
            for tag, n in lst:
                e += 1
                f.write("%d %d 2 %d %d %s\n" % (e, tmap[len(n)], tag, tag,
                                                " ".join(str(v + 1) for v in n)))
        f.write("$EndElements\n")


# 2 次要素の辺の並び (Gmsh の tet10 / tri6)
TET_EDGE = ((0, 1), (1, 2), (2, 0), (3, 0), (3, 2), (3, 1))
TRI_EDGE = ((0, 1), (1, 2), (2, 0))


def to_order2(nodes, tets, tris, snap=None):
    """1 次の格子を 2 次 (tet10 / tri6) に上げる

    共有される辺には同じ中間節点を割り当てる (辺を節点番号の組で識別する)。
    既定では辺の中点に置くので要素は直線のままになる。snap(pa, pb, mid) を
    渡すとその戻り値を中間節点の座標に使えるので、境界を曲面に載せられる
    (等パラメトリック要素)。
    """
    nodes = list(nodes)
    mid = {}

    def midnode(ia, ib):
        key = (ia, ib) if ia < ib else (ib, ia)
        if key not in mid:
            pa, pb = nodes[ia], nodes[ib]
            p = tuple((pa[d] + pb[d]) / 2 for d in range(3))
            if snap is not None:
                p = snap(pa, pb, p)
            mid[key] = len(nodes)
            nodes.append(p)
        return mid[key]

    tets2 = [(tag, list(n) + [midnode(n[a], n[b]) for a, b in TET_EDGE])
             for tag, n in tets]
    tris2 = [(tag, list(n) + [midnode(n[a], n[b]) for a, b in TRI_EDGE])
             for tag, n in tris]
    return nodes, tets2, tris2


def snap_cylinder(pa, pb, p):
    """円柱面に載せる : 両端の半径が等しい辺だけ、中点をその半径まで押し出す

    半径方向の辺 (両端の半径が違う) と軸方向の辺は中点のままにする。
    """
    ra = math.hypot(pa[0], pa[1])
    rb = math.hypot(pb[0], pb[1])
    if (ra <= 0) or (abs(ra - rb) > 1e-12 * ra):
        return p
    rm = math.hypot(p[0], p[1])
    if rm <= 0:
        return p
    return (p[0] * ra / rm, p[1] * ra / rm, p[2])


# 六面体 (節点 8 個) を 6 四面体に分割する定型分割
HEX2TET = ((0, 1, 3, 7), (0, 1, 7, 5), (0, 5, 7, 4),
           (0, 3, 2, 7), (0, 6, 4, 7), (0, 2, 6, 7))


def make_box(nx=8, ny=8, nz=6, lx=1e-3, ly=1e-3, lz=0.2e-3):
    """平行平板コンデンサ : C = eps0 epsr A / d"""
    nodes = []
    idx = {}
    for i in range(nx + 1):
        for j in range(ny + 1):
            for k in range(nz + 1):
                idx[(i, j, k)] = len(nodes)
                nodes.append((lx * i / nx, ly * j / ny, lz * k / nz))

    tets, tris = [], []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                # 局所節点 : ビット (i, j, k)
                h = [idx[(i + ((b >> 2) & 1), j + ((b >> 1) & 1), k + (b & 1))]
                     for b in range(8)]
                for t in HEX2TET:
                    tets.append((1, [h[t[0]], h[t[1]], h[t[2]], h[t[3]]]))

    # 電極面 (z = 0 と z = lz) を三角形で覆う
    for i in range(nx):
        for j in range(ny):
            for k, tag in ((0, 10), (nz, 11)):
                a = idx[(i, j, k)]
                b = idx[(i + 1, j, k)]
                c = idx[(i + 1, j + 1, k)]
                d = idx[(i, j + 1, k)]
                tris.append((tag, [a, b, c]))
                tris.append((tag, [a, c, d]))
    return nodes, tets, tris


def make_coax(nr=16, nt=48, ra=0.5e-3, rb=1.5e-3, lz=0.1e-3):
    """同軸線路の円環 : C' = 2 pi eps / ln(b/a)、円形境界に適合する"""
    nodes = []
    idx = {}
    for i in range(nr + 1):
        # 半径方向は対数等分にすると場の変化に合う
        r = ra * ((rb / ra) ** (i / nr))
        for j in range(nt):
            th = 2 * math.pi * j / nt
            for k in range(2):
                idx[(i, j, k)] = len(nodes)
                nodes.append((r * math.cos(th), r * math.sin(th), lz * k))

    tets, tris = [], []
    for i in range(nr):
        for j in range(nt):
            j2 = (j + 1) % nt
            # 四角形柱を 2 つの三角形柱に割り、各三角形柱を 3 四面体にする
            quad = ((i, j), (i + 1, j), (i + 1, j2), (i, j2))
            for tri in (((i, j), (i + 1, j), (i + 1, j2)),
                        ((i, j), (i + 1, j2), (i, j2))):
                b0 = [idx[(p[0], p[1], 0)] for p in tri]
                b1 = [idx[(p[0], p[1], 1)] for p in tri]
                # 三角形柱 (b0, b1) -> 3 四面体 (節点番号順で整合する定型分割)
                tets.append((1, [b0[0], b0[1], b0[2], b1[2]]))
                tets.append((1, [b0[0], b0[1], b1[2], b1[1]]))
                tets.append((1, [b0[0], b1[1], b1[2], b1[0]]))
            del quad

    # 電極面 : 内側 r=a (タグ 11)、外側 r=b (タグ 10)
    for j in range(nt):
        j2 = (j + 1) % nt
        for i, tag in ((0, 11), (nr, 10)):
            a = idx[(i, j, 0)]
            b = idx[(i, j2, 0)]
            c = idx[(i, j2, 1)]
            d = idx[(i, j, 1)]
            tris.append((tag, [a, b, c]))
            tris.append((tag, [a, c, d]))
    return nodes, tets, tris


def make_bar(nx=24, ny=2, nz=24, lx=2e-3, ly=0.25e-3, lz=1e-3, grade=1.0):
    """導体棒 (3 次元渦電流 A-φ の検証)

    x 方向に電流を流す。1 次元の厳密解が成り立つよう、A の接線成分を 0 にする
    面を z=0, z=lz, x=0, x=lx にとる。y 面は自然境界条件で、1 次元解が厳密に満たす。

    **要素は等方的にする。** 1 次 Nedelec 要素の誤差は「場が変化する方向の刻み」
    ではなく要素の最大寸法で決まるので、場が x, y に一様でも dx, dy を粗くすると
    誤差が出る (実測 : dz を 1/40 まで細かくしても dx = 1.7mm では R が 13% ずれ、
    dx を 0.42mm にすると 1%、0.083mm で 0.06% になる)。
    grade で z を両表面に向けて等比に細かくできるが、その分 dx も詰める必要が
    あるので既定は等間隔 (grade = 1.0)。

    物理タグ : 1 = 体積、10 = x=0 面 (電極 0)、11 = x=lx 面 (電極 1)、
               20 = A_t = 0 の面 (z=0, z=lz, x=0, x=lx)
    電極面は 10/11 と 20 の 2 つのタグで二重に出力する (Gmsh でも同じ扱い)。
    """
    if nz % 2:
        raise ValueError("nz must be even (symmetric grading)")

    # z の分割 : 両側から等比、中央で対称
    half = nz // 2
    w = [grade ** i for i in range(half)]
    s = sum(w)
    zs = [0.0]
    for i in range(half):
        zs.append(zs[-1] + (lz / 2) * w[i] / s)
    for i in range(half - 1, -1, -1):
        zs.append(zs[-1] + (lz / 2) * w[i] / s)
    zs[-1] = lz

    nodes = []
    idx = {}
    for i in range(nx + 1):
        for j in range(ny + 1):
            for k in range(nz + 1):
                idx[(i, j, k)] = len(nodes)
                nodes.append((lx * i / nx, ly * j / ny, zs[k]))

    tets, tris = [], []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                h = [idx[(i + ((b >> 2) & 1), j + ((b >> 1) & 1), k + (b & 1))]
                     for b in range(8)]
                for t in HEX2TET:
                    tets.append((1, [h[t[0]], h[t[1]], h[t[2]], h[t[3]]]))

    def quad(a, b, c, d, tag):
        tris.append((tag, [a, b, c]))
        tris.append((tag, [a, c, d]))

    # z = 0 と z = lz : A_t = 0
    for i in range(nx):
        for j in range(ny):
            for k in (0, nz):
                quad(idx[(i, j, k)], idx[(i + 1, j, k)],
                     idx[(i + 1, j + 1, k)], idx[(i, j + 1, k)], 20)
    # x = 0 (電極 0) と x = lx (電極 1) : 電極かつ A_t = 0
    for j in range(ny):
        for k in range(nz):
            for i, tag in ((0, 10), (nx, 11)):
                a = idx[(i, j, k)]
                b = idx[(i, j + 1, k)]
                c = idx[(i, j + 1, k + 1)]
                d = idx[(i, j, k + 1)]
                quad(a, b, c, d, tag)
                quad(a, b, c, d, 20)
    return nodes, tets, tris


def make_bar_air(nx=24, ny=2, nz=24, nza=4,
                 lx=2e-3, ly=0.25e-3, lz=1e-3, gz=0.5e-3):
    """導体棒 + 非導電層 (3 次元渦電流 A-φ で空気を含む系の検証)

    make_bar の上に厚さ gz の空気層を載せる。解は依然 1 次元で、空気層は
    界面で Robin 条件 A(t) + (g/mur) A'(t) = 0 に潰れるため閉形式が残る:

        Z = γ ℓ / (σ W [sinh(γt) − X (cosh(γt) − 1)])
        X = (cosh(γt) − 1 + (g/mur) γ sinh(γt)) / (sinh(γt) + (g/mur) γ cosh(γt))

    g → 0 で make_bar の 2 tanh(γt/2) に厳密に戻る。

    **A_t = 0 の面は空気側も覆う。** 導体側だけタグを付けると、空気部分の
    x 面が自然境界条件になって接線 H が 0 に強制され、1 次元解が崩れる。
    そのため x 面はタグ 20 を全高さに、電極タグ 10/11 は導体部分だけに付ける。

    物理タグ : 1 = 導体、2 = 空気、10 = x=0 の導体断面 (電極 0)、
               11 = x=lx の導体断面 (電極 1)、20 = A_t = 0 の面
    """
    zs = [lz * k / nz for k in range(nz + 1)] \
       + [lz + gz * k / nza for k in range(1, nza + 1)]
    nzt = nz + nza

    nodes = []
    idx = {}
    for i in range(nx + 1):
        for j in range(ny + 1):
            for k in range(nzt + 1):
                idx[(i, j, k)] = len(nodes)
                nodes.append((lx * i / nx, ly * j / ny, zs[k]))

    tets, tris = [], []
    for i in range(nx):
        for j in range(ny):
            for k in range(nzt):
                tag = 1 if k < nz else 2
                h = [idx[(i + ((b >> 2) & 1), j + ((b >> 1) & 1), k + (b & 1))]
                     for b in range(8)]
                for t in HEX2TET:
                    tets.append((tag, [h[t[0]], h[t[1]], h[t[2]], h[t[3]]]))

    def quad(a, b, c, d, tag):
        tris.append((tag, [a, b, c]))
        tris.append((tag, [a, c, d]))

    # z = 0 と z = lz + gz : A_t = 0
    for i in range(nx):
        for j in range(ny):
            for k in (0, nzt):
                quad(idx[(i, j, k)], idx[(i + 1, j, k)],
                     idx[(i + 1, j + 1, k)], idx[(i, j + 1, k)], 20)
    # x = 0 / x = lx : 全高さに A_t = 0、電極は導体部分だけ
    for j in range(ny):
        for k in range(nzt):
            for i, tag in ((0, 10), (nx, 11)):
                a = idx[(i, j, k)]
                b = idx[(i, j + 1, k)]
                c = idx[(i, j + 1, k + 1)]
                d = idx[(i, j, k + 1)]
                quad(a, b, c, d, 20)
                if k < nz:
                    quad(a, b, c, d, tag)
    return nodes, tets, tris


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 1
    kind, path = sys.argv[1], sys.argv[2]
    # 追加引数 : -order 2 で 2 次要素、-nt / -nr で同軸の分割数
    opt = {}
    a = 3
    while a + 1 < len(sys.argv):
        opt[sys.argv[a].lstrip("-")] = int(sys.argv[a + 1])
        a += 2
    order = opt.pop("order", 1)

    snap = None
    if kind == "box":
        nodes, tets, tris = make_box(**opt)
    elif kind == "coax":
        nodes, tets, tris = make_coax(**opt)
        snap = snap_cylinder			# 2 次にするとき円筒面に載せる
    elif kind == "bar":
        nodes, tets, tris = make_bar(**opt)
    elif kind == "bar_air":
        nodes, tets, tris = make_bar_air(**opt)
    else:
        print("unknown mesh kind: %s" % kind)
        return 1
    if order == 2:
        nodes, tets, tris = to_order2(nodes, tets, tris, snap)
    write_msh(path, nodes, tets, tris)
    print("%s : order %d, %d nodes, %d tets, %d tris"
          % (path, order, len(nodes), len(tets), len(tris)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
