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

使い方:
  python3 mkmesh.py box  box_tet.msh
  python3 mkmesh.py coax coax_tet.msh
"""

import math
import sys


def write_msh(path, nodes, tets, tris):
    with open(path, "w") as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        f.write("$Nodes\n%d\n" % len(nodes))
        for i, (x, y, z) in enumerate(nodes):
            f.write("%d %.10g %.10g %.10g\n" % (i + 1, x, y, z))
        f.write("$EndNodes\n")
        f.write("$Elements\n%d\n" % (len(tets) + len(tris)))
        e = 0
        for tag, n in tris:
            e += 1
            f.write("%d 2 2 %d %d %d %d %d\n" % (e, tag, tag,
                                                 n[0] + 1, n[1] + 1, n[2] + 1))
        for tag, n in tets:
            e += 1
            f.write("%d 4 2 %d %d %d %d %d %d\n" % (e, tag, tag,
                                                    n[0] + 1, n[1] + 1, n[2] + 1, n[3] + 1))
        f.write("$EndElements\n")


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


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 1
    kind, path = sys.argv[1], sys.argv[2]
    if kind == "box":
        nodes, tets, tris = make_box()
    elif kind == "coax":
        nodes, tets, tris = make_coax()
    else:
        print("unknown mesh kind: %s" % kind)
        return 1
    write_msh(path, nodes, tets, tris)
    print("%s : %d nodes, %d tets, %d tris" % (path, len(nodes), len(tets), len(tris)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
