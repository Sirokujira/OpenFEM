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

  plate2d : 平行平板線路の**断面 2 次元**格子 (三角形が体積要素、M / F 用)。
         -swap 1 で面内 2 軸を入れ替え、-rot <度> で面内に任意角だけ回す。
         伝送線路軸は x。物理タグ 1 = 真空、2 = 導体の材料、
         10 = 導体 0 の断面、11 = 導体 1 の断面

  bar_air : bar に非導電層 (空気) を載せたもの。解は 1 次元のままで
         閉形式が残るので、空気を含む系の検証に使える。
         物理タグ 1 = 導体、2 = 空気、10/11 = 導体断面の電極、20 = A_t = 0

  box_pyr : 純ピラミッド格子 (各六面体セルを中心点で 6 分割)。-warp 1 で
         内部節点を乱数で動かす (底面が平行四辺形でないピラミッドを作る)。
         物理タグ 1 = 体積、10 = z 下面 (電極 0)、11 = z 上面 (電極 1)

  box_hexpyrtet : 六面体 -> ピラミッド -> 四面体の遷移格子 (下 nzh 層が六面体、
         その上の遷移層が底面ピラミッド + 四面体 10 個、残りが四面体)。
         物理タグ 1 = 六面体、2 = ピラミッド、3 = 四面体、
         10 = z 下面 (電極 0、四角形)、11 = z 上面 (電極 1、三角形)

使い方:
  python3 mkmesh.py box  box_tet.msh
  python3 mkmesh.py coax coax_tet.msh
  python3 mkmesh.py bar  bar_tet.msh
  python3 mkmesh.py bar_air bar_air.msh

  -order 2 を付けると 2 次要素 (tet10 / tri6) にする。coax では中間節点を
  円筒面に載せるので、境界が折れ線でなく円になる (等パラメトリック要素):
  python3 mkmesh.py coax coax_p2.msh -order 2 -nr 4 -nt 12

  分割数は -<キーワード引数> <値> で上書きできる (例 -nt 12 -nr 4)。

  -v41 1 を付けると Gmsh ASCII **4.1** 形式で書く (読み込みの検証用)。
  同じ形状を 2.2 と 4.1 で書いて結果が完全に一致することを rlc_check.sh で見る。

バイナリ形式の検証用ファイル (box_bin*.msh / plate2d_bin*.msh / box_p2_bin*.msh /
box_hexpyrtet_bin*.msh / box_hexpyrtet_41.msh) だけは、**このスクリプトでは
書きません**。自作の書き手と読み手が同じ誤解を
共有しているとテストが素通りするので、本物の gmsh (4.12.1) に変換させています:

  python3 mkmesh.py box     small.msh  -nx 3 -ny 3 -nz 2
  python3 mkmesh.py plate2d p2d.msh    -nw 6 -nt 2 -ng 2
  python3 mkmesh.py box     p2ord.msh  -nx 2 -ny 2 -nz 2 -order 2
  for src base in (small box_bin) (p2d plate2d_bin) (p2ord box_p2_bin):
      gmsh -0 $src.msh -o $base.msh     -format msh22
      gmsh -0 $src.msh -o ${base}_22.msh -format msh22 -bin
      gmsh -0 $src.msh -o ${base}_41.msh -format msh41 -bin

**gmsh は変換のたびに節点を振り直す**ので、比較する ASCII 側 (base.msh) も
同じ変換で出したものを置いてあります (元の出力と直接比べると丸めの順序が
変わって一致しません)。
"""

import math
import random
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


def write_msh41(path, nodes, tets, tris):
    """同じ格子を Gmsh ASCII 4.1 形式で書く (読み込みの検証用)

    4.1 では要素の行に物理タグが無く、**エンティティに付く**。そこで物理タグ
    ごとに 1 つのエンティティを作り、そのエンティティのブロックに要素を入れる。
    節点は 1 ブロックにまとめる (4.1 はブロック内で「番号がまとめて、そのあと
    座標がまとめて」の並びになる点が 2.2 と違う)。
    """
    ttype = {4: 4, 10: 11}
    stype = {3: 2, 6: 9}
    # 物理タグ -> (次元, エンティティ番号)。次元は三角形が 2、四面体が 3
    ent = {}
    for tag, n in tris:
        ent.setdefault((2, tag), len(ent) + 1)
    for tag, n in tets:
        ent.setdefault((3, tag), len(ent) + 1)

    with open(path, "w") as f:
        f.write("$MeshFormat\n4.1 0 8\n$EndMeshFormat\n")
        # $Entities : 点 0、曲線 0、面は三角形のタグ、体積は四面体のタグ
        surf = [(t, e) for (d, t), e in ent.items() if d == 2]
        vol = [(t, e) for (d, t), e in ent.items() if d == 3]
        # 点エンティティも書く。実際の Gmsh 出力には必ずあり、**点は座標 3 個、
        # それ以外は外接直方体 6 個**と項目数が違うので、読み手の分岐がここで
        # しか実行されない (点を省くとその分岐が死んだコードになる)
        f.write("$Entities\n2 0 %d %d\n" % (len(surf), len(vol)))
        f.write("1 0 0 0 0\n")
        f.write("2 0 0 0 0\n")
        for tag, e in surf:
            f.write("%d 0 0 0 0 0 0 1 %d 0\n" % (e, tag))
        for tag, e in vol:
            f.write("%d 0 0 0 0 0 0 1 %d 0\n" % (e, tag))
        f.write("$EndEntities\n")

        f.write("$Nodes\n1 %d 1 %d\n" % (len(nodes), len(nodes)))
        f.write("0 1 0 %d\n" % len(nodes))
        for i in range(len(nodes)):
            f.write("%d\n" % (i + 1))
        for (x, y, z) in nodes:
            f.write("%.16g %.16g %.16g\n" % (x, y, z))
        f.write("$EndNodes\n")

        # 要素は (次元, タグ, 型) ごとにブロックにまとめる
        blocks = {}
        for lst, dim, tmap in ((tris, 2, stype), (tets, 3, ttype)):
            for tag, n in lst:
                blocks.setdefault((dim, tag, tmap[len(n)]), []).append(n)
        ne = sum(len(v) for v in blocks.values())
        f.write("$Elements\n%d %d 1 %d\n" % (len(blocks), ne, ne))
        eid = 0
        for (dim, tag, ty), lst in blocks.items():
            f.write("%d %d %d %d\n" % (dim, ent[(dim, tag)], ty, len(lst)))
            for n in lst:
                eid += 1
                f.write("%d %s\n" % (eid, " ".join(str(v + 1) for v in n)))
        f.write("$EndElements\n")


def write_msh_hex(path, nodes, hexes, quads):
    """六面体格子を Gmsh ASCII 2.2 で書く (型 5 = hex8、型 3 = quad4)

    hex8 の局所節点の並びは Gmsh の規約 (**下面を反時計回り、その上に上面**)。
    VTK_HEXAHEDRON も同じ並びなので、書き出しでの並べ替えは要らない。
    向きが逆でもソルバー側は |det J| を使い符号の一定性だけを見るので解ける。
    """
    with open(path, "w") as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        f.write("$Nodes\n%d\n" % len(nodes))
        for i, (x, y, z) in enumerate(nodes):
            f.write("%d %.16g %.16g %.16g\n" % (i + 1, x, y, z))
        f.write("$EndNodes\n")
        f.write("$Elements\n%d\n" % (len(hexes) + len(quads)))
        e = 0
        for lst, ty in ((quads, 3), (hexes, 5)):
            for tag, n in lst:
                e += 1
                f.write("%d %d 2 %d %d %s\n" % (e, ty, tag, tag,
                                                 " ".join(str(v + 1) for v in n)))
        f.write("$EndElements\n")


def write_msh_cells(path, nodes, cells):
    """任意の要素型を Gmsh ASCII 2.2 で書く。cells = [(物理タグ, Gmsh 型, 節点)]

    型は 2 = 三角形、3 = 四角形、5 = 六面体、6 = 角柱、7 = ピラミッド。
    局所節点の並びは
    Gmsh の規約に合わせること (VTK も同じ並びなので書き出しでの入れ替えは不要)。
    """
    with open(path, "w") as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        f.write("$Nodes\n%d\n" % len(nodes))
        for i, (x, y, z) in enumerate(nodes):
            f.write("%d %.16g %.16g %.16g\n" % (i + 1, x, y, z))
        f.write("$EndNodes\n")
        f.write("$Elements\n%d\n" % len(cells))
        for e, (tag, ty, n) in enumerate(cells):
            f.write("%d %d 2 %d %d %s\n" % (e + 1, ty, tag, tag,
                                             " ".join(str(v + 1) for v in n)))
        f.write("$EndElements\n")


def make_box_mixed(nx=3, ny=3, nz=4, lx=1e-3, ly=1e-3, lz=0.2e-3, warp=0, seed=20260803):
    """**角柱 + 四面体**の混在格子 (境界層を角柱、内部を四面体にした形)

    下 2 層を角柱、上の層を四面体にする。境界面は三角形どうしで接するので
    **適合する** (どちらも面上で 1 次)。六面体と四面体は四角形面と三角形面で
    面上の解が食い違うため直接は隣り合わせにできない (ソルバー側で弾く)。

    warp = 1 で内部節点を乱数で動かす (等パラメトリック写像の検証用)。

    **角柱と四面体で物理タグを分ける** (1 / 2)。同じタグにすると要素番号から
    材料を引く経路 (種別ごとに配列が別で添字のずれ方も違う) が単一材料になり、
    添字を取り違えても答えが変わらない。

    物理タグ : 1 = 角柱の領域、2 = 四面体の領域、
               10 = z 下面 (電極 0)、11 = z 上面 (電極 1)
    """
    rnd = random.Random(seed)
    nodes = []
    idx = {}
    for i in range(nx + 1):
        for j in range(ny + 1):
            for k in range(nz + 1):
                x, y, z = lx * i / nx, ly * j / ny, lz * k / nz
                if warp:
                    if 0 < i < nx: x += 0.30 * (rnd.random() - 0.5) * lx / nx
                    if 0 < j < ny: y += 0.30 * (rnd.random() - 0.5) * ly / ny
                    if 0 < k < nz: z += 0.30 * (rnd.random() - 0.5) * lz / nz
                idx[(i, j, k)] = len(nodes)
                nodes.append((x, y, z))

    # 四角柱を対角線で 2 つの三角柱に割る (角柱と四面体で同じ三角形分割を使う)
    tris = (((0, 0), (1, 0), (1, 1)), ((0, 0), (1, 1), (0, 1)))
    nprism = 2                          # 下から 2 層を角柱に
    cells = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for tri in tris:
                    b = [idx[(i + a, j + c, k)] for a, c in tri]
                    t = [idx[(i + a, j + c, k + 1)] for a, c in tri]
                    if k < nprism:
                        cells.append((1, 6, b + t))
                    else:
                        # 三角柱を 3 四面体に割る (面が三角形で揃うので適合)
                        p = b + t
                        for q in ((0, 1, 2, 5), (0, 1, 5, 4), (0, 4, 5, 3)):
                            cells.append((2, 4, [p[q[0]], p[q[1]], p[q[2]], p[q[3]]]))
    for i in range(nx):
        for j in range(ny):
            for k, tag in ((0, 10), (nz, 11)):
                for tri in tris:
                    cells.append((tag, 2, [idx[(i + a, j + c, k)] for a, c in tri]))

    return nodes, cells


def make_box_prism_warp(nx=3, ny=3, nz=3, lx=1e-3, ly=1e-3, lz=0.2e-3, seed=20260803):
    """**ゆがんだ**角柱格子 (等パラメトリック写像の自己検証用)

    四角形を対角線で 2 つの三角形に割って z に押し出す。内部節点だけを乱数で
    動かすので、**上下の三角形が平行移動の関係でなくなる** — これがないと
    ヤコビアンが ζ に依らなくなり、「J を 1 回だけ評価する」誤りが検出できない。

    物理タグ : 1 = 体積、10 = z 下面 (電極 0)、11 = z 上面 (電極 1)
    """
    rnd = random.Random(seed)
    nodes = []
    idx = {}
    for i in range(nx + 1):
        for j in range(ny + 1):
            for k in range(nz + 1):
                x, y, z = lx * i / nx, ly * j / ny, lz * k / nz
                if 0 < i < nx: x += 0.30 * (rnd.random() - 0.5) * lx / nx
                if 0 < j < ny: y += 0.30 * (rnd.random() - 0.5) * ly / ny
                if 0 < k < nz: z += 0.30 * (rnd.random() - 0.5) * lz / nz
                idx[(i, j, k)] = len(nodes)
                nodes.append((x, y, z))

    cells = []
    # 下面の三角形の分割 : (i,j)-(i+1,j)-(i+1,j+1) と (i,j)-(i+1,j+1)-(i,j+1)
    tris = (((0, 0), (1, 0), (1, 1)), ((0, 0), (1, 1), (0, 1)))
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for tri in tris:
                    b = [idx[(i + a, j + c, k)] for a, c in tri]
                    t = [idx[(i + a, j + c, k + 1)] for a, c in tri]
                    cells.append((1, 6, b + t))
    for i in range(nx):
        for j in range(ny):
            for k, tag in ((0, 10), (nz, 11)):
                for tri in tris:
                    cells.append((tag, 2, [idx[(i + a, j + c, k)] for a, c in tri]))

    return nodes, cells


def make_box_pyr(nx=3, ny=3, nz=3, lx=1e-3, ly=1e-3, lz=0.2e-3, warp=0, seed=20260803):
    """純**ピラミッド**格子: 各六面体セルを中心点で 6 つのピラミッドに割る

    どのピラミッドも底面 (セルの面) は隣のセルのピラミッドと四角形面どうしで、
    側面 (三角形) は同じセルの隣のピラミッドと接するので**適合する**。

    warp = 1 で内部の格子節点とセル中心を乱数で動かす (等パラメトリック写像の
    検証用)。底面が平行四辺形のピラミッドでは剛性が厳密になるので、
    ゆがみが無いと有理式の積分近似が一度も実行されない。

    物理タグ : 1 = 体積、10 = z 下面 (電極 0)、11 = z 上面 (電極 1)
    """
    rnd = random.Random(seed)
    nodes = []
    idx = {}
    for i in range(nx + 1):
        for j in range(ny + 1):
            for k in range(nz + 1):
                x, y, z = lx * i / nx, ly * j / ny, lz * k / nz
                if warp:
                    if 0 < i < nx: x += 0.30 * (rnd.random() - 0.5) * lx / nx
                    if 0 < j < ny: y += 0.30 * (rnd.random() - 0.5) * ly / ny
                    if 0 < k < nz: z += 0.30 * (rnd.random() - 0.5) * lz / nz
                idx[(i, j, k)] = len(nodes)
                nodes.append((x, y, z))

    cells = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                # 頂点 (セル中心)。ゆがんだセルでも内部に入るよう角の平均にとり、
                # warp 時はさらに乱数で動かす
                cs = [nodes[idx[(i + a, j + b, k + c)]]
                      for a in (0, 1) for b in (0, 1) for c in (0, 1)]
                cx = [sum(p[d] for p in cs) / 8 for d in range(3)]
                if warp:
                    cx[0] += 0.15 * (rnd.random() - 0.5) * lx / nx
                    cx[1] += 0.15 * (rnd.random() - 0.5) * ly / ny
                    cx[2] += 0.15 * (rnd.random() - 0.5) * lz / nz
                c0 = len(nodes)
                nodes.append(tuple(cx))

                def v(a, b, c):
                    return idx[(i + a, j + b, k + c)]
                # 6 つの面 (巡回順) を底面に、中心を頂点にする
                faces = ((v(0, 0, 0), v(1, 0, 0), v(1, 1, 0), v(0, 1, 0)),   # z-
                         (v(0, 0, 1), v(1, 0, 1), v(1, 1, 1), v(0, 1, 1)),   # z+
                         (v(0, 0, 0), v(0, 1, 0), v(0, 1, 1), v(0, 0, 1)),   # x-
                         (v(1, 0, 0), v(1, 1, 0), v(1, 1, 1), v(1, 0, 1)),   # x+
                         (v(0, 0, 0), v(1, 0, 0), v(1, 0, 1), v(0, 0, 1)),   # y-
                         (v(0, 1, 0), v(1, 1, 0), v(1, 1, 1), v(0, 1, 1)))   # y+
                for f in faces:
                    cells.append((1, 7, list(f) + [c0]))
    for i in range(nx):
        for j in range(ny):
            for k, tag in ((0, 10), (nz, 11)):
                cells.append((tag, 3, [idx[(i, j, k)], idx[(i + 1, j, k)],
                                       idx[(i + 1, j + 1, k)], idx[(i, j + 1, k)]]))

    return nodes, cells


def make_box_hexpyrtet(nx=3, ny=3, nz=4, nzh=2, lx=1e-3, ly=1e-3, lz=0.2e-3,
                       warp=0, seed=20260803):
    """**六面体 -> ピラミッド -> 四面体**の遷移格子 (hex-to-tet transition)

    下 nzh 層は六面体。その上の 1 層が遷移層で、各セルを
      「底面のピラミッド (頂点 = セル中心) 1 個 + 側面・上面の三角形と
       セル中心を結ぶ四面体 10 個」
    に割る。六面体とは底面の四角形面どうしで、上の四面体層とは上面の
    三角形どうしで接するので**適合する**。残りの層は角柱分割由来の四面体。

    セル間で共有される面の対角線は格子座標の規則で固定する (両側のセルが
    同じ割り方をしないと適合しない)。

    **物理タグを種別ごとに分ける** (1 = 六面体、2 = ピラミッド、3 = 四面体)。
    要素番号 -> 材料の対応が種別の境目でずれると必ず落ちる (box_mixed と同じ)。
    ピラミッドと四面体は同じ高さ範囲に混ざるので、.ofe 側でタグ 2 と 3 に
    同じ材料を割り当てること (場を 1 次元に保つため)。

    物理タグ : 10 = z 下面 (電極 0、四角形)、11 = z 上面 (電極 1、三角形)
    """
    rnd = random.Random(seed)
    nodes = []
    idx = {}
    for i in range(nx + 1):
        for j in range(ny + 1):
            for k in range(nz + 1):
                x, y, z = lx * i / nx, ly * j / ny, lz * k / nz
                if warp:
                    if 0 < i < nx: x += 0.30 * (rnd.random() - 0.5) * lx / nx
                    if 0 < j < ny: y += 0.30 * (rnd.random() - 0.5) * ly / ny
                    if 0 < k < nz: z += 0.30 * (rnd.random() - 0.5) * lz / nz
                idx[(i, j, k)] = len(nodes)
                nodes.append((x, y, z))

    cells = []
    # 六面体の層 (k = 0 .. nzh-1)
    for i in range(nx):
        for j in range(ny):
            for k in range(nzh):
                b = [idx[(i, j, k)], idx[(i + 1, j, k)],
                     idx[(i + 1, j + 1, k)], idx[(i, j + 1, k)]]
                t = [idx[(i, j, k + 1)], idx[(i + 1, j, k + 1)],
                     idx[(i + 1, j + 1, k + 1)], idx[(i, j + 1, k + 1)]]
                cells.append((1, 5, b + t))

    # 遷移層 (k = nzh)。対角線の規則 : 面の格子座標 (s, t) で
    # (s,t)-(s+1,t+1) を結ぶ (隣のセルから見ても同じ線になる)
    tris = (((0, 0), (1, 0), (1, 1)), ((0, 0), (1, 1), (0, 1)))
    kt = nzh
    for i in range(nx):
        for j in range(ny):
            cs = [nodes[idx[(i + a, j + b, kt + c)]]
                  for a in (0, 1) for b in (0, 1) for c in (0, 1)]
            cx = [sum(p[d] for p in cs) / 8 for d in range(3)]
            if warp:
                cx[0] += 0.15 * (rnd.random() - 0.5) * lx / nx
                cx[1] += 0.15 * (rnd.random() - 0.5) * ly / ny
                cx[2] += 0.15 * (rnd.random() - 0.5) * lz / nz
            c0 = len(nodes)
            nodes.append(tuple(cx))

            # 底面のピラミッド (六面体の上面と四角形面どうしで接する)
            cells.append((2, 7, [idx[(i, j, kt)], idx[(i + 1, j, kt)],
                                 idx[(i + 1, j + 1, kt)], idx[(i, j + 1, kt)], c0]))
            # 側面 4 枚 (各 2 三角形) と上面 (2 三角形) + セル中心 = 四面体 10 個
            side = []
            for ii in (i, i + 1):			# x 一定の面 : (s, t) = (j, k)
                q = [(ii, j, kt), (ii, j + 1, kt), (ii, j + 1, kt + 1), (ii, j, kt + 1)]
                side += [(q[0], q[1], q[2]), (q[0], q[2], q[3])]
            for jj in (j, j + 1):			# y 一定の面 : (s, t) = (i, k)
                q = [(i, jj, kt), (i + 1, jj, kt), (i + 1, jj, kt + 1), (i, jj, kt + 1)]
                side += [(q[0], q[1], q[2]), (q[0], q[2], q[3])]
            for tri in tris:				# 上面 (上の四面体層と接する)
                side.append(tuple((i + a, j + b, kt + 1) for a, b in tri))
            for tri in side:
                cells.append((3, 4, [idx[t] for t in tri] + [c0]))

    # 四面体の層 (k = nzh+1 .. nz-1、角柱を 3 四面体に割る)
    for i in range(nx):
        for j in range(ny):
            for k in range(nzh + 1, nz):
                for tri in tris:
                    b = [idx[(i + a, j + c, k)] for a, c in tri]
                    t = [idx[(i + a, j + c, k + 1)] for a, c in tri]
                    p = b + t
                    for q in ((0, 1, 2, 5), (0, 1, 5, 4), (0, 4, 5, 3)):
                        cells.append((3, 4, [p[q[0]], p[q[1]], p[q[2]], p[q[3]]]))

    # 電極 : 下面は四角形、上面は三角形 (対角線は四面体側の分割と同じ)
    for i in range(nx):
        for j in range(ny):
            cells.append((10, 3, [idx[(i, j, 0)], idx[(i + 1, j, 0)],
                                  idx[(i + 1, j + 1, 0)], idx[(i, j + 1, 0)]]))
            for tri in tris:
                cells.append((11, 2, [idx[(i + a, j + b, nz)] for a, b in tri]))

    return nodes, cells


# Gmsh hex27 の局所節点 -> (ξ, η, ζ) ∈ {-1,0,1} (gmsh 4.12.1 の出力を分類して
# 実測した並び。C 側の HEX2_XI と同じ表)
HEX27_XI = ((-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1),
            (-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1),
            (0, -1, -1), (-1, 0, -1), (-1, -1, 0), (1, 0, -1), (1, -1, 0),
            (0, 1, -1), (1, 1, 0), (-1, 1, 0), (0, -1, 1), (-1, 0, 1),
            (1, 0, 1), (0, 1, 1),
            (0, 0, -1), (0, -1, 0), (-1, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1),
            (0, 0, 0))

# Gmsh prism18 の局所節点 -> (2 倍した三角形格子座標のオフセット, ζ レベル 0/1/2)。
# 三角形の頂点を A, B, C とすると (実測: 辺 (0,1)(0,2)(0,3)(1,2)(1,4)(2,5)
# (3,4)(3,5)(4,5)、四角形面 (0,1,4,3)(0,2,5,3)(1,2,5,4))
PRISM18_TZ = (("A", 0), ("B", 0), ("C", 0), ("A", 2), ("B", 2), ("C", 2),
              ("AB", 0), ("CA", 0), ("A", 1), ("BC", 0), ("B", 1), ("C", 1),
              ("AB", 2), ("CA", 2), ("BC", 2),
              ("AB", 1), ("CA", 1), ("BC", 1))


def make_coax_hex2(nr=4, nt=12, ra=0.5e-3, rb=1.5e-3, lz=0.1e-3):
    """同軸線路を**曲がった 2 次六面体** (hex27) で切る : C' = 2 pi eps / ln(b/a)

    (r, theta) を半分の刻みまで含めて全節点を**厳密に円筒面に載せる**ので、
    等パラメトリック写像が境界を円弧として表す (coax_p2 の六面体版)。
    半径方向は対数等分。同じ nr, nt の 1 次六面体 (coax_hex 相当) より
    誤差が大きく減ることを rlc_check.sh で対にして見る。

    物理タグ : 1 = 体積、10 = 外側 r=b (電極 0)、11 = 内側 r=a (電極 1)
    """
    nodes = []
    idx = {}

    def node(i, j, k):
        # i : 半径 (0..2nr)、j : 周方向 (mod 2nt)、k : z (0..2)
        j = j % (2 * nt)
        if (i, j, k) not in idx:
            r = ra * ((rb / ra) ** (i / (2.0 * nr)))
            th = math.pi * j / nt
            idx[(i, j, k)] = len(nodes)
            nodes.append((r * math.cos(th), r * math.sin(th), lz * k / 2))
        return idx[(i, j, k)]

    cells = []
    for i in range(nr):
        for j in range(nt):
            nd = [node((2 * i) + 1 + sx, (2 * j) + 1 + sy, 1 + sz)
                  for sx, sy, sz in HEX27_XI]
            cells.append((1, 12, nd))
    # 電極面 (内外の円筒面、9 節点四角形)。並びは角 4 + 辺 4 + 中心
    for j in range(nt):
        for i2, tag in ((0, 11), (2 * nr, 10)):
            q = [(2 * j, 0), ((2 * j) + 2, 0), ((2 * j) + 2, 2), (2 * j, 2),
                 ((2 * j) + 1, 0), ((2 * j) + 2, 1), ((2 * j) + 1, 2), (2 * j, 1),
                 ((2 * j) + 1, 1)]
            cells.append((tag, 10, [node(i2, a, k) for a, k in q]))

    return nodes, cells


def make_box_hex2_warp(nx=3, ny=3, nz=3, lx=1e-3, ly=1e-3, lz=0.2e-3, seed=20260803):
    """**ゆがんだ** 2 次六面体格子 (hex27、analysis = P の自己検証用)

    角の格子節点だけを乱数で動かし、中間節点 (辺・面・体心) は**ゆがんだ角の
    三重線形の位置** (対応する角の平均) に置く。要素は「まっすぐだが平行六面体
    でない」形になり、線形場の恒等式が有理式の積分近似ごしに検査される。
    中間節点が 1 次の写像の位置にあるので det J は 1 次と同じ多項式のままで、
    3 点則と 4 点則の体積は機械精度で一致しなければならない。

    物理タグ : 1 = 体積、10 = z 下面 (電極 0)、11 = z 上面 (電極 1)
    """
    rnd = random.Random(seed)
    corner = {}
    for i in range(nx + 1):
        for j in range(ny + 1):
            for k in range(nz + 1):
                x, y, z = lx * i / nx, ly * j / ny, lz * k / nz
                if 0 < i < nx: x += 0.30 * (rnd.random() - 0.5) * lx / nx
                if 0 < j < ny: y += 0.30 * (rnd.random() - 0.5) * ly / ny
                if 0 < k < nz: z += 0.30 * (rnd.random() - 0.5) * lz / nz
                corner[(2 * i, 2 * j, 2 * k)] = (x, y, z)

    nodes = []
    idx = {}

    def node(gi, gj, gk):
        # 偶数座標は角。奇数を含む座標は「関係する角の平均」(三重線形の位置)
        if (gi, gj, gk) not in idx:
            ii = ((gi - 1, gi + 1) if gi % 2 else (gi,))
            jj = ((gj - 1, gj + 1) if gj % 2 else (gj,))
            kk = ((gk - 1, gk + 1) if gk % 2 else (gk,))
            cs = [corner[(a, b, c)] for a in ii for b in jj for c in kk]
            idx[(gi, gj, gk)] = len(nodes)
            nodes.append(tuple(sum(p[d] for p in cs) / len(cs) for d in range(3)))
        return idx[(gi, gj, gk)]

    cells = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                nd = [node((2 * i) + 1 + sx, (2 * j) + 1 + sy, (2 * k) + 1 + sz)
                      for sx, sy, sz in HEX27_XI]
                cells.append((1, 12, nd))
    for i in range(nx):
        for j in range(ny):
            for gk, tag in ((0, 10), (2 * nz, 11)):
                gi, gj = 2 * i, 2 * j
                q = [(gi, gj), (gi + 2, gj), (gi + 2, gj + 2), (gi, gj + 2),
                     (gi + 1, gj), (gi + 2, gj + 1), (gi + 1, gj + 2), (gi, gj + 1),
                     (gi + 1, gj + 1)]
                cells.append((tag, 10, [node(a, b, gk) for a, b in q]))

    return nodes, cells


def make_box_prism2_warp(nx=3, ny=3, nz=3, lx=1e-3, ly=1e-3, lz=0.2e-3, seed=20260803):
    """**ゆがんだ** 2 次角柱格子 (prism18、analysis = P の自己検証用)

    四角形を対角線で 2 つの三角形に割って z に押し出す。角だけを乱数で動かし、
    中間節点は対応する角の平均 (1 次の写像の位置) に置く。

    物理タグ : 1 = 体積、10 = z 下面 (電極 0)、11 = z 上面 (電極 1)
    """
    rnd = random.Random(seed)
    corner = {}
    for i in range(nx + 1):
        for j in range(ny + 1):
            for k in range(nz + 1):
                x, y, z = lx * i / nx, ly * j / ny, lz * k / nz
                if 0 < i < nx: x += 0.30 * (rnd.random() - 0.5) * lx / nx
                if 0 < j < ny: y += 0.30 * (rnd.random() - 0.5) * ly / ny
                if 0 < k < nz: z += 0.30 * (rnd.random() - 0.5) * lz / nz
                corner[(2 * i, 2 * j, 2 * k)] = (x, y, z)

    nodes = []
    idx = {}

    def node(gi, gj, gk):
        if (gi, gj, gk) not in idx:
            # 奇数座標は関係する角の平均。面内は三角形の辺の中点 (角 2 個) か
            # 対角線の中点 (これも角 2 個)、z 方向は上下の平均
            pts = [(gi, gj)]
            if (gi % 2) and (gj % 2):
                # 対角線 (2i,2j)-(2i+2,2j+2) の中点
                pts = [(gi - 1, gj - 1), (gi + 1, gj + 1)]
            elif gi % 2:
                pts = [(gi - 1, gj), (gi + 1, gj)]
            elif gj % 2:
                pts = [(gi, gj - 1), (gi, gj + 1)]
            kk = ((gk - 1, gk + 1) if gk % 2 else (gk,))
            cs = [corner[(a, b, c)] for a, b in pts for c in kk]
            idx[(gi, gj, gk)] = len(nodes)
            nodes.append(tuple(sum(p[d] for p in cs) / len(cs) for d in range(3)))
        return idx[(gi, gj, gk)]

    tris = (((0, 0), (1, 0), (1, 1)), ((0, 0), (1, 1), (0, 1)))
    cells = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for tri in tris:
                    tp = {"A": (2 * (i + tri[0][0]), 2 * (j + tri[0][1])),
                          "B": (2 * (i + tri[1][0]), 2 * (j + tri[1][1])),
                          "C": (2 * (i + tri[2][0]), 2 * (j + tri[2][1]))}
                    for e in (("AB", "A", "B"), ("BC", "B", "C"), ("CA", "C", "A")):
                        tp[e[0]] = ((tp[e[1]][0] + tp[e[2]][0]) // 2,
                                    (tp[e[1]][1] + tp[e[2]][1]) // 2)
                    nd = [node(tp[t][0], tp[t][1], (2 * k) + z)
                          for t, z in PRISM18_TZ]
                    cells.append((1, 13, nd))
    for i in range(nx):
        for j in range(ny):
            for gk, tag in ((0, 10), (2 * nz, 11)):
                for tri in tris:
                    v = [(2 * (i + a), 2 * (j + b)) for a, b in tri]
                    m = [((v[0][0] + v[1][0]) // 2, (v[0][1] + v[1][1]) // 2),
                         ((v[1][0] + v[2][0]) // 2, (v[1][1] + v[2][1]) // 2),
                         ((v[2][0] + v[0][0]) // 2, (v[2][1] + v[0][1]) // 2)]
                    cells.append((tag, 9, [node(a, b, gk) for a, b in v + m]))

    return nodes, cells


def make_coax_hex(nr=8, nt=32, ra=0.5e-3, rb=1.5e-3, lz=0.1e-3):
    """同軸線路の円環を**六面体**で切る : C' = 2 pi eps / ln(b/a)

    (r, theta) の構造格子を z に 1 層押し出す。要素は台形 (r 方向に幅が変わる)
    なので**アフィンでない** — 要素内でヤコビアンが変化する。等パラメトリック
    写像を Gauss 点毎に評価していないと、この形でだけ誤差が出る
    (剛体回転した直方体はアフィンなので、そちらでは検出できない)。

    物理タグ : 1 = 体積、10 = 外側 r=b (電極 0)、11 = 内側 r=a (電極 1)
    """
    nodes = []
    idx = {}
    for i in range(nr + 1):
        # 半径方向は対数等分にすると場の変化に合う (四面体版と同じ)
        r = ra * ((rb / ra) ** (i / nr))
        for j in range(nt):
            th = 2 * math.pi * j / nt
            for k in range(2):
                idx[(i, j, k)] = len(nodes)
                nodes.append((r * math.cos(th), r * math.sin(th), lz * k))

    hexes, quads = [], []
    for i in range(nr):
        for j in range(nt):
            jn = (j + 1) % nt			# 周方向は閉じる
            # 下面 (k=0) を反時計回りに、その上に上面 (k=1)
            b0 = [idx[(i, j, 0)], idx[(i + 1, j, 0)], idx[(i + 1, jn, 0)], idx[(i, jn, 0)]]
            t0 = [idx[(i, j, 1)], idx[(i + 1, j, 1)], idx[(i + 1, jn, 1)], idx[(i, jn, 1)]]
            hexes.append((1, b0 + t0))
    # 電極面 (内外の円筒面)
    for j in range(nt):
        jn = (j + 1) % nt
        quads.append((11, [idx[(0, j, 0)], idx[(0, jn, 0)],
                           idx[(0, jn, 1)], idx[(0, j, 1)]]))
        quads.append((10, [idx[(nr, j, 0)], idx[(nr, jn, 0)],
                           idx[(nr, jn, 1)], idx[(nr, j, 1)]]))

    return nodes, hexes, quads


def make_box_hex_warp(nx=3, ny=3, nz=3, lx=1e-3, ly=1e-3, lz=0.2e-3, seed=20260803):
    """**ゆがんだ**六面体格子 (等パラメトリック写像の自己検証用)

    内部節点だけを乱数で動かす (境界面は平面のまま = 電極が平面を保つ)。
    平行六面体では要素内でヤコビアンが一定になり、「J を要素中心で 1 回だけ
    評価する」という手抜きが厳密に正しくなってしまうため、写像の検証には
    **アフィンでない要素**が要る。乱数の種は固定して結果を再現可能にする。

    物理タグ : 1 = 体積、10 = z 下面 (電極 0)、11 = z 上面 (電極 1)
    """
    rnd = random.Random(seed)
    nodes = []
    idx = {}
    for i in range(nx + 1):
        for j in range(ny + 1):
            for k in range(nz + 1):
                x, y, z = lx * i / nx, ly * j / ny, lz * k / nz
                if 0 < i < nx: x += 0.30 * (rnd.random() - 0.5) * lx / nx
                if 0 < j < ny: y += 0.30 * (rnd.random() - 0.5) * ly / ny
                if 0 < k < nz: z += 0.30 * (rnd.random() - 0.5) * lz / nz
                idx[(i, j, k)] = len(nodes)
                nodes.append((x, y, z))

    hexes, quads = [], []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                b = [idx[(i, j, k)], idx[(i + 1, j, k)],
                     idx[(i + 1, j + 1, k)], idx[(i, j + 1, k)]]
                t = [idx[(i, j, k + 1)], idx[(i + 1, j, k + 1)],
                     idx[(i + 1, j + 1, k + 1)], idx[(i, j + 1, k + 1)]]
                hexes.append((1, b + t))
    for i in range(nx):
        for j in range(ny):
            for k, tag in ((0, 10), (nz, 11)):
                quads.append((tag, [idx[(i, j, k)], idx[(i + 1, j, k)],
                                    idx[(i + 1, j + 1, k)], idx[(i, j + 1, k)]]))

    return nodes, hexes, quads


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


def make_plate2d(nw=40, nt=12, ng=8, w=1e-3, t=0.05e-3, d=0.2e-3, swap=0, rot=0):
    """平行平板線路の**断面 2 次元**格子 (三角形が体積要素、M / F 用)

    構造格子版 plate_line_dc / plate_line_ac と同じ形状を y-z 断面で切る
    (伝送線路軸は x)。1 次元厳密解がそのまま使えるので、非構造格子の
    M / F をそこに突き合わせられる。

      導体 1 : z = d .. d+t     (幅 W)
      導体 0 : z = -t .. 0      (幅 W、帰路)
      その間 (z = 0 .. d) は真空

    z 方向は導体内を nt 分割、間隙を ng 分割。表皮効果を刻めるよう
    導体内は両面に向けて等比に細かくする。

    swap = 1 で面内の 2 軸 (y, z) を入れ替える。形状を面内で 90 度回すだけなので
    答えは変わらないが、Az が変化する向きが p 軸から q 軸に移る。1 次元解では
    面内勾配の片方しか立たないので、**両方の格子を回さないと面内テンソルの
    片側の成分が一度も検査されない** (実測: c_pp を壊す変異が素通りした)。

    rot = <度> で面内に任意角だけ回す。**格子の位相は変えないので、回した格子は
    元の格子と合同な離散問題そのもの**になる。材料テンソルも同じ角で回せば
    答えは変わらないはずで、これが非対角成分を検査する唯一の手掛かりになる
    (1 次元解では ∂Az/∂p = 0 なので非対角項が効かず、同軸は回転対称なので
    テンソルを回しても答えが変わらない。斜めに回した非対称断面だけが効く)。

    物理タグ : 1 = 真空、2 = 導体の材料、10 = 導体 0 の断面、11 = 導体 1 の断面
    """
    def graded(z0, z1, n, ratio=1.35):
        # 両端に向けて細かくする対称分割
        h = n // 2
        wts = [ratio ** i for i in range(h)]
        wts = wts + wts[::-1] if (n % 2 == 0) else wts + [ratio ** h] + wts[::-1]
        s = sum(wts)
        out, acc = [z0], 0.0
        for wt in wts:
            acc += wt
            out.append(z0 + (z1 - z0) * acc / s)
        out[-1] = z1
        return out

    # z の分割点と、各層の (材料タグ, 導体タグ)
    zs, lay = [], []
    for (a0, b0, n, tag, cnd) in ((-t, 0.0, nt, 2, 10),
                                  (0.0, d, ng, 1, None),
                                  (d, d + t, nt, 2, 11)):
        g = graded(a0, b0, n)
        if zs:
            g = g[1:]
        zs += g
        lay += [(tag, cnd)] * (len(g) if not lay else len(g))
    # lay は層 (区間) 毎なので長さを合わせ直す
    lay = ([ (2, 10) ] * nt) + ([ (1, None) ] * ng) + ([ (2, 11) ] * nt)

    ys = [(-w / 2) + (w * j / nw) for j in range(nw + 1)]
    nzt = len(zs) - 1

    nodes, idx = [], {}
    th = math.radians(rot)
    cs, sn = math.cos(th), math.sin(th)
    for j in range(nw + 1):
        for k in range(nzt + 1):
            idx[(j, k)] = len(nodes)
            a1, a2 = (zs[k], ys[j]) if swap else (ys[j], zs[k])
            nodes.append((0.0, (a1 * cs) - (a2 * sn), (a1 * sn) + (a2 * cs)))

    tris = []
    for j in range(nw):
        for k in range(nzt):
            mtag, ctag = lay[k]
            tag = ctag if ctag is not None else mtag
            a = idx[(j, k)]
            b = idx[(j + 1, k)]
            c = idx[(j + 1, k + 1)]
            e = idx[(j, k + 1)]
            tris.append((tag, [a, b, c]))
            tris.append((tag, [a, c, e]))
    return nodes, [], tris


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
    opt41 = opt.pop("v41", 0)			# -v41 1 で Gmsh 4.1 形式で書く

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
    elif kind == "plate2d":
        nodes, tets, tris = make_plate2d(**opt)
    elif kind in ("box_prism_warp", "box_mixed"):
        nodes, cells = (make_box_prism_warp(**opt) if kind == "box_prism_warp"
                        else make_box_mixed(**opt))
        write_msh_cells(path, nodes, cells)
        print("%s : prisms, %d nodes, %d cells" % (path, len(nodes), len(cells)))
        return 0
    elif kind in ("box_pyr", "box_hexpyrtet"):
        nodes, cells = (make_box_pyr(**opt) if kind == "box_pyr"
                        else make_box_hexpyrtet(**opt))
        write_msh_cells(path, nodes, cells)
        print("%s : pyramids, %d nodes, %d cells" % (path, len(nodes), len(cells)))
        return 0
    elif kind in ("coax_hex2", "box_hex2_warp", "box_prism2_warp"):
        nodes, cells = (make_coax_hex2(**opt) if kind == "coax_hex2"
                        else make_box_hex2_warp(**opt) if kind == "box_hex2_warp"
                        else make_box_prism2_warp(**opt))
        write_msh_cells(path, nodes, cells)
        print("%s : order 2, %d nodes, %d cells" % (path, len(nodes), len(cells)))
        return 0
    elif kind in ("coax_hex", "box_hex_warp"):
        nodes, hexes, quads = (make_coax_hex(**opt) if kind == "coax_hex"
                               else make_box_hex_warp(**opt))
        write_msh_hex(path, nodes, hexes, quads)
        print("%s : hexahedra, %d nodes, %d hexes, %d quads"
              % (path, len(nodes), len(hexes), len(quads)))
        return 0
    else:
        print("unknown mesh kind: %s" % kind)
        return 1
    if order == 2:
        nodes, tets, tris = to_order2(nodes, tets, tris, snap)
    if opt41:
        write_msh41(path, nodes, tets, tris)
    else:
        write_msh(path, nodes, tets, tris)
    print("%s : order %d, %d nodes, %d tets, %d tris"
          % (path, order, len(nodes), len(tets), len(tris)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
