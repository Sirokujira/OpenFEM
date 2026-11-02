#!/usr/bin/env python3
"""Gmsh のバイナリ格子 (2.2 / 3.0 / 4.1) のバイト順を入れ替える (別のエンディアンの計算機で gmsh が書いたファイルを作る)。

    python3 gmsh_swap_endian.py <入力 .msh> <出力 .msh>

手元に別のエンディアンの gmsh が無いので、gmsh が書いたバイナリの数値をすべて自然な幅 (int32 / size_t / double) の単位で
並べ替えて作る。並びは gmsh のマニュアル (MSH file format) の各版の節の定義から書いた (ソルバーの読み手とは独立)。
ASCII の行 ($MeshFormat の版の行・2.2 / 3.0 の個数の行・3.0 の $Entities・$PhysicalNames) はそのまま。

2 回かけると元に戻ること (バイト単位で一致) を検証に使う。並びを読み違えるとその時点で節の終わりの行が見つからずに止まる。
"""
import struct
import sys

# 要素型 -> 節点数 (gmsh の要素型の表。使う型と、境界の点・線分)
NNODE = {1: 2, 2: 3, 3: 4, 4: 4, 5: 8, 6: 6, 7: 5, 8: 3, 9: 6, 10: 9, 11: 10, 12: 27, 13: 18, 14: 14, 15: 1,
         16: 8, 17: 20, 18: 15, 19: 13, 20: 9, 21: 10, 22: 12, 23: 15, 24: 15, 25: 21, 26: 4, 27: 5, 28: 6,
         29: 20, 30: 35, 31: 56, 36: 16, 37: 25, 38: 36, 42: 28, 43: 36, 71: 84, 72: 120,
         90: 40, 91: 75, 92: 64, 93: 125, 94: 216, 106: 126}


class Swapper:
    def __init__(self, data):
        self.d = data
        self.p = 0
        self.out = bytearray()
        self.src = '<'			# 入力のエンディアン (マーカーで決める)

    def line(self):
        q = self.d.index(b'\n', self.p) + 1
        s = self.d[self.p:q]
        self.out += s
        self.p = q
        return s.decode('latin-1').strip()

    def nums(self, fmt, count):
        """count 個の数を読んで、逆のエンディアンで書く。値を返す"""
        size = struct.calcsize('<' + fmt)
        vals = struct.unpack_from(self.src + fmt * count, self.d, self.p)
        self.p += size * count
        dst = '>' if self.src == '<' else '<'
        self.out += struct.pack(dst + fmt * count, *vals)
        return vals

    def rest_newline(self):
        # バイナリの塊の後ろの改行 (無いこともある)
        while self.p < len(self.d) and self.d[self.p:self.p + 1] in (b'\n', b'\r'):
            self.out += self.d[self.p:self.p + 1]
            self.p += 1

    def expect(self, tag):
        self.rest_newline()
        s = self.line()
        if s != tag:
            raise SystemExit('layout error: expected %s, got %r at byte %d' % (tag, s[:40], self.p))


def run(data):
    w = Swapper(data)
    ver = None
    while w.p < len(data):
        s = w.line()
        if s == '$MeshFormat':
            hdr = w.line().split()
            ver = hdr[0]
            if hdr[1] != '1':
                raise SystemExit('not a binary file')
            if hdr[2] != '8':
                raise SystemExit('only 8-byte floats')
            raw = struct.unpack_from('<i', data, w.p)[0]
            w.src = '<' if raw == 1 else '>'
            w.nums('i', 1)
            w.expect('$EndMeshFormat')
        elif s == '$Entities' and ver.startswith('4'):
            npt, ncu, nsu, nvo = w.nums('Q', 4)
            for _ in range(npt):
                w.nums('i', 1)
                w.nums('d', 3)
                (nph,) = w.nums('Q', 1)
                w.nums('i', nph)
            for n in (ncu, nsu, nvo):
                for _ in range(n):
                    w.nums('i', 1)
                    w.nums('d', 6)
                    (nph,) = w.nums('Q', 1)
                    w.nums('i', nph)
                    (nb,) = w.nums('Q', 1)
                    w.nums('i', nb)
            w.expect('$EndEntities')
        elif s == '$Nodes':
            if ver.startswith('4'):
                nblk, nn, _, _ = w.nums('Q', 4)
                for _ in range(nblk):
                    _, _, par = w.nums('i', 3)
                    if par != 0:
                        raise SystemExit('parametric nodes are not handled')
                    (cnt,) = w.nums('Q', 1)
                    w.nums('Q', cnt)
                    w.nums('d', 3 * cnt)
            elif ver.startswith('3'):
                nn = int(w.line())
                for _ in range(nn):
                    w.nums('i', 1)
                    w.nums('d', 3)
                    w.nums('i', 1)		# エンティティ (3.0 の節点は 1 行に「番号 x y z エンティティ」)
            else:
                nn = int(w.line())
                for _ in range(nn):
                    w.nums('i', 1)
                    w.nums('d', 3)
            w.expect('$EndNodes')
        elif s == '$Elements':
            if ver.startswith('4'):
                nblk, ne, _, _ = w.nums('Q', 4)
                for _ in range(nblk):
                    _, _, typ = w.nums('i', 3)
                    (cnt,) = w.nums('Q', 1)
                    w.nums('Q', cnt * (1 + NNODE[typ]))
            elif ver.startswith('3'):
                ne = int(w.line())
                for _ in range(ne):
                    _, typ, _, nv = w.nums('i', 4)
                    w.nums('i', nv)
            else:
                ne = int(w.line())
                done = 0
                while done < ne:
                    typ, cnt, ntag = w.nums('i', 3)
                    w.nums('i', cnt * (1 + ntag + NNODE[typ]))
                    done += cnt
            w.expect('$EndElements')
        # 他の節 ($PhysicalNames、3.0 の $Entities) は ASCII なので行ごとにそのまま写す
    return bytes(w.out)


if __name__ == '__main__':
    with open(sys.argv[1], 'rb') as f:
        data = f.read()
    out = run(data)
    if run(out) != data:
        raise SystemExit('swapping twice does not give the original back')
    with open(sys.argv[2], 'wb') as f:
        f.write(out)
