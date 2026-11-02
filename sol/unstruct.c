/*
unstruct.c

非構造格子 (4 節点 / 10 節点四面体) の読み込み・CRS 構築・要素行列。

格子は Gmsh ASCII 2.2 形式 (.msh) で与える。物理タグで領域と電極を指定し、
`region` / `electrode` キーで材料番号・導体番号に対応づける。

構造格子との違いは「節点の並びと隣接関係」だけなので、Dirichlet の扱い・
反作用からの電荷抽出・反復解法は構造格子版をそのまま使える。
現時点では静電系 (C / L / R) のみ対応する (M / F は断面 2 次元の定式化なので
構造格子専用)。

要素次数は格子ファイルから決める (Gmsh の要素型 4 / 2 なら 1 次、
11 / 9 なら 2 次)。混在は受け付けない。2 次要素は等パラメトリックなので、
辺上の中間節点を境界形状に載せれば曲面をそのまま表せる。
*/

#include "fem.h"
#include "fem_prototype.h"

// 格子ファイルの 1 行の最大長 (BUFSIZ は MSVC で 512 しかなく、64 節点の要素の行が入らない)
#define MESH_LINE (8192)

#define ARRAY_INC (100000)


// ---- Gmsh ASCII 2.2 の読み込み ----

static int read_nodes(FILE *fp, int32_t **idmap, int32_t *maxid)
{
	char line[MESH_LINE];

	if (fgets(line, sizeof(line), fp) == NULL) return 1;
	const int nn = atoi(line);
	if (nn < 4) {
		printf("*** mesh : too few nodes (%d)\n", nn);
		return 1;
	}

	NNode = nn;
	Xp = (double *)malloc((size_t)nn * sizeof(double));
	Yp = (double *)malloc((size_t)nn * sizeof(double));
	Zp = (double *)malloc((size_t)nn * sizeof(double));

	// Gmsh の節点番号は 1 始まりだが連続とは限らないので対応表を作る
	int32_t *id = (int32_t *)malloc((size_t)nn * sizeof(int32_t));
	int32_t mx = 0;
	for (int i = 0; i < nn; i++) {
		if (fgets(line, sizeof(line), fp) == NULL) return 1;
		long gid = 0;
		double x = 0, y = 0, z = 0;
		if (sscanf(line, "%ld %lf %lf %lf", &gid, &x, &y, &z) != 4) {
			printf("*** mesh : invalid node line %d\n", i + 1);
			free(id);
			return 1;
		}
		id[i] = (int32_t)gid;
		if (id[i] > mx) mx = id[i];
		Xp[i] = x;
		Yp[i] = y;
		Zp[i] = z;
	}

	int32_t *map = (int32_t *)malloc(((size_t)mx + 1) * sizeof(int32_t));
	for (int32_t i = 0; i <= mx; i++) map[i] = -1;
	for (int i = 0; i < nn; i++) map[id[i]] = i;
	free(id);

	*idmap = map;
	*maxid = mx;

	return 0;
}


// 要素の格納は 2.2 と 4.1 で共通 (読み方だけが違う)

// 最初の同じ種別の要素と節点数の違う四面体・六面体・角柱 (elem_finish の ho_raise で最高次数に上げる)
static int PendN = 0, PendCap = 0;
static int *PendKind = NULL, *PendNen = NULL, *PendTag = NULL;	// 種別 (MESHELEM_TET / HEX / PRISM)、節点数
static int32_t *PendNd = NULL;		// [ELEM_NMAX * PendN] gmsh の局所の並び
static int ho_raise(void);
static void ho_pend(int kind, int nen, int tag, const int32_t *nd)
{
	if (PendN == PendCap) {
		PendCap = ((PendCap > 0) ? (2 * PendCap) : 256);
		PendKind = (int *)realloc(PendKind, (size_t)PendCap * sizeof(int));
		PendNen = (int *)realloc(PendNen, (size_t)PendCap * sizeof(int));
		PendTag = (int *)realloc(PendTag, (size_t)PendCap * sizeof(int));
		PendNd = (int32_t *)realloc(PendNd, (size_t)PendCap * ELEM_NMAX * sizeof(int32_t));
	}
	PendKind[PendN] = kind;
	PendNen[PendN] = nen;
	PendTag[PendN] = tag;
	for (int l = 0; l < nen; l++) PendNd[((size_t)PendN * ELEM_NMAX) + l] = nd[l];
	PendN++;
}

static void elem_reset(void)
{
	PendN = 0;
	MeshRaisedElem = 0;
	MeshRaisedNode = 0;
	MeshRaisedFrom = 0;
	MeshRaisedKinds = 0;
	NTet = 0;
	NTri = 0;
	Tet = NULL;
	TetTag = NULL;
	Tri = NULL;
	TriTag = NULL;
	Tet2 = NULL;
	Tri2 = NULL;
	Tet3 = NULL;
	Tri3 = NULL;
	TetOrder = 0;			// 最初に出た四面体で決まる
	NHex = 0;
	Hex = NULL;
	HexTag = NULL;
	HexNen = 0;			// 最初に出た六面体で決まる (8 / 20 / 27)
	Hex2 = NULL;
	NQuad = 0;
	Quad = NULL;
	QuadTag = NULL;
	Quad2 = NULL;
	NPrism = 0;
	Prism = NULL;
	PrismTag = NULL;
	PrismNen = 0;		// 最初に出た角柱で決まる (6 / 15 / 18)
	Prism2 = NULL;
	NPyr = 0;
	Pyr = NULL;
	PyrTag = NULL;
	MeshElem = MESHELEM_TET;
}


static int ho_tables_init(int p);
static int tri_rule_check(void);

// 三角形の Tri3 の欄 (頂点以外の節点、使わない欄は -1) から次数 (3〜7、それ以外は 0) を決める
static int tri3_order(const int32_t *t3)
{
	int nho = 0;
	while ((nho < TRI_HO_MAX) && (t3[nho] >= 0)) nho++;
	for (int p = 3; p <= 7; p++) {
		if (((((p + 1) * (p + 2)) / 2) - 3) == nho) return p;
	}
	return 0;
}


// 要素 1 個を格納する。type は Gmsh の要素型
// (4/11/29 = 四面体、2/9/21 = 三角形、5 = 六面体、6 = 角柱、7 = ピラミッド、3 = 四角形)
static int elem_store(int type, int tag, const int32_t *nd)
{
	const int order = (((type == 11) || (type == 9)) ? 2
	                 : ((type == 29) || (type == 21)) ? 3
	                 : ((type == 30) || (type == 23)) ? 4
	                 : ((type == 31) || (type == 25)) ? 5
	                 : ((type == 71) || (type == 42)) ? 6
	                 : ((type == 72) || (type == 43)) ? 7 : 1);

	if ((type == 5) || (type == 12) || (type == 17) || (type == 92) || (type == 93) || (type == 94)) {
		// 六面体 (8 / 27 / 20 節点)。局所の並びは Gmsh のまま使う。
		// 頂点 8 個を Hex に、2 次の残り (辺 12 + 面 6 + 体心 1) を Hex2 に
		// 入れる。次数は最初の 1 個で決め、以後は混在を許さない
		const int hn = ((type == 5) ? 8 : (type == 12) ? 27 : (type == 92) ? 64 : (type == 93) ? 125 : (type == 94) ? 216 : 20);
		if (HexNen == 0) HexNen = hn;
		if (HexNen != hn) {
			// 節点数の違う六面体 : 取っておき、読み終えてから最高次数に上げる (ho_raise)
			ho_pend(MESHELEM_HEX, hn, tag, nd);
			return 0;
		}
		if (NHex % ARRAY_INC == 0) {
			Hex = (int32_t *)realloc(Hex, (size_t)(NHex + ARRAY_INC) * 8 * sizeof(int32_t));
			HexTag = (int *)realloc(HexTag, (size_t)(NHex + ARRAY_INC) * sizeof(int));
			if (hn > 8) {
				Hex2 = (int32_t *)realloc(Hex2, (size_t)(NHex + ARRAY_INC) * HEX2_W * sizeof(int32_t));
			}
		}
		for (int l = 0; l < 8; l++) Hex[(NHex * 8) + l] = nd[l];
		for (int l = 8; l < hn; l++) Hex2[(NHex * HEX2_W) + (l - 8)] = nd[l];
		HexTag[NHex] = tag;
		NHex++;

		return 0;
	}
	if ((type == 6) || (type == 13) || (type == 18) || (type == 90) || (type == 91) || (type == 106)) {
		// 角柱 (6 / 18 / 15 節点)。局所の並びは Gmsh のまま使う。
		// 頂点 6 個を Prism に、2 次の残り (辺 9 + 四角形面 3) を Prism2 に入れる
		const int pn = ((type == 6) ? 6 : (type == 13) ? 18 : (type == 90) ? 40 : (type == 91) ? 75 : (type == 106) ? 126 : 15);
		if (PrismNen == 0) PrismNen = pn;
		if (PrismNen != pn) {
			// 節点数の違う角柱 : 取っておき、読み終えてから最高次数に上げる (ho_raise)
			ho_pend(MESHELEM_PRISM, pn, tag, nd);
			return 0;
		}
		if (NPrism % ARRAY_INC == 0) {
			Prism = (int32_t *)realloc(Prism, (size_t)(NPrism + ARRAY_INC) * 6 * sizeof(int32_t));
			PrismTag = (int *)realloc(PrismTag, (size_t)(NPrism + ARRAY_INC) * sizeof(int));
			if (pn > 6) {
				Prism2 = (int32_t *)realloc(Prism2, (size_t)(NPrism + ARRAY_INC) * PRISM2_W * sizeof(int32_t));
			}
		}
		for (int l = 0; l < 6; l++) Prism[(NPrism * 6) + l] = nd[l];
		for (int l = 6; l < pn; l++) Prism2[(NPrism * PRISM2_W) + (l - 6)] = nd[l];
		PrismTag[NPrism] = tag;
		NPrism++;

		return 0;
	}
	if (type == 7) {
		// ピラミッド (5 節点)。局所の並びは Gmsh のまま使う
		if (NPyr % ARRAY_INC == 0) {
			Pyr = (int32_t *)realloc(Pyr, (size_t)(NPyr + ARRAY_INC) * 5 * sizeof(int32_t));
			PyrTag = (int *)realloc(PyrTag, (size_t)(NPyr + ARRAY_INC) * sizeof(int));
		}
		for (int l = 0; l < 5; l++) Pyr[(NPyr * 5) + l] = nd[l];
		PyrTag[NPyr] = tag;
		NPyr++;

		return 0;
	}
	if ((type == 3) || (type == 16) || (type == 10) || (type == 36) || (type == 37) || (type == 38)) {
		// 四角形 (境界面、4 / 8 / 9 節点)。六面体格子の電極面はこれになる。
		// 2 次の中間節点 (辺 4 + 中心) は Quad2 に入れる。無い位置は -1 で
		// 埋める (realloc した領域は未初期化なので、必ず全 5 要素を書くこと)
		const int qn = ((type == 3) ? 4 : (type == 16) ? 8 : (type == 36) ? 16 : (type == 37) ? 25 : (type == 38) ? 36 : 9);
		if (NQuad % ARRAY_INC == 0) {
			Quad = (int32_t *)realloc(Quad, (size_t)(NQuad + ARRAY_INC) * 4 * sizeof(int32_t));
			QuadTag = (int *)realloc(QuadTag, (size_t)(NQuad + ARRAY_INC) * sizeof(int));
			Quad2 = (int32_t *)realloc(Quad2, (size_t)(NQuad + ARRAY_INC) * QUAD2_W * sizeof(int32_t));
		}
		for (int l = 0; l < 4; l++) Quad[(NQuad * 4) + l] = nd[l];
		for (int l = 0; l < QUAD2_W; l++) {
			Quad2[(NQuad * QUAD2_W) + l] = (((4 + l) < qn) ? nd[4 + l] : -1);
		}
		QuadTag[NQuad] = tag;
		NQuad++;

		return 0;
	}

	if ((type == 4) || (type == 11) || (type == 29) || (type == 30) || (type == 31) || (type == 71) || (type == 72)) {
		// 四面体。次数は最初の 1 個で決め、以後は混在を許さない
		if (TetOrder == 0) TetOrder = order;
		if (TetOrder != order) {
			// 次数の違う四面体 : 取っておき、読み終えてから全部を最高次数に上げる (ho_raise)
			ho_pend(MESHELEM_TET, (((order + 1) * (order + 2) * (order + 3)) / 6), tag, nd);
			return 0;
		}
		if (NTet % ARRAY_INC == 0) {
			Tet = (int32_t *)realloc(Tet, (size_t)(NTet + ARRAY_INC) * 4 * sizeof(int32_t));
			TetTag = (int *)realloc(TetTag, (size_t)(NTet + ARRAY_INC) * sizeof(int));
			if (order == 2) {
				Tet2 = (int32_t *)realloc(Tet2, (size_t)(NTet + ARRAY_INC) * 6 * sizeof(int32_t));
			}
			if (order >= 3) {
				Tet3 = (int32_t *)realloc(Tet3, (size_t)(NTet + ARRAY_INC) * TET_HO_MAX * sizeof(int32_t));
			}
		}
		for (int l = 0; l < 4; l++) Tet[(NTet * 4) + l] = nd[l];
		if (order == 2) {
			for (int l = 0; l < 6; l++) Tet2[(NTet * 6) + l] = nd[4 + l];
		}
		if (order >= 3) {
			// 頂点以外 (3 次 16、4 次 31、5 次 52、6 次 80、7 次 116)。使わない欄は -1
			const int nho = ((((order + 1) * (order + 2) * (order + 3)) / 6) - 4);
			for (int l = 0; l < TET_HO_MAX; l++) Tet3[(NTet * TET_HO_MAX) + l] = ((l < nho) ? nd[4 + l] : -1);
		}
		TetTag[NTet] = tag;
		NTet++;
	}
	else {
		if (NTri % ARRAY_INC == 0) {
			Tri = (int32_t *)realloc(Tri, (size_t)(NTri + ARRAY_INC) * 3 * sizeof(int32_t));
			TriTag = (int *)realloc(TriTag, (size_t)(NTri + ARRAY_INC) * sizeof(int));
			Tri2 = (int32_t *)realloc(Tri2, (size_t)(NTri + ARRAY_INC) * 3 * sizeof(int32_t));
			Tri3 = (int32_t *)realloc(Tri3, (size_t)(NTri + ARRAY_INC) * TRI_HO_MAX * sizeof(int32_t));
		}
		for (int l = 0; l < 3; l++) Tri[(NTri * 3) + l] = nd[l];
		// 3 次 (10 節点) / 4 次 (15 節点) の頂点以外の 7 / 12 節点。使わない欄は -1 で埋める
		// (realloc した領域は未初期化なので必ず TRI_HO_MAX 要素とも書くこと)
		{
			const int nho = ((order >= 3) ? ((((order + 1) * (order + 2)) / 2) - 3) : 0);
			for (int l = 0; l < TRI_HO_MAX; l++) {
				Tri3[(NTri * TRI_HO_MAX) + l] = ((l < nho) ? nd[3 + l] : -1);
			}
		}
		// 1 次の三角形では中間節点が無いので頂点で埋める (Dirichlet の
		// 塗り分けは重複しても同じ値になるので無害)
		for (int l = 0; l < 3; l++) {
			Tri2[(NTri * 3) + l] = ((order == 2) ? nd[3 + l] : nd[l]);
		}
		TriTag[NTri] = tag;
		NTri++;
	}

	return 0;
}


static int cmp_face3(const void *a, const void *b)
{
	const int32_t *x = (const int32_t *)a;
	const int32_t *y = (const int32_t *)b;

	for (int i = 0; i < 3; i++) {
		if (x[i] != y[i]) return ((x[i] < y[i]) ? -1 : 1);
	}

	return 0;
}


static int cmp_face4(const void *a, const void *b)
{
	const int32_t *x = (const int32_t *)a;
	const int32_t *y = (const int32_t *)b;

	for (int i = 0; i < 4; i++) {
		if (x[i] != y[i]) return ((x[i] < y[i]) ? -1 : 1);
	}

	return 0;
}


// 面の節点を昇順に並べる (表裏・巡回によらず同じ面を同じ鍵にする)
static void face_sort(int32_t *v, int n)
{
	for (int i = 1; i < n; i++) {
		const int32_t t = v[i];
		int j = i;
		for (; (j > 0) && (v[j - 1] > t); j--) v[j] = v[j - 1];
		v[j] = t;
	}
}


/*
「見かけだけ接している」四角形面 - 三角形面の接続の検出。

六面体・角柱・ピラミッドの四角形面を四面体 (やピラミッド) の三角形 2 枚で
覆うと、**節点はすべて共有されるのに面上の解が食い違う** (四角形面の
トレースは双 1 次で交差項を持ち、三角形 2 枚の 1 次のトレースでは表せない)。
節点の共有しか見ない検査 (線形場の恒等式) はこれを検出できない — 1 次の場は
どちらの側でも厳密に表せるため。そこで位相で捕まえる:

  適合した内部の四角形面は表裏 2 回現れる。**1 回しか現れない四角形面**の
  4 節点のうち 3 節点が、**1 回しか現れない三角形面**と一致していれば、
  その四角形は対角線で 2 枚の三角形に割られている = 非適合として弾く。

三角形側も 1 回しか現れないものに限るのは、内部の三角形面 (2 回現れる =
適合) との偶然の一致を除くため。境界の折り目が四角形面と 3 節点を共有する
形は偽陽性になり得るが、それは面どうしが幾何的に交差する病的な格子に限る。
正しいケースで 1 件も出ないことは rlc_check.sh の全ケースが兼ねる。
*/
static int elem_face_check(void)
{
	const int64_t nq = ((int64_t)6 * NHex) + ((int64_t)3 * NPrism) + NPyr;
	const int64_t nt = ((int64_t)4 * NTet) + ((int64_t)2 * NPrism) + ((int64_t)4 * NPyr);

	if ((nq < 1) || (nt < 1)) return 0;

	// 局所節点 -> 面 (並びは Gmsh の hex8 / prism6 / pyramid5)
	static const int hexf[6][4] = {{0, 1, 2, 3}, {4, 5, 6, 7}, {0, 1, 5, 4},
	                               {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}};
	static const int prif[3][4] = {{0, 1, 4, 3}, {1, 2, 5, 4}, {2, 0, 3, 5}};
	static const int tetf[4][3] = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};
	static const int pyrf[4][3] = {{0, 1, 4}, {1, 2, 4}, {2, 3, 4}, {3, 0, 4}};

	int32_t *fq = (int32_t *)malloc((size_t)nq * 4 * sizeof(int32_t));
	int64_t k = 0;
	for (int e = 0; e < NHex; e++) {
		for (int f = 0; f < 6; f++) {
			for (int i = 0; i < 4; i++) fq[(k * 4) + i] = Hex[(e * 8) + hexf[f][i]];
			face_sort(&fq[k * 4], 4);
			k++;
		}
	}
	for (int e = 0; e < NPrism; e++) {
		for (int f = 0; f < 3; f++) {
			for (int i = 0; i < 4; i++) fq[(k * 4) + i] = Prism[(e * 6) + prif[f][i]];
			face_sort(&fq[k * 4], 4);
			k++;
		}
	}
	for (int e = 0; e < NPyr; e++) {
		for (int i = 0; i < 4; i++) fq[(k * 4) + i] = Pyr[(e * 5) + i];
		face_sort(&fq[k * 4], 4);
		k++;
	}
	qsort(fq, (size_t)nq, 4 * sizeof(int32_t), cmp_face4);

	int32_t *ft = (int32_t *)malloc((size_t)nt * 3 * sizeof(int32_t));
	k = 0;
	for (int e = 0; e < NTet; e++) {
		for (int f = 0; f < 4; f++) {
			for (int i = 0; i < 3; i++) ft[(k * 3) + i] = Tet[(e * 4) + tetf[f][i]];
			face_sort(&ft[k * 3], 3);
			k++;
		}
	}
	for (int e = 0; e < NPrism; e++) {
		for (int f = 0; f < 2; f++) {
			for (int i = 0; i < 3; i++) ft[(k * 3) + i] = Prism[(e * 6) + (f * 3) + i];
			face_sort(&ft[k * 3], 3);
			k++;
		}
	}
	for (int e = 0; e < NPyr; e++) {
		for (int f = 0; f < 4; f++) {
			for (int i = 0; i < 3; i++) ft[(k * 3) + i] = Pyr[(e * 5) + pyrf[f][i]];
			face_sort(&ft[k * 3], 3);
			k++;
		}
	}
	qsort(ft, (size_t)nt, 3 * sizeof(int32_t), cmp_face3);

	// 1 回しか現れない三角形面だけ残す (詰め直し)
	int64_t nt1 = 0;
	for (int64_t p = 0; p < nt; ) {
		int64_t q = p + 1;
		while ((q < nt) && (cmp_face3(&ft[p * 3], &ft[q * 3]) == 0)) q++;
		if ((q - p) == 1) {
			for (int i = 0; i < 3; i++) ft[(nt1 * 3) + i] = ft[(p * 3) + i];
			nt1++;
		}
		p = q;
	}

	// 1 回しか現れない四角形面の 3 節点部分集合を照合する
	int bad = 0;
	for (int64_t p = 0; (p < nq) && !bad; ) {
		int64_t q = p + 1;
		while ((q < nq) && (cmp_face4(&fq[p * 4], &fq[q * 4]) == 0)) q++;
		if ((q - p) == 1) {
			for (int drop = 0; (drop < 4) && !bad; drop++) {
				int32_t key[3];
				int m = 0;
				for (int i = 0; i < 4; i++) {
					if (i != drop) key[m++] = fq[(p * 4) + i];
				}
				if ((nt1 > 0) && (bsearch(key, ft, (size_t)nt1,
						3 * sizeof(int32_t), cmp_face3) != NULL)) {
					bad = 1;
				}
			}
		}
		p = q;
	}

	free(fq);
	free(ft);

	if (bad) {
		printf("%s\n", "*** mesh : a quadrilateral face is covered by triangular "
			"faces, so the mesh is not conforming (the bilinear trace on the "
			"quadrilateral cannot match the linear traces on the triangles); "
			"put pyramids on such faces to make the transition");
		return 1;
	}

	return 0;
}


// node_merge_check() の比較関数 (qsort に渡すので座標へのポインタは静的に持つ)
static const double *NodeSortX, *NodeSortY, *NodeSortZ;

static int cmp_node_xyz(const void *va, const void *vb)
{
	const int32_t a = *(const int32_t *)va;
	const int32_t b = *(const int32_t *)vb;

	if (NodeSortX[a] != NodeSortX[b]) return ((NodeSortX[a] < NodeSortX[b]) ? -1 : 1);
	if (NodeSortY[a] != NodeSortY[b]) return ((NodeSortY[a] < NodeSortY[b]) ? -1 : 1);
	if (NodeSortZ[a] != NodeSortZ[b]) return ((NodeSortZ[a] < NodeSortZ[b]) ? -1 : 1);

	return 0;
}


/*
**未マージの格子 (同じ座標に別番号の節点がある) を弾く。**

複数の部品を別々に切って貼り合わせた格子で節点をマージし忘れると、界面で
要素どうしが自由度を共有しなくなり、有限要素空間が「裂けた」ものになる。
**これは自己検証のどの恒等式にも引っかからない** (実測: 界面の 21 節点を
複製した格子で、辺要素の (a)〜(d) と接線連続性 (g) がすべて機械精度で通り、
連結成分が 1 -> 2 に増えたことだけが痕跡として残った。裂けた面はもはや
「共有面」ではないので接線連続性の検査すら素通りする)。要素毎の量の和で
書ける恒等式は裂け目を見ないので、格子そのものを見るしかない。

比較は**厳密な等値**で行う (許容値を持ち込まない)。同じ点から複製された節点は
書き出しの丸めまで一致するので、これで十分に捕まる。
種別の混在を許すと部品ごとに切る使い方が増えるので、この検査が要る。

**見るのは要素の頂点だけ** (2 次要素の中間節点は除く)。曲面に載せた 2 次格子
では、**別々の辺**の中点が丸めまで同じ点に落ちることがある (実測: coax_p2 で
48 組。どれも異なる頂点対の中点で、同じ辺の中点が 2 つある「本物の裂け目」は
1 組も無かった)。中間節点まで見ると、この正しい格子を弾く誤検知になる。
裂けた格子は必ず**頂点**を複製するので、頂点だけで検出できる。
*/
static int node_merge_check(void)
{
	const int n = NNode;
	if (n < 2) return 0;

	// 要素の頂点 (1 次の節点) だけを対象にする
	unsigned char *cor = (unsigned char *)malloc((size_t)n * sizeof(unsigned char));
	memset(cor, 0, (size_t)n * sizeof(unsigned char));
	for (int e = 0; e < NTet; e++) {
		for (int l = 0; l < 4; l++) cor[Tet[(e * 4) + l]] = 1;
	}
	for (int e = 0; e < NHex; e++) {
		for (int l = 0; l < 8; l++) cor[Hex[(e * 8) + l]] = 1;
	}
	for (int e = 0; e < NPrism; e++) {
		for (int l = 0; l < 6; l++) cor[Prism[(e * 6) + l]] = 1;
	}
	for (int e = 0; e < NPyr; e++) {
		for (int l = 0; l < 5; l++) cor[Pyr[(e * 5) + l]] = 1;
	}
	for (int e = 0; e < NTri; e++) {
		for (int l = 0; l < 3; l++) cor[Tri[(e * 3) + l]] = 1;
	}
	for (int e = 0; e < NQuad; e++) {
		for (int l = 0; l < 4; l++) cor[Quad[(e * 4) + l]] = 1;
	}

	int32_t *ord = (int32_t *)malloc((size_t)n * sizeof(int32_t));
	int nc = 0;
	for (int i = 0; i < n; i++) {
		if (cor[i]) ord[nc++] = (int32_t)i;
	}
	free(cor);
	if (nc < 2) {
		free(ord);
		return 0;
	}
	NodeSortX = Xp;
	NodeSortY = Yp;
	NodeSortZ = Zp;
	qsort(ord, (size_t)nc, sizeof(int32_t), cmp_node_xyz);

	int ndup = 0;
	int32_t a = -1, b = -1;
	for (int i = 1; i < nc; i++) {
		const int32_t p = ord[i - 1], q = ord[i];
		if ((Xp[p] == Xp[q]) && (Yp[p] == Yp[q]) && (Zp[p] == Zp[q])) {
			if (!ndup) { a = p; b = q; }
			ndup++;
		}
	}
	free(ord);

	if (ndup > 0) {
		printf("*** mesh : %d node(s) share the coordinates of another node "
			"(for instance %d and %d at %.6e %.6e %.6e); the parts of the mesh "
			"were not merged, so the elements on the two sides do not share "
			"degrees of freedom and the finite element space is torn\n",
			ndup, (int)a + 1, (int)b + 1, Xp[a], Yp[a], Zp[a]);
		return 1;
	}

	return 0;
}


// 全要素を読んだあとの整合性検査 (格子の次元と次数)
static int hexp_order(int nen);
static int prismp_order(int nen);

static int elem_finish(void)
{
	if (node_merge_check()) return 1;
	if (ho_raise()) return 1;

	// 3 次元の格子では 3 次 (10 節点) の三角形は 20 節点の四面体の面としてだけ使う。
	// 四面体が 3 次でないのに 10 節点の三角形があると、その面内の 7 節点がどの要素にも
	// 属さず行列が特異になる。断面 2 次元 (四面体が無い) の 3 次三角形は体積要素として
	// 下で次数を決める
	for (int t = 0; (NTet > 0) && (t < NTri); t++) {
		const int to = tri3_order(&Tri3[t * TRI_HO_MAX]);
		if ((to >= 3) && (TetOrder != to)) {
			printf("*** mesh : %d-node (order %d) triangles are supported only as the faces of %d-node "
				"tetrahedra in a 3-D mesh (no mixing with other orders; the tetrahedra are order %d)\n",
				((to + 1) * (to + 2)) / 2, to, ((to + 1) * (to + 2) * (to + 3)) / 6, TetOrder);
			return 1;
		}
		if ((TetOrder >= 3) && (to != TetOrder)) {
			printf("%s\n", "*** mesh : the tetrahedra are order 3 / 4 but a triangle is not "
				"(the electrode faces must be 10-node / 15-node triangles; regenerate the mesh "
				"with the same -order)");
			return 1;
		}
	}

	/*
	3 次元の要素種別を決める。**混在 (四面体 + 六面体 + 角柱 + ピラミッド) を
	許す。**

	要素は「四面体 -> 六面体 -> 角柱 -> ピラミッド」の連番で扱う (`elem3d_*`)。
	純粋な四面体格子ではこの順序が従来と完全に一致するので、既存の答えは
	ビット単位で変わらない。2 次要素は四面体にしかないので、混ざるときは
	1 次に限る。
	*/
	{
		const int nk = ((NTet > 0) ? 1 : 0) + ((NHex > 0) ? 1 : 0)
		             + ((NPrism > 0) ? 1 : 0) + ((NPyr > 0) ? 1 : 0);
		if (nk > 1) {
			/*
			**六面体と四面体は直接隣り合わせにできない。**
			六面体の面は四角形で、その上の解は双 1 次 (xy の項を持つ)。
			四面体の面は三角形で解は 1 次なので、六面体の四角形面を三角形 2 枚で
			覆っても**面上の解が一致しない** (節点は合っていても要素間で場が
			食い違う)。これを埋めるのがピラミッド (四角形面 1 + 三角形面 4 の
			遷移要素) で、格子にピラミッドが 1 つも無ければその意図も無いと
			みなして弾く。ピラミッドがあっても四角形面を三角形で覆った
			非適合の面が残っていないかは elem_face_check() が別に見る。

			角柱を挟むのは問題ない:
			  四面体 - 角柱 : 三角形の面どうし (どちらも 1 次) -> 適合
			  六面体 - 角柱 : 四角形の面どうし (どちらも双 1 次) -> 適合
			境界層を角柱で切って内部を四面体にする、という一番よくある形は通る。
			角柱があれば六面体 - 角柱 (側面) - 四面体 (上下の面) の 3 種別もつながるので弾かない
			(六面体と四面体が直接接していないかは同じく elem_face_check() が位相で見る)。
			*/
			if ((NTet > 0) && (NHex > 0) && (NPyr < 1) && (NPrism < 1)) {
				printf("%s\n", "*** mesh : hexahedra and tetrahedra cannot share a "
					"face conformingly (the quadrilateral face carries a bilinear "
					"trace, the triangular face a linear one); put prisms or "
					"pyramids between them, or use one element type");
				return 1;
			}
			/*
			**種別の混在は次数が揃っていれば高次でもよい** (四面体・六面体・角柱が同じ次数 p。ピラミッドは
			1 次だけ)。面のトレースは三角形なら P_p、四角形なら Q_p (完全 Lagrange) で、隣の種別と同じ空間
			になるので適合する。2 次の六面体・角柱は serendipity (20 / 15 節点) と完全 (27 / 18) を混ぜない
			(四角形面が 8 節点と 9 節点で食い違う)。次数の違う要素が面を共有する格子は遷移要素が要るので弾く
			*/
			const int ot = ((NTet > 0) ? TetOrder : 0);
			const int oh = ((NHex > 0) ? ((HexNen == 8) ? 1 : (HexNen <= 27) ? 2 : hexp_order(HexNen)) : 0);
			const int op = ((NPrism > 0) ? ((PrismNen == 6) ? 1 : (PrismNen <= 18) ? 2 : prismp_order(PrismNen)) : 0);
			const int oy = ((NPyr > 0) ? 1 : 0);
			const int pmix = ((ot > 0) ? ot : (oh > 0) ? oh : (op > 0) ? op : 1);
			if (((ot > 0) && (ot != pmix)) || ((oh > 0) && (oh != pmix)) || ((op > 0) && (op != pmix)) || ((oy > 0) && (oy != pmix))) {
				printf("*** mesh : mixed element kinds must have the same order (tetrahedra %d, hexahedra %d, "
					"prisms %d, pyramids %d; 0 = none). Elements of different orders cannot share a face "
					"without transition elements (pyramids are first order only)\n", ot, oh, op, oy);
				return 1;
			}
			if ((NHex > 0) && (NPrism > 0) && (((HexNen == 20) && (PrismNen != 15)) || ((HexNen == 27) && (PrismNen != 18))
			                                   || ((PrismNen == 15) && (HexNen != 20)) || ((PrismNen == 18) && (HexNen != 27)))) {
				printf("*** mesh : %d-node hexahedra and %d-node prisms do not share conforming quadrilateral faces "
					"(serendipity 20 / 15 and complete 27 / 18 must not be mixed)\n", HexNen, PrismNen);
				return 1;
			}
			// 高次の混在格子の電極面は同じ次数 (面の内部の節点まで固定する)
			if (pmix >= 2) {
				for (int t = 0; t < NTri; t++) {
					const int t3 = tri3_order(&Tri3[t * TRI_HO_MAX]);
					const int o = ((t3 >= 3) ? t3 : (Tri2[(t * 3)] != Tri[(t * 3)]) ? 2 : 1);
					if (o != pmix) {
						printf("*** mesh : the mixed mesh is order %d but a triangle is order %d\n", pmix, o);
						return 1;
					}
				}
				const int nqx = ((pmix == 2) ? (((HexNen == 20) || (PrismNen == 15)) ? 4 : 5) : (((pmix + 1) * (pmix + 1)) - 4));
				for (int t = 0; t < NQuad; t++) {
					if ((Quad2[(t * QUAD2_W) + nqx - 1] < 0) || ((nqx < QUAD2_W) && (Quad2[(t * QUAD2_W) + nqx] >= 0))) {
						printf("*** mesh : the mixed mesh is order %d but a quadrilateral does not have %d nodes\n",
							pmix, nqx + 4);
						return 1;
					}
				}
			}
			MeshElem = MESHELEM_MIXED;
			MeshDim = 3;
			TetOrder = pmix;
			if ((NTri < 1) && (NQuad < 1)) {
				printf("%s\n", "*** mesh : a mixed mesh needs triangular or "
					"quadrilateral boundary faces for the electrodes");
				return 1;
			}

			return elem_face_check();
		}
	}
	if (NPyr > 0) {
		// 純ピラミッド格子 (六面体を中心点で 6 分割した形など)。
		// 電極面は底面の四角形でも側面の三角形でもよい
		MeshElem = MESHELEM_PYR;
		MeshDim = 3;
		TetOrder = 1;
		if ((NTri < 1) && (NQuad < 1)) {
			printf("%s\n", "*** mesh : a pyramid mesh needs triangular or "
				"quadrilateral boundary faces for the electrodes");
			return 1;
		}

		return elem_face_check();
	}
	if (NPrism > 0) {
		// 角柱格子。電極面は上下の三角形でも側面の四角形でもよい
		MeshElem = MESHELEM_PRISM;
		MeshDim = 3;
		TetOrder = (prismp_order(PrismNen) ? prismp_order(PrismNen) : (PrismNen > 6) ? 2 : 1);
		if ((NTri < 1) && (NQuad < 1)) {
			printf("%s\n", "*** mesh : a prism mesh needs triangular or "
				"quadrilateral boundary faces for the electrodes");
			return 1;
		}
		// 3 次 / 4 次の角柱 (40 / 75 節点) の境界面は 10 / 15 節点三角形と 16 / 25 節点四角形
		// (面の内部の節点まで固定する。次数の違う面は電極が「穴だらけ」になるので弾く)
		if (TetOrder >= 3) {
			const int p = TetOrder, nqx = ((p + 1) * (p + 1)) - 4;
			for (int t = 0; t < NTri; t++) {
				const int o = tri3_order(&Tri3[t * TRI_HO_MAX]);
				if (o != p) {
					printf("*** mesh : the prisms are order %d (%d-node) but a triangle is not a %d-node "
						"triangle (regenerate the mesh with -order %d)\n", p, PrismNen, (p + 1) * (p + 2) / 2, p);
					return 1;
				}
			}
			for (int t = 0; t < NQuad; t++) {
				if ((Quad2[(t * QUAD2_W) + nqx - 1] < 0) || ((nqx < QUAD2_W) && (Quad2[(t * QUAD2_W) + nqx] >= 0))) {
					printf("*** mesh : the prisms are order %d (%d-node) but a quadrilateral is not a %d-node "
						"quadrilateral (regenerate the mesh with -order %d)\n", p, PrismNen, (p + 1) * (p + 1), p);
					return 1;
				}
			}
			return elem_face_check();
		}
		// 2 次の角柱に 1 次の境界面が混ざると電極の中間節点が固定されず、
		// 電極が「穴だらけ」になる (収束はするが Q が合わない)。tet10 と同じ規則
		if (TetOrder >= 2) {
			for (int t = 0; t < NTri; t++) {
				if (Tri2[(t * 3)] == Tri[(t * 3)]) {
					printf("%s\n", "*** mesh : the prisms are order 2 but a "
						"triangle is order 1 (regenerate the mesh with -order 2)");
					return 1;
				}
			}
			for (int t = 0; t < NQuad; t++) {
				if (Quad2[(t * QUAD2_W)] < 0) {
					printf("%s\n", "*** mesh : the prisms are order 2 but a "
						"quadrilateral is order 1 (regenerate the mesh with -order 2)");
					return 1;
				}
				// 18 節点の角柱の側面は面心つきの 9 節点四角形 (六面体と同じ穴)
				if ((PrismNen == 18) && (Quad2[(t * QUAD2_W) + 4] < 0)) {
					printf("%s\n", "*** mesh : the prisms are 18-node but a "
						"quadrilateral has no centre node (9-node quadrilaterals "
						"are required; do not mix Mesh.SecondOrderIncomplete)");
					return 1;
				}
			}
		}

		return elem_face_check();
	}
	if (NHex > 0) {
		// 六面体格子 (1 次の 8 節点、または 2 次の 20 / 27 節点)
		MeshElem = MESHELEM_HEX;
		MeshDim = 3;
		TetOrder = (hexp_order(HexNen) ? hexp_order(HexNen) : (HexNen > 8) ? 2 : 1);
		// 3 次 / 4 次の六面体 (64 / 125 節点) の電極面は 16 / 25 節点四角形 (面の内部の節点まで固定する)
		if (TetOrder >= 3) {
			const int p = TetOrder, nqx = ((p + 1) * (p + 1)) - 4;
			for (int t = 0; t < NQuad; t++) {
				if ((Quad2[(t * QUAD2_W) + nqx - 1] < 0) || ((nqx < QUAD2_W) && (Quad2[(t * QUAD2_W) + nqx] >= 0))) {
					printf("*** mesh : the hexahedra are order %d (%d-node) but a quadrilateral is not a %d-node "
						"quadrilateral (regenerate the mesh with -order %d)\n", p, HexNen, (p + 1) * (p + 1), p);
					return 1;
				}
			}
		}
		if (NQuad < 1) {
			printf("%s\n", "*** mesh : a hexahedral mesh needs quadrilateral "
				"boundary faces for the electrodes");
			return 1;
		}
		// 2 次の六面体に 1 次の四角形が混ざると電極の中間節点が固定されない。
		// **27 節点の六面体は面心節点まで要る** (8 節点四角形では辺の中間節点が
		// 埋まるので、先頭だけ見る検査を素通りして C が -41% になる。実測)
		if (TetOrder == 2) {
			for (int t = 0; t < NQuad; t++) {
				if (Quad2[(t * QUAD2_W)] < 0) {
					printf("%s\n", "*** mesh : the hexahedra are order 2 but a "
						"quadrilateral is order 1 (regenerate the mesh with -order 2)");
					return 1;
				}
				if ((HexNen == 27) && (Quad2[(t * QUAD2_W) + 4] < 0)) {
					printf("%s\n", "*** mesh : the hexahedra are 27-node but a "
						"quadrilateral has no centre node (9-node quadrilaterals "
						"are required; do not mix Mesh.SecondOrderIncomplete)");
					return 1;
				}
			}
		}

		return 0;
	}

	// 四面体が 1 つも無ければ断面 2 次元の格子として扱う (三角形が体積要素)。
	// M / F は断面 2 次元の定式化なので、この形でしか非構造格子に載らない
	MeshDim = ((NTet > 0) ? 3 : 2);
	if (MeshDim == 2) {
		if (NTri < 1) {
			printf("%s\n", "*** mesh : no tetrahedron and no triangle found");
			return 1;
		}
		// 三角形の次数は「中間節点が頂点と違うか」で決まる (1 次では頂点で埋めてある)。
		// 3 次 (10 節点) は Tri3 に辺の 6 節点と面心が入る (3 次でなければ -1)
		TetOrder = ((tri3_order(Tri3) >= 3) ? tri3_order(Tri3) : (Tri2[0] != Tri[0]) ? 2 : 1);
		for (int t = 0; t < NTri; t++) {
			const int t3 = tri3_order(&Tri3[t * TRI_HO_MAX]);
			const int o = ((t3 >= 3) ? t3 : (Tri2[(t * 3)] != Tri[(t * 3)]) ? 2 : 1);
			if (o != TetOrder) {
				printf("%s\n", "*** mesh : mixed triangle orders in a 2-D mesh");
				return 1;
			}
		}
	}

	// 三角形の次数が四面体と食い違うと電極面の中間節点が固定されず、
	// 電極が「穴だらけ」になる (収束はするが Q が合わない) ので弾く
	if ((MeshDim == 3) && (TetOrder == 2) && (NTri > 0)) {
		for (int t = 0; t < NTri; t++) {
			if (Tri2[(t * 3)] == Tri[(t * 3)]) {
				printf("%s\n", "*** mesh : the tetrahedra are order 2 but the "
					"triangles are order 1 (regenerate the mesh with -order 2)");
				return 1;
			}
		}
	}
	if (ho_tables_init(TetOrder)) return 1;

	return 0;
}


static int read_elements(FILE *fp, const int32_t *idmap, int32_t maxid)
{
	char line[MESH_LINE];

	if (fgets(line, sizeof(line), fp) == NULL) return 1;
	const int ne = atoi(line);

	elem_reset();

	for (int e = 0; e < ne; e++) {
		if (fgets(line, sizeof(line), fp) == NULL) return 1;

		// id type ntags tag1 ... tagn node1 ...
		int nv = 0;
		long v[ELEM_NMAX + 16];		// id type ntags tags... + 節点 (64 節点の六面体まで)
		char *p = line;
		while ((nv < (ELEM_NMAX + 16)) && (*p != '\0')) {
			while ((*p == ' ') || (*p == '\t')) p++;
			if ((*p == '\0') || (*p == '\n') || (*p == '\r')) break;
			v[nv++] = strtol(p, &p, 10);
		}
		if (nv < 4) continue;

		const int type = (int)v[1];
		const int ntag = (int)v[2];
		if (nv < 3 + ntag) continue;
		const int tag = ((ntag > 0) ? (int)v[3] : 0);
		const int off = 3 + ntag;

		// 節点数 : 四面体 (型 4 / 11) と三角形 (型 2 / 9)
		const int nn = ((type == 4) ? 4 : (type == 11) ? 10 : (type == 29) ? 20 : (type == 30) ? 35 : (type == 31) ? 56 : (type == 71) ? 84 : (type == 72) ? 120
		              : (type == 2) ? 3 : (type == 9) ? 6 : (type == 21) ? 10 : (type == 23) ? 15 : (type == 25) ? 21
		              : (type == 42) ? 28 : (type == 43) ? 36
		              : (type == 5) ? 8 : (type == 3) ? 4 : (type == 6) ? 6
		              : (type == 7) ? 5
		              : (type == 12) ? 27 : (type == 17) ? 20
		              : (type == 13) ? 18 : (type == 18) ? 15
		              : (type == 10) ? 9 : (type == 16) ? 8
		              : (type == 92) ? 64 : (type == 90) ? 40 : (type == 36) ? 16
		              : (type == 93) ? 125 : (type == 91) ? 75 : (type == 37) ? 25
		              : (type == 94) ? 216 : (type == 106) ? 126 : (type == 38) ? 36 : 0);
		if (nn == 0) continue;			// 点・線分など、使わない要素型
		if (nv < off + nn) continue;

		// 節点番号を先に解決する (未解決なら打ち切り。ここで continue すると
		// 配列の該当要素が未初期化のまま確定し、setup_unstruct() が
		// それを添字に使って領域外書き込みになる)
		int32_t nd[ELEM_NMAX];
		for (int l = 0; l < nn; l++) {
			const long g = v[off + l];
			if ((g < 0) || (g > maxid) || (idmap[g] < 0)) {
				printf("*** mesh : element %d refers to an unknown node %ld\n", e + 1, g);
				return 1;
			}
			nd[l] = idmap[g];
		}

		if (elem_store(type, tag, nd)) return 1;
	}

	return elem_finish();
}


// ---- Gmsh ASCII 4.1 の読み込み ----
//
// 2.2 との違いは 2 つ:
//   ・節点も要素も「エンティティのブロック」に分かれる (ブロック内は同じ型)
//   ・**要素の行に物理タグが無い**。物理タグはエンティティに付いているので、
//     $Entities を読んで (次元, エンティティ番号) -> 物理タグ の表を作る
// この 2 つ目が要点で、$Entities を無視すると全要素の物理タグが 0 になり、
// 材料も電極も割り当たらないまま「電極に節点が無い」で落ちる。

#define MAXENT (4096)
#define ENT_PMAX (8)

typedef struct {
	int dim;
	long tag;
	int phys;			// 物理タグ (無ければ 0)
	int np;				// 物理タグの数 (ENT_PMAX まで) と全部の並び
	int pt[ENT_PMAX];
} entity_t;

static entity_t Ent[MAXENT];
static int NEnt;

static int ent_phys(int dim, long tag)
{
	for (int i = 0; i < NEnt; i++) {
		if ((Ent[i].dim == dim) && (Ent[i].tag == tag)) return Ent[i].phys;
	}

	return 0;
}


static void ent_add(int dim, long tag, int np, const int *pt)
{
	if (NEnt >= MAXENT) return;
	Ent[NEnt].dim = dim;
	Ent[NEnt].tag = tag;
	Ent[NEnt].phys = ((np > 0) ? pt[0] : 0);
	Ent[NEnt].np = ((np < ENT_PMAX) ? np : ENT_PMAX);
	for (int i = 0; i < Ent[NEnt].np; i++) Ent[NEnt].pt[i] = pt[i];
	NEnt++;
}


/*
エンティティから物理タグを引いて要素を格納する (3.0 / 4.1)。**面 (電極・awall) が複数の物理グループに属するときは
グループごとに 1 回ずつ格納する** — 2.2 は同じ面をグループごとに別の要素として書くので、それと同じになる。
最初のタグだけにすると、電極の面を awall にも入れた格子で awall が黙って落ちる (gmsh 製の 3 次の六面体の棒の 100 kHz で
R +27%、実測)。体積要素は最初のタグだけ (2 回格納すると剛性が 2 重になる)
*/
static int elem_store_ent(int type, int dim, long tag, const int32_t *nd)
{
	if (dim == 2) {
		for (int i = 0; i < NEnt; i++) {
			if ((Ent[i].dim != dim) || (Ent[i].tag != tag) || (Ent[i].np < 2)) continue;
			for (int k = 0; k < Ent[i].np; k++) {
				if (elem_store(type, Ent[i].pt[k], nd)) return 1;
			}
			return 0;
		}
	}

	return elem_store(type, ent_phys(dim, tag), nd);
}


static int read_entities_v41(FILE *fp)
{
	long n[4];

	NEnt = 0;
	for (int d = 0; d < 4; d++) {
		if (fscanf(fp, "%ld", &n[d]) != 1) return 1;
	}
	for (int d = 0; d < 4; d++) {
		for (long e = 0; e < n[d]; e++) {
			long tag = 0;
			double b[6];
			if (fscanf(fp, "%ld", &tag) != 1) return 1;
			// 点は座標 3 個、それ以外は外接直方体 6 個
			const int nb = ((d == 0) ? 3 : 6);
			for (int i = 0; i < nb; i++) {
				if (fscanf(fp, "%lf", &b[i]) != 1) return 1;
			}
			long np = 0;
			if (fscanf(fp, "%ld", &np) != 1) return 1;
			int pts[ENT_PMAX], k = 0;
			for (long i = 0; i < np; i++) {
				long pt = 0;
				if (fscanf(fp, "%ld", &pt) != 1) return 1;
				if (k < ENT_PMAX) pts[k++] = (int)pt;
			}
			// 境界エンティティの並び (点エンティティには無い)
			if (d > 0) {
				long nbd = 0;
				if (fscanf(fp, "%ld", &nbd) != 1) return 1;
				for (long i = 0; i < nbd; i++) {
					long dummy = 0;
					if (fscanf(fp, "%ld", &dummy) != 1) return 1;
				}
			}
			ent_add(d, tag, k, pts);
		}
	}

	return 0;
}


static int read_nodes_v41(FILE *fp, int32_t **idmap, int32_t *maxid)
{
	long nblk = 0, nn = 0, mn = 0, mx = 0;

	if (fscanf(fp, "%ld %ld %ld %ld", &nblk, &nn, &mn, &mx) != 4) return 1;
	if ((nn < 4) || (mx < 1)) {
		printf("*** mesh : too few nodes (%ld)\n", nn);
		return 1;
	}

	NNode = (int)nn;
	Xp = (double *)malloc((size_t)nn * sizeof(double));
	Yp = (double *)malloc((size_t)nn * sizeof(double));
	Zp = (double *)malloc((size_t)nn * sizeof(double));
	int32_t *map = (int32_t *)malloc(((size_t)mx + 1) * sizeof(int32_t));
	for (int32_t i = 0; i <= (int32_t)mx; i++) map[i] = -1;

	long *tags = (long *)malloc((size_t)nn * sizeof(long));
	int k = 0;
	for (long b = 0; b < nblk; b++) {
		long dim = 0, tag = 0, par = 0, cnt = 0;
		if (fscanf(fp, "%ld %ld %ld %ld", &dim, &tag, &par, &cnt) != 4) {
			free(map); free(tags); return 1;
		}
		// ブロック内は「節点番号がまとめて、そのあと座標がまとめて」
		if (k + cnt > nn) { free(map); free(tags); return 1; }
		for (long i = 0; i < cnt; i++) {
			if (fscanf(fp, "%ld", &tags[k + i]) != 1) { free(map); free(tags); return 1; }
		}
		for (long i = 0; i < cnt; i++) {
			if (fscanf(fp, "%lf %lf %lf", &Xp[k + i], &Yp[k + i], &Zp[k + i]) != 3) {
				free(map); free(tags); return 1;
			}
		}
		k += (int)cnt;
	}
	if (k != nn) { free(map); free(tags); return 1; }

	for (int i = 0; i < NNode; i++) {
		const long g = tags[i];
		if ((g < 1) || (g > mx)) { free(map); free(tags); return 1; }
		map[g] = i;
	}
	free(tags);

	*idmap = map;
	*maxid = (int32_t)mx;

	return 0;
}


static int read_elements_v41(FILE *fp, const int32_t *idmap, int32_t maxid)
{
	long nblk = 0, ne = 0, mn = 0, mx = 0;

	if (fscanf(fp, "%ld %ld %ld %ld", &nblk, &ne, &mn, &mx) != 4) return 1;

	elem_reset();

	for (long b = 0; b < nblk; b++) {
		long dim = 0, tag = 0, type = 0, cnt = 0;
		if (fscanf(fp, "%ld %ld %ld %ld", &dim, &tag, &type, &cnt) != 4) return 1;

		const int nn = ((type == 4) ? 4 : (type == 11) ? 10 : (type == 29) ? 20 : (type == 30) ? 35 : (type == 31) ? 56 : (type == 71) ? 84 : (type == 72) ? 120
		              : (type == 2) ? 3 : (type == 9) ? 6 : (type == 21) ? 10 : (type == 23) ? 15 : (type == 25) ? 21
		              : (type == 42) ? 28 : (type == 43) ? 36
		              : (type == 5) ? 8 : (type == 3) ? 4 : (type == 6) ? 6
		              : (type == 7) ? 5
		              : (type == 12) ? 27 : (type == 17) ? 20
		              : (type == 13) ? 18 : (type == 18) ? 15
		              : (type == 10) ? 9 : (type == 16) ? 8
		              : (type == 92) ? 64 : (type == 90) ? 40 : (type == 36) ? 16
		              : (type == 93) ? 125 : (type == 91) ? 75 : (type == 37) ? 25
		              : (type == 94) ? 216 : (type == 106) ? 126 : (type == 38) ? 36 : 0);
		// 物理タグはエンティティ側にある (2.2 と違い要素の行には無い)

		for (long e = 0; e < cnt; e++) {
			long etag = 0;
			if (fscanf(fp, "%ld", &etag) != 1) return 1;
			if (nn == 0) {
				// 使わない要素型 (点・線分など)。節点数が分からないので
				// 行末まで読み飛ばす
				int c;
				while (((c = fgetc(fp)) != EOF) && (c != '\n')) ;
				continue;
			}
			int32_t nd[ELEM_NMAX];
			for (int l = 0; l < nn; l++) {
				long g = 0;
				if (fscanf(fp, "%ld", &g) != 1) return 1;
				if ((g < 0) || (g > maxid) || (idmap[g] < 0)) {
					printf("*** mesh : element %ld refers to an unknown node %ld\n", etag, g);
					return 1;
				}
				nd[l] = idmap[g];
			}
			if (elem_store_ent((int)type, (int)dim, tag, nd)) return 1;
		}
	}

	if ((NTet < 1) && (NTri < 1) && (NHex < 1) && (NPrism < 1) && (NPyr < 1)) {
		printf("%s\n", "*** mesh : no tetrahedron, hexahedron, prism, pyramid "
			"or triangle found");
		return 1;
	}

	return elem_finish();
}


// ---- Gmsh バイナリの読み込み (2.2 / 4.1) ----
/*
ASCII との違いは「データの並べ方」だけなので、格納 (elem_store / elem_finish)
と物理タグの解決 (ent_phys) は共通のものをそのまま使う。

**エンディアンは変換しない。** $MeshFormat の直後に書かれている int32 が 1 と
読めるかどうかで判定し、違えば「別エンディアンの格子は非対応」と言って落とす。
黙って読み違えるより落ちる方がよく、手元で検証できない変換を書いても
「動くつもりのコード」が増えるだけになる。

  2.2 : $Nodes / $Elements の**個数だけが ASCII 行**で、以降はバイナリ。
        節点は (int32 tag, double x, y, z) の並び。要素は**型ごとのブロック**で、
        ブロック頭が (型, 個数, タグ数) の int32 3 個。ASCII と違って
        1 要素ごとに型が書かれていない。
  4.1 : 個数も含めて全部バイナリ。**size_t (ヘッダ 3 番目 = 8) と int32 が
        混在する**ので、どちらかを取り違えると即座にずれる。

読み飛ばしのために全要素型の節点数表を持つ (Gmsh のマニュアルの表)。
表に無い型が出たら、何バイト進めばよいか分からないので**落とす**。
ASCII では 1 要素 1 行なので行末まで飛ばせば済むが、バイナリでは
「知らない型は読み飛ばせない」。
*/

// Gmsh の要素型 -> 節点数 (1..31)。0 は「未知」
static const int GMSH_NNODE[32] = {
	0,  2,  3,  4,  4,  8,  6,  5,  3,  6,  9,		//  0..10
	10, 27, 18, 14,  1,  8, 20, 15, 13,  9,			// 11..20
	10, 12, 15, 15, 21,  4,  5,  6, 20, 35,			// 21..30
	56												// 31
};

static int gmsh_nnode(long type)
{
	// 3 次の六面体・角柱・四角形 (表の外の番号)
	if (type == 92) return 64;
	if (type == 90) return 40;
	if (type == 36) return 16;
	if (type == 93) return 125;
	if (type == 91) return 75;
	if (type == 37) return 25;
	// 5 次の六面体・角柱・四角形
	if (type == 94) return 216;
	if (type == 106) return 126;
	if (type == 38) return 36;
	// 6 次 / 7 次の四面体・三角形
	if (type == 71) return 84;
	if (type == 72) return 120;
	if (type == 42) return 28;
	if (type == 43) return 36;
	if ((type < 1) || (type > 31)) return 0;

	return GMSH_NNODE[type];
}


static int rd_raw(FILE *fp, void *p, size_t n)
{
	return (fread(p, 1, n, fp) != n);
}

static int rd_i32(FILE *fp, int32_t *v)
{
	return rd_raw(fp, v, 4);
}

static int rd_i64(FILE *fp, int64_t *v)
{
	return rd_raw(fp, v, 8);
}


// 節点 (2.2 バイナリ)
static int read_nodes_bin22(FILE *fp, int32_t **idmap, int32_t *maxid)
{
	char line[MESH_LINE];

	if (fgets(line, sizeof(line), fp) == NULL) return 1;
	const int nn = atoi(line);
	if (nn < 4) {
		printf("*** mesh : too few nodes (%d)\n", nn);
		return 1;
	}

	NNode = nn;
	Xp = (double *)malloc((size_t)nn * sizeof(double));
	Yp = (double *)malloc((size_t)nn * sizeof(double));
	Zp = (double *)malloc((size_t)nn * sizeof(double));

	int32_t *id = (int32_t *)malloc((size_t)nn * sizeof(int32_t));
	int32_t mx = 0;
	for (int i = 0; i < nn; i++) {
		double xyz[3];
		if (rd_i32(fp, &id[i]) || rd_raw(fp, xyz, 3 * sizeof(double))) {
			printf("*** mesh : truncated binary node %d\n", i + 1);
			free(id);
			return 1;
		}
		if (id[i] > mx) mx = id[i];
		Xp[i] = xyz[0];
		Yp[i] = xyz[1];
		Zp[i] = xyz[2];
	}
	if (mx < 1) {
		free(id);
		return 1;
	}

	int32_t *map = (int32_t *)malloc(((size_t)mx + 1) * sizeof(int32_t));
	for (int32_t i = 0; i <= mx; i++) map[i] = -1;
	for (int i = 0; i < nn; i++) map[id[i]] = i;
	free(id);

	*idmap = map;
	*maxid = mx;

	return 0;
}


// 要素 (2.2 バイナリ)。型ごとのブロックで並ぶ
static int read_elements_bin22(FILE *fp, const int32_t *idmap, int32_t maxid)
{
	char line[MESH_LINE];

	if (fgets(line, sizeof(line), fp) == NULL) return 1;
	const long ne = atol(line);

	elem_reset();

	long got = 0;
	while (got < ne) {
		int32_t hdr[3];
		if (rd_i32(fp, &hdr[0]) || rd_i32(fp, &hdr[1]) || rd_i32(fp, &hdr[2])) {
			printf("%s\n", "*** mesh : truncated binary element block");
			return 1;
		}
		const long type = hdr[0], cnt = hdr[1], ntag = hdr[2];
		const int nn = gmsh_nnode(type);
		if ((nn < 1) || (cnt < 0) || (ntag < 0)) {
			printf("*** mesh : unknown element type %ld in a binary file "
				"(cannot skip it)\n", type);
			return 1;
		}
		if ((got + cnt) > ne) {
			printf("%s\n", "*** mesh : binary element blocks overrun the count");
			return 1;
		}
		const int use = ((type == 4) || (type == 11) || (type == 2) || (type == 9)
		              || (type == 29) || (type == 21) || (type == 30) || (type == 23) || (type == 31) || (type == 25)
		              || (type == 71) || (type == 42) || (type == 72) || (type == 43)
		              || (type == 5) || (type == 3) || (type == 6) || (type == 7)
		              || (type == 12) || (type == 17) || (type == 13) || (type == 18)
		              || (type == 10) || (type == 16) || (type == 92) || (type == 90) || (type == 36)
		              || (type == 93) || (type == 91) || (type == 37) || (type == 94) || (type == 106) || (type == 38));

		for (long e = 0; e < cnt; e++) {
			int32_t etag = 0;
			if (rd_i32(fp, &etag)) return 1;
			int tag = 0;
			for (long t = 0; t < ntag; t++) {
				int32_t v = 0;
				if (rd_i32(fp, &v)) return 1;
				if (t == 0) tag = (int)v;		// 1 つ目が物理タグ
			}
			int32_t nd[ELEM_NMAX];
			for (int l = 0; l < nn; l++) {
				int32_t g = 0;
				if (rd_i32(fp, &g)) return 1;
				if (!use) continue;				// 使わない型は読み捨てる
				if ((g < 0) || (g > maxid) || (idmap[g] < 0)) {
					printf("*** mesh : element %d refers to an unknown node %d\n",
						(int)etag, (int)g);
					return 1;
				}
				nd[l] = idmap[g];
			}
			if (use && elem_store((int)type, tag, nd)) return 1;
		}
		got += cnt;
	}

	return elem_finish();
}


// エンティティ (4.1 バイナリ)。物理タグはここにしか無い
static int read_entities_bin41(FILE *fp)
{
	int64_t n[4];

	NEnt = 0;
	for (int d = 0; d < 4; d++) {
		if (rd_i64(fp, &n[d])) return 1;
	}
	for (int d = 0; d < 4; d++) {
		for (int64_t e = 0; e < n[d]; e++) {
			int32_t tag = 0;
			double b[6];
			// 点は座標 3 個、それ以外は外接直方体 6 個
			const int nb = ((d == 0) ? 3 : 6);
			if (rd_i32(fp, &tag) || rd_raw(fp, b, (size_t)nb * sizeof(double))) return 1;
			int64_t np = 0;
			if (rd_i64(fp, &np)) return 1;
			int pts[ENT_PMAX], k = 0;
			for (int64_t i = 0; i < np; i++) {
				int32_t pt = 0;
				if (rd_i32(fp, &pt)) return 1;
				if (k < ENT_PMAX) pts[k++] = (int)pt;
			}
			if (d > 0) {
				int64_t nbd = 0;
				if (rd_i64(fp, &nbd)) return 1;
				for (int64_t i = 0; i < nbd; i++) {
					int32_t dummy = 0;
					if (rd_i32(fp, &dummy)) return 1;
				}
			}
			ent_add(d, tag, k, pts);
		}
	}

	return 0;
}


// 節点 (4.1 バイナリ)
static int read_nodes_bin41(FILE *fp, int32_t **idmap, int32_t *maxid)
{
	int64_t nblk = 0, nn = 0, mn = 0, mx = 0;

	if (rd_i64(fp, &nblk) || rd_i64(fp, &nn) || rd_i64(fp, &mn) || rd_i64(fp, &mx)) return 1;
	if ((nn < 4) || (mx < 1) || (mx > INT32_MAX)) {
		printf("*** mesh : too few nodes (%lld)\n", (long long)nn);
		return 1;
	}

	NNode = (int)nn;
	Xp = (double *)malloc((size_t)nn * sizeof(double));
	Yp = (double *)malloc((size_t)nn * sizeof(double));
	Zp = (double *)malloc((size_t)nn * sizeof(double));
	int32_t *map = (int32_t *)malloc(((size_t)mx + 1) * sizeof(int32_t));
	for (int32_t i = 0; i <= (int32_t)mx; i++) map[i] = -1;

	int64_t *tags = (int64_t *)malloc((size_t)nn * sizeof(int64_t));
	int k = 0;
	for (int64_t b = 0; b < nblk; b++) {
		int32_t dim = 0, tag = 0, par = 0;
		int64_t cnt = 0;
		if (rd_i32(fp, &dim) || rd_i32(fp, &tag) || rd_i32(fp, &par) || rd_i64(fp, &cnt)
		 || (cnt < 0) || ((k + cnt) > nn)) {
			free(map); free(tags);
			return 1;
		}
		// ブロック内は「節点番号がまとめて、そのあと座標がまとめて」
		if (rd_raw(fp, &tags[k], (size_t)cnt * sizeof(int64_t))) {
			free(map); free(tags);
			return 1;
		}
		for (int64_t i = 0; i < cnt; i++) {
			double xyz[3];
			if (rd_raw(fp, xyz, 3 * sizeof(double))) {
				free(map); free(tags);
				return 1;
			}
			Xp[k + i] = xyz[0];
			Yp[k + i] = xyz[1];
			Zp[k + i] = xyz[2];
		}
		k += (int)cnt;
	}
	if (k != nn) { free(map); free(tags); return 1; }

	for (int i = 0; i < NNode; i++) {
		const int64_t g = tags[i];
		if ((g < 1) || (g > mx)) { free(map); free(tags); return 1; }
		map[g] = i;
	}
	free(tags);

	*idmap = map;
	*maxid = (int32_t)mx;

	return 0;
}


// 要素 (4.1 バイナリ)
static int read_elements_bin41(FILE *fp, const int32_t *idmap, int32_t maxid)
{
	int64_t nblk = 0, ne = 0, mn = 0, mx = 0;

	if (rd_i64(fp, &nblk) || rd_i64(fp, &ne) || rd_i64(fp, &mn) || rd_i64(fp, &mx)) return 1;

	elem_reset();

	for (int64_t b = 0; b < nblk; b++) {
		int32_t dim = 0, tag = 0, type = 0;
		int64_t cnt = 0;
		if (rd_i32(fp, &dim) || rd_i32(fp, &tag) || rd_i32(fp, &type) || rd_i64(fp, &cnt)
		 || (cnt < 0)) {
			printf("%s\n", "*** mesh : truncated binary element block");
			return 1;
		}
		const int nn = gmsh_nnode(type);
		if (nn < 1) {
			printf("*** mesh : unknown element type %d in a binary file "
				"(cannot skip it)\n", (int)type);
			return 1;
		}
		// 物理タグはエンティティ側にある (2.2 と違い要素には無い)
		const int use = ((type == 4) || (type == 11) || (type == 2) || (type == 9)
		              || (type == 29) || (type == 21) || (type == 30) || (type == 23) || (type == 31) || (type == 25)
		              || (type == 71) || (type == 42) || (type == 72) || (type == 43)
		              || (type == 5) || (type == 3) || (type == 6) || (type == 7)
		              || (type == 12) || (type == 17) || (type == 13) || (type == 18)
		              || (type == 10) || (type == 16) || (type == 92) || (type == 90) || (type == 36)
		              || (type == 93) || (type == 91) || (type == 37) || (type == 94) || (type == 106) || (type == 38));

		for (int64_t e = 0; e < cnt; e++) {
			int64_t etag = 0;
			if (rd_i64(fp, &etag)) return 1;
			int32_t nd[ELEM_NMAX];
			for (int l = 0; l < nn; l++) {
				int64_t g = 0;
				if (rd_i64(fp, &g)) return 1;
				if (!use) continue;
				if ((g < 0) || (g > maxid) || (idmap[g] < 0)) {
					printf("*** mesh : element %lld refers to an unknown node %lld\n",
						(long long)etag, (long long)g);
					return 1;
				}
				nd[l] = idmap[g];
			}
			if (use && elem_store_ent((int)type, (int)dim, tag, nd)) return 1;
		}
	}

	if ((NTet < 1) && (NTri < 1) && (NHex < 1) && (NPrism < 1) && (NPyr < 1)) {
		printf("%s\n", "*** mesh : no tetrahedron, hexahedron, prism, pyramid "
			"or triangle found");
		return 1;
	}

	return elem_finish();
}


/*
バイナリ格子の読み込み本体。

**セクションの探索は行単位だが、バイナリ領域を行として読み飛ばしてはいけない**
(データに 0x0A が現れる)。各セクションの読み手がちょうどの長さを消費するので、
戻ってきた位置から次のセクション見出しを探す形にしてある。
*/
// ---- Gmsh 3.0 (ASCII / バイナリ) の読み込み ----
/*
gmsh 4.15.2 に `Mesh.MshFileVersion = 3` で書かせたファイルを分類して決めた並び (2.2 と 4.1 の間の形):
  $Entities : 4 つの個数のあと、点は「番号 物理タグ数 物理タグ…」、曲線・面・体積は
              「番号 境界の数 境界… 物理タグ数 物理タグ…」(**4.1 と違い外接直方体が無く、境界が物理タグより前**)。
              バイナリのファイルでもこの節は ASCII
  $Nodes    : 個数のあと 1 行 1 節点「番号 x y z エンティティ」(バイナリは int32, double x 3, int32 の 28 バイト)
  $Elements : 個数のあと 1 行 1 要素「番号 型 エンティティ 節点数 節点…」(バイナリは同じ並びの int32)。
              **物理タグは要素の行に無い** (4.1 と同じくエンティティから引く)
節点数を行が持つので、知らない型 (点・線分) はバイナリでも読み飛ばせる。
*/

// 使う要素型 (3 次元要素と境界面)。2.2 のバイナリの読み手と同じ集合
static int gmsh_type_used(long type)
{
	return ((type == 4) || (type == 11) || (type == 2) || (type == 9)
	     || (type == 29) || (type == 21) || (type == 30) || (type == 23) || (type == 31) || (type == 25)
	     || (type == 71) || (type == 42) || (type == 72) || (type == 43)
	     || (type == 5) || (type == 3) || (type == 6) || (type == 7)
	     || (type == 12) || (type == 17) || (type == 13) || (type == 18)
	     || (type == 10) || (type == 16) || (type == 92) || (type == 90) || (type == 36)
	     || (type == 93) || (type == 91) || (type == 37) || (type == 94) || (type == 106) || (type == 38));
}

// 使う要素型の次元 (物理タグをエンティティから引くのに要る)
static int gmsh_type_dim(long type)
{
	return (((type == 2) || (type == 9) || (type == 21) || (type == 23) || (type == 25) || (type == 42) || (type == 43)
	      || (type == 3) || (type == 10) || (type == 16) || (type == 36) || (type == 37) || (type == 38)) ? 2 : 3);
}


static int read_entities_v3(FILE *fp)
{
	long n[4];

	NEnt = 0;
	for (int d = 0; d < 4; d++) {
		if (fscanf(fp, "%ld", &n[d]) != 1) return 1;
	}
	for (int d = 0; d < 4; d++) {
		for (long e = 0; e < n[d]; e++) {
			long tag = 0;
			if (fscanf(fp, "%ld", &tag) != 1) return 1;
			// 境界エンティティの並び (点エンティティには無い)。物理タグより前
			if (d > 0) {
				long nbd = 0;
				if (fscanf(fp, "%ld", &nbd) != 1) return 1;
				for (long i = 0; i < nbd; i++) {
					long dummy = 0;
					if (fscanf(fp, "%ld", &dummy) != 1) return 1;
				}
			}
			long np = 0;
			if (fscanf(fp, "%ld", &np) != 1) return 1;
			int pts[ENT_PMAX], k = 0;
			for (long i = 0; i < np; i++) {
				long pt = 0;
				if (fscanf(fp, "%ld", &pt) != 1) return 1;
				if (k < ENT_PMAX) pts[k++] = (int)pt;
			}
			ent_add(d, tag, k, pts);
		}
	}

	return 0;
}


// 節点の番号表を作る (2.2 / 3.0 共通の後始末)
static int nodes_idmap(int32_t *id, int nn, int32_t **idmap, int32_t *maxid)
{
	int32_t mx = 0;
	for (int i = 0; i < nn; i++) {
		if (id[i] < 1) return 1;
		if (id[i] > mx) mx = id[i];
	}
	int32_t *map = (int32_t *)malloc(((size_t)mx + 1) * sizeof(int32_t));
	for (int32_t i = 0; i <= mx; i++) map[i] = -1;
	for (int i = 0; i < nn; i++) map[id[i]] = i;
	*idmap = map;
	*maxid = mx;

	return 0;
}


// 節点 (3.0)。bin = 1 ならバイナリ
static int read_nodes_v3(FILE *fp, int bin, int32_t **idmap, int32_t *maxid)
{
	char line[MESH_LINE];

	if (fgets(line, sizeof(line), fp) == NULL) return 1;
	const long nn = atol(line);
	if ((nn < 4) || (nn >= INT32_MAX)) {
		printf("*** mesh : too few nodes (%ld)\n", nn);
		return 1;
	}

	NNode = (int)nn;
	Xp = (double *)malloc((size_t)nn * sizeof(double));
	Yp = (double *)malloc((size_t)nn * sizeof(double));
	Zp = (double *)malloc((size_t)nn * sizeof(double));
	int32_t *id = (int32_t *)malloc((size_t)nn * sizeof(int32_t));
	for (long i = 0; i < nn; i++) {
		if (bin) {
			double xyz[3];
			int32_t ent = 0;
			if (rd_i32(fp, &id[i]) || rd_raw(fp, xyz, 3 * sizeof(double)) || rd_i32(fp, &ent)) {
				printf("*** mesh : truncated binary node %ld\n", i + 1);
				free(id);
				return 1;
			}
			Xp[i] = xyz[0];
			Yp[i] = xyz[1];
			Zp[i] = xyz[2];
		}
		else {
			long g = 0;
			if ((fgets(line, sizeof(line), fp) == NULL)
			 || (sscanf(line, "%ld %lf %lf %lf", &g, &Xp[i], &Yp[i], &Zp[i]) != 4) || (g < 1) || (g >= INT32_MAX)) {
				printf("*** mesh : cannot parse node %ld\n", i + 1);
				free(id);
				return 1;
			}
			id[i] = (int32_t)g;
		}
	}
	const int ierr = nodes_idmap(id, NNode, idmap, maxid);
	free(id);

	return ierr;
}


// バイナリの節の直後が終わりの行であること (並びの読み違いを黙って通さない)
static int bin_section_end(FILE *fp, const char *tag)
{
	int c;
	while (((c = fgetc(fp)) == '\n') || (c == '\r')) ;
	if (c == EOF) return 1;
	ungetc(c, fp);
	char line[MESH_LINE];
	if ((fgets(line, sizeof(line), fp) == NULL) || strncmp(line, tag, strlen(tag))) {
		printf("*** mesh : unexpected binary layout before %s (only gmsh's MSH 3.0 without parametric "
			"coordinates is supported)\n", tag);
		return 1;
	}

	return 0;
}


// 要素 (3.0)。bin = 1 ならバイナリ
static int read_elements_v3(FILE *fp, int bin, const int32_t *idmap, int32_t maxid)
{
	long ne = 0;

	if (bin) {
		char line[MESH_LINE];
		if (fgets(line, sizeof(line), fp) == NULL) return 1;
		ne = atol(line);
	}
	else if (fscanf(fp, "%ld", &ne) != 1) return 1;

	elem_reset();

	for (long e = 0; e < ne; e++) {
		long hdr[4];
		for (int k = 0; k < 4; k++) {
			if (bin) {
				int32_t v = 0;
				if (rd_i32(fp, &v)) {
					printf("%s\n", "*** mesh : truncated binary element");
					return 1;
				}
				hdr[k] = v;
			}
			else if (fscanf(fp, "%ld", &hdr[k]) != 1) return 1;
		}
		const long etag = hdr[0], type = hdr[1], ent = hdr[2], nv = hdr[3];
		const int use = gmsh_type_used(type);
		if ((nv < 1) || (nv > ELEM_NMAX) || (use && (nv != gmsh_nnode(type)))) {
			printf("*** mesh : element %ld of type %ld has %ld nodes\n", etag, type, nv);
			return 1;
		}
		int32_t nd[ELEM_NMAX];
		for (long l = 0; l < nv; l++) {
			long g = 0;
			if (bin) {
				int32_t v = 0;
				if (rd_i32(fp, &v)) return 1;
				g = v;
			}
			else if (fscanf(fp, "%ld", &g) != 1) return 1;
			if (!use) continue;				// 使わない型 (点・線分) は読み捨てる
			if ((g < 0) || (g > maxid) || (idmap[g] < 0)) {
				printf("*** mesh : element %ld refers to an unknown node %ld\n", etag, g);
				return 1;
			}
			nd[l] = idmap[g];
		}
		// 物理タグはエンティティ側にある
		if (use && elem_store_ent((int)type, gmsh_type_dim(type), ent, nd)) return 1;
	}

	if ((NTet < 1) && (NTri < 1) && (NHex < 1) && (NPrism < 1) && (NPyr < 1)) {
		printf("%s\n", "*** mesh : no tetrahedron, hexahedron, prism, pyramid "
			"or triangle found");
		return 1;
	}

	return elem_finish();
}


static int mesh_read_binary(const char *fname, int v41, int v3)
{
	char line[MESH_LINE];
	int32_t *idmap = NULL;
	int32_t maxid = 0;
	int ierr = 0;
	int have_nodes = 0, have_elements = 0;

	NEnt = 0;

	FILE *fp = fopen(fname, "rb");
	if (fp == NULL) {
		printf("*** mesh file %s open error.\n", fname);
		return 1;
	}

	while (fgets(line, sizeof(line), fp) != NULL) {
		if      (!strncmp(line, "$MeshFormat", 11)) {
			if (fgets(line, sizeof(line), fp) == NULL) { ierr = 1; break; }
			double ver = 0;
			long bin = 0, sz = 0;
			if (sscanf(line, "%lf %ld %ld", &ver, &bin, &sz) != 3) { ierr = 1; break; }
			if (sz != 8) {
				printf("*** mesh : the binary file uses %ld-byte floats "
					"(only 8 is supported)\n", sz);
				ierr = 1;
				break;
			}
			// 直後の int32 が 1 なら同じエンディアン
			int32_t one = 0;
			if (rd_i32(fp, &one)) { ierr = 1; break; }
			if (one != 1) {
				// マーカーの 4 バイトが全部印字可能なら、ヘッダだけ
				// バイナリを名乗って中身が ASCII のファイル (手で書き換えた
				// もの) である可能性が高い。原因の違う 2 つを同じ文言で
				// 報告すると、直しようのない診断になる
				const unsigned char *b = (const unsigned char *)&one;
				int text = 1;
				for (int i = 0; i < 4; i++) {
					if ((b[i] < 0x20) || (b[i] > 0x7e)) text = 0;
				}
				if (text) {
					printf("%s\n", "*** mesh : the header says binary but the "
						"content is ASCII (the second field of $MeshFormat "
						"must be 0 for ASCII files)");
				}
				else {
					printf("*** mesh : the binary file has the opposite byte order "
						"(endianness marker = %d); convert it with gmsh on this "
						"machine\n", (int)one);
				}
				ierr = 1;
				break;
			}
			// マーカーの後ろの改行を読み捨てる
			int c = fgetc(fp);
			if (c != '\n') ungetc(c, fp);
		}
		else if (!strncmp(line, "$Entities", 9)) {
			// 3.0 のエンティティはバイナリのファイルでも ASCII
			ierr = (v3 ? read_entities_v3(fp) : v41 ? read_entities_bin41(fp) : 0);
			if (ierr) {
				printf("%s\n", "*** mesh : cannot parse the binary $Entities");
				break;
			}
		}
		else if (!strncmp(line, "$Nodes", 6)) {
			ierr = (v3 ? (read_nodes_v3(fp, 1, &idmap, &maxid) || bin_section_end(fp, "$EndNodes"))
			      : v41 ? read_nodes_bin41(fp, &idmap, &maxid)
			            : read_nodes_bin22(fp, &idmap, &maxid));
			if (ierr) break;
			have_nodes = 1;
		}
		else if (!strncmp(line, "$Elements", 9)) {
			if (!have_nodes) {
				printf("%s\n", "*** mesh : $Elements before $Nodes");
				ierr = 1;
				break;
			}
			ierr = (v3 ? (read_elements_v3(fp, 1, idmap, maxid) || bin_section_end(fp, "$EndElements"))
			      : v41 ? read_elements_bin41(fp, idmap, maxid)
			            : read_elements_bin22(fp, idmap, maxid));
			if (ierr) break;
			have_elements = 1;
		}
	}

	fclose(fp);
	free(idmap);

	if (!ierr && (!have_nodes || !have_elements)) {
		printf("%s\n", "*** mesh : $Nodes or $Elements is missing");
		ierr = 1;
	}

	return ierr;
}


int mesh_read(const char *fname)
{
	char line[MESH_LINE];
	int32_t *idmap = NULL;
	int32_t maxid = 0;
	int ierr = 0;
	int have_nodes = 0, have_elements = 0;
	int v41 = 0, v3 = 0;

	NEnt = 0;

	FILE *fp = fopen(fname, "r");
	if (fp == NULL) {
		printf("*** mesh file %s open error.\n", fname);
		return 1;
	}

	while (fgets(line, sizeof(line), fp) != NULL) {
		if      (!strncmp(line, "$MeshFormat", 11)) {
			if (fgets(line, sizeof(line), fp) == NULL) break;
			if      (line[0] == '2') v41 = 0;
			else if (line[0] == '3') v3 = 1;
			else if (line[0] == '4') v41 = 1;
			else {
				printf("*** mesh : only Gmsh format 2.x / 3.0 / 4.x is supported (got %s)", line);
				ierr = 1;
				break;
			}
			// バイナリは 2 番目の数字が 1 (ASCII は 0)。**別の読み手に渡す**。
			// ASCII 用のこの経路はテキストモードで開いた FILE* を使っており、
			// Windows では改行が変換されてバイナリが壊れるので、
			// 開き直して最初から読む (ASCII の経路には一切手を入れない)
			{
				double ver = 0;
				long bin = 0, sz = 0;
				if (sscanf(line, "%lf %ld %ld", &ver, &bin, &sz) >= 2) {
					if (bin != 0) {
						fclose(fp);
						free(idmap);
						return mesh_read_binary(fname, v41, v3);
					}
				}
			}
		}
		else if (!strncmp(line, "$Entities", 9)) {
			ierr = (v3 ? read_entities_v3(fp) : read_entities_v41(fp));
			if (ierr) {
				printf("%s\n", "*** mesh : cannot parse $Entities");
				break;
			}
		}
		else if (!strncmp(line, "$Nodes", 6)) {
			ierr = (v3 ? read_nodes_v3(fp, 0, &idmap, &maxid)
			      : v41 ? read_nodes_v41(fp, &idmap, &maxid)
			            : read_nodes(fp, &idmap, &maxid));
			if (ierr) break;
			have_nodes = 1;
		}
		else if (!strncmp(line, "$Elements", 9)) {
			if (!have_nodes) {
				printf("%s\n", "*** mesh : $Elements before $Nodes");
				ierr = 1;
				break;
			}
			ierr = (v3 ? read_elements_v3(fp, 0, idmap, maxid)
			      : v41 ? read_elements_v41(fp, idmap, maxid)
			            : read_elements(fp, idmap, maxid));
			if (ierr) break;
			have_elements = 1;
		}
	}

	fclose(fp);
	free(idmap);

	if (!ierr && (!have_nodes || !have_elements)) {
		printf("%s\n", "*** mesh : $Nodes or $Elements is missing");
		ierr = 1;
	}

	return ierr;
}


// ---- 一般 CRS の構築 ----

static int cmp_int32(const void *a, const void *b)
{
	const int32_t x = *(const int32_t *)a;
	const int32_t y = *(const int32_t *)b;

	return ((x < y) ? -1 : ((x > y) ? 1 : 0));
}


// 要素 e の局所節点をまとめて取り出す。戻り値は節点数 (1 次 4、2 次 10、3 次 20、4 次 35)。
// 並びは Gmsh の tet10 / tet20 / tet35 と同じ (頂点 4 個のあと辺・面・内部の節点)
int tet_nodes(int e, int32_t nd[TET_NMAX])
{
	for (int l = 0; l < 4; l++) nd[l] = Tet[(e * 4) + l];
	if (TetOrder < 2) return 4;
	if (TetOrder >= 3) {
		const int nho = tet_nen() - 4;
		for (int l = 0; l < nho; l++) nd[4 + l] = Tet3[(e * TET_HO_MAX) + l];
		return 4 + nho;
	}
	for (int l = 0; l < 6; l++) nd[4 + l] = Tet2[(e * 6) + l];

	return 10;
}


// 四面体あたりの節点数 (行列のパターンの見積りに使う)
int tet_nen(void)
{
	return ((TetOrder >= 2) ? (((TetOrder + 1) * (TetOrder + 2) * (TetOrder + 3)) / 6) : 4);
}


/*
要素の連結から節点の隣接関係を作る (対角成分を含む)。

四面体でも六面体でも手順は同じ (要素あたりの節点数と取り出し方だけが違う)
ので、`get` に要素の節点を書き出す関数を渡して共通化する。
**四面体の経路は結果も並び順もこれまでと同一** (既存の回帰がバイト単位で
一致することで確かめている)。
*/
static void crs_alloc_conn(crs_t *A, int nelem, int nenmax, int (*get)(int, int32_t *))
{
	const int n = NNode;

	// 節点毎の要素数を数える
	int *cnt = (int *)malloc((size_t)n * sizeof(int));
	memset(cnt, 0, (size_t)n * sizeof(int));
	for (int e = 0; e < nelem; e++) {
		int32_t nd[ELEM_NMAX];
		const int nen = get(e, nd);
		for (int l = 0; l < nen; l++) cnt[nd[l]]++;
	}
	int64_t *nptr = (int64_t *)malloc(((size_t)n + 1) * sizeof(int64_t));
	nptr[0] = 0;
	for (int i = 0; i < n; i++) nptr[i + 1] = nptr[i] + cnt[i];
	int32_t *nlist = (int32_t *)malloc((size_t)nptr[n] * sizeof(int32_t));
	memset(cnt, 0, (size_t)n * sizeof(int));
	for (int e = 0; e < nelem; e++) {
		int32_t nd[ELEM_NMAX];
		const int nen = get(e, nd);
		for (int l = 0; l < nen; l++) {
			const int32_t i = nd[l];
			nlist[nptr[i] + cnt[i]] = (int32_t)e;
			cnt[i]++;
		}
	}

	// 行毎に隣接節点を集めて整列・重複除去する
	A->n = n;
	A->rowptr = (int64_t *)malloc(((size_t)n + 1) * sizeof(int64_t));

	int cap = 64;
	int32_t *work = (int32_t *)malloc((size_t)cap * sizeof(int32_t));
	int *rown = (int *)malloc((size_t)n * sizeof(int));

	for (int i = 0; i < n; i++) {
		const int64_t p0 = nptr[i], p1 = nptr[i + 1];
		const int need = (int)(p1 - p0) * nenmax;
		if (need > cap) {
			cap = need;
			work = (int32_t *)realloc(work, (size_t)cap * sizeof(int32_t));
		}
		int m = 0;
		for (int64_t p = p0; p < p1; p++) {
			int32_t nd[ELEM_NMAX];
			const int nen = get(nlist[p], nd);
			for (int l = 0; l < nen; l++) work[m++] = nd[l];
		}
		qsort(work, (size_t)m, sizeof(int32_t), cmp_int32);
		int u = 0;
		for (int q = 0; q < m; q++) {
			if ((q == 0) || (work[q] != work[q - 1])) u++;
		}
		rown[i] = u;
	}

	A->rowptr[0] = 0;
	for (int i = 0; i < n; i++) A->rowptr[i + 1] = A->rowptr[i] + rown[i];
	A->nnz = A->rowptr[n];
	A->col = (int32_t *)malloc((size_t)A->nnz * sizeof(int32_t));
	A->val = (double *)malloc((size_t)A->nnz * sizeof(double));

	for (int i = 0; i < n; i++) {
		const int64_t p0 = nptr[i], p1 = nptr[i + 1];
		int m = 0;
		for (int64_t p = p0; p < p1; p++) {
			int32_t nd[ELEM_NMAX];
			const int nen = get(nlist[p], nd);
			for (int l = 0; l < nen; l++) work[m++] = nd[l];
		}
		qsort(work, (size_t)m, sizeof(int32_t), cmp_int32);
		int64_t w = A->rowptr[i];
		for (int q = 0; q < m; q++) {
			if ((q == 0) || (work[q] != work[q - 1])) A->col[w++] = work[q];
		}
	}

	free(cnt);
	free(nptr);
	free(nlist);
	free(work);
	free(rown);

	crs_zero(A);
}


/*
3 次元要素の統一した見方 (種別の混在を許すため)。

番号は **「四面体 -> 六面体 -> 角柱 -> ピラミッド」の連番**。この順序は
変えないこと: 材料・場の出力・HDF5 がこの順に並ぶ。純粋な四面体格子では
この並びが従来と完全に一致するので、既存の答えがビット単位で変わらない。
*/
int elem3d_count(void)
{
	return (NTet + NHex + NPrism + NPyr);
}


int elem3d_kind(int e)
{
	if (e < NTet) return MESHELEM_TET;
	if (e < (NTet + NHex)) return MESHELEM_HEX;
	if (e < (NTet + NHex + NPrism)) return MESHELEM_PRISM;

	return MESHELEM_PYR;
}


// 六面体 e の局所節点をまとめて取り出す。戻り値は節点数 (8 / 20 / 27)。
// 並びは Gmsh の hex8 / hex20 / hex27 と同じ (頂点 8 個のあと辺・面・体心)
int hex_nodes(int e, int32_t nd[ELEM_NMAX])
{
	for (int l = 0; l < 8; l++) nd[l] = Hex[(e * 8) + l];
	for (int l = 8; l < HexNen; l++) nd[l] = Hex2[(e * HEX2_W) + (l - 8)];

	return HexNen;
}


// 角柱 e の局所節点をまとめて取り出す。戻り値は節点数 (6 / 15 / 18)
int prism_nodes(int e, int32_t nd[PRISMP_NMAX])
{
	for (int l = 0; l < 6; l++) nd[l] = Prism[(e * 6) + l];
	for (int l = 6; l < PrismNen; l++) nd[l] = Prism2[(e * PRISM2_W) + (l - 6)];

	return PrismNen;
}


// 要素 e の節点。戻り値は節点数
// (四面体 4 or 10、六面体 8/20/27、角柱 6/15/18、ピラミッド 5)
int elem3d_nodes(int e, int32_t nd[ELEM_NMAX])
{
	if (e < NTet) return tet_nodes(e, nd);
	if (e < (NTet + NHex)) return hex_nodes(e - NTet, nd);
	if (e < (NTet + NHex + NPrism)) return prism_nodes(e - NTet - NHex, nd);
	{
		const int i = e - NTet - NHex - NPrism;
		for (int l = 0; l < 5; l++) nd[l] = Pyr[(i * 5) + l];
	}

	return 5;
}


int elem3d_mat(int e)
{
	if (e < NTet) return ((TetMat != NULL) ? TetMat[e] : 0);
	if (e < (NTet + NHex)) return ((HexMat != NULL) ? HexMat[e - NTet] : 0);
	if (e < (NTet + NHex + NPrism)) {
		return ((PrismMat != NULL) ? PrismMat[e - NTet - NHex] : 0);
	}

	return ((PyrMat != NULL) ? PyrMat[e - NTet - NHex - NPrism] : 0);
}


int elem3d_tag(int e)
{
	if (e < NTet) return TetTag[e];
	if (e < (NTet + NHex)) return HexTag[e - NTet];
	if (e < (NTet + NHex + NPrism)) return PrismTag[e - NTet - NHex];

	return PyrTag[e - NTet - NHex - NPrism];
}


static int get_tet(int e, int32_t *nd)
{
	return tet_nodes(e, nd);
}


void crs_alloc_tet(crs_t *A)
{
	crs_alloc_conn(A, NTet, tet_nen(), get_tet);
}


// 3 次元の非構造格子 (種別の混在を含む)
void crs_alloc_elem3d(crs_t *A)
{
	// work バッファの見積りに使う要素あたりの最大節点数 (大きい分には無害)
	int nenmax = ((TetOrder >= 2) ? tet_nen() : 8);
	if (HexNen > nenmax) nenmax = HexNen;
	if (PrismNen > nenmax) nenmax = PrismNen;
	crs_alloc_conn(A, elem3d_count(), nenmax, elem3d_nodes);
}


// 行 row の中で列 col の位置を二分探索で求める (列は昇順)
static int64_t crs_find(const crs_t *A, int32_t row, int32_t col)
{
	int64_t lo = A->rowptr[row];
	int64_t hi = A->rowptr[row + 1] - 1;

	while (lo <= hi) {
		const int64_t mid = (lo + hi) / 2;
		if      (A->col[mid] < col) lo = mid + 1;
		else if (A->col[mid] > col) hi = mid - 1;
		else return mid;
	}

	return -1;
}


// ---- 四面体の要素行列 ----

// 1 次四面体の形状関数勾配 (要素内で一定) と体積
// 戻り値 : 0 = 正常、1 = 退化要素
int tet_grad_pub(const int32_t nd[4], double g[4][3], double *vol)
{
	const double x0 = Xp[nd[0]], y0 = Yp[nd[0]], z0 = Zp[nd[0]];
	const double j00 = Xp[nd[1]] - x0, j01 = Xp[nd[2]] - x0, j02 = Xp[nd[3]] - x0;
	const double j10 = Yp[nd[1]] - y0, j11 = Yp[nd[2]] - y0, j12 = Yp[nd[3]] - y0;
	const double j20 = Zp[nd[1]] - z0, j21 = Zp[nd[2]] - z0, j22 = Zp[nd[3]] - z0;

	const double det = (j00 * ((j11 * j22) - (j12 * j21)))
	                 - (j01 * ((j10 * j22) - (j12 * j20)))
	                 + (j02 * ((j10 * j21) - (j11 * j20)));
	if (fabs(det) <= 0) return 1;

	// ∇ξ_i は J^-1 の第 i 行
	const double a[3][3] = {
		{ ((j11 * j22) - (j12 * j21)) / det,
		 -((j01 * j22) - (j02 * j21)) / det,
		  ((j01 * j12) - (j02 * j11)) / det},
		{-((j10 * j22) - (j12 * j20)) / det,
		  ((j00 * j22) - (j02 * j20)) / det,
		 -((j00 * j12) - (j02 * j10)) / det},
		{ ((j10 * j21) - (j11 * j20)) / det,
		 -((j00 * j21) - (j01 * j20)) / det,
		  ((j00 * j11) - (j01 * j10)) / det}
	};

	for (int d = 0; d < 3; d++) {
		g[1][d] = a[0][d];
		g[2][d] = a[1][d];
		g[3][d] = a[2][d];
		g[0][d] = -(a[0][d] + a[1][d] + a[2][d]);
	}
	*vol = fabs(det) / 6;

	return 0;
}


// ---- 2 次四面体 (10 節点、等パラメトリック) ----
//
// 積分則は Duffy 変換 + 各方向 3 点 Gauss-Legendre (27 点)。
//
//   λ1 = u, λ2 = v(1-u), λ3 = w(1-u)(1-v),  λ0 = 1 - λ1 - λ2 - λ3
//   dλ1 dλ2 dλ3 = (1-u)^2 (1-v) du dv dw
//
// これを選んだ理由は「重みが全部正で、正しさを手で追える」から。四面体の
// 少点数則 (Keast 等) は次数が上がると負の重みが出るうえ、係数を暗記に頼ると
// 検算できない。Duffy は 1 次元 Gauss の積で書けるので導出が閉じている。
// 被積分関数 ∇N・∇N は λ について 2 次なので (1-u) の因子を含めても u で 4 次、
// 3 点則 (5 次まで厳密) で直線要素なら厳密。曲がった要素では有理式になるので
// 厳密ではないが、27 点あれば形状誤差より十分小さい。
#define NQTET (27)

// 積分点 (バリセントリック λ[4]) と重み (∫ f dλ1dλ2dλ3 の重み、合計 1/6)
static void tet_quad(double lam[NQTET][4], double wq[NQTET])
{
	// [0,1] の 3 点 Gauss-Legendre
	const double gp[3] = {0.5 - (0.5 * 0.77459666924148337704),
	                      0.5,
	                      0.5 + (0.5 * 0.77459666924148337704)};
	const double gw[3] = {5.0 / 18, 8.0 / 18, 5.0 / 18};

	int q = 0;
	for (int a = 0; a < 3; a++) {
	for (int b = 0; b < 3; b++) {
	for (int c = 0; c < 3; c++) {
		const double u = gp[a], v = gp[b], w = gp[c];
		lam[q][1] = u;
		lam[q][2] = v * (1 - u);
		lam[q][3] = w * (1 - u) * (1 - v);
		lam[q][0] = 1 - lam[q][1] - lam[q][2] - lam[q][3];
		wq[q] = gw[a] * gw[b] * gw[c] * (1 - u) * (1 - u) * (1 - v);
		q++;
	}
	}
	}
}


// 10 節点四面体の形状関数の λ 微分 dl[i][a] = ∂N_i/∂λ_a
static void tet10_dlam(const double lam[4], double dl[10][4])
{
	// 辺 (Gmsh の tet10 の並び)
	static const int ed[6][2] = {{0, 1}, {1, 2}, {2, 0}, {3, 0}, {3, 2}, {3, 1}};

	memset(dl, 0, sizeof(double) * 10 * 4);
	for (int a = 0; a < 4; a++) {
		dl[a][a] = (4 * lam[a]) - 1;			// N_a = λ_a (2λ_a - 1)
	}
	for (int l = 0; l < 6; l++) {
		const int a = ed[l][0], b = ed[l][1];	// N = 4 λ_a λ_b
		dl[4 + l][a] = 4 * lam[b];
		dl[4 + l][b] = 4 * lam[a];
	}
}


// 積分点 1 点での物理座標勾配 ∇N_i とヤコビアン行列式。
// 戻り値 : 0 = 正常、1 = 退化 (det <= 0)
static int tet10_grad(const int32_t nd[10], const double lam[4],
	double gn[10][3], double *det)
{
	double dl[10][4];
	tet10_dlam(lam, dl);

	// 参照座標 ξ = (λ1, λ2, λ3) についての微分 (λ0 = 1 - λ1 - λ2 - λ3)
	double dx[10][3];
	for (int i = 0; i < 10; i++) {
		for (int k = 0; k < 3; k++) dx[i][k] = dl[i][k + 1] - dl[i][0];
	}

	// J[r][k] = Σ_i x_i[r] dN_i/dξ_k
	double j[3][3];
	for (int r = 0; r < 3; r++) {
		const double *p = ((r == 0) ? Xp : (r == 1) ? Yp : Zp);
		for (int k = 0; k < 3; k++) {
			double s = 0;
			for (int i = 0; i < 10; i++) s += p[nd[i]] * dx[i][k];
			j[r][k] = s;
		}
	}

	const double d = (j[0][0] * ((j[1][1] * j[2][2]) - (j[1][2] * j[2][1])))
	               - (j[0][1] * ((j[1][0] * j[2][2]) - (j[1][2] * j[2][0])))
	               + (j[0][2] * ((j[1][0] * j[2][1]) - (j[1][1] * j[2][0])));
	// 符号は節点の並び順で決まるので絶対値で潰さず、そのまま返して
	// 呼び出し側で「要素内で符号が一定か」を見る (曲がった要素の裏返り検出)
	if (d == 0) return 1;

	// Jinv[k][r] = ∂ξ_k/∂x_r
	const double ji[3][3] = {
		{ ((j[1][1] * j[2][2]) - (j[1][2] * j[2][1])) / d,
		 -((j[0][1] * j[2][2]) - (j[0][2] * j[2][1])) / d,
		  ((j[0][1] * j[1][2]) - (j[0][2] * j[1][1])) / d},
		{-((j[1][0] * j[2][2]) - (j[1][2] * j[2][0])) / d,
		  ((j[0][0] * j[2][2]) - (j[0][2] * j[2][0])) / d,
		 -((j[0][0] * j[1][2]) - (j[0][2] * j[1][0])) / d},
		{ ((j[1][0] * j[2][1]) - (j[1][1] * j[2][0])) / d,
		 -((j[0][0] * j[2][1]) - (j[0][1] * j[2][0])) / d,
		  ((j[0][0] * j[1][1]) - (j[0][1] * j[1][0])) / d}
	};

	for (int i = 0; i < 10; i++) {
		for (int r = 0; r < 3; r++) {
			gn[i][r] = (dx[i][0] * ji[0][r]) + (dx[i][1] * ji[1][r]) + (dx[i][2] * ji[2][r]);
		}
	}
	*det = d;

	return 0;
}


// 要素 e の 10 節点要素行列 ke[10][10] と体積。戻り値 : 0 = 正常、1 = 退化
int tet10_element(int e, const double c[6], double ke[10][10], double *vol)
{
	int32_t nd[TET_NMAX];
	if (tet_nodes(e, nd) != 10) return 1;

	double lam[NQTET][4], wq[NQTET];
	tet_quad(lam, wq);

	memset(ke, 0, sizeof(double) * 10 * 10);
	double v = 0, sgn = 0;
	for (int q = 0; q < NQTET; q++) {
		double gn[10][3], det;
		if (tet10_grad(nd, lam[q], gn, &det)) return 1;
		// 節点の並びで det の符号は変わる。要素内で一定なら向きが揃っている
		if (sgn == 0) sgn = ((det > 0) ? 1 : -1);
		else if (det * sgn <= 0) return 1;		// 曲がりすぎて裏返っている
		const double dw = wq[q] * fabs(det);
		v += dw;
		for (int l = 0; l < 10; l++) {
			for (int m = 0; m < 10; m++) {
				ke[l][m] += dw * ((c[0] * gn[l][0] * gn[m][0])
				                + (c[1] * gn[l][1] * gn[m][1])
				                + (c[2] * gn[l][2] * gn[m][2])
				                + (c[3] * ((gn[l][0] * gn[m][1]) + (gn[l][1] * gn[m][0])))
				                + (c[4] * ((gn[l][1] * gn[m][2]) + (gn[l][2] * gn[m][1])))
				                + (c[5] * ((gn[l][2] * gn[m][0]) + (gn[l][0] * gn[m][2]))));
			}
		}
	}
	*vol = v;

	return 0;
}


// ---- 3 次四面体 (20 節点、等パラメトリック) ----
//
// 局所節点の並びは **gmsh 4.15.2 の出力を分類した実測の表** (多重指数 3λ)。
// 暗記や文献の表は書き手と読み手が同じ誤解を共有し得るので使わない
// (2 次の六面体・角柱と同じ方針)。表があれば形状関数は次数によらない式
//   N_a(λ) = Π_v Π_{m=0}^{a_v - 1} (p λ_v - m) / (m + 1)       (Σ a_v = p)
// で書ける (頂点 ½λ(3λ-1)(3λ-2)、辺 (9/2)λ_iλ_j(3λ_i-1)、面 27λ_iλ_jλ_k に一致)。
// 積分は Duffy 変換 + 各方向 4 点 Gauss-Legendre (64 点)。∇N・∇N は λ の 4 次で、
// (1-u)^2 の因子を含めて u の 6 次 -> 4 点則 (7 次まで厳密) で直線要素なら厳密。
static const int TET3_A[20][4] = {
	{3, 0, 0, 0}, {0, 3, 0, 0}, {0, 0, 3, 0}, {0, 0, 0, 3},
	{2, 1, 0, 0}, {1, 2, 0, 0},		// 辺 (0,1)
	{0, 2, 1, 0}, {0, 1, 2, 0},		// 辺 (1,2)
	{1, 0, 2, 0}, {2, 0, 1, 0},		// 辺 (2,0)
	{1, 0, 0, 2}, {2, 0, 0, 1},		// 辺 (3,0)
	{0, 0, 1, 2}, {0, 0, 2, 1},		// 辺 (3,2)
	{0, 1, 0, 2}, {0, 2, 0, 1},		// 辺 (3,1)
	{1, 1, 1, 0}, {1, 1, 0, 1}, {1, 0, 1, 1}, {0, 1, 1, 1}	// 面
};
#define NQTET4 (64)

// 4 次 (35 節点、gmsh の tet35 = 型 30) の多重指数 4λ。**gmsh 4.15.2 の出力 (box_p4.msh の
// 108 要素がすべて同じ並び) を重心座標で分類した実測の表**。頂点、辺 (0,1) (1,2) (2,0) (3,0)
// (3,2) (3,1) に 3 個ずつ (始点の側から)、面 (0,1,2) (0,1,3) (0,2,3) (1,2,3) に 3 個ずつ、内部 1 個
static const int TET4_A[35][4] = {
	{4, 0, 0, 0}, {0, 4, 0, 0}, {0, 0, 4, 0}, {0, 0, 0, 4},
	{3, 1, 0, 0}, {2, 2, 0, 0}, {1, 3, 0, 0},
	{0, 3, 1, 0}, {0, 2, 2, 0}, {0, 1, 3, 0},
	{1, 0, 3, 0}, {2, 0, 2, 0}, {3, 0, 1, 0},
	{1, 0, 0, 3}, {2, 0, 0, 2}, {3, 0, 0, 1},
	{0, 0, 1, 3}, {0, 0, 2, 2}, {0, 0, 3, 1},
	{0, 1, 0, 3}, {0, 2, 0, 2}, {0, 3, 0, 1},
	{2, 1, 1, 0}, {1, 1, 2, 0}, {1, 2, 1, 0},
	{2, 1, 0, 1}, {1, 2, 0, 1}, {1, 1, 0, 2},
	{2, 0, 1, 1}, {1, 0, 1, 2}, {1, 0, 2, 1},
	{0, 1, 1, 2}, {0, 2, 1, 1}, {0, 1, 2, 1},
	{1, 1, 1, 1}
};
// 5 次 (56 節点、gmsh の tet56 = 型 31) の多重指数 5λ。**gmsh 4.15.2 の出力 (立方体の 1140 要素がすべて同じ並び)
// を重心座標で分類した実測の表** (同じ分類で TET4_A / TRI4_A が全項一致することを確認してから作った)。
// 頂点、辺 (0,1) (1,2) (2,0) (3,0) (3,2) (3,1) に 4 個ずつ、面 (0,1,2) (0,1,3) (0,2,3) (1,2,3) に 6 個ずつ、内部 4 個
static const int TET5_A[56][4] = {
	{5, 0, 0, 0}, {0, 5, 0, 0}, {0, 0, 5, 0}, {0, 0, 0, 5},
	{4, 1, 0, 0}, {3, 2, 0, 0}, {2, 3, 0, 0}, {1, 4, 0, 0},
	{0, 4, 1, 0}, {0, 3, 2, 0}, {0, 2, 3, 0}, {0, 1, 4, 0},
	{1, 0, 4, 0}, {2, 0, 3, 0}, {3, 0, 2, 0}, {4, 0, 1, 0},
	{1, 0, 0, 4}, {2, 0, 0, 3}, {3, 0, 0, 2}, {4, 0, 0, 1},
	{0, 0, 1, 4}, {0, 0, 2, 3}, {0, 0, 3, 2}, {0, 0, 4, 1},
	{0, 1, 0, 4}, {0, 2, 0, 3}, {0, 3, 0, 2}, {0, 4, 0, 1},
	{3, 1, 1, 0}, {1, 1, 3, 0}, {1, 3, 1, 0}, {2, 1, 2, 0},
	{1, 2, 2, 0}, {2, 2, 1, 0}, {3, 1, 0, 1}, {1, 3, 0, 1},
	{1, 1, 0, 3}, {2, 2, 0, 1}, {1, 2, 0, 2}, {2, 1, 0, 2},
	{3, 0, 1, 1}, {1, 0, 1, 3}, {1, 0, 3, 1}, {2, 0, 1, 2},
	{1, 0, 2, 2}, {2, 0, 2, 1}, {0, 1, 1, 3}, {0, 3, 1, 1},
	{0, 1, 3, 1}, {0, 2, 1, 2}, {0, 2, 2, 1}, {0, 1, 2, 2},
	{2, 1, 1, 1}, {1, 2, 1, 1}, {1, 1, 2, 1}, {1, 1, 1, 2}
};
#define NQTETMAX (729)		// 9 x 9 x 9 (7 次の自己検証)

// 2 次 (10 節点) の多重指数 2λ (Gmsh の tet10 : 辺 (0,1) (1,2) (2,0) (3,0) (3,2) (3,1))。
// 熱連成の積分点 (elem3d_ip) で次数によらない Lagrange の式を使うためのもの
// (電気の 2 次の要素行列は従来の tet10_grad のまま)
static const int TET2_A[10][4] = {
	{2, 0, 0, 0}, {0, 2, 0, 0}, {0, 0, 2, 0}, {0, 0, 0, 2},
	{1, 1, 0, 0}, {0, 1, 1, 0}, {1, 0, 1, 0}, {1, 0, 0, 1}, {0, 0, 1, 1}, {0, 1, 0, 1}
};

// 6 次以上 (6 次 84 節点 = gmsh の型 71、7 次 120 節点 = 型 72) の多重指数は gmsh の並びの**再帰の規則**で作る:
//   頂点 0..3、辺 (0,1) (1,2) (2,0) (3,0) (3,2) (3,1) に p−1 個ずつ (始点の側から)、
//   面 (0,2,1) (0,1,3) (0,3,2) (3,1,2) に (1,1,1) + 次数 p−3 の三角形 (同じ規則、面の頂点の順に写す)、
//   内部に (1,1,1,1) + 次数 p−4 の四面体 (同じ規則)。三角形は頂点、辺 (0,1) (1,2) (2,0)、内部 (1,1,1) + 次数 p−3。
// この規則は gmsh 4.15.2 の出力 (立方体の 1140 要素 / 540 面) を重心座標で分類したものと 2〜7 次の全項で一致し
// (data/sample/gmsh_ho_order.py)、実測の表 TET3_A / TET4_A / TET5_A / TRI3_A / TRI4_A / TRI5_A とも全項一致する
// (ho_tables_init がそれを毎回確かめる。一致しなければ入力エラー)
static int tri_rule(int p, int (*out)[3])
{
	if (p == 0) { out[0][0] = out[0][1] = out[0][2] = 0; return 1; }
	int n = 0;
	for (int v = 0; v < 3; v++) {
		for (int c = 0; c < 3; c++) out[n][c] = ((c == v) ? p : 0);
		n++;
	}
	static const int ed[3][2] = {{0, 1}, {1, 2}, {2, 0}};
	for (int e = 0; e < 3; e++) {
		for (int k = 1; k < p; k++) {
			for (int c = 0; c < 3; c++) out[n][c] = 0;
			out[n][ed[e][0]] = p - k;
			out[n][ed[e][1]] = k;
			n++;
		}
	}
	if (p >= 3) {
		const int m = tri_rule(p - 3, &out[n]);
		for (int s = 0; s < m; s++) {
			for (int c = 0; c < 3; c++) out[n + s][c] += 1;
		}
		n += m;
	}
	return n;
}

static int tet_rule(int p, int (*out)[4])
{
	if (p == 0) { out[0][0] = out[0][1] = out[0][2] = out[0][3] = 0; return 1; }
	int n = 0;
	for (int v = 0; v < 4; v++) {
		for (int c = 0; c < 4; c++) out[n][c] = ((c == v) ? p : 0);
		n++;
	}
	static const int ed[6][2] = {{0, 1}, {1, 2}, {2, 0}, {3, 0}, {3, 2}, {3, 1}};
	for (int e = 0; e < 6; e++) {
		for (int k = 1; k < p; k++) {
			for (int c = 0; c < 4; c++) out[n][c] = 0;
			out[n][ed[e][0]] = p - k;
			out[n][ed[e][1]] = k;
			n++;
		}
	}
	if (p >= 3) {
		static const int fc[4][3] = {{0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {3, 1, 2}};
		int sub[TRI_NMAX][3];
		const int m = tri_rule(p - 3, sub);
		for (int f = 0; f < 4; f++) {
			for (int s = 0; s < m; s++) {
				for (int c = 0; c < 4; c++) out[n][c] = 0;
				for (int c = 0; c < 3; c++) out[n][fc[f][c]] = sub[s][c] + 1;
				n++;
			}
		}
	}
	if (p >= 4) {
		const int m = tet_rule(p - 4, &out[n]);
		for (int s = 0; s < m; s++) {
			for (int c = 0; c < 4; c++) out[n + s][c] += 1;
		}
		n += m;
	}
	return n;
}

static int TETG_A[TET_NMAX][4], TRIG_A[TRI_NMAX][3], HoGenOrder = 0;

// 6 次以上の表を作り、規則が 3〜5 次の実測の表と全項一致することを確かめる (格子を読んだ直後、直列で呼ぶ)
static int ho_tables_init(int p)
{
	int t4[TET_NMAX][4];
	for (int q = 3; q <= 5; q++) {
		const int nt = tet_rule(q, t4);
		const int (*mt)[4] = ((q == 5) ? TET5_A : (q == 4) ? TET4_A : TET3_A);
		int bad = (nt != (((q + 1) * (q + 2) * (q + 3)) / 6));
		for (int i = 0; (i < nt) && !bad; i++) {
			for (int c = 0; c < 4; c++) bad |= (t4[i][c] != mt[i][c]);
		}
		if (!bad) bad = (tri_rule_check() != 0) * (q == 3);
		if (bad) {
			printf("*** internal error : the node-ordering rule disagrees with the measured order-%d tables "
				"(tetrahedra, or triangles of orders 3-5)\n", q);
			return 1;
		}
	}
	if ((p >= 6) && (HoGenOrder != p)) {
		tet_rule(p, TETG_A);
		tri_rule(p, TRIG_A);
		HoGenOrder = p;
	}
	return 0;
}

// 次数 TetOrder (2〜7) の多重指数の表
static const int *tet_ho_a(int i)
{
	return ((TetOrder >= 6) ? TETG_A[i] : (TetOrder == 5) ? TET5_A[i] : (TetOrder == 4) ? TET4_A[i]
	        : (TetOrder == 3) ? TET3_A[i] : TET2_A[i]);
}


// 公開 : 局所節点 i の多重指数 (重心座標 x 次数。3 次以上の Nedelec の eddy3d_pr.c が λ の引き戻しに使う)
const int *tet_ho_index(int i)
{
	return tet_ho_a(i);
}


// [0,1] の n 点 Gauss-Legendre (n = 3〜7 は表、8 以上 (GL_NMAX まで) は Legendre 多項式の根を Newton で求める。
// 3 点は 3 次の剛性には足りない : 変異解析用)。**表に無い n を黙って 5 点に落とさないこと** (以前は 8 点を
// 頼むと 5 点が返った)
static void gauss01(int n, double *x, double *w)
{
	if (n >= 8) {
		for (int i = 0; i < n; i++) {
			double t = cos(PI * (i + 0.75) / (n + 0.5)), dp = 1;
			for (int it = 0; it < 100; it++) {
				double p0 = 1, p1 = t;
				for (int k = 2; k <= n; k++) {
					const double p2 = (((2 * k - 1) * t * p1) - ((k - 1) * p0)) / k;
					p0 = p1; p1 = p2;
				}
				dp = n * ((t * p1) - p0) / ((t * t) - 1);
				const double dt = p1 / dp;
				t -= dt;
				if (fabs(dt) < 1e-16) break;
			}
			x[i] = 0.5 * (1 - t);		// t は降順なので x は昇順
			w[i] = 1 / ((1 - (t * t)) * dp * dp);
		}
		return;
	}
	if (n == 3) {
		x[0] = 0.5 - (0.5 * 0.77459666924148337704);
		x[1] = 0.5;
		x[2] = 0.5 + (0.5 * 0.77459666924148337704);
		w[0] = w[2] = 5.0 / 18;
		w[1] = 8.0 / 18;
		return;
	}
	if (n == 4) {
		const double t[2] = {0.33998104358485626480, 0.86113631159405257522};
		const double a[2] = {0.65214515486254614263, 0.34785484513745385737};
		for (int k = 0; k < 2; k++) {
			x[k] = 0.5 - (0.5 * t[1 - k]);
			w[k] = 0.5 * a[1 - k];
			x[3 - k] = 0.5 + (0.5 * t[1 - k]);
			w[3 - k] = 0.5 * a[1 - k];
		}
		return;
	}
	if (n == 7) {
		const double t[4] = {0, 0.40584515137739716691, 0.74153118559939443986, 0.94910791234275852453};
		const double a[4] = {0.41795918367346938776, 0.38183005050511894495,
		                     0.27970539148927666790, 0.12948496616886969327};
		x[3] = 0.5;
		w[3] = 0.5 * a[0];
		for (int k = 1; k < 4; k++) {
			x[3 - k] = 0.5 - (0.5 * t[k]);
			x[3 + k] = 0.5 + (0.5 * t[k]);
			w[3 - k] = w[3 + k] = 0.5 * a[k];
		}
		return;
	}
	if (n == 6) {
		const double t[3] = {0.23861918608319690863, 0.66120938646626451366, 0.93246951420315202781};
		const double a[3] = {0.46791393457269104739, 0.36076157304813860757, 0.17132449237917034504};
		for (int k = 0; k < 3; k++) {
			x[2 - k] = 0.5 - (0.5 * t[k]);
			x[3 + k] = 0.5 + (0.5 * t[k]);
			w[2 - k] = w[3 + k] = 0.5 * a[k];
		}
		return;
	}
	{
		const double t[3] = {0, 0.53846931010568309104, 0.90617984593866399280};
		const double a[3] = {0.56888888888888888889, 0.47862867049936646804,
		                     0.23692688505618908751};
		x[2] = 0.5;
		w[2] = 0.5 * a[0];
		for (int k = 1; k < 3; k++) {
			x[2 - k] = 0.5 - (0.5 * t[k]);
			x[2 + k] = 0.5 + (0.5 * t[k]);
			w[2 - k] = w[2 + k] = 0.5 * a[k];
		}
	}
}


// Duffy 変換 + 各方向 n 点 (n = 4 : 64 点、n = 5 : 125 点)。重みの合計は 1/6
static int tet_quad_n(int n, double lam[][4], double *wq)
{
	double gp[GL_NMAX], gw[GL_NMAX];
	gauss01(n, gp, gw);
	int q = 0;
	for (int a = 0; a < n; a++) {
	for (int b = 0; b < n; b++) {
	for (int c = 0; c < n; c++) {
		const double u = gp[a], v = gp[b], w = gp[c];
		lam[q][1] = u;
		lam[q][2] = v * (1 - u);
		lam[q][3] = w * (1 - u) * (1 - v);
		lam[q][0] = 1 - lam[q][1] - lam[q][2] - lam[q][3];
		wq[q] = gw[a] * gw[b] * gw[c] * (1 - u) * (1 - u) * (1 - v);
		q++;
	}
	}
	}
	return q;
}


// 次数 p の Lagrange 基底 (多重指数 a) の値と ∂N/∂λ_v
static double lagr_tet(int p, const int a[4], const double lam[4], double dn[4])
{
	double f[4], df[4];
	for (int v = 0; v < 4; v++) {
		double val = 1, der = 0;
		for (int m = 0; m < a[v]; m++) {
			const double t = ((p * lam[v]) - m) / (m + 1);
			der = (der * t) + (val * p / (m + 1));
			val *= t;
		}
		f[v] = val;
		df[v] = der;
	}
	for (int v = 0; v < 4; v++) {
		double d = df[v];
		for (int w = 0; w < 4; w++) {
			if (w != v) d *= f[w];
		}
		dn[v] = d;
	}
	return f[0] * f[1] * f[2] * f[3];
}




// 20 節点四面体の積分点 1 点での ∇N_i と det J (tet10_grad の 3 次版)。
// ヤコビアンは 20 節点すべてから作る (等パラメトリック)
static int tet20_grad(const int32_t *nd, const double lam[4], double gn[TET_NMAX][3], double *det)
{
	const int nen = tet_nen();
	double dx[TET_NMAX][3];
	for (int i = 0; i < nen; i++) {
		double dl[4];
		lagr_tet(TetOrder, tet_ho_a(i), lam, dl);
		for (int k = 0; k < 3; k++) dx[i][k] = dl[k + 1] - dl[0];
	}
	double j[3][3];
	for (int r = 0; r < 3; r++) {
		const double *p = ((r == 0) ? Xp : (r == 1) ? Yp : Zp);
		for (int k = 0; k < 3; k++) {
			double s = 0;
			for (int i = 0; i < nen; i++) s += p[nd[i]] * dx[i][k];
			j[r][k] = s;
		}
	}
	const double d = (j[0][0] * ((j[1][1] * j[2][2]) - (j[1][2] * j[2][1])))
	               - (j[0][1] * ((j[1][0] * j[2][2]) - (j[1][2] * j[2][0])))
	               + (j[0][2] * ((j[1][0] * j[2][1]) - (j[1][1] * j[2][0])));
	if (d == 0) return 1;
	const double ji[3][3] = {
		{ ((j[1][1] * j[2][2]) - (j[1][2] * j[2][1])) / d,
		 -((j[0][1] * j[2][2]) - (j[0][2] * j[2][1])) / d,
		  ((j[0][1] * j[1][2]) - (j[0][2] * j[1][1])) / d},
		{-((j[1][0] * j[2][2]) - (j[1][2] * j[2][0])) / d,
		  ((j[0][0] * j[2][2]) - (j[0][2] * j[2][0])) / d,
		 -((j[0][0] * j[1][2]) - (j[0][2] * j[1][0])) / d},
		{ ((j[1][0] * j[2][1]) - (j[1][1] * j[2][0])) / d,
		 -((j[0][0] * j[2][1]) - (j[0][1] * j[2][0])) / d,
		  ((j[0][0] * j[1][1]) - (j[0][1] * j[1][0])) / d}
	};
	for (int i = 0; i < nen; i++) {
		for (int r = 0; r < 3; r++) {
			gn[i][r] = (dx[i][0] * ji[0][r]) + (dx[i][1] * ji[1][r]) + (dx[i][2] * ji[2][r]);
		}
	}
	*det = d;

	return 0;
}


// 要素 e の 20 / 35 節点要素行列 ke と体積。戻り値 : 0 = 正常、1 = 退化・裏返り。
// 積分は Duffy + 各方向 4 点 (3 次 : ∇N・∇N は λ の 4 次、u の 6 次) / 5 点 (4 次 : 6 次、u の 8 次)
int tet20_element(int e, const double c[6], double ke[TET_NMAX][TET_NMAX], double *vol)
{
	int32_t nd[TET_NMAX];
	const int nen = tet_nodes(e, nd);
	if (nen < 20) return 1;

	double lam[NQTETMAX][4], wq[NQTETMAX];
	const int nq = tet_quad_n(TetOrder + 1, lam, wq);

	memset(ke, 0, sizeof(double) * TET_NMAX * TET_NMAX);
	double v = 0, sgn = 0;
	for (int q = 0; q < nq; q++) {
		double gn[TET_NMAX][3], det;
		if (tet20_grad(nd, lam[q], gn, &det)) return 1;
		if (sgn == 0) sgn = ((det > 0) ? 1 : -1);
		else if (det * sgn <= 0) return 1;
		const double dw = wq[q] * fabs(det);
		v += dw;
		for (int l = 0; l < nen; l++) {
			for (int m = 0; m < nen; m++) {
				ke[l][m] += dw * ((c[0] * gn[l][0] * gn[m][0])
				                + (c[1] * gn[l][1] * gn[m][1])
				                + (c[2] * gn[l][2] * gn[m][2])
				                + (c[3] * ((gn[l][0] * gn[m][1]) + (gn[l][1] * gn[m][0])))
				                + (c[4] * ((gn[l][1] * gn[m][2]) + (gn[l][2] * gn[m][1])))
				                + (c[5] * ((gn[l][2] * gn[m][0]) + (gn[l][0] * gn[m][2]))));
			}
		}
	}
	*vol = v;

	return 0;
}


// 要素の重心での ∇N_i (場の出力に使う。1 次では要素内一定なので同じ値)
int tet_grad_center(int e, double gn[TET_NMAX][3], int *nen)
{
	int32_t nd[TET_NMAX];
	const int n = tet_nodes(e, nd);

	*nen = n;
	if (n == 4) {
		double g[4][3], vol;
		if (tet_grad_pub(nd, g, &vol)) return 1;
		for (int i = 0; i < 4; i++) {
			for (int r = 0; r < 3; r++) gn[i][r] = g[i][r];
		}
		return 0;
	}

	const double lam[4] = {0.25, 0.25, 0.25, 0.25};
	double det;

	if (n >= 20) return tet20_grad(nd, lam, gn, &det);
	double g10[10][3];
	const int r = tet10_grad(nd, lam, g10, &det);
	for (int i = 0; i < 10; i++) {
		for (int k = 0; k < 3; k++) gn[i][k] = g10[i][k];
	}
	return r;
}


// 全体行列の作成 (非構造格子)
//   K_ij = ∫ (∇N_i)^T C (∇N_j) dV
// 1 次四面体では被積分関数が一定なので体積を掛けるだけで厳密。
// 2 次四面体は Gauss 積分 (tet10_element)。
/*
六面体 (8 節点、三重線形の等パラメトリック要素)。

局所節点の並びは **Gmsh / VTK と同じ**「下面を反時計回り、その上に上面」で、
局所座標 (ξ, η, ζ) ∈ [-1,1]^3 の符号は下の HEX_SGN のとおり。構造格子側の
局所番号 (i,j,k のビット) とは並びが違うので、表を共有してはいけない。

  N_a = (1 + ξ_a ξ)(1 + η_a η)(1 + ζ_a ζ) / 8

積分は各方向 2 点の Gauss-Legendre (計 8 点)。**直方体では厳密**で、
一般の六面体でも三重線形要素の標準の選択。ヤコビアンは 8 節点すべてから
作る (等パラメトリック) ので、平行六面体でない要素も正しく積分できる。

det の符号は節点の並び順で決まるので絶対値で潰さず、**要素内で符号が一定か
だけを見る** (裏返った要素の検出。四面体の 2 次要素と同じ規約)。
*/
static const signed char HEX_SGN[8][3] = {
	{-1, -1, -1}, {+1, -1, -1}, {+1, +1, -1}, {-1, +1, -1},
	{-1, -1, +1}, {+1, -1, +1}, {+1, +1, +1}, {-1, +1, +1}
};


/*
Gauss 点 q での物理座標の勾配 g[a][i] = ∂N_a/∂x_i と |det J| を求める。
q < 0 のときは要素中心 (ξ = η = ζ = 0) で評価する (場の出力で使う)。
戻り値 : 0 = 正常、1 = 退化 (det = 0)
*/
static int hex_shape(const int32_t *nd, int q, double g[8][3], double *det)
{
	const double gp = 1.0 / sqrt(3.0);
	const double xi  = ((q < 0) ? 0 : (gp * HEX_SGN[q][0]));
	const double eta = ((q < 0) ? 0 : (gp * HEX_SGN[q][1]));
	const double ze  = ((q < 0) ? 0 : (gp * HEX_SGN[q][2]));

	// 局所座標での微分
	double dn[8][3];
	for (int a = 0; a < 8; a++) {
		const double sx = HEX_SGN[a][0], sy = HEX_SGN[a][1], sz = HEX_SGN[a][2];
		dn[a][0] = sx * (1 + (sy * eta)) * (1 + (sz * ze)) / 8;
		dn[a][1] = sy * (1 + (sx * xi))  * (1 + (sz * ze)) / 8;
		dn[a][2] = sz * (1 + (sx * xi))  * (1 + (sy * eta)) / 8;
	}

	// ヤコビアン J[i][j] = ∂x_i/∂ξ_j = Σ_a x_a[i] ∂N_a/∂ξ_j
	double jm[3][3];
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) jm[i][j] = 0;
	}
	for (int a = 0; a < 8; a++) {
		const int32_t v = nd[a];
		const double p[3] = {Xp[v], Yp[v], Zp[v]};
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) jm[i][j] += p[i] * dn[a][j];
		}
	}
	const double d = (jm[0][0] * ((jm[1][1] * jm[2][2]) - (jm[1][2] * jm[2][1])))
	               - (jm[0][1] * ((jm[1][0] * jm[2][2]) - (jm[1][2] * jm[2][0])))
	               + (jm[0][2] * ((jm[1][0] * jm[2][1]) - (jm[1][1] * jm[2][0])));
	if (d == 0) return 1;

	// 逆行列 ji[j][i] = ∂ξ_j/∂x_i
	double ji[3][3];
	ji[0][0] = ((jm[1][1] * jm[2][2]) - (jm[1][2] * jm[2][1])) / d;
	ji[0][1] = ((jm[0][2] * jm[2][1]) - (jm[0][1] * jm[2][2])) / d;
	ji[0][2] = ((jm[0][1] * jm[1][2]) - (jm[0][2] * jm[1][1])) / d;
	ji[1][0] = ((jm[1][2] * jm[2][0]) - (jm[1][0] * jm[2][2])) / d;
	ji[1][1] = ((jm[0][0] * jm[2][2]) - (jm[0][2] * jm[2][0])) / d;
	ji[1][2] = ((jm[0][2] * jm[1][0]) - (jm[0][0] * jm[1][2])) / d;
	ji[2][0] = ((jm[1][0] * jm[2][1]) - (jm[1][1] * jm[2][0])) / d;
	ji[2][1] = ((jm[0][1] * jm[2][0]) - (jm[0][0] * jm[2][1])) / d;
	ji[2][2] = ((jm[0][0] * jm[1][1]) - (jm[0][1] * jm[1][0])) / d;

	// ∂N_a/∂x_i = Σ_j (∂N_a/∂ξ_j)(∂ξ_j/∂x_i)
	// ji は転置の形 (ji[j][i]) で持っているので、添字を取り違えないこと。
	// 直交格子では ji が対角になって取り違えても答えが同じになるため、
	// 検証には**軸に平行でない六面体**が要る
	for (int a = 0; a < 8; a++) {
		for (int i = 0; i < 3; i++) {
			g[a][i] = (dn[a][0] * ji[0][i]) + (dn[a][1] * ji[1][i]) + (dn[a][2] * ji[2][i]);
		}
	}
	*det = d;

	return 0;
}


// 3 点 Gauss-Legendre で積分した体積 (2 点則の検算に使う独立な計算)
static double hex_volume3(const int32_t *nd)
{
	const double gx[3] = {-0.7745966692414834, 0.0, 0.7745966692414834};
	const double gw[3] = {5.0 / 9, 8.0 / 9, 5.0 / 9};
	double v = 0;

	for (int i = 0; i < 3; i++) {
	for (int j = 0; j < 3; j++) {
	for (int k = 0; k < 3; k++) {
		double dn[8][3];
		for (int b = 0; b < 8; b++) {
			const double sx = HEX_SGN[b][0], sy = HEX_SGN[b][1], sz = HEX_SGN[b][2];
			dn[b][0] = sx * (1 + (sy * gx[j])) * (1 + (sz * gx[k])) / 8;
			dn[b][1] = sy * (1 + (sx * gx[i])) * (1 + (sz * gx[k])) / 8;
			dn[b][2] = sz * (1 + (sx * gx[i])) * (1 + (sy * gx[j])) / 8;
		}
		double jm[3][3];
		for (int p = 0; p < 3; p++) {
			for (int q = 0; q < 3; q++) jm[p][q] = 0;
		}
		for (int b = 0; b < 8; b++) {
			const int32_t v2 = nd[b];
			const double pt[3] = {Xp[v2], Yp[v2], Zp[v2]};
			for (int p = 0; p < 3; p++) {
				for (int q = 0; q < 3; q++) jm[p][q] += pt[p] * dn[b][q];
			}
		}
		const double d = (jm[0][0] * ((jm[1][1] * jm[2][2]) - (jm[1][2] * jm[2][1])))
		               - (jm[0][1] * ((jm[1][0] * jm[2][2]) - (jm[1][2] * jm[2][0])))
		               + (jm[0][2] * ((jm[1][0] * jm[2][1]) - (jm[1][1] * jm[2][0])));
		v += gw[i] * gw[j] * gw[k] * ((d > 0) ? d : -d);
	}
	}
	}

	return v;
}


// 六面体の要素行列。vol には体積 Σ w |det J| を返す
static int hex_element(const int32_t *nd, const double c[6], double ke[8][8], double *vol)
{
	double v = 0;
	int sgn = 0;

	for (int l = 0; l < 8; l++) {
		for (int m = 0; m < 8; m++) ke[l][m] = 0;
	}

	for (int q = 0; q < 8; q++) {
		double g[8][3], det;
		if (hex_shape(nd, q, g, &det)) return 1;
		const int sq = ((det > 0) ? 1 : -1);
		if (q == 0) sgn = sq;
		else if (sq != sgn) return 1;			// 要素内で符号が変わる = 裏返り
		const double w = ((det > 0) ? det : -det);	// 重みは 1
		v += w;
		for (int l = 0; l < 8; l++) {
			for (int m = 0; m < 8; m++) {
				ke[l][m] += w * ((c[0] * g[l][0] * g[m][0])
				               + (c[1] * g[l][1] * g[m][1])
				               + (c[2] * g[l][2] * g[m][2])
				               + (c[3] * ((g[l][0] * g[m][1]) + (g[l][1] * g[m][0])))
				               + (c[4] * ((g[l][1] * g[m][2]) + (g[l][2] * g[m][1])))
				               + (c[5] * ((g[l][2] * g[m][0]) + (g[l][0] * g[m][2]))));
			}
		}
	}
	*vol = v;

	return 0;
}


// 要素中心での勾配 (場の出力用)。戻り値 0 = 正常
int hex_grad_center(int e, double g[8][3])
{
	double det;

	return hex_shape(&Hex[e * 8], -1, g, &det);
}


/*
角柱 (6 節点、三角形を押し出した等パラメトリック要素)。

局所節点の並びは **Gmsh / VTK と同じ**「下面の三角形 (0,1,2)、その真上に
上面 (3,4,5)」。形状関数は三角形の面積座標と 1 次元の線形関数の積:

	N_a = λ_{i(a)}(u, v) (1 + s_a ζ) / 2,   s_a = -1 (下面) / +1 (上面)
	λ = (1 - u - v, u, v)

積分は**三角形 3 点則 × ζ 方向 2 点 Gauss** の計 6 点。
det J は (u,v) について 2 次、ζ について 2 次までなので、この組で体積は厳密。
六面体と同じく det の符号は絶対値で潰さず要素内の一定性だけを見る。

なぜ「六面体の 1 種」として書けないか: 角柱は三角形方向が面積座標 (λ の和が 1)
で、局所座標の取り方も積分則も違う。無理に同じ表にすると重みが合わない。
*/

// 局所節点 -> (三角形の頂点番号, ζ の符号)
static const signed char PRISM_TRI[6] = {0, 1, 2, 0, 1, 2};
static const signed char PRISM_SGN[6] = {-1, -1, -1, +1, +1, +1};

// 三角形 3 点則 (次数 2 まで厳密)。重みは面積 1/2 を 3 等分した 1/6
static const double PRISM_TU[3] = {1.0 / 6, 2.0 / 3, 1.0 / 6};
static const double PRISM_TV[3] = {1.0 / 6, 1.0 / 6, 2.0 / 3};


/*
局所座標 (u, v, ζ) での物理勾配 g[a][i] = ∂N_a/∂x_i と det J を求める。
q < 0 のときは要素中心 (u = v = 1/3, ζ = 0) で評価する。
*/
static int prism_shape_at(const int32_t *nd, double u, double v, double ze,
	double g[6][3], double *det)
{
	// λ とその (u, v) 微分
	const double lam[3] = {1 - u - v, u, v};
	static const double dlu[3] = {-1, 1, 0};
	static const double dlv[3] = {-1, 0, 1};

	double dn[6][3];
	for (int a = 0; a < 6; a++) {
		const int i = PRISM_TRI[a];
		const double sz = PRISM_SGN[a];
		const double h = (1 + (sz * ze)) / 2;
		dn[a][0] = dlu[i] * h;
		dn[a][1] = dlv[i] * h;
		dn[a][2] = lam[i] * sz / 2;
	}

	double jm[3][3];
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) jm[i][j] = 0;
	}
	for (int a = 0; a < 6; a++) {
		const int32_t w = nd[a];
		const double p[3] = {Xp[w], Yp[w], Zp[w]};
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) jm[i][j] += p[i] * dn[a][j];
		}
	}
	const double d = (jm[0][0] * ((jm[1][1] * jm[2][2]) - (jm[1][2] * jm[2][1])))
	               - (jm[0][1] * ((jm[1][0] * jm[2][2]) - (jm[1][2] * jm[2][0])))
	               + (jm[0][2] * ((jm[1][0] * jm[2][1]) - (jm[1][1] * jm[2][0])));
	if (d == 0) return 1;

	double ji[3][3];
	ji[0][0] = ((jm[1][1] * jm[2][2]) - (jm[1][2] * jm[2][1])) / d;
	ji[0][1] = ((jm[0][2] * jm[2][1]) - (jm[0][1] * jm[2][2])) / d;
	ji[0][2] = ((jm[0][1] * jm[1][2]) - (jm[0][2] * jm[1][1])) / d;
	ji[1][0] = ((jm[1][2] * jm[2][0]) - (jm[1][0] * jm[2][2])) / d;
	ji[1][1] = ((jm[0][0] * jm[2][2]) - (jm[0][2] * jm[2][0])) / d;
	ji[1][2] = ((jm[0][2] * jm[1][0]) - (jm[0][0] * jm[1][2])) / d;
	ji[2][0] = ((jm[1][0] * jm[2][1]) - (jm[1][1] * jm[2][0])) / d;
	ji[2][1] = ((jm[0][1] * jm[2][0]) - (jm[0][0] * jm[2][1])) / d;
	ji[2][2] = ((jm[0][0] * jm[1][1]) - (jm[0][1] * jm[1][0])) / d;

	for (int a = 0; a < 6; a++) {
		for (int i = 0; i < 3; i++) {
			g[a][i] = (dn[a][0] * ji[0][i]) + (dn[a][1] * ji[1][i]) + (dn[a][2] * ji[2][i]);
		}
	}
	*det = d;

	return 0;
}


// 積分点 q (0..5) での評価。q < 0 なら要素中心
static int prism_shape(const int32_t *nd, int q, double g[6][3], double *det)
{
	const double gp = 1.0 / sqrt(3.0);

	if (q < 0) return prism_shape_at(nd, 1.0 / 3, 1.0 / 3, 0, g, det);

	return prism_shape_at(nd, PRISM_TU[q % 3], PRISM_TV[q % 3],
		((q < 3) ? -gp : gp), g, det);
}


// 高次の積分則で求めた体積 (6 点則の検算に使う独立な計算)。
// 三角形 6 点則 (次数 4) x ζ 方向 3 点 Gauss
static double prism_volume_hi(const int32_t *nd)
{
	// Dunavant の 6 点則 (次数 4)。重みの和は 1 (面積 1/2 は最後に掛ける)
	static const double tu[6] = {0.44594849091597, 0.44594849091597, 0.10810301816807,
	                             0.09157621350977, 0.09157621350977, 0.81684757298046};
	static const double tv[6] = {0.44594849091597, 0.10810301816807, 0.44594849091597,
	                             0.09157621350977, 0.81684757298046, 0.09157621350977};
	static const double tw[6] = {0.22338158967801, 0.22338158967801, 0.22338158967801,
	                             0.10995174365532, 0.10995174365532, 0.10995174365532};
	static const double gx[3] = {-0.7745966692414834, 0.0, 0.7745966692414834};
	static const double gw[3] = {5.0 / 9, 8.0 / 9, 5.0 / 9};
	double v = 0;

	for (int t = 0; t < 6; t++) {
		for (int k = 0; k < 3; k++) {
			double g[6][3], det;
			if (prism_shape_at(nd, tu[t], tv[t], gx[k], g, &det)) continue;
			v += 0.5 * tw[t] * gw[k] * ((det > 0) ? det : -det);
		}
	}

	return v;
}


// 角柱の要素行列。vol には体積 Σ w |det J| を返す
static int prism_element(const int32_t *nd, const double c[6], double ke[6][6], double *vol)
{
	double v = 0;
	int sgn = 0;

	for (int l = 0; l < 6; l++) {
		for (int m = 0; m < 6; m++) ke[l][m] = 0;
	}

	for (int q = 0; q < 6; q++) {
		double g[6][3], det;
		if (prism_shape(nd, q, g, &det)) return 1;
		const int sq = ((det > 0) ? 1 : -1);
		if (q == 0) sgn = sq;
		else if (sq != sgn) return 1;
		// 重み : 三角形 (1/2 x 1/3) x ζ 方向 1
		const double w = (1.0 / 6) * ((det > 0) ? det : -det);
		v += w;
		for (int l = 0; l < 6; l++) {
			for (int m = 0; m < 6; m++) {
				ke[l][m] += w * ((c[0] * g[l][0] * g[m][0])
				               + (c[1] * g[l][1] * g[m][1])
				               + (c[2] * g[l][2] * g[m][2])
				               + (c[3] * ((g[l][0] * g[m][1]) + (g[l][1] * g[m][0])))
				               + (c[4] * ((g[l][1] * g[m][2]) + (g[l][2] * g[m][1])))
				               + (c[5] * ((g[l][2] * g[m][0]) + (g[l][0] * g[m][2]))));
			}
		}
	}
	*vol = v;

	return 0;
}


int prism_grad_center(int e, double g[6][3])
{
	double det;

	return prism_shape(&Prism[e * 6], -1, g, &det);
}


void assemble_prism(crs_t *A, int mode)
{
	crs_zero(A);

	for (int e = 0; e < NPrism; e++) {
		double c[6];
		material_coef_pub(PrismMat[e], mode, c);
		if ((c[0] <= 0) && (c[1] <= 0) && (c[2] <= 0)) continue;

		const int32_t *nd = &Prism[e * 6];
		double ke[6][6], vol;
		if (prism_element(nd, c, ke, &vol)) continue;

		for (int l = 0; l < 6; l++) {
			for (int m = 0; m < 6; m++) {
				const int64_t p = crs_find(A, nd[l], nd[m]);
				if (p >= 0) A->val[p] += ke[l][m];
			}
		}
	}
}


// ---- 2 次の六面体 (20 / 27 節点、等パラメトリック) ----
//
// 局所節点 -> (ξ, η, ζ) ∈ {-1, 0, +1}^3。並びは Gmsh の実測 (fem.h 参照):
// 頂点 8、辺 12、面 6 (z- y- x- x+ y+ z+)、体心。20 節点 (serendipity) は
// この表の先頭 20 個。
//
// 27 節点は 1 次元 2 次 Lagrange の 3 重積 N = l(ξ) l(η) l(ζ)。積分は各方向
// 3 点 Gauss-Legendre (27 点) で、**平行六面体では厳密** (被積分関数が各方向
// 4 次まで、3 点則は 5 次まで厳密)。20 節点は標準の serendipity 基底
// (角: (1+ξaξ)(1+ηaη)(1+ζaζ)(ξaξ+ηaη+ζaζ-2)/8、辺: (1-ξ²)(1+ηaη)(1+ζaζ)/4)。
// どちらも P2 (完全 2 次多項式) を含むので、直線要素なら 2 次の場を厳密に表す。
//
// ヤコビアンは 1 次と同じく**積分点毎に全節点から**作る (等パラメトリック)。
static const signed char HEX2_XI[27][3] = {
	{-1, -1, -1}, {+1, -1, -1}, {+1, +1, -1}, {-1, +1, -1},
	{-1, -1, +1}, {+1, -1, +1}, {+1, +1, +1}, {-1, +1, +1},
	// 辺 (0,1)(0,3)(0,4)(1,2)(1,5)(2,3)(2,6)(3,7)(4,5)(4,7)(5,6)(6,7)
	{0, -1, -1}, {-1, 0, -1}, {-1, -1, 0}, {+1, 0, -1}, {+1, -1, 0}, {0, +1, -1},
	{+1, +1, 0}, {-1, +1, 0}, {0, -1, +1}, {-1, 0, +1}, {+1, 0, +1}, {0, +1, +1},
	// 面 (z- y- x- x+ y+ z+) と体心
	{0, 0, -1}, {0, -1, 0}, {-1, 0, 0}, {+1, 0, 0}, {0, +1, 0}, {0, 0, +1},
	{0, 0, 0}
};

// 1 次元 2 次 Lagrange (節点 x = -1, 0, +1) とその微分
static double lag3(int s, double x)
{
	if (s < 0) return x * (x - 1) / 2;
	if (s > 0) return x * (x + 1) / 2;

	return 1 - (x * x);
}

static double lag3d(int s, double x)
{
	if (s < 0) return x - 0.5;
	if (s > 0) return x + 0.5;

	return -2 * x;
}


/*
3 次の六面体 (hexa64、Gmsh の型 92) と角柱 (prism40、型 90)。

局所節点の並びは **gmsh 4.15.2 の出力を分類した実測の表** (g3/cls.py で、要素の角の節点から
作った局所座標を 3 倍して丸めたもの)。六面体は各方向の格子位置 (0..3、ξ = -1 + 2i/3)、角柱は
三角形の重心座標の多重指数 (a, b, c) (a + b + c = 3、λ0 が角 0) と ζ の格子位置 (0..3)。
形状関数は次数によらない Lagrange の式 (1 次元の等間隔の Lagrange の積 / 三角形の多重指数の式)
なので、表さえ正しければ式を暗記に頼らない。
*/
static const signed char HEX3_I[64][3] = {
	{0, 0, 0}, {3, 0, 0}, {3, 3, 0}, {0, 3, 0}, {0, 0, 3}, {3, 0, 3}, {3, 3, 3}, {0, 3, 3},
	{1, 0, 0}, {2, 0, 0}, {0, 1, 0}, {0, 2, 0}, {0, 0, 1}, {0, 0, 2}, {3, 1, 0}, {3, 2, 0},
	{3, 0, 1}, {3, 0, 2}, {2, 3, 0}, {1, 3, 0}, {3, 3, 1}, {3, 3, 2}, {0, 3, 1}, {0, 3, 2},
	{1, 0, 3}, {2, 0, 3}, {0, 1, 3}, {0, 2, 3}, {3, 1, 3}, {3, 2, 3}, {2, 3, 3}, {1, 3, 3},
	{1, 1, 0}, {1, 2, 0}, {2, 2, 0}, {2, 1, 0}, {1, 0, 1}, {2, 0, 1}, {2, 0, 2}, {1, 0, 2},
	{0, 1, 1}, {0, 1, 2}, {0, 2, 2}, {0, 2, 1}, {3, 1, 1}, {3, 2, 1}, {3, 2, 2}, {3, 1, 2},
	{2, 3, 1}, {1, 3, 1}, {1, 3, 2}, {2, 3, 2}, {1, 1, 3}, {2, 1, 3}, {2, 2, 3}, {1, 2, 3},
	{1, 1, 1}, {2, 1, 1}, {2, 2, 1}, {1, 2, 1}, {1, 1, 2}, {2, 1, 2}, {2, 2, 2}, {1, 2, 2},
};
static const signed char PRISM3_A[40][4] = {
	{3, 0, 0, 0}, {0, 3, 0, 0}, {0, 0, 3, 0}, {3, 0, 0, 3}, {0, 3, 0, 3}, {0, 0, 3, 3},
	{2, 1, 0, 0}, {1, 2, 0, 0}, {2, 0, 1, 0}, {1, 0, 2, 0}, {3, 0, 0, 1}, {3, 0, 0, 2},
	{0, 2, 1, 0}, {0, 1, 2, 0}, {0, 3, 0, 1}, {0, 3, 0, 2}, {0, 0, 3, 1}, {0, 0, 3, 2},
	{2, 1, 0, 3}, {1, 2, 0, 3}, {2, 0, 1, 3}, {1, 0, 2, 3}, {0, 2, 1, 3}, {0, 1, 2, 3},
	{1, 1, 1, 0}, {1, 1, 1, 3}, {2, 1, 0, 1}, {1, 2, 0, 1}, {1, 2, 0, 2}, {2, 1, 0, 2},
	{2, 0, 1, 1}, {2, 0, 1, 2}, {1, 0, 2, 2}, {1, 0, 2, 1}, {0, 2, 1, 1}, {0, 1, 2, 1},
	{0, 1, 2, 2}, {0, 2, 1, 2}, {1, 1, 1, 1}, {1, 1, 1, 2},
};
// 4 次の六面体 (hexa125、型 93) / 角柱 (prism75、型 91)。同じく gmsh 4.15.2 の出力を分類した表
static const signed char HEX4_I[125][3] = {
	{0, 0, 0}, {4, 0, 0}, {4, 4, 0}, {0, 4, 0}, {0, 0, 4}, {4, 0, 4}, {4, 4, 4}, {0, 4, 4},
	{1, 0, 0}, {2, 0, 0}, {3, 0, 0}, {0, 1, 0}, {0, 2, 0}, {0, 3, 0}, {0, 0, 1}, {0, 0, 2},
	{0, 0, 3}, {4, 1, 0}, {4, 2, 0}, {4, 3, 0}, {4, 0, 1}, {4, 0, 2}, {4, 0, 3}, {3, 4, 0},
	{2, 4, 0}, {1, 4, 0}, {4, 4, 1}, {4, 4, 2}, {4, 4, 3}, {0, 4, 1}, {0, 4, 2}, {0, 4, 3},
	{1, 0, 4}, {2, 0, 4}, {3, 0, 4}, {0, 1, 4}, {0, 2, 4}, {0, 3, 4}, {4, 1, 4}, {4, 2, 4},
	{4, 3, 4}, {3, 4, 4}, {2, 4, 4}, {1, 4, 4}, {1, 1, 0}, {1, 3, 0}, {3, 3, 0}, {3, 1, 0},
	{1, 2, 0}, {2, 3, 0}, {3, 2, 0}, {2, 1, 0}, {2, 2, 0}, {1, 0, 1}, {3, 0, 1}, {3, 0, 3},
	{1, 0, 3}, {2, 0, 1}, {3, 0, 2}, {2, 0, 3}, {1, 0, 2}, {2, 0, 2}, {0, 1, 1}, {0, 1, 3},
	{0, 3, 3}, {0, 3, 1}, {0, 1, 2}, {0, 2, 3}, {0, 3, 2}, {0, 2, 1}, {0, 2, 2}, {4, 1, 1},
	{4, 3, 1}, {4, 3, 3}, {4, 1, 3}, {4, 2, 1}, {4, 3, 2}, {4, 2, 3}, {4, 1, 2}, {4, 2, 2},
	{3, 4, 1}, {1, 4, 1}, {1, 4, 3}, {3, 4, 3}, {2, 4, 1}, {1, 4, 2}, {2, 4, 3}, {3, 4, 2},
	{2, 4, 2}, {1, 1, 4}, {3, 1, 4}, {3, 3, 4}, {1, 3, 4}, {2, 1, 4}, {3, 2, 4}, {2, 3, 4},
	{1, 2, 4}, {2, 2, 4}, {1, 1, 1}, {3, 1, 1}, {3, 3, 1}, {1, 3, 1}, {1, 1, 3}, {3, 1, 3},
	{3, 3, 3}, {1, 3, 3}, {2, 1, 1}, {1, 2, 1}, {1, 1, 2}, {3, 2, 1}, {3, 1, 2}, {2, 3, 1},
	{3, 3, 2}, {1, 3, 2}, {2, 1, 3}, {1, 2, 3}, {3, 2, 3}, {2, 3, 3}, {2, 2, 1}, {2, 1, 2},
	{1, 2, 2}, {3, 2, 2}, {2, 3, 2}, {2, 2, 3}, {2, 2, 2},
};
static const signed char PRISM4_A[75][4] = {
	{4, 0, 0, 0}, {0, 4, 0, 0}, {0, 0, 4, 0}, {4, 0, 0, 4}, {0, 4, 0, 4}, {0, 0, 4, 4},
	{3, 1, 0, 0}, {2, 2, 0, 0}, {1, 3, 0, 0}, {3, 0, 1, 0}, {2, 0, 2, 0}, {1, 0, 3, 0},
	{4, 0, 0, 1}, {4, 0, 0, 2}, {4, 0, 0, 3}, {0, 3, 1, 0}, {0, 2, 2, 0}, {0, 1, 3, 0},
	{0, 4, 0, 1}, {0, 4, 0, 2}, {0, 4, 0, 3}, {0, 0, 4, 1}, {0, 0, 4, 2}, {0, 0, 4, 3},
	{3, 1, 0, 4}, {2, 2, 0, 4}, {1, 3, 0, 4}, {3, 0, 1, 4}, {2, 0, 2, 4}, {1, 0, 3, 4},
	{0, 3, 1, 4}, {0, 2, 2, 4}, {0, 1, 3, 4}, {2, 1, 1, 0}, {1, 1, 2, 0}, {1, 2, 1, 0},
	{2, 1, 1, 4}, {1, 2, 1, 4}, {1, 1, 2, 4}, {3, 1, 0, 1}, {1, 3, 0, 1}, {1, 3, 0, 3},
	{3, 1, 0, 3}, {2, 2, 0, 1}, {1, 3, 0, 2}, {2, 2, 0, 3}, {3, 1, 0, 2}, {2, 2, 0, 2},
	{3, 0, 1, 1}, {3, 0, 1, 3}, {1, 0, 3, 3}, {1, 0, 3, 1}, {3, 0, 1, 2}, {2, 0, 2, 3},
	{1, 0, 3, 2}, {2, 0, 2, 1}, {2, 0, 2, 2}, {0, 3, 1, 1}, {0, 1, 3, 1}, {0, 1, 3, 3},
	{0, 3, 1, 3}, {0, 2, 2, 1}, {0, 1, 3, 2}, {0, 2, 2, 3}, {0, 3, 1, 2}, {0, 2, 2, 2},
	{2, 1, 1, 1}, {2, 1, 1, 3}, {2, 1, 1, 2}, {1, 2, 1, 1}, {1, 2, 1, 3}, {1, 2, 1, 2},
	{1, 1, 2, 1}, {1, 1, 2, 3}, {1, 1, 2, 2},
};

// 次数 p (3 / 4) の表の行
/*
5 次 (hexa216 = 型 94 / prism126 = 型 106、境界面の 36 節点四角形 = 型 38) も同じ方法の**実測の表**。gmsh 4.15.2 の
出力を分類する同じ手順 (data/sample/gmsh_hp_order.py) が 3 次・4 次の表 (HEX3_I〜PRISM4_A) を全項再現することを
確かめてから作った (六面体 8 要素・角柱 16 要素がすべて同じ並び)
*/
static const signed char HEX5_I[216][3] = {
	{0, 0, 0}, {5, 0, 0}, {5, 5, 0}, {0, 5, 0}, {0, 0, 5}, {5, 0, 5}, {5, 5, 5}, {0, 5, 5}, {1, 0, 0},
	{2, 0, 0}, {3, 0, 0}, {4, 0, 0}, {0, 1, 0}, {0, 2, 0}, {0, 3, 0}, {0, 4, 0}, {0, 0, 1}, {0, 0, 2},
	{0, 0, 3}, {0, 0, 4}, {5, 1, 0}, {5, 2, 0}, {5, 3, 0}, {5, 4, 0}, {5, 0, 1}, {5, 0, 2}, {5, 0, 3},
	{5, 0, 4}, {4, 5, 0}, {3, 5, 0}, {2, 5, 0}, {1, 5, 0}, {5, 5, 1}, {5, 5, 2}, {5, 5, 3}, {5, 5, 4},
	{0, 5, 1}, {0, 5, 2}, {0, 5, 3}, {0, 5, 4}, {1, 0, 5}, {2, 0, 5}, {3, 0, 5}, {4, 0, 5}, {0, 1, 5},
	{0, 2, 5}, {0, 3, 5}, {0, 4, 5}, {5, 1, 5}, {5, 2, 5}, {5, 3, 5}, {5, 4, 5}, {4, 5, 5}, {3, 5, 5},
	{2, 5, 5}, {1, 5, 5}, {1, 1, 0}, {1, 4, 0}, {4, 4, 0}, {4, 1, 0}, {1, 2, 0}, {1, 3, 0}, {2, 4, 0},
	{3, 4, 0}, {4, 3, 0}, {4, 2, 0}, {3, 1, 0}, {2, 1, 0}, {2, 2, 0}, {2, 3, 0}, {3, 3, 0}, {3, 2, 0},
	{1, 0, 1}, {4, 0, 1}, {4, 0, 4}, {1, 0, 4}, {2, 0, 1}, {3, 0, 1}, {4, 0, 2}, {4, 0, 3}, {3, 0, 4},
	{2, 0, 4}, {1, 0, 3}, {1, 0, 2}, {2, 0, 2}, {3, 0, 2}, {3, 0, 3}, {2, 0, 3}, {0, 1, 1}, {0, 1, 4},
	{0, 4, 4}, {0, 4, 1}, {0, 1, 2}, {0, 1, 3}, {0, 2, 4}, {0, 3, 4}, {0, 4, 3}, {0, 4, 2}, {0, 3, 1},
	{0, 2, 1}, {0, 2, 2}, {0, 2, 3}, {0, 3, 3}, {0, 3, 2}, {5, 1, 1}, {5, 4, 1}, {5, 4, 4}, {5, 1, 4},
	{5, 2, 1}, {5, 3, 1}, {5, 4, 2}, {5, 4, 3}, {5, 3, 4}, {5, 2, 4}, {5, 1, 3}, {5, 1, 2}, {5, 2, 2},
	{5, 3, 2}, {5, 3, 3}, {5, 2, 3}, {4, 5, 1}, {1, 5, 1}, {1, 5, 4}, {4, 5, 4}, {3, 5, 1}, {2, 5, 1},
	{1, 5, 2}, {1, 5, 3}, {2, 5, 4}, {3, 5, 4}, {4, 5, 3}, {4, 5, 2}, {3, 5, 2}, {2, 5, 2}, {2, 5, 3},
	{3, 5, 3}, {1, 1, 5}, {4, 1, 5}, {4, 4, 5}, {1, 4, 5}, {2, 1, 5}, {3, 1, 5}, {4, 2, 5}, {4, 3, 5},
	{3, 4, 5}, {2, 4, 5}, {1, 3, 5}, {1, 2, 5}, {2, 2, 5}, {3, 2, 5}, {3, 3, 5}, {2, 3, 5}, {1, 1, 1},
	{4, 1, 1}, {4, 4, 1}, {1, 4, 1}, {1, 1, 4}, {4, 1, 4}, {4, 4, 4}, {1, 4, 4}, {2, 1, 1}, {3, 1, 1},
	{1, 2, 1}, {1, 3, 1}, {1, 1, 2}, {1, 1, 3}, {4, 2, 1}, {4, 3, 1}, {4, 1, 2}, {4, 1, 3}, {3, 4, 1},
	{2, 4, 1}, {4, 4, 2}, {4, 4, 3}, {1, 4, 2}, {1, 4, 3}, {2, 1, 4}, {3, 1, 4}, {1, 2, 4}, {1, 3, 4},
	{4, 2, 4}, {4, 3, 4}, {3, 4, 4}, {2, 4, 4}, {2, 2, 1}, {2, 3, 1}, {3, 3, 1}, {3, 2, 1}, {2, 1, 2},
	{3, 1, 2}, {3, 1, 3}, {2, 1, 3}, {1, 2, 2}, {1, 2, 3}, {1, 3, 3}, {1, 3, 2}, {4, 2, 2}, {4, 3, 2},
	{4, 3, 3}, {4, 2, 3}, {3, 4, 2}, {2, 4, 2}, {2, 4, 3}, {3, 4, 3}, {2, 2, 4}, {3, 2, 4}, {3, 3, 4},
	{2, 3, 4}, {2, 2, 2}, {3, 2, 2}, {3, 3, 2}, {2, 3, 2}, {2, 2, 3}, {3, 2, 3}, {3, 3, 3}, {2, 3, 3}
};
static const signed char PRISM5_A[126][4] = {
	{5, 0, 0, 0}, {0, 5, 0, 0}, {0, 0, 5, 0}, {5, 0, 0, 5}, {0, 5, 0, 5}, {0, 0, 5, 5}, {4, 1, 0, 0}, {3, 2, 0, 0},
	{2, 3, 0, 0}, {1, 4, 0, 0}, {4, 0, 1, 0}, {3, 0, 2, 0}, {2, 0, 3, 0}, {1, 0, 4, 0}, {5, 0, 0, 1}, {5, 0, 0, 2},
	{5, 0, 0, 3}, {5, 0, 0, 4}, {0, 4, 1, 0}, {0, 3, 2, 0}, {0, 2, 3, 0}, {0, 1, 4, 0}, {0, 5, 0, 1}, {0, 5, 0, 2},
	{0, 5, 0, 3}, {0, 5, 0, 4}, {0, 0, 5, 1}, {0, 0, 5, 2}, {0, 0, 5, 3}, {0, 0, 5, 4}, {4, 1, 0, 5}, {3, 2, 0, 5},
	{2, 3, 0, 5}, {1, 4, 0, 5}, {4, 0, 1, 5}, {3, 0, 2, 5}, {2, 0, 3, 5}, {1, 0, 4, 5}, {0, 4, 1, 5}, {0, 3, 2, 5},
	{0, 2, 3, 5}, {0, 1, 4, 5}, {3, 1, 1, 0}, {1, 1, 3, 0}, {1, 3, 1, 0}, {2, 1, 2, 0}, {1, 2, 2, 0}, {2, 2, 1, 0},
	{3, 1, 1, 5}, {1, 3, 1, 5}, {1, 1, 3, 5}, {2, 2, 1, 5}, {1, 2, 2, 5}, {2, 1, 2, 5}, {4, 1, 0, 1}, {1, 4, 0, 1},
	{1, 4, 0, 4}, {4, 1, 0, 4}, {3, 2, 0, 1}, {2, 3, 0, 1}, {1, 4, 0, 2}, {1, 4, 0, 3}, {2, 3, 0, 4}, {3, 2, 0, 4},
	{4, 1, 0, 3}, {4, 1, 0, 2}, {3, 2, 0, 2}, {2, 3, 0, 2}, {2, 3, 0, 3}, {3, 2, 0, 3}, {4, 0, 1, 1}, {4, 0, 1, 4},
	{1, 0, 4, 4}, {1, 0, 4, 1}, {4, 0, 1, 2}, {4, 0, 1, 3}, {3, 0, 2, 4}, {2, 0, 3, 4}, {1, 0, 4, 3}, {1, 0, 4, 2},
	{2, 0, 3, 1}, {3, 0, 2, 1}, {3, 0, 2, 2}, {3, 0, 2, 3}, {2, 0, 3, 3}, {2, 0, 3, 2}, {0, 4, 1, 1}, {0, 1, 4, 1},
	{0, 1, 4, 4}, {0, 4, 1, 4}, {0, 3, 2, 1}, {0, 2, 3, 1}, {0, 1, 4, 2}, {0, 1, 4, 3}, {0, 2, 3, 4}, {0, 3, 2, 4},
	{0, 4, 1, 3}, {0, 4, 1, 2}, {0, 3, 2, 2}, {0, 2, 3, 2}, {0, 2, 3, 3}, {0, 3, 2, 3}, {3, 1, 1, 1}, {3, 1, 1, 4},
	{3, 1, 1, 2}, {3, 1, 1, 3}, {1, 3, 1, 1}, {1, 3, 1, 4}, {1, 3, 1, 2}, {1, 3, 1, 3}, {1, 1, 3, 1}, {1, 1, 3, 4},
	{1, 1, 3, 2}, {1, 1, 3, 3}, {2, 2, 1, 1}, {2, 2, 1, 4}, {2, 2, 1, 2}, {2, 2, 1, 3}, {1, 2, 2, 1}, {1, 2, 2, 4},
	{1, 2, 2, 2}, {1, 2, 2, 3}, {2, 1, 2, 1}, {2, 1, 2, 4}, {2, 1, 2, 2}, {2, 1, 2, 3}
};
static const signed char *hexp_i(int p, int a) { return ((p == 5) ? HEX5_I[a] : (p == 4) ? HEX4_I[a] : HEX3_I[a]); }
static const signed char *prismp_a(int p, int a) { return ((p == 5) ? PRISM5_A[a] : (p == 4) ? PRISM4_A[a] : PRISM3_A[a]); }
// 節点数から次数 (3 次以上の六面体 / 角柱、それ以外は 0)
static int hexp_order(int nen) { return ((nen == 216) ? 5 : (nen == 125) ? 4 : (nen == 64) ? 3 : 0); }
static int prismp_order(int nen) { return ((nen == 126) ? 5 : (nen == 75) ? 4 : (nen == 40) ? 3 : 0); }

// 六面体の局所節点 a の格子位置 ix (各方向 0..p) と次数 p (完全 Lagrange の 8 / 27 / 64 / 125 / 216 節点。serendipity の
// 20 節点は格子に乗らないので 0)。高次の六面体の辺要素 (eddy3d_pr.c) が部分辺の番号付けに使う
int hex_lattice(int a, int ix[3])
{
	const int p = hexp_order(HexNen);
	if (p) {
		const signed char *s = hexp_i(p, a);
		for (int c = 0; c < 3; c++) ix[c] = s[c];
		return p;
	}
	if ((HexNen != 8) && (HexNen != 27)) return 0;
	for (int c = 0; c < 3; c++) ix[c] = ((HexNen == 27) ? (HEX2_XI[a][c] + 1) : ((HEX2_XI[a][c] + 1) / 2));
	return ((HexNen == 27) ? 2 : 1);
}




// 1 次元の p 次 Lagrange (等間隔の節点 x_i = -1 + 2i/p) の値 v と微分 d
static void lagp(int p, int i, double x, double *v, double *d)
{
	double val = 1, der = 0;
	const double xi = -1 + ((2.0 * i) / p);
	for (int j = 0; j <= p; j++) {
		if (j == i) continue;
		const double xj = -1 + ((2.0 * j) / p);
		const double f = (x - xj) / (xi - xj);
		der = (der * f) + (val / (xi - xj));
		val *= f;
	}
	*v = val;
	*d = der;
}


// 三角形の p 次 Lagrange (多重指数 al、重心座標 lam) の値と dT/dλ_m
static double lagp_tri(int p, const signed char al[3], const double lam[3], double dl[3])
{
	double g[3], gd[3];
	for (int m = 0; m < 3; m++) {
		double val = 1, der = 0;
		for (int s = 0; s < al[m]; s++) {
			const double f = ((p * lam[m]) - s) / (s + 1);
			der = (der * f) + (val * p / (s + 1));
			val *= f;
		}
		g[m] = val;
		gd[m] = der;
	}
	dl[0] = gd[0] * g[1] * g[2];
	dl[1] = g[0] * gd[1] * g[2];
	dl[2] = g[0] * g[1] * gd[2];
	return g[0] * g[1] * g[2];
}


// p 次 (3〜5) の六面体の値と局所微分
static void hexp_shape(int p, double x, double y, double z, double nv[HEXP_NMAX], double dn[HEXP_NMAX][3])
{
	double lx[GL_NMAX], dx[GL_NMAX], ly[GL_NMAX], dy[GL_NMAX], lz[GL_NMAX], dz[GL_NMAX];		// p + 1 個
	for (int i = 0; i <= p; i++) {
		lagp(p, i, x, &lx[i], &dx[i]);
		lagp(p, i, y, &ly[i], &dy[i]);
		lagp(p, i, z, &lz[i], &dz[i]);
	}
	const int nen = (p + 1) * (p + 1) * (p + 1);
	for (int a = 0; a < nen; a++) {
		const signed char *ix = hexp_i(p, a);
		const int i = ix[0], j = ix[1], k = ix[2];
		if (nv) nv[a] = lx[i] * ly[j] * lz[k];
		if (dn) {
			dn[a][0] = dx[i] * ly[j] * lz[k];
			dn[a][1] = lx[i] * dy[j] * lz[k];
			dn[a][2] = lx[i] * ly[j] * dz[k];
		}
	}
}


// p 次 (3〜5) の角柱の値と局所微分 (u, v は λ1, λ2、λ0 = 1 - u - v)
static void prismp_shape(int p, double u, double v, double ze, double nv[PRISMP_NMAX], double dn[PRISMP_NMAX][3])
{
	const double lam[3] = {1 - u - v, u, v};
	double lz[GL_NMAX], dz[GL_NMAX];		// p + 1 個
	for (int k = 0; k <= p; k++) lagp(p, k, ze, &lz[k], &dz[k]);
	const int nen = ((p + 1) * (p + 2) / 2) * (p + 1);
	for (int a = 0; a < nen; a++) {
		const signed char *al = prismp_a(p, a);
		double dl[3];
		const double t = lagp_tri(p, al, lam, dl);
		const int k = al[3];
		if (nv) nv[a] = t * lz[k];
		if (dn) {
			dn[a][0] = (dl[1] - dl[0]) * lz[k];
			dn[a][1] = (dl[2] - dl[0]) * lz[k];
			dn[a][2] = t * dz[k];
		}
	}
}


// 局所座標での微分 dn[a][k] = ∂N_a/∂(ξ,η,ζ)
static void hex2_dshape(int nen, double x, double y, double z, double dn[HEXP_NMAX][3])
{
	if (hexp_order(nen)) { hexp_shape(hexp_order(nen), x, y, z, NULL, dn); return; }
	if (nen == 27) {
		for (int a = 0; a < 27; a++) {
			const int sx = HEX2_XI[a][0], sy = HEX2_XI[a][1], sz = HEX2_XI[a][2];
			dn[a][0] = lag3d(sx, x) * lag3(sy, y) * lag3(sz, z);
			dn[a][1] = lag3(sx, x) * lag3d(sy, y) * lag3(sz, z);
			dn[a][2] = lag3(sx, x) * lag3(sy, y) * lag3d(sz, z);
		}
		return;
	}

	// 20 節点 serendipity
	for (int a = 0; a < 8; a++) {
		const double sx = HEX2_XI[a][0], sy = HEX2_XI[a][1], sz = HEX2_XI[a][2];
		const double fx = 1 + (sx * x), fy = 1 + (sy * y), fz = 1 + (sz * z);
		const double r = (sx * x) + (sy * y) + (sz * z) - 2;
		dn[a][0] = sx * fy * fz * (r + fx) / 8;
		dn[a][1] = sy * fx * fz * (r + fy) / 8;
		dn[a][2] = sz * fx * fy * (r + fz) / 8;
	}
	for (int a = 8; a < 20; a++) {
		const double sx = HEX2_XI[a][0], sy = HEX2_XI[a][1], sz = HEX2_XI[a][2];
		if (sx == 0) {
			dn[a][0] = -2 * x * (1 + (sy * y)) * (1 + (sz * z)) / 4;
			dn[a][1] = sy * (1 - (x * x)) * (1 + (sz * z)) / 4;
			dn[a][2] = sz * (1 - (x * x)) * (1 + (sy * y)) / 4;
		}
		else if (sy == 0) {
			dn[a][0] = sx * (1 - (y * y)) * (1 + (sz * z)) / 4;
			dn[a][1] = -2 * y * (1 + (sx * x)) * (1 + (sz * z)) / 4;
			dn[a][2] = sz * (1 - (y * y)) * (1 + (sx * x)) / 4;
		}
		else {
			dn[a][0] = sx * (1 - (z * z)) * (1 + (sy * y)) / 4;
			dn[a][1] = sy * (1 - (z * z)) * (1 + (sx * x)) / 4;
			dn[a][2] = -2 * z * (1 + (sx * x)) * (1 + (sy * y)) / 4;
		}
	}
}


// 局所座標 (x, y, z) での物理勾配 g[a][i] = ∂N_a/∂x_i と det J
static int hex2_grad(const int32_t *nd, int nen, double x, double y, double z,
	double g[HEXP_NMAX][3], double *det)
{
	double dn[HEXP_NMAX][3];
	hex2_dshape(nen, x, y, z, dn);

	double jm[3][3];
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) jm[i][j] = 0;
	}
	for (int a = 0; a < nen; a++) {
		const int32_t w = nd[a];
		const double p[3] = {Xp[w], Yp[w], Zp[w]};
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) jm[i][j] += p[i] * dn[a][j];
		}
	}
	const double d = (jm[0][0] * ((jm[1][1] * jm[2][2]) - (jm[1][2] * jm[2][1])))
	               - (jm[0][1] * ((jm[1][0] * jm[2][2]) - (jm[1][2] * jm[2][0])))
	               + (jm[0][2] * ((jm[1][0] * jm[2][1]) - (jm[1][1] * jm[2][0])));
	if (d == 0) return 1;

	double ji[3][3];
	ji[0][0] = ((jm[1][1] * jm[2][2]) - (jm[1][2] * jm[2][1])) / d;
	ji[0][1] = ((jm[0][2] * jm[2][1]) - (jm[0][1] * jm[2][2])) / d;
	ji[0][2] = ((jm[0][1] * jm[1][2]) - (jm[0][2] * jm[1][1])) / d;
	ji[1][0] = ((jm[1][2] * jm[2][0]) - (jm[1][0] * jm[2][2])) / d;
	ji[1][1] = ((jm[0][0] * jm[2][2]) - (jm[0][2] * jm[2][0])) / d;
	ji[1][2] = ((jm[0][2] * jm[1][0]) - (jm[0][0] * jm[1][2])) / d;
	ji[2][0] = ((jm[1][0] * jm[2][1]) - (jm[1][1] * jm[2][0])) / d;
	ji[2][1] = ((jm[0][1] * jm[2][0]) - (jm[0][0] * jm[2][1])) / d;
	ji[2][2] = ((jm[0][0] * jm[1][1]) - (jm[0][1] * jm[1][0])) / d;

	for (int a = 0; a < nen; a++) {
		for (int i = 0; i < 3; i++) {
			g[a][i] = (dn[a][0] * ji[0][i]) + (dn[a][1] * ji[1][i]) + (dn[a][2] * ji[2][i]);
		}
	}
	*det = d;

	return 0;
}


// 局所座標での形状関数の値 (hex2_dshape と同じ基底。熱連成の積分点 elem3d_ip で使う)
static void hex2_nval(int nen, double x, double y, double z, double nv[HEXP_NMAX])
{
	if (hexp_order(nen)) { hexp_shape(hexp_order(nen), x, y, z, nv, NULL); return; }
	if (nen == 27) {
		for (int a = 0; a < 27; a++) {
			nv[a] = lag3(HEX2_XI[a][0], x) * lag3(HEX2_XI[a][1], y) * lag3(HEX2_XI[a][2], z);
		}
		return;
	}
	for (int a = 0; a < 8; a++) {
		const double sx = HEX2_XI[a][0], sy = HEX2_XI[a][1], sz = HEX2_XI[a][2];
		nv[a] = (1 + (sx * x)) * (1 + (sy * y)) * (1 + (sz * z)) * ((sx * x) + (sy * y) + (sz * z) - 2) / 8;
	}
	for (int a = 8; a < 20; a++) {
		const double sx = HEX2_XI[a][0], sy = HEX2_XI[a][1], sz = HEX2_XI[a][2];
		if (sx == 0) nv[a] = (1 - (x * x)) * (1 + (sy * y)) * (1 + (sz * z)) / 4;
		else if (sy == 0) nv[a] = (1 - (y * y)) * (1 + (sx * x)) * (1 + (sz * z)) / 4;
		else nv[a] = (1 - (z * z)) * (1 + (sx * x)) * (1 + (sy * y)) / 4;
	}
}


// 3 点 Gauss-Legendre (各方向) の位置と重み
static const double GX3[3] = {-0.7745966692414834, 0.0, 0.7745966692414834};
static const double GW3[3] = {5.0 / 9, 8.0 / 9, 5.0 / 9};
// 4 点 Gauss-Legendre (3 次の要素の各方向、2 次の体積の検算)
static const double GX4[4] = {-0.8611363115940526, -0.3399810435848563,
                               0.3399810435848563,  0.8611363115940526};
static const double GW4[4] = {0.3478548451374538, 0.6521451548625461,
                              0.6521451548625461, 0.3478548451374538};


// 各方向 n 点 (3 / 4 / 5) の Gauss-Legendre ([-1, 1])。3 / 4 点は従来の定数 (2 次・3 次の答えを変えない)
static void ho_gauss(int n, const double **x, const double **w)
{
	// [-1,1] の n 点 (3 / 4 は表、5 以上は gauss01 から作って n ごとに保持する。**表に無い n を 5 点に落とさない**
	// こと : 以前は 5 以外の n でも 5 点が返った)。最初の呼び出しは格子の組み立て (直列) で起きる
	static double xs[GL_NMAX + 1][GL_NMAX], ws[GL_NMAX + 1][GL_NMAX];
	static int init[GL_NMAX + 1];
	if (n == 3) { *x = GX3; *w = GW3; return; }
	if (n == 4) { *x = GX4; *w = GW4; return; }
	if (n > GL_NMAX) n = GL_NMAX;
	if (!init[n]) {
		double t[GL_NMAX], u[GL_NMAX];
		gauss01(n, t, u);
		for (int k = 0; k < n; k++) { xs[n][k] = (2 * t[k]) - 1; ws[n][k] = 2 * u[k]; }
		init[n] = 1;
	}
	*x = xs[n];
	*w = ws[n];
}


// 2 次六面体の要素行列 (3x3x3 Gauss)。vol には体積 Σ w |det J| を返す
static int hex2_element(int e, const double c[6], double ke[HEXP_NMAX][HEXP_NMAX], double *vol)
{
	int32_t nd[ELEM_NMAX];
	const int nen = hex_nodes(e, nd);

	memset(ke, 0, sizeof(double) * HEXP_NMAX * HEXP_NMAX);
	// 2 次は各方向 3 点、3 次は 4 点 (直方体の剛性・質量が厳密)
	const int ng = (hexp_order(nen) ? (hexp_order(nen) + 1) : 3);
	const double *gx, *gw;
	ho_gauss(ng, &gx, &gw);
	double v = 0, sgn = 0;
	for (int i = 0; i < ng; i++) {
	for (int j = 0; j < ng; j++) {
	for (int k = 0; k < ng; k++) {
		double g[HEXP_NMAX][3], det;
		if (hex2_grad(nd, nen, gx[i], gx[j], gx[k], g, &det)) return 1;
		if (sgn == 0) sgn = ((det > 0) ? 1 : -1);
		else if (det * sgn <= 0) return 1;		// 要素内で符号が変わる = 裏返り
		const double dw = gw[i] * gw[j] * gw[k] * fabs(det);
		v += dw;
		for (int l = 0; l < nen; l++) {
			for (int m = 0; m < nen; m++) {
				ke[l][m] += dw * ((c[0] * g[l][0] * g[m][0])
				                + (c[1] * g[l][1] * g[m][1])
				                + (c[2] * g[l][2] * g[m][2])
				                + (c[3] * ((g[l][0] * g[m][1]) + (g[l][1] * g[m][0])))
				                + (c[4] * ((g[l][1] * g[m][2]) + (g[l][2] * g[m][1])))
				                + (c[5] * ((g[l][2] * g[m][0]) + (g[l][0] * g[m][2]))));
			}
		}
	}
	}
	}
	*vol = v;

	return 0;
}


// 要素中心での勾配 (場の出力用)
int hex2_grad_center(int e, double g[HEXP_NMAX][3], int *nen)
{
	int32_t nd[ELEM_NMAX];
	double det;

	*nen = hex_nodes(e, nd);

	return hex2_grad(nd, *nen, 0, 0, 0, g, &det);
}


// ---- 2 次の角柱 (15 / 18 節点、等パラメトリック) ----
//
// 局所節点 -> (6 節点三角形の基底番号 t, ζ の格子位置 z ∈ {-1,0,+1})。
// 並びは Gmsh の実測 (fem.h 参照)。三角形の基底は tri6 と同じ
// (頂点 0,1,2、辺 (0,1)(1,2)(2,0) = 3,4,5)。
//
// 18 節点は積 N = T_t(λ) l_z(ζ) (完全 2 次)。15 節点は serendipity:
//   角 (t<3, z=±1)     : N = λ(2λ-1)(1+sζ)/2 - λ(1-ζ²)/2
//   三角形の辺 (z=±1)  : N = 2λiλj(1+sζ)
//   鉛直の辺 (t<3, z=0): N = λ(1-ζ²)
// (和が 1 になることは代数的に確認済み。線形場の恒等式でも毎回検査される)
//
// 積分は面内 Duffy 9 点 (tri_quad) x ζ 方向 3 点 Gauss の 27 点。面内は
// λ の 4 次まで、ζ は 5 次まで厳密なので、**まっすぐな角柱 (上下が平行移動)
// では剛性・体積とも厳密**。
static const signed char PRISM2_T[18] = {0, 1, 2, 0, 1, 2,
                                         3, 5, 0, 4, 1, 2, 3, 5, 4,
                                         3, 5, 4};
static const signed char PRISM2_Z[18] = {-1, -1, -1, +1, +1, +1,
                                         -1, -1, 0, -1, 0, 0, +1, +1, +1,
                                         0, 0, 0};
// tri6 の辺基底番号 -> 頂点対
static const signed char TRI6_ED[3][2] = {{0, 1}, {1, 2}, {2, 0}};

#define NQTRI (9)
static void tri_quad(double lam[NQTRI][3], double wq[NQTRI]);


// 局所座標 (u, v, ζ) での微分 dn[a][k] (λ = (1-u-v, u, v))
static void prism2_dshape(int nen, double u, double v, double ze, double dn[PRISMP_NMAX][3])
{
	if (prismp_order(nen)) { prismp_shape(prismp_order(nen), u, v, ze, NULL, dn); return; }
	const double lam[3] = {1 - u - v, u, v};
	static const signed char dlu[3] = {-1, 1, 0};
	static const signed char dlv[3] = {-1, 0, 1};

	// tri6 の基底 T_t(λ) とその (u, v) 微分
	double tn[6], tu[6], tv[6];
	for (int t = 0; t < 3; t++) {
		tn[t] = lam[t] * ((2 * lam[t]) - 1);
		tu[t] = ((4 * lam[t]) - 1) * dlu[t];
		tv[t] = ((4 * lam[t]) - 1) * dlv[t];
	}
	for (int t = 0; t < 3; t++) {
		const int a = TRI6_ED[t][0], b = TRI6_ED[t][1];
		tn[3 + t] = 4 * lam[a] * lam[b];
		tu[3 + t] = 4 * ((dlu[a] * lam[b]) + (lam[a] * dlu[b]));
		tv[3 + t] = 4 * ((dlv[a] * lam[b]) + (lam[a] * dlv[b]));
	}

	if (nen == 18) {
		for (int a = 0; a < 18; a++) {
			const int t = PRISM2_T[a], z = PRISM2_Z[a];
			dn[a][0] = tu[t] * lag3(z, ze);
			dn[a][1] = tv[t] * lag3(z, ze);
			dn[a][2] = tn[t] * lag3d(z, ze);
		}
		return;
	}

	// 15 節点 serendipity
	for (int a = 0; a < 15; a++) {
		const int t = PRISM2_T[a], z = PRISM2_Z[a];
		if (z == 0) {
			// 鉛直の辺 : N = λt (1 - ζ²)
			dn[a][0] = dlu[t] * (1 - (ze * ze));
			dn[a][1] = dlv[t] * (1 - (ze * ze));
			dn[a][2] = lam[t] * (-2 * ze);
		}
		else if (t < 3) {
			// 角 : N = λ(2λ-1)(1+sζ)/2 - λ(1-ζ²)/2
			const double f = (1 + (z * ze)) / 2;
			const double dl = (((4 * lam[t]) - 1) * f) - ((1 - (ze * ze)) / 2);
			dn[a][0] = dlu[t] * dl;
			dn[a][1] = dlv[t] * dl;
			dn[a][2] = (lam[t] * ((2 * lam[t]) - 1) * z / 2) + (lam[t] * ze);
		}
		else {
			// 三角形の辺 : N = 2 λi λj (1+sζ)
			const int i = TRI6_ED[t - 3][0], j = TRI6_ED[t - 3][1];
			const double f = 1 + (z * ze);
			dn[a][0] = 2 * f * ((dlu[i] * lam[j]) + (lam[i] * dlu[j]));
			dn[a][1] = 2 * f * ((dlv[i] * lam[j]) + (lam[i] * dlv[j]));
			dn[a][2] = 2 * lam[i] * lam[j] * z;
		}
	}
}


// 局所座標での形状関数の値 (prism2_dshape と同じ基底)
static void prism2_nval(int nen, double u, double v, double ze, double nv[PRISMP_NMAX])
{
	if (prismp_order(nen)) { prismp_shape(prismp_order(nen), u, v, ze, nv, NULL); return; }
	const double lam[3] = {1 - u - v, u, v};
	double tn[6];
	for (int t = 0; t < 3; t++) tn[t] = lam[t] * ((2 * lam[t]) - 1);
	for (int t = 0; t < 3; t++) tn[3 + t] = 4 * lam[TRI6_ED[t][0]] * lam[TRI6_ED[t][1]];
	if (nen == 18) {
		for (int a = 0; a < 18; a++) nv[a] = tn[PRISM2_T[a]] * lag3(PRISM2_Z[a], ze);
		return;
	}
	for (int a = 0; a < 15; a++) {
		const int t = PRISM2_T[a], z = PRISM2_Z[a];
		if (z == 0) nv[a] = lam[t] * (1 - (ze * ze));
		else if (t < 3) nv[a] = (lam[t] * ((2 * lam[t]) - 1) * (1 + (z * ze)) / 2) - (lam[t] * (1 - (ze * ze)) / 2);
		else nv[a] = 2 * lam[TRI6_ED[t - 3][0]] * lam[TRI6_ED[t - 3][1]] * (1 + (z * ze));
	}
}


// ---- 次数の違う要素が面を共有する格子 (四面体・六面体・角柱・三角形) ----
//
// 読み込んだ後で、低い次数の要素を格子の最高次数 P に上げる (P = 2 の六面体・角柱は完全 Lagrange の 27 / 18 節点。
// serendipity と完全の混在もこれで適合する)。P 次の節点は参照座標の格子位置で表し、辺・面の上の節点は全体の頂点番号で
// 決まる鍵で隣の要素と共有する (隣が元から P 次ならその節点をそのまま使う):
//   辺 : 両端の (頂点番号, 重み) の組 (頂点番号の昇順。重みの和が P)
//   三角形の面 : 3 頂点の (頂点番号, 重心座標 x P) の組 (昇順)
//   四角形の面 : 4 頂点の番号 (昇順) と、最小番号の角を原点・その隣で番号の小さい方を第 1 軸にした格子位置 (s, t)
// 新しい節点の位置は元の要素の q 次の等パラメトリック写像の上 (低い次数の要素の形は変わらない)。元の q 次の節点で
// 位置の一致するもの (参照座標が同じ) は使い回す。使われなくなった節点は行列の行が空になるだけ。
typedef struct { int32_t k[9]; int32_t id; } hokey_t;
static hokey_t *HoTab = NULL;
static size_t HoCap = 0;

static size_t ho_hash(const int32_t k[9])
{
	uint64_t h = 1469598103934665603ULL;
	for (int i = 0; i < 9; i++) { h ^= (uint64_t)(uint32_t)k[i]; h *= 1099511628211ULL; }
	return (size_t)(h & (uint64_t)(HoCap - 1));
}

// 頂点 g[nv] と重み a[nv] から鍵 (a > 0 の組を頂点番号の昇順、k[0] = 組の数) を作り、組の数を返す
static int ho_key(const int32_t *g, const int *a, int nv, int32_t k[9])
{
	int32_t pg[4];
	int pa[4], n = 0;
	for (int v = 0; v < nv; v++) {
		if (a[v] > 0) { pg[n] = g[v]; pa[n] = a[v]; n++; }
	}
	for (int x = 1; x < n; x++) {
		for (int y = x; (y > 0) && (pg[y] < pg[y - 1]); y--) {
			const int32_t tg = pg[y]; pg[y] = pg[y - 1]; pg[y - 1] = tg;
			const int ta = pa[y]; pa[y] = pa[y - 1]; pa[y - 1] = ta;
		}
	}
	for (int i = 0; i < 9; i++) k[i] = -1;
	k[0] = n;
	for (int m = 0; (m < n) && (m < 3); m++) { k[1 + (2 * m)] = pg[m]; k[2 + (2 * m)] = pa[m]; }
	return n;
}

// 四角形の面の内部の節点の鍵。角 c[4] は面の格子 (u, v) の (0,0) (P,0) (P,P) (0,P) の頂点番号
static void ho_qkey(const int32_t c[4], int u, int v, int P, int32_t k[9])
{
	static const int cu[4] = {0, 1, 1, 0}, cv[4] = {0, 0, 1, 1};
	int o = 0;
	for (int m = 1; m < 4; m++) {
		if (c[m] < c[o]) o = m;
	}
	const int ou = cu[o] * P, ov = cv[o] * P;
	int nu = -1, nv = -1;
	for (int m = 0; m < 4; m++) {
		if ((cu[m] != cu[o]) && (cv[m] == cv[o])) nu = m;
		if ((cu[m] == cu[o]) && (cv[m] != cv[o])) nv = m;
	}
	const int du = ((u > ou) ? (u - ou) : (ou - u)), dv = ((v > ov) ? (v - ov) : (ov - v));
	int32_t sg[4] = {c[0], c[1], c[2], c[3]};
	for (int x = 1; x < 4; x++) {
		for (int y = x; (y > 0) && (sg[y] < sg[y - 1]); y--) { const int32_t t = sg[y]; sg[y] = sg[y - 1]; sg[y - 1] = t; }
	}
	k[0] = 9;
	for (int m = 0; m < 4; m++) k[1 + m] = sg[m];
	k[5] = ((c[nu] < c[nv]) ? du : dv);
	k[6] = ((c[nu] < c[nv]) ? dv : du);
	k[7] = k[8] = -1;
}

// 鍵を引く (無ければ id を登録して id を返す。id < 0 なら引くだけで、無ければ -1)
static int32_t ho_find(const int32_t k[9], int32_t id)
{
	size_t h = ho_hash(k);
	for (;;) {
		hokey_t *e = &HoTab[h];
		if (e->id < 0) {
			if (id >= 0) {
				for (int i = 0; i < 9; i++) e->k[i] = k[i];
				e->id = id;
			}
			return id;
		}
		int same = 1;
		for (int i = 0; i < 9; i++) same &= (e->k[i] == k[i]);
		if (same) return e->id;
		h = (h + 1) & (HoCap - 1);
	}
}

static int HoNodeCap = 0, HoMiss = 0;
static int32_t ho_new_node(const double x[3])
{
	if (NNode >= HoNodeCap) {
		HoNodeCap = ((HoNodeCap > NNode) ? HoNodeCap : NNode) * 2 + 64;
		Xp = (double *)realloc(Xp, (size_t)HoNodeCap * sizeof(double));
		Yp = (double *)realloc(Yp, (size_t)HoNodeCap * sizeof(double));
		Zp = (double *)realloc(Zp, (size_t)HoNodeCap * sizeof(double));
	}
	Xp[NNode] = x[0]; Yp[NNode] = x[1]; Zp[NNode] = x[2];
	MeshRaisedNode++;
	return NNode++;
}

// 要素の種別と節点数から次数 (六面体・角柱の 20 / 15 節点は 2)
static int ho_kind_order(int kind, int nen)
{
	if (kind == MESHELEM_HEX) return ((nen == 8) ? 1 : (nen <= 27) ? 2 : hexp_order(nen));
	if (kind == MESHELEM_PRISM) return ((nen == 6) ? 1 : (nen <= 18) ? 2 : prismp_order(nen));
	for (int p = 1; p <= 7; p++) {
		if ((((p + 1) * (p + 2) * (p + 3)) / 6) == nen) return p;
	}
	return 0;
}

// 種別 kind・節点数 nen の要素の局所節点 a の参照座標の格子位置 d[4] (次数 q = ho_kind_order):
// 四面体 = 重心座標 x q、六面体 = (i, j, k) ∈ [0, q]、角柱 = (三角形の重心座標 x q, ζ の格子位置)
static void ho_desc(int kind, int nen, int a, int d[4])
{
	const int q = ho_kind_order(kind, nen);
	d[0] = d[1] = d[2] = d[3] = 0;
	if (kind == MESHELEM_HEX) {
		if (q >= 3) {
			const signed char *ix = hexp_i(q, a);
			for (int c = 0; c < 3; c++) d[c] = ix[c];
		}
		else {
			for (int c = 0; c < 3; c++) d[c] = ((q == 2) ? (HEX2_XI[a][c] + 1) : ((HEX2_XI[a][c] + 1) / 2));
		}
		return;
	}
	if (kind == MESHELEM_PRISM) {
		if (q >= 3) {
			const signed char *al = prismp_a(q, a);
			for (int c = 0; c < 4; c++) d[c] = al[c];
			return;
		}
		const int t = PRISM2_T[a], z = PRISM2_Z[a];
		if (q == 2) {
			if (t < 3) d[t] = 2;
			else { d[TRI6_ED[t - 3][0]] = 1; d[TRI6_ED[t - 3][1]] = 1; }
			d[3] = z + 1;
		}
		else {
			d[t] = 1;
			d[3] = (z + 1) / 2;
		}
		return;
	}
	int ta[TET_NMAX][4];
	tet_rule(q, ta);
	for (int c = 0; c < 4; c++) d[c] = ta[a][c];
}

// 次数 P の目標の要素の節点数 (六面体・角柱の 2 次は完全 Lagrange)
static int ho_target_nen(int kind, int P)
{
	if (kind == MESHELEM_HEX) return ((P == 2) ? 27 : ((P + 1) * (P + 1) * (P + 1)));
	if (kind == MESHELEM_PRISM) return ((P == 2) ? 18 : ((((P + 1) * (P + 2)) / 2) * (P + 1)));
	return (((P + 1) * (P + 2) * (P + 3)) / 6);
}

// 頂点 nd[] の要素の P 次の格子位置 d の節点の分類 : 戻り値 0 = 頂点 (*vid)、1 = 鍵つき (辺・面、k)、2 = 内部
static int ho_classify(int kind, const int32_t *nd, const int d[4], int P, int32_t *vid, int32_t k[9])
{
	if (kind == MESHELEM_HEX) {
		int nb = 0, fr[3], nf = 0;
		for (int c = 0; c < 3; c++) {
			if ((d[c] == 0) || (d[c] == P)) nb++;
			else fr[nf++] = c;
		}
		// 格子位置 g (各成分 0 / 1) の角の頂点番号
#define HO_HC(g0, g1, g2) nd[((g0) == 0) ? (((g1) == 0) ? (((g2) == 0) ? 0 : 4) : (((g2) == 0) ? 3 : 7)) : (((g1) == 0) ? (((g2) == 0) ? 1 : 5) : (((g2) == 0) ? 2 : 6))]
		if (nb == 3) { *vid = HO_HC(d[0] / P, d[1] / P, d[2] / P); return 0; }
		if (nb == 2) {
			int ca[3], cb[3];
			for (int c = 0; c < 3; c++) { ca[c] = d[c] / P; cb[c] = d[c] / P; }
			ca[fr[0]] = 0; cb[fr[0]] = 1;
			const int32_t g[2] = {HO_HC(ca[0], ca[1], ca[2]), HO_HC(cb[0], cb[1], cb[2])};
			const int w[2] = {P - d[fr[0]], d[fr[0]]};
			ho_key(g, w, 2, k);
			return 1;
		}
		if (nb == 1) {
			int fx = 0;
			for (int c = 0; c < 3; c++) {
				if ((c != fr[0]) && (c != fr[1])) fx = c;
			}
			int g4[4][3];
			static const int qu[4] = {0, 1, 1, 0}, qv[4] = {0, 0, 1, 1};
			int32_t cg[4];
			for (int m = 0; m < 4; m++) {
				g4[m][fx] = d[fx] / P;
				g4[m][fr[0]] = qu[m];
				g4[m][fr[1]] = qv[m];
				cg[m] = HO_HC(g4[m][0], g4[m][1], g4[m][2]);
			}
#undef HO_HC
			ho_qkey(cg, d[fr[0]], d[fr[1]], P, k);
			return 1;
		}
		return 2;
	}
	if (kind == MESHELEM_PRISM) {
		int nz = 0, iz[3];
		for (int c = 0; c < 3; c++) {
			if (d[c] > 0) iz[nz++] = c;
		}
		const int kb = ((d[3] == 0) || (d[3] == P));
		const int lv = ((d[3] == P) ? 3 : 0);
		if ((nz == 1) && kb) { *vid = nd[lv + iz[0]]; return 0; }
		if (nz == 1) {
			const int32_t g[2] = {nd[iz[0]], nd[3 + iz[0]]};
			const int w[2] = {P - d[3], d[3]};
			ho_key(g, w, 2, k);
			return 1;
		}
		if ((nz == 2) && kb) {
			const int32_t g[2] = {nd[lv + iz[0]], nd[lv + iz[1]]};
			const int w[2] = {d[iz[0]], d[iz[1]]};
			ho_key(g, w, 2, k);
			return 1;
		}
		if (nz == 2) {
			// 側面 : 三角形の辺 (i, j) x ζ。u = j の重み (i -> j)、v = ζ の格子位置
			const int i = iz[0], j = iz[1];
			const int32_t cg[4] = {nd[i], nd[j], nd[3 + j], nd[3 + i]};
			ho_qkey(cg, d[j], d[3], P, k);
			return 1;
		}
		if (kb) {
			const int32_t g[3] = {nd[lv], nd[lv + 1], nd[lv + 2]};
			ho_key(g, d, 3, k);
			return 1;
		}
		return 2;
	}
	const int cnt = ho_key(nd, d, 4, k);
	if (cnt == 1) {
		for (int v = 0; v < 4; v++) {
			if (d[v] > 0) *vid = nd[v];
		}
		return 0;
	}
	return ((cnt == 4) ? 2 : 1);
}

// 要素 (種別 kind、節点数 nen、節点 nd は gmsh の並び) の参照座標 d (次数 P の格子位置) の物理座標 (元の写像)
static void ho_position(int kind, int nen, const int32_t *nd, const int d[4], int P, double x[3])
{
	double nv[ELEM_NMAX];
	if (kind == MESHELEM_HEX) {
		const double xi[3] = {-1 + ((2.0 * d[0]) / P), -1 + ((2.0 * d[1]) / P), -1 + ((2.0 * d[2]) / P)};
		if (nen == 8) {
			for (int a = 0; a < 8; a++) {
				nv[a] = (1 + (HEX2_XI[a][0] * xi[0])) * (1 + (HEX2_XI[a][1] * xi[1])) * (1 + (HEX2_XI[a][2] * xi[2])) / 8;
			}
		}
		else hex2_nval(nen, xi[0], xi[1], xi[2], nv);
	}
	else if (kind == MESHELEM_PRISM) {
		const double u = (double)d[1] / P, v = (double)d[2] / P, ze = -1 + ((2.0 * d[3]) / P);
		if (nen == 6) {
			const double lam[3] = {1 - u - v, u, v};
			for (int a = 0; a < 6; a++) nv[a] = lam[PRISM2_T[a]] * (1 + (PRISM2_Z[a] * ze)) / 2;
		}
		else prism2_nval(nen, u, v, ze, nv);
	}
	else {
		const int q = ho_kind_order(kind, nen);
		double lam[4];
		for (int v = 0; v < 4; v++) lam[v] = (double)d[v] / P;
		int aq[TET_NMAX][4];
		tet_rule(q, aq);		// 新しい節点を作るときだけ呼ばれる
		for (int b = 0; b < nen; b++) {
			double dn[4];
			nv[b] = ((q == 1) ? lam[b] : lagr_tet(q, aq[b], lam, dn));
		}
	}
	x[0] = x[1] = x[2] = 0;
	for (int b = 0; b < nen; b++) {
		x[0] += nv[b] * Xp[nd[b]]; x[1] += nv[b] * Yp[nd[b]]; x[2] += nv[b] * Zp[nd[b]];
	}
}

// 要素を P 次に上げた節点 out[] (目標の並び)。reg = 1 なら目標と同じ形の要素の節点を登録するだけ
static void ho_desc_all(int kind, int nen, int (*d)[4])
{
	if (kind == MESHELEM_TET) {
		tet_rule(ho_kind_order(kind, nen), d);
		return;
	}
	for (int a = 0; a < nen; a++) ho_desc(kind, nen, a, d[a]);
}

static void ho_elem3(int kind, int nen, const int32_t *nd, int P, int32_t *out, int reg)
{
	const int nt = ho_target_nen(kind, P), q = ho_kind_order(kind, nen);
	int dp[ELEM_NMAX][4], dq[ELEM_NMAX][4];
	ho_desc_all(kind, nt, dp);
	if (!reg) ho_desc_all(kind, nen, dq);
	for (int i = 0; i < nt; i++) {
		const int *d = dp[i];
		int32_t vid = -1, k[9];
		const int cl = ho_classify(kind, nd, d, P, &vid, k);
		if (reg) {
			if (cl == 1) (void)ho_find(k, nd[i]);
			continue;
		}
		if (cl == 0) { out[i] = vid; continue; }
		int32_t id = ((cl == 1) ? ho_find(k, -1) : -1);
		if (id < 0) {
			// 元の q 次の節点で同じ参照座標のもの
			for (int j = 0; (j < nen) && (id < 0); j++) {
				int eq = 1;
				for (int c = 0; c < 4; c++) eq &= ((dq[j][c] * P) == (d[c] * q));
				if (eq) id = nd[j];
			}
		}
		if (id < 0) {
			double x[3];
			ho_position(kind, nen, nd, d, P, x);
			id = ho_new_node(x);
		}
		if (cl == 1) (void)ho_find(k, id);
		out[i] = id;
	}
}

// 三角形 t の次数 (1〜7)
static int ho_tri_order(int t)
{
	const int t3 = tri3_order(&Tri3[t * TRI_HO_MAX]);
	return ((t3 >= 3) ? t3 : (Tri2[t * 3] != Tri[t * 3]) ? 2 : 1);
}

// 三角形 t の節点 (gmsh の並び)
static int ho_tri_nodes(int t, int32_t *nd)
{
	const int q = ho_tri_order(t);
	for (int l = 0; l < 3; l++) nd[l] = Tri[(t * 3) + l];
	if (q == 2) { for (int l = 0; l < 3; l++) nd[3 + l] = Tri2[(t * 3) + l]; }
	if (q >= 3) {
		const int nho = (((q + 1) * (q + 2)) / 2) - 3;
		for (int l = 0; l < nho; l++) nd[3 + l] = Tri3[(t * TRI_HO_MAX) + l];
	}
	return q;
}

// 三角形 (q 次、節点 nd) を P 次に上げた節点 out[] (gmsh の並び)。reg = 1 なら登録するだけ。
// 断面 2 次元では三角形が要素 (内部の節点は自分だけのもの)、3 次元では要素の面 (内部も面の鍵で引く)
static void ho_tri(int q, const int32_t *nd, int P, int32_t *out, int reg, int dim2)
{
	int aq[TRI_NMAX][3], ap[TRI_NMAX][3];
	const int nq = tri_rule(q, aq), np = tri_rule(P, ap);
	for (int i = 0; i < np; i++) {
		int32_t k[9];
		const int cnt = ho_key(nd, ap[i], 3, k);
		if (cnt == 1) {
			if (!reg) {
				for (int v = 0; v < 3; v++) {
					if (ap[i][v] > 0) out[i] = nd[v];
				}
			}
			continue;
		}
		const int own = (dim2 && (cnt == 3));		// 断面 2 次元の三角形の内部
		if (reg) {
			if (!own) (void)ho_find(k, nd[i]);
			continue;
		}
		int32_t id = (own ? -1 : ho_find(k, -1));
		if ((id < 0) && !dim2) { HoMiss = 1; out[i] = -1; continue; }
		for (int j = 0; (j < nq) && (id < 0); j++) {
			if (((aq[j][0] * P) == (ap[i][0] * q)) && ((aq[j][1] * P) == (ap[i][1] * q)) && ((aq[j][2] * P) == (ap[i][2] * q))) id = nd[j];
		}
		if (id < 0) {
			double lam[4] = {(double)ap[i][0] / P, (double)ap[i][1] / P, (double)ap[i][2] / P, 0}, x[3] = {0, 0, 0};
			for (int b = 0; b < nq; b++) {
				double dn[4];
				const int a4[4] = {aq[b][0], aq[b][1], aq[b][2], 0};
				const double w = ((q == 1) ? lam[b] : lagr_tet(q, a4, lam, dn));
				x[0] += w * Xp[nd[b]]; x[1] += w * Yp[nd[b]]; x[2] += w * Zp[nd[b]];
			}
			id = ho_new_node(x);
		}
		if (!own) (void)ho_find(k, id);
		out[i] = id;
	}
}

// 三角形 t を P 次の節点 nd で書き直す
static void ho_tri_store(int t, int P, const int32_t *nd)
{
	const int nho = (P >= 3) ? ((((P + 1) * (P + 2)) / 2) - 3) : 0;
	for (int l = 0; l < 3; l++) Tri2[(t * 3) + l] = ((P == 2) ? nd[3 + l] : nd[l]);
	for (int l = 0; l < TRI_HO_MAX; l++) Tri3[(t * TRI_HO_MAX) + l] = ((l < nho) ? nd[3 + l] : -1);
}

// 格納済みの 3 次元要素 (種別 kind の e 番目) の節点数と節点
static int ho_stored(int kind, int e, int32_t *nd)
{
	if (kind == MESHELEM_HEX) {
		for (int l = 0; l < 8; l++) nd[l] = Hex[(e * 8) + l];
		for (int l = 8; l < HexNen; l++) nd[l] = Hex2[((size_t)e * HEX2_W) + (l - 8)];
		return HexNen;
	}
	if (kind == MESHELEM_PRISM) {
		for (int l = 0; l < 6; l++) nd[l] = Prism[(e * 6) + l];
		for (int l = 6; l < PrismNen; l++) nd[l] = Prism2[((size_t)e * PRISM2_W) + (l - 6)];
		return PrismNen;
	}
	const int nen = (((TetOrder + 1) * (TetOrder + 2) * (TetOrder + 3)) / 6);
	for (int l = 0; l < 4; l++) nd[l] = Tet[(e * 4) + l];
	if (TetOrder == 2) { for (int l = 0; l < 6; l++) nd[4 + l] = Tet2[(e * 6) + l]; }
	if (TetOrder >= 3) { for (int l = 0; l < (nen - 4); l++) nd[4 + l] = Tet3[((size_t)e * TET_HO_MAX) + l]; }
	return ((TetOrder >= 1) ? nen : 4);
}

static int ho_raise(void)
{
	static const int KINDS[3] = {MESHELEM_TET, MESHELEM_HEX, MESHELEM_PRISM};
	const int nst[3] = {NTet, NHex, NPrism};
	const int dim3 = ((NTet + NHex + NPrism + PendN) > 0);
	int32_t nd[ELEM_NMAX], out[ELEM_NMAX];
	int P = 0, qmin = 99, ser = 0, cpl = 0, onlytet = 1;
	int nall[3] = {0, 0, 0};
	for (int s = 0; s < 3; s++) {
		const int nen0 = ((s == 0) ? (((TetOrder + 1) * (TetOrder + 2) * (TetOrder + 3)) / 6) : (s == 1) ? HexNen : PrismNen);
		nall[s] = nst[s];
		if (nst[s] > 0) {
			const int q = ho_kind_order(KINDS[s], nen0);
			if (q > P) P = q;
			if (q < qmin) qmin = q;
			if ((s > 0) && (nen0 == ((s == 1) ? 20 : 15))) ser = 1;
			if ((s > 0) && (nen0 == ((s == 1) ? 27 : 18))) cpl = 1;
			if (s > 0) onlytet = 0;
		}
	}
	for (int i = 0; i < PendN; i++) {
		const int q = ho_kind_order(PendKind[i], PendNen[i]);
		const int s = ((PendKind[i] == MESHELEM_HEX) ? 1 : (PendKind[i] == MESHELEM_PRISM) ? 2 : 0);
		nall[s]++;
		if (q > P) P = q;
		if (q < qmin) qmin = q;
		if ((s > 0) && (PendNen[i] == ((s == 1) ? 20 : 15))) ser = 1;
		if ((s > 0) && (PendNen[i] == ((s == 1) ? 27 : 18))) cpl = 1;
		if (s > 0) onlytet = 0;
	}
	int trimix = 0;
	if (!dim3) {
		// 断面 2 次元 : 三角形が要素
		if (NTri < 1) return 0;
		P = qmin = ho_tri_order(0);
		for (int t = 1; t < NTri; t++) {
			const int q = ho_tri_order(t);
			if (q > P) P = q;
			if (q < qmin) qmin = q;
		}
		if (P == qmin) return 0;
		trimix = 1;
	}
	else if ((PendN == 0) && (P == qmin) && !(ser && cpl)) return 0;
	if (dim3 && (NPyr > 0)) {
		printf("%s\n", "*** mesh : elements of different orders cannot be raised in a mesh with pyramids "
			"(pyramids are first order only)");
		return 1;
	}
	if (dim3) {
		for (int t = 0; t < NTri; t++) {
			if (ho_tri_order(t) > P) {
				printf("*** mesh : a triangle is order %d but the highest %s order is %d\n", ho_tri_order(t),
					(onlytet ? "tetrahedron" : "element"), P);
				return 1;
			}
		}
	}
	size_t want = (size_t)NTri * TRI_NMAX * 2 + 64;
	for (int s = 0; s < 3; s++) want += (size_t)nall[s] * (size_t)ho_target_nen(KINDS[s], P) * 2;
	HoCap = 64;
	while (HoCap < want) HoCap *= 2;
	HoTab = (hokey_t *)malloc(HoCap * sizeof(hokey_t));
	if (HoTab == NULL) { printf("%s\n", "*** mesh : out of memory raising the element orders"); return 1; }
	for (size_t h = 0; h < HoCap; h++) HoTab[h].id = -1;
	HoNodeCap = NNode;

	// (1) 元から目標の形 (P 次、六面体・角柱は 2 次なら完全) の要素と三角形の節点を登録する
	if (dim3) {
		for (int s = 0; s < 3; s++) {
			const int kind = KINDS[s], nt = ho_target_nen(kind, P);
			for (int e = 0; e < nst[s]; e++) {
				if (ho_stored(kind, e, nd) == nt) ho_elem3(kind, nt, nd, P, NULL, 1);
			}
		}
		for (int i = 0; i < PendN; i++) {
			if (PendNen[i] == ho_target_nen(PendKind[i], P)) ho_elem3(PendKind[i], PendNen[i], &PendNd[(size_t)i * ELEM_NMAX], P, NULL, 1);
		}
	}
	// 断面 2 次元の三角形は要素なので、元から P 次のものを登録する。3 次元の三角形・四角形は電極の面で、(3) (4) で要素の
	// 面の節点に置き換える (**次数が目標と同じ面も**。面の下の要素が上げられたとき、元の面の節点は要素のものでない。
	// 実測 : 曲がった同軸の 9 節点四角形の電極を残したら内側の電極が「穴だらけ」になり C' −16%)
	if (!dim3) {
		for (int t = 0; t < NTri; t++) {
			if (ho_tri_nodes(t, nd) == P) ho_tri(P, nd, P, NULL, 1, 1);
		}
	}

	// (2) 3 次元の要素を目標の形に上げて格納し直す (種別ごとに、格納済み -> 取っておいたもの の順)
	int nraise = 0;
	if (dim3) {
		for (int s = 0; s < 3; s++) {
			const int kind = KINDS[s], nt = ho_target_nen(kind, P);
			if (nall[s] == 0) continue;
			const int nv = ((kind == MESHELEM_HEX) ? 8 : (kind == MESHELEM_PRISM) ? 6 : 4);
			const int w = ((kind == MESHELEM_HEX) ? HEX2_W : (kind == MESHELEM_PRISM) ? PRISM2_W : ((P == 2) ? 6 : TET_HO_MAX));
			int32_t *nvx = (int32_t *)malloc((size_t)nall[s] * nv * sizeof(int32_t));
			int *ntag = (int *)malloc((size_t)nall[s] * sizeof(int));
			int32_t *nh = (int32_t *)malloc((size_t)nall[s] * w * sizeof(int32_t));
			int m = 0;
			for (int pass = 0; pass < 2; pass++) {
				const int cnt = ((pass == 0) ? nst[s] : PendN);
				for (int e = 0; e < cnt; e++) {
					int nen, tag;
					const int32_t *src;
					if (pass == 0) {
						nen = ho_stored(kind, e, nd);
						src = nd;
						tag = ((kind == MESHELEM_HEX) ? HexTag[e] : (kind == MESHELEM_PRISM) ? PrismTag[e] : TetTag[e]);
					}
					else {
						if (PendKind[e] != kind) continue;
						nen = PendNen[e];
						src = &PendNd[(size_t)e * ELEM_NMAX];
						tag = PendTag[e];
					}
					if (nen == nt) { for (int l = 0; l < nt; l++) out[l] = src[l]; }
					else { ho_elem3(kind, nen, src, P, out, 0); nraise++; }
					for (int l = 0; l < nv; l++) nvx[(m * nv) + l] = out[l];
					for (int l = 0; l < w; l++) nh[((size_t)m * w) + l] = (((nv + l) < nt) ? out[nv + l] : -1);
					ntag[m] = tag;
					m++;
				}
			}
			if (kind == MESHELEM_HEX) {
				free(Hex); free(HexTag); free(Hex2);
				Hex = nvx; HexTag = ntag; Hex2 = nh; NHex = m; HexNen = nt;
			}
			else if (kind == MESHELEM_PRISM) {
				free(Prism); free(PrismTag); free(Prism2);
				Prism = nvx; PrismTag = ntag; Prism2 = nh; NPrism = m; PrismNen = nt;
			}
			else {
				free(Tet); free(TetTag); free(Tet2); free(Tet3);
				Tet = nvx; TetTag = ntag; Tet2 = NULL; Tet3 = NULL;
				if (P == 2) Tet2 = nh;
				else Tet3 = nh;
				NTet = m;
			}
		}
		TetOrder = P;
		free(PendKind); free(PendNen); free(PendTag); free(PendNd);
		PendKind = PendNen = PendTag = NULL; PendNd = NULL; PendN = PendCap = 0;
	}

	// (3) 三角形を P 次に (断面 2 次元では低い次数の要素を上げる。3 次元では電極の面をすべて要素の面の節点にする)
	HoMiss = 0;
	for (int t = 0; t < NTri; t++) {
		const int q = ho_tri_nodes(t, nd);
		if ((q == P) && !dim3) continue;
		ho_tri(q, nd, P, out, 0, !dim3);
		if (HoMiss) {
			printf("*** mesh : triangle %d is not a face of any element (cannot raise it to order %d)\n", t + 1, P);
			free(HoTab); HoTab = NULL; HoCap = 0;
			return 1;
		}
		ho_tri_store(t, P, out);
		if (trimix) nraise++;
	}
	if (trimix) TetOrder = P;

	// (4) 四角形 (電極の面) の節点をすべて要素の面の節点で埋め直す (固定する節点の集合なので並びは問わない)
	if (dim3 && (NQuad > 0)) {
		for (int t = 0; t < NQuad; t++) {
			const int32_t c[4] = {Quad[t * 4], Quad[(t * 4) + 1], Quad[(t * 4) + 2], Quad[(t * 4) + 3]};
			int l = 0;
			for (int u = 0; u <= P; u++) {
				for (int v = 0; v <= P; v++) {
					const int bu = ((u == 0) || (u == P)), bv = ((v == 0) || (v == P));
					if (bu && bv) continue;
					int32_t k[9];
					if (bu || bv) {
						// 辺 : (u, v) の辺の両端
						int32_t g[2];
						int wv[2];
						if (bv) { g[0] = c[(v == 0) ? 0 : 3]; g[1] = c[(v == 0) ? 1 : 2]; wv[0] = P - u; wv[1] = u; }
						else { g[0] = c[(u == 0) ? 0 : 1]; g[1] = c[(u == 0) ? 3 : 2]; wv[0] = P - v; wv[1] = v; }
						ho_key(g, wv, 2, k);
					}
					else ho_qkey(c, u, v, P, k);
					const int32_t id = ho_find(k, -1);
					if (id < 0) {
						printf("*** mesh : quadrilateral %d is not a face of any element (cannot raise it to order %d)\n", t + 1, P);
						free(HoTab); HoTab = NULL; HoCap = 0;
						return 1;
					}
					Quad2[(t * QUAD2_W) + (l++)] = id;
				}
			}
			while (l < QUAD2_W) Quad2[(t * QUAD2_W) + (l++)] = -1;
		}
	}
	free(HoTab); HoTab = NULL; HoCap = 0;
	MeshRaisedElem = nraise;
	MeshRaisedFrom = qmin;
	MeshRaisedKinds = !onlytet;
	return 0;
}


static int prism2_grad(const int32_t *nd, int nen, double u, double v, double ze,
	double g[PRISMP_NMAX][3], double *det)
{
	double dn[PRISMP_NMAX][3];
	prism2_dshape(nen, u, v, ze, dn);

	double jm[3][3];
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) jm[i][j] = 0;
	}
	for (int a = 0; a < nen; a++) {
		const int32_t w = nd[a];
		const double p[3] = {Xp[w], Yp[w], Zp[w]};
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) jm[i][j] += p[i] * dn[a][j];
		}
	}
	const double d = (jm[0][0] * ((jm[1][1] * jm[2][2]) - (jm[1][2] * jm[2][1])))
	               - (jm[0][1] * ((jm[1][0] * jm[2][2]) - (jm[1][2] * jm[2][0])))
	               + (jm[0][2] * ((jm[1][0] * jm[2][1]) - (jm[1][1] * jm[2][0])));
	if (d == 0) return 1;

	double ji[3][3];
	ji[0][0] = ((jm[1][1] * jm[2][2]) - (jm[1][2] * jm[2][1])) / d;
	ji[0][1] = ((jm[0][2] * jm[2][1]) - (jm[0][1] * jm[2][2])) / d;
	ji[0][2] = ((jm[0][1] * jm[1][2]) - (jm[0][2] * jm[1][1])) / d;
	ji[1][0] = ((jm[1][2] * jm[2][0]) - (jm[1][0] * jm[2][2])) / d;
	ji[1][1] = ((jm[0][0] * jm[2][2]) - (jm[0][2] * jm[2][0])) / d;
	ji[1][2] = ((jm[0][2] * jm[1][0]) - (jm[0][0] * jm[1][2])) / d;
	ji[2][0] = ((jm[1][0] * jm[2][1]) - (jm[1][1] * jm[2][0])) / d;
	ji[2][1] = ((jm[0][1] * jm[2][0]) - (jm[0][0] * jm[2][1])) / d;
	ji[2][2] = ((jm[0][0] * jm[1][1]) - (jm[0][1] * jm[1][0])) / d;

	for (int a = 0; a < nen; a++) {
		for (int i = 0; i < 3; i++) {
			g[a][i] = (dn[a][0] * ji[0][i]) + (dn[a][1] * ji[1][i]) + (dn[a][2] * ji[2][i]);
		}
	}
	*det = d;

	return 0;
}


// 2 次角柱の要素行列 (面内 Duffy 9 点 x ζ 3 点 Gauss)。vol には体積を返す
static int prism2_element(int e, const double c[6], double ke[PRISMP_NMAX][PRISMP_NMAX], double *vol)
{
	int32_t nd[PRISMP_NMAX];
	const int nen = prism_nodes(e, nd);

	// 2 次は面内 Duffy 9 点 x ζ 3 点、3 次は 16 点 x 4 点
	double lam[QTRI_MAX][3], wq[QTRI_MAX];
	int ntq = NQTRI, ngz = 3;
	if (prismp_order(nen)) { ngz = prismp_order(nen) + 1; ntq = tri_quad_pub(ngz, lam, wq); }
	else tri_quad(lam, wq);
	const double *gz, *gwz;
	ho_gauss(ngz, &gz, &gwz);

	memset(ke, 0, sizeof(double) * PRISMP_NMAX * PRISMP_NMAX);
	double v = 0, sgn = 0;
	for (int q = 0; q < ntq; q++) {
	for (int k = 0; k < ngz; k++) {
		double g[PRISMP_NMAX][3], det;
		if (prism2_grad(nd, nen, lam[q][1], lam[q][2], gz[k], g, &det)) return 1;
		if (sgn == 0) sgn = ((det > 0) ? 1 : -1);
		else if (det * sgn <= 0) return 1;
		const double dw = wq[q] * gwz[k] * fabs(det);
		v += dw;
		for (int l = 0; l < nen; l++) {
			for (int m = 0; m < nen; m++) {
				ke[l][m] += dw * ((c[0] * g[l][0] * g[m][0])
				                + (c[1] * g[l][1] * g[m][1])
				                + (c[2] * g[l][2] * g[m][2])
				                + (c[3] * ((g[l][0] * g[m][1]) + (g[l][1] * g[m][0])))
				                + (c[4] * ((g[l][1] * g[m][2]) + (g[l][2] * g[m][1])))
				                + (c[5] * ((g[l][2] * g[m][0]) + (g[l][0] * g[m][2]))));
			}
		}
	}
	}
	*vol = v;

	return 0;
}


// 要素中心 (u = v = 1/3, ζ = 0) での勾配 (場の出力用)
// 角柱の局所節点 a の格子位置 : 三角形の多重指数 al[0..2] (和 p、λ0 = 1 − u − v, λ1 = u, λ2 = v) と ζ の位置 al[3] (0..p)。
// 戻り値は次数 p (1 次の 6 節点と完全 Lagrange の 18 / 40 / 75 / 126 節点。serendipity の 15 節点は 0)
int prism_lattice(int a, int al[4])
{
	const int p = prismp_order(PrismNen);
	if (p) {
		const signed char *s = prismp_a(p, a);
		for (int c = 0; c < 4; c++) al[c] = s[c];
		return p;
	}
	if (PrismNen == 6) {
		// 1 次 : 頂点 0, 1, 2 が下 (ζ = −1)、3, 4, 5 が上
		for (int c = 0; c < 3; c++) al[c] = ((c == (a % 3)) ? 1 : 0);
		al[3] = a / 3;
		return 1;
	}
	if (PrismNen != 18) return 0;
	// tri6 の基底番号 (頂点 0, 1, 2、辺 (0,1) (1,2) (2,0)) -> 多重指数
	static const signed char T6[6][3] = {{2, 0, 0}, {0, 2, 0}, {0, 0, 2}, {1, 1, 0}, {0, 1, 1}, {1, 0, 1}};
	for (int c = 0; c < 3; c++) al[c] = T6[PRISM2_T[a]][c];
	al[3] = PRISM2_Z[a] + 1;
	return 2;
}


int prism2_grad_center(int e, double g[PRISMP_NMAX][3], int *nen)
{
	int32_t nd[PRISMP_NMAX];
	double det;

	*nen = prism_nodes(e, nd);

	return prism2_grad(nd, *nen, 1.0 / 3, 1.0 / 3, 0, g, &det);
}


/*
ピラミッド (5 節点、四角形の底面 + 頂点)。

局所節点の並びは **Gmsh / VTK と同じ**「底面を反時計回り (0..3)、最後に
頂点 (4)」。形状関数は底面座標を (1-ζ) で縮めた**立方体座標** (u, v, ζ) で書く:

    ξ = u (1-ζ),  η = v (1-ζ)      (u, v ∈ [-1,1], ζ ∈ [0,1])
    N_a = (1 + u_a u)(1 + v_a v)(1-ζ)/4  (底面 a = 0..3),   N_4 = ζ

これはよく知られた有理形状関数 [(1-ζ) + u_a ξ + v_a η + u_a v_a ξη/(1-ζ)]/4
と同じものだが、立方体座標では**多項式**になる。頂点 (ζ = 1) の特異性は
積分点がそこを踏まないので現れない。トレースは底面の四角形で双 1 次
(六面体と適合)、4 枚の三角形面で 1 次 (四面体と適合) — 例えば u = +1 の面では
β = (1-ζ)(1-v)/2, γ = (1-ζ)(1+v)/2 とおくと N_1 = β, N_2 = γ, N_4 = ζ で
β + γ + ζ = 1、つまり N がその三角形の**面積座標に厳密に退化**する。
これが「六面体 - ピラミッド - 四面体」の遷移が適合する理由の全て。

積分は**立方体座標のまま** 2x2x2 Gauss (ζ は [0,1] の 2 点)。ピラミッド座標に
戻さないのが要点で、dV = |det(∂x/∂(u,v,ζ))| du dv dζ の中に錐の縮み (1-ζ)² が
自動的に入るため、特異な重み関数も専用の積分則も要らない。det J は
どの 5 節点ピラミッドでも各方向 2 次までの多項式なので**体積は常に厳密**。
剛性は「基準ピラミッド (底面 [-1,1]²、頂点が底面中心の真上) のアフィン像」で
厳密 (六面体の「直方体で厳密」と同格の主張)。底面が平行四辺形でない・頂点が
ずれた一般の形では有理式になるので厳密ではない (六面体と同じ)。

六面体・角柱と同じく、ヤコビアンは**積分点毎に** 5 節点すべてから作る。
ピラミッドでは (1-ζ)² の因子のため J が要素内で一定になることは**決してない**
ので、「中心で 1 回だけ評価する」手抜きはどんな形でも答えを壊す
(六面体と違い、この誤りを隠す特別な形が存在しない)。
*/
static const signed char PYR_SGN[4][2] = {
	{-1, -1}, {+1, -1}, {+1, +1}, {-1, +1}
};


// 立方体座標 (u, v, ζ) での物理勾配 g[a][i] = ∂N_a/∂x_i と det J
static int pyr_shape_at(const int32_t *nd, double u, double v, double ze,
	double g[5][3], double *det)
{
	// 立方体座標での微分
	double dn[5][3];
	for (int a = 0; a < 4; a++) {
		const double su = PYR_SGN[a][0], sv = PYR_SGN[a][1];
		dn[a][0] = su * (1 + (sv * v)) * (1 - ze) / 4;
		dn[a][1] = sv * (1 + (su * u)) * (1 - ze) / 4;
		dn[a][2] = -(1 + (su * u)) * (1 + (sv * v)) / 4;
	}
	dn[4][0] = 0;
	dn[4][1] = 0;
	dn[4][2] = 1;

	double jm[3][3];
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) jm[i][j] = 0;
	}
	for (int a = 0; a < 5; a++) {
		const int32_t w = nd[a];
		const double p[3] = {Xp[w], Yp[w], Zp[w]};
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) jm[i][j] += p[i] * dn[a][j];
		}
	}
	const double d = (jm[0][0] * ((jm[1][1] * jm[2][2]) - (jm[1][2] * jm[2][1])))
	               - (jm[0][1] * ((jm[1][0] * jm[2][2]) - (jm[1][2] * jm[2][0])))
	               + (jm[0][2] * ((jm[1][0] * jm[2][1]) - (jm[1][1] * jm[2][0])));
	if (d == 0) return 1;

	double ji[3][3];
	ji[0][0] = ((jm[1][1] * jm[2][2]) - (jm[1][2] * jm[2][1])) / d;
	ji[0][1] = ((jm[0][2] * jm[2][1]) - (jm[0][1] * jm[2][2])) / d;
	ji[0][2] = ((jm[0][1] * jm[1][2]) - (jm[0][2] * jm[1][1])) / d;
	ji[1][0] = ((jm[1][2] * jm[2][0]) - (jm[1][0] * jm[2][2])) / d;
	ji[1][1] = ((jm[0][0] * jm[2][2]) - (jm[0][2] * jm[2][0])) / d;
	ji[1][2] = ((jm[0][2] * jm[1][0]) - (jm[0][0] * jm[1][2])) / d;
	ji[2][0] = ((jm[1][0] * jm[2][1]) - (jm[1][1] * jm[2][0])) / d;
	ji[2][1] = ((jm[0][1] * jm[2][0]) - (jm[0][0] * jm[2][1])) / d;
	ji[2][2] = ((jm[0][0] * jm[1][1]) - (jm[0][1] * jm[1][0])) / d;

	for (int a = 0; a < 5; a++) {
		for (int i = 0; i < 3; i++) {
			g[a][i] = (dn[a][0] * ji[0][i]) + (dn[a][1] * ji[1][i]) + (dn[a][2] * ji[2][i]);
		}
	}
	*det = d;

	return 0;
}


// 積分点 q (0..7) での評価。q < 0 なら要素中心 (u = v = 0, ζ = 1/4 = 重心の高さ)
static int pyr_shape(const int32_t *nd, int q, double g[5][3], double *det)
{
	const double gp = 1.0 / sqrt(3.0);

	if (q < 0) return pyr_shape_at(nd, 0, 0, 0.25, g, det);

	return pyr_shape_at(nd, gp * PYR_SGN[q % 4][0], gp * PYR_SGN[q % 4][1],
		((q < 4) ? (0.5 - (0.5 * gp)) : (0.5 + (0.5 * gp))), g, det);
}


// ---- 1 次の 3 次元要素の積分点 (熱連成の発熱・質量・要素の温度で使う) ----
//
// 剛性と**同じ積分点・同じ重み**で N・∇N・重み w を返す (そうすると Σ_i ∫σ|∇φ|²N_i =
// φᵀKφ が丸めを除いて厳密になる)。四面体 1 点 (重心、∇N は一定)、六面体 2x2x2、
// 角柱 6 点、ピラミッド 2x2x2 (立方体座標)。2 次の六面体 3x3x3、2 次の角柱 9 x 3 点。
// 高次の四面体 e の局所の面 (頂点 v_opp の反対側) 上の点 (面の重心座標 mu、面の頂点の
// 局所番号 fa[3]) での形状関数 nv[] と、面の 2 つの接線方向の ∂x/∂μ1, ∂x/∂μ2 (面積要素 =
// |t1 × t2|)。熱連成の露出面 (曲がった面) の積分に使う。戻り値は節点数
int tet_face_shape(int e, const int fa[3], const double mu[3], double nv[TET_NMAX],
	double t1[3], double t2[3])
{
	int32_t nd[TET_NMAX];
	const int nen = tet_nodes(e, nd);
	double lam[4] = {0, 0, 0, 0};
	for (int k = 0; k < 3; k++) lam[fa[k]] = mu[k];
	for (int c = 0; c < 3; c++) { t1[c] = 0; t2[c] = 0; }
	for (int a = 0; a < nen; a++) {
		double dl[4];
		nv[a] = lagr_tet(TetOrder, tet_ho_a(a), lam, dl);
		const double d1 = dl[fa[1]] - dl[fa[0]], d2 = dl[fa[2]] - dl[fa[0]];
		const double p[3] = {Xp[nd[a]], Yp[nd[a]], Zp[nd[a]]};
		for (int c = 0; c < 3; c++) {
			t1[c] += d1 * p[c];
			t2[c] += d2 * p[c];
		}
	}
	return nen;
}


// 三角形の Duffy 積分点 (各方向 n 点、n = 3..6)。mu は重心座標、重みの合計は 1/2
int tri_quad_pub(int n, double mu[][3], double *wq)
{
	double gp[GL_NMAX], gw[GL_NMAX];
	gauss01(n, gp, gw);
	int q = 0;
	for (int a = 0; a < n; a++) {
	for (int b = 0; b < n; b++) {
		const double u = gp[a], v = gp[b];
		mu[q][1] = u;
		mu[q][2] = v * (1 - u);
		mu[q][0] = 1 - mu[q][1] - mu[q][2];
		wq[q] = gw[a] * gw[b] * (1 - u);
		q++;
	}
	}
	return q;
}


// 2 次の六面体・角柱 e の面 (角の節点 cn[0..k-1]、k = 3 / 4) の積分点 q (0..8) での親の要素の
// 形状関数 nv[] と面の重み w (曲がった面の上の面積要素 |∂x/∂s × ∂x/∂t| を含む)。熱連成の
// 露出面の積分に使う。面は角の局所座標から決める: 六面体は 1 つの座標が ±1 の面 (残り 2 方向
// 3 点 Gauss)、角柱は ζ = ±1 の三角形 (Duffy 9 点) か三角形の辺 x ζ の四角形 (3 x 3 Gauss)。
// 戻り値は節点数 (失敗は 0)
// elem3d_face_ip の積分点の数 (2 次 9、3 次 16)
int elem3d_face_nip(int e)
{
	int32_t nd[ELEM_NMAX];
	const int nen = elem3d_nodes(e, nd);
	const int p = (hexp_order(nen) ? hexp_order(nen) : prismp_order(nen));
	return (p ? ((p + 1) * (p + 1)) : 9);
}


int elem3d_face_ip(int e, const int32_t *cn, int k, int q, double nv[ELEM_NMAX], double *w)
{
	const int kind = elem3d_kind(e);
	int32_t nd[ELEM_NMAX];
	const int nen = elem3d_nodes(e, nd);
	const int ncor = ((kind == MESHELEM_HEX) ? 8 : 6);
	// 2 次は 9 点 (各方向 3 点)、3 次は 16 点 (各方向 4 点)
	const int pho = ((kind == MESHELEM_HEX) ? hexp_order(nen) : prismp_order(nen));
	const int o3 = (pho >= 3);
	const int ng = (o3 ? (pho + 1) : 3);
	const double *gx, *gw;
	ho_gauss(ng, &gx, &gw);
	if ((q < 0) || (q >= (ng * ng))) return 0;
	if (((kind == MESHELEM_HEX) && ((HexNen <= 8) || (k != 4)))
	 || ((kind == MESHELEM_PRISM) && (PrismNen <= 6))
	 || ((kind != MESHELEM_HEX) && (kind != MESHELEM_PRISM))) return 0;
	int loc[4];
	for (int m = 0; m < k; m++) {
		loc[m] = -1;
		for (int a = 0; a < ncor; a++) {
			if (nd[a] == cn[m]) loc[m] = a;
		}
		if (loc[m] < 0) return 0;
	}
	double dn[HEXP_NMAX][3], d1[3] = {0, 0, 0}, d2[3] = {0, 0, 0}, wt;
	if (kind == MESHELEM_HEX) {
		int d = -1;
		for (int dd = 0; dd < 3; dd++) {
			int same = 1;
			for (int m = 1; m < 4; m++) same &= (HEX2_XI[loc[m]][dd] == HEX2_XI[loc[0]][dd]);
			if (same) d = dd;
		}
		if (d < 0) return 0;
		const int da = (d + 1) % 3, db = (d + 2) % 3;
		double xi[3];
		xi[d] = HEX2_XI[loc[0]][d];
		xi[da] = gx[q / ng];
		xi[db] = gx[q % ng];
		wt = gw[q / ng] * gw[q % ng];
		d1[da] = 1;
		d2[db] = 1;
		hex2_dshape(nen, xi[0], xi[1], xi[2], dn);
		hex2_nval(nen, xi[0], xi[1], xi[2], nv);
	}
	else {
		double u, v, ze, dpr[PRISMP_NMAX][3], vv[PRISMP_NMAX];
		if (k == 3) {
			// ζ = ±1 の三角形
			for (int m = 1; m < 3; m++) {
				if (PRISM2_Z[loc[m]] != PRISM2_Z[loc[0]]) return 0;
			}
			double lam[QTRI_MAX][3], wq[QTRI_MAX];
			if (o3) tri_quad_pub(ng, lam, wq);
			else tri_quad(lam, wq);
			u = lam[q][1];
			v = lam[q][2];
			ze = PRISM2_Z[loc[0]];
			wt = wq[q];
			d1[0] = 1;
			d2[1] = 1;
		}
		else {
			// 三角形の辺 (λ_i → λ_j) x ζ の四角形
			const int i = PRISM2_T[loc[0]];
			int j = -1;
			for (int m = 1; m < 4; m++) {
				if (PRISM2_T[loc[m]] != i) j = PRISM2_T[loc[m]];
			}
			if (j < 0) return 0;
			const double sp = (gx[q / ng] + 1) / 2;
			double lam[3] = {0, 0, 0};
			lam[i] = 1 - sp;
			lam[j] = sp;
			u = lam[1];
			v = lam[2];
			ze = gx[q % ng];
			wt = gw[q / ng] * gw[q % ng] / 2;
			d1[0] = (double)((j == 1) - (i == 1));
			d1[1] = (double)((j == 2) - (i == 2));
			d2[2] = 1;
		}
		prism2_dshape(nen, u, v, ze, dpr);
		prism2_nval(nen, u, v, ze, vv);
		for (int a = 0; a < nen; a++) {
			nv[a] = vv[a];
			for (int c = 0; c < 3; c++) dn[a][c] = dpr[a][c];
		}
	}
	double t1[3] = {0, 0, 0}, t2[3] = {0, 0, 0};
	for (int a = 0; a < nen; a++) {
		const double s1 = (dn[a][0] * d1[0]) + (dn[a][1] * d1[1]) + (dn[a][2] * d1[2]);
		const double s2 = (dn[a][0] * d2[0]) + (dn[a][1] * d2[1]) + (dn[a][2] * d2[2]);
		const double p[3] = {Xp[nd[a]], Yp[nd[a]], Zp[nd[a]]};
		for (int c = 0; c < 3; c++) {
			t1[c] += s1 * p[c];
			t2[c] += s2 * p[c];
		}
	}
	const double cx = (t1[1] * t2[2]) - (t1[2] * t2[1]);
	const double cy = (t1[2] * t2[0]) - (t1[0] * t2[2]);
	const double cz = (t1[0] * t2[1]) - (t1[1] * t2[0]);
	*w = wt * sqrt((cx * cx) + (cy * cy) + (cz * cz));
	return nen;
}


// 3 次元要素 e の境界面 (角の全体節点 cn[k]、三角形 k = 3 / 四角形 k = 4 は巡回順) の媒介変数 (s, t) の点。
// 参照空間で三角形は 角0 + s (角1 − 角0) + t (角2 − 角0)、四角形は角の双 1 次 (面は参照空間で平らなので厳密)。
// 形状関数の値 nv、物理座標 x、接ベクトル t1 = ∂x/∂s、t2 = ∂x/∂t を返す (曲がった面の相互放射が使う)。
// 戻り値は節点数 (失敗は 0)
int elem3d_face_point(int e, const int32_t *cn, int k, double s, double t, double nv[ELEM_NMAX],
	double x[3], double t1[3], double t2[3])
{
	const int kind = elem3d_kind(e);
	int32_t nd[ELEM_NMAX];
	const int nen = elem3d_nodes(e, nd);
	const int ncor = ((kind == MESHELEM_HEX) ? 8 : (kind == MESHELEM_PRISM) ? 6 : 4);
	int loc[4];
	for (int m = 0; m < k; m++) {
		loc[m] = -1;
		for (int a = 0; a < ncor; a++) {
			if (nd[a] == cn[m]) loc[m] = a;
		}
		if (loc[m] < 0) return 0;
	}
	for (int c = 0; c < 3; c++) x[c] = t1[c] = t2[c] = 0;
	if (kind == MESHELEM_TET) {
		if (k != 3) return 0;
		const double mu[3] = {1 - s - t, s, t};
		if (TetOrder >= 2) {
			tet_face_shape(e, loc, mu, nv, t1, t2);
		}
		else {
			for (int a = 0; a < 4; a++) nv[a] = 0;
			for (int m = 0; m < 3; m++) nv[loc[m]] = mu[m];
			const double p0[3] = {Xp[nd[loc[0]]], Yp[nd[loc[0]]], Zp[nd[loc[0]]]};
			t1[0] = Xp[nd[loc[1]]] - p0[0]; t1[1] = Yp[nd[loc[1]]] - p0[1]; t1[2] = Zp[nd[loc[1]]] - p0[2];
			t2[0] = Xp[nd[loc[2]]] - p0[0]; t2[1] = Yp[nd[loc[2]]] - p0[1]; t2[2] = Zp[nd[loc[2]]] - p0[2];
		}
		for (int a = 0; a < nen; a++) {
			x[0] += nv[a] * Xp[nd[a]]; x[1] += nv[a] * Yp[nd[a]]; x[2] += nv[a] * Zp[nd[a]];
		}
		return nen;
	}
	if (((kind == MESHELEM_HEX) && ((nen <= 8) || (k != 4))) || ((kind == MESHELEM_PRISM) && (nen <= 6))
	 || ((kind != MESHELEM_HEX) && (kind != MESHELEM_PRISM))) return 0;
	// 角の参照座標 (六面体 ξ、角柱 (u, v, ζ))
	double rc[4][3];
	for (int m = 0; m < k; m++) {
		if (kind == MESHELEM_HEX) {
			for (int c = 0; c < 3; c++) rc[m][c] = HEX2_XI[loc[m]][c];
		}
		else {
			const int iv = PRISM2_T[loc[m]];
			rc[m][0] = (iv == 1); rc[m][1] = (iv == 2); rc[m][2] = PRISM2_Z[loc[m]];
		}
	}
	double r[3], dr1[3], dr2[3];
	for (int c = 0; c < 3; c++) {
		if (k == 3) {
			r[c] = rc[0][c] + (s * (rc[1][c] - rc[0][c])) + (t * (rc[2][c] - rc[0][c]));
			dr1[c] = rc[1][c] - rc[0][c];
			dr2[c] = rc[2][c] - rc[0][c];
		}
		else {
			r[c] = ((1 - s) * (1 - t) * rc[0][c]) + (s * (1 - t) * rc[1][c]) + (s * t * rc[2][c]) + ((1 - s) * t * rc[3][c]);
			dr1[c] = ((1 - t) * (rc[1][c] - rc[0][c])) + (t * (rc[2][c] - rc[3][c]));
			dr2[c] = ((1 - s) * (rc[3][c] - rc[0][c])) + (s * (rc[2][c] - rc[1][c]));
		}
	}
	double dn[HEXP_NMAX][3];
	if (kind == MESHELEM_HEX) {
		hex2_nval(nen, r[0], r[1], r[2], nv);
		hex2_dshape(nen, r[0], r[1], r[2], dn);
	}
	else {
		double vv[PRISMP_NMAX], dpr[PRISMP_NMAX][3];
		prism2_nval(nen, r[0], r[1], r[2], vv);
		prism2_dshape(nen, r[0], r[1], r[2], dpr);
		for (int a = 0; a < nen; a++) {
			nv[a] = vv[a];
			for (int c = 0; c < 3; c++) dn[a][c] = dpr[a][c];
		}
	}
	for (int a = 0; a < nen; a++) {
		const double p[3] = {Xp[nd[a]], Yp[nd[a]], Zp[nd[a]]};
		const double s1 = (dn[a][0] * dr1[0]) + (dn[a][1] * dr1[1]) + (dn[a][2] * dr1[2]);
		const double s2 = (dn[a][0] * dr2[0]) + (dn[a][1] * dr2[1]) + (dn[a][2] * dr2[2]);
		for (int c = 0; c < 3; c++) {
			x[c] += nv[a] * p[c];
			t1[c] += s1 * p[c];
			t2[c] += s2 * p[c];
		}
	}
	return nen;
}


// 面の媒介変数 (s, t) の積分点 : 三角形は Duffy (n x n、重みの和 1/2)、四角形は [0,1]² の Gauss (n x n、和 1)。
// n は 3..7。戻り値は点の数
int face_quad_pub(int k, int n, double st[][2], double *w)
{
	if (k == 3) {
		double mu[QTRI_MAX][3], wq[QTRI_MAX];
		const int nq = tri_quad_pub(n, mu, wq);
		for (int q = 0; q < nq; q++) { st[q][0] = mu[q][1]; st[q][1] = mu[q][2]; w[q] = wq[q]; }
		return nq;
	}
	double gp[GL_NMAX], gw[GL_NMAX];
	gauss01(n, gp, gw);
	int q = 0;
	for (int a = 0; a < n; a++) {
		for (int b = 0; b < n; b++) {
			st[q][0] = gp[a]; st[q][1] = gp[b]; w[q] = gw[a] * gw[b];
			q++;
		}
	}
	return q;
}


// 高次の四面体の任意の参照点 (重心座標 lam) での形状関数の値と物理空間の勾配、det J
// (elem3d_ip の高次の四面体の分岐を積分点でない点で使う。2 次の A の場の出力が重心で使う)
int tet_ho_point(int e, const double lam[4], double nv[ELEM_NMAX], double g[ELEM_NMAX][3], double *det)
{
	if (TetOrder < 2) return 1;
	int32_t nd[TET_NMAX];
	const int nen = tet_nodes(e, nd);
	double gn[TET_NMAX][3];
	if (tet20_grad(nd, lam, gn, det)) return 1;
	for (int a = 0; a < nen; a++) {
		double dl[4];
		nv[a] = lagr_tet(TetOrder, tet_ho_a(a), lam, dl);
		for (int c = 0; c < 3; c++) g[a][c] = gn[a][c];
	}
	return 0;
}


int elem3d_nip(int e)
{
	const int kind = elem3d_kind(e);
	// 高次の四面体 : Duffy + 各方向 p + 2 点 (質量 N_i N_j の 2p 次が直線要素で厳密)
	if (kind == MESHELEM_TET) return ((TetOrder >= 2) ? ((TetOrder + 2) * (TetOrder + 2) * (TetOrder + 2)) : 1);
	// 2 次の六面体 3x3x3 Gauss、2 次の角柱 面内 Duffy 9 点 x ζ 3 点 (剛性の組み立てと同じ則)
	if (kind == MESHELEM_HEX) {
		const int p = hexp_order(HexNen);
		return (p ? ((p + 1) * (p + 1) * (p + 1)) : (HexNen > 8) ? 27 : 8);
	}
	if (kind == MESHELEM_PRISM) {
		const int p = prismp_order(PrismNen);
		return (p ? ((p + 1) * (p + 1) * (p + 1)) : (PrismNen > 6) ? (NQTRI * 3) : 6);
	}
	return 8;
}


int elem3d_ip(int e, int q, double nv[ELEM_NMAX], double g[ELEM_NMAX][3], double *w)
{
	const int kind = elem3d_kind(e);
	if ((kind == MESHELEM_TET) && (TetOrder >= 2)) {
		// 高次の四面体 (等パラメトリック) : 積分点を次数ごとに一度だけ作っておく
		static double lq[NQTETMAX][4], wq[NQTETMAX];
		static int nqo = 0, ordo = 0;
		if (ordo != TetOrder) {
			nqo = tet_quad_n(TetOrder + 2, lq, wq);
			ordo = TetOrder;
		}
		if ((q < 0) || (q >= nqo)) return 1;
		int32_t nd[TET_NMAX];
		const int nen = tet_nodes(e, nd);
		double gn[TET_NMAX][3], det;
		if (tet20_grad(nd, lq[q], gn, &det)) return 1;
		for (int a = 0; a < nen; a++) {
			double dl[4];
			nv[a] = lagr_tet(TetOrder, tet_ho_a(a), lq[q], dl);
			for (int c = 0; c < 3; c++) g[a][c] = gn[a][c];
		}
		*w = wq[q] * fabs(det);
		return 0;
	}
	if (kind == MESHELEM_TET) {
		double gg[4][3], vol;
		if (tet_grad_pub(&Tet[e * 4], gg, &vol)) return 1;
		for (int a = 0; a < 4; a++) {
			nv[a] = 0.25;
			for (int c = 0; c < 3; c++) g[a][c] = gg[a][c];
		}
		*w = vol;
		return 0;
	}
	if ((kind == MESHELEM_HEX) && (HexNen > 8)) {
		const int ng = (hexp_order(HexNen) ? (hexp_order(HexNen) + 1) : 3);
		const double *gx, *gw;
		ho_gauss(ng, &gx, &gw);
		if ((q < 0) || (q >= (ng * ng * ng))) return 1;
		int32_t nd[ELEM_NMAX];
		const int nen = hex_nodes(e - NTet, nd);
		const int qi = q / (ng * ng), qj = (q / ng) % ng, qk = q % ng;
		const double x = gx[qi], y = gx[qj], z = gx[qk];
		double gg[HEXP_NMAX][3], vv[HEXP_NMAX], det;
		if (hex2_grad(nd, nen, x, y, z, gg, &det)) return 1;
		hex2_nval(nen, x, y, z, vv);
		for (int a = 0; a < nen; a++) {
			nv[a] = vv[a];
			for (int c = 0; c < 3; c++) g[a][c] = gg[a][c];
		}
		*w = gw[qi] * gw[qj] * gw[qk] * fabs(det);
		return 0;
	}
	if ((kind == MESHELEM_PRISM) && (PrismNen > 6)) {
		const int pp = prismp_order(PrismNen);
		const int ngz = (pp ? (pp + 1) : 3), ntq = (pp ? (ngz * ngz) : NQTRI);
		const double *gz, *gwz;
		ho_gauss(ngz, &gz, &gwz);
		if ((q < 0) || (q >= (ntq * ngz))) return 1;
		double lam[QTRI_MAX][3], wq[QTRI_MAX];
		if (pp) tri_quad_pub(ngz, lam, wq);
		else tri_quad(lam, wq);
		int32_t nd[PRISMP_NMAX];
		const int nen = prism_nodes(e - NTet - NHex, nd);
		const int iq = q / ngz, k = q % ngz;
		double gg[PRISMP_NMAX][3], vv[PRISMP_NMAX], det;
		if (prism2_grad(nd, nen, lam[iq][1], lam[iq][2], gz[k], gg, &det)) return 1;
		prism2_nval(nen, lam[iq][1], lam[iq][2], gz[k], vv);
		for (int a = 0; a < nen; a++) {
			nv[a] = vv[a];
			for (int c = 0; c < 3; c++) g[a][c] = gg[a][c];
		}
		*w = wq[iq] * gwz[k] * fabs(det);
		return 0;
	}
	if (kind == MESHELEM_HEX) {
		const int32_t *nd = &Hex[(e - NTet) * 8];
		double det;
		if (hex_shape(nd, q, g, &det)) return 1;
		const double gp = 1.0 / sqrt(3.0);
		const double xi = gp * HEX_SGN[q][0], eta = gp * HEX_SGN[q][1], ze = gp * HEX_SGN[q][2];
		for (int a = 0; a < 8; a++) {
			nv[a] = (1 + (HEX_SGN[a][0] * xi)) * (1 + (HEX_SGN[a][1] * eta))
			      * (1 + (HEX_SGN[a][2] * ze)) / 8;
		}
		*w = fabs(det);
		return 0;
	}
	if (kind == MESHELEM_PRISM) {
		const int32_t *nd = &Prism[(e - NTet - NHex) * 6];
		double gg[6][3], det;
		if (prism_shape(nd, q, gg, &det)) return 1;
		const double gp = 1.0 / sqrt(3.0);
		const double u = PRISM_TU[q % 3], v = PRISM_TV[q % 3], ze = ((q < 3) ? -gp : gp);
		const double lam[3] = {1 - u - v, u, v};
		for (int a = 0; a < 6; a++) {
			nv[a] = lam[PRISM_TRI[a]] * (1 + (PRISM_SGN[a] * ze)) / 2;
			for (int c = 0; c < 3; c++) g[a][c] = gg[a][c];
		}
		*w = fabs(det) / 6;
		return 0;
	}
	{
		const int32_t *nd = &Pyr[(e - NTet - NHex - NPrism) * 5];
		double gg[5][3], det;
		if (pyr_shape(nd, q, gg, &det)) return 1;
		const double gp = 1.0 / sqrt(3.0);
		const double u = gp * PYR_SGN[q % 4][0], v = gp * PYR_SGN[q % 4][1];
		const double ze = ((q < 4) ? (0.5 - (0.5 * gp)) : (0.5 + (0.5 * gp)));
		for (int a = 0; a < 4; a++) {
			nv[a] = (1 + (PYR_SGN[a][0] * u)) * (1 + (PYR_SGN[a][1] * v)) * (1 - ze) / 4;
		}
		nv[4] = ze;
		for (int a = 0; a < 5; a++) {
			for (int c = 0; c < 3; c++) g[a][c] = gg[a][c];
		}
		*w = 0.5 * fabs(det);
		return 0;
	}
}


// 3x3x3 Gauss で積分した体積 (2x2x2 則の検算に使う独立な計算)
static double pyr_volume_hi(const int32_t *nd)
{
	static const double gx[3] = {-0.7745966692414834, 0.0, 0.7745966692414834};
	static const double gw[3] = {5.0 / 9, 8.0 / 9, 5.0 / 9};
	double v = 0;

	for (int i = 0; i < 3; i++) {
	for (int j = 0; j < 3; j++) {
	for (int k = 0; k < 3; k++) {
		double g[5][3], det;
		if (pyr_shape_at(nd, gx[i], gx[j], (1 + gx[k]) / 2, g, &det)) continue;
		// ζ を [0,1] に縮めるので ζ 方向の重みは gw/2
		v += gw[i] * gw[j] * (gw[k] / 2) * ((det > 0) ? det : -det);
	}
	}
	}

	return v;
}


// ピラミッドの要素行列。vol には体積 Σ w |det J| を返す
static int pyr_element(const int32_t *nd, const double c[6], double ke[5][5], double *vol)
{
	double v = 0;
	int sgn = 0;

	for (int l = 0; l < 5; l++) {
		for (int m = 0; m < 5; m++) ke[l][m] = 0;
	}

	for (int q = 0; q < 8; q++) {
		double g[5][3], det;
		if (pyr_shape(nd, q, g, &det)) return 1;
		const int sq = ((det > 0) ? 1 : -1);
		if (q == 0) sgn = sq;
		else if (sq != sgn) return 1;			// 要素内で符号が変わる = 裏返り
		// 重み : (u, v) 方向 1 x 1、ζ 方向 1/2 ([0,1] の 2 点 Gauss)
		const double w = 0.5 * ((det > 0) ? det : -det);
		v += w;
		for (int l = 0; l < 5; l++) {
			for (int m = 0; m < 5; m++) {
				ke[l][m] += w * ((c[0] * g[l][0] * g[m][0])
				               + (c[1] * g[l][1] * g[m][1])
				               + (c[2] * g[l][2] * g[m][2])
				               + (c[3] * ((g[l][0] * g[m][1]) + (g[l][1] * g[m][0])))
				               + (c[4] * ((g[l][1] * g[m][2]) + (g[l][2] * g[m][1])))
				               + (c[5] * ((g[l][2] * g[m][0]) + (g[l][0] * g[m][2]))));
			}
		}
	}
	*vol = v;

	return 0;
}


int pyr_grad_center(int e, double g[5][3])
{
	double det;

	return pyr_shape(&Pyr[e * 5], -1, g, &det);
}


/*
混在格子の自己検証 (線形場の恒等式)。

φ = a・r は四面体・六面体・角柱のいずれでも厳密に補間できるので、種別が
混ざっていても φᵀKφ = (aᵀCa) V が厳密に成り立つ。**種別をまたぐ面で節点が
共有されていない (格子が割れている) と、この恒等式が破れる**ので、混在で
一番危ない「面の不整合」がここで落ちる。
*/
/*
高次の混在格子 (同じ次数 p の四面体・六面体・角柱) の自己検証 : 種別ごとの要素の行列 (各種別の純粋な格子で
p 次の場の恒等式を通したもの) を組み、線形場 φ = a・r で φᵀKφ = aᵀCa V (どんな形でも成立)。
組み立ての振り分け (種別ごとの添字のずらし方・次数ごとの要素の関数) の検査。界面の適合は位相
(elem_face_check) と曲がった同軸の収束で見る
*/
static int nodal_test_mixed_ho(FILE *fp_log)
{
	int ierr = 0;
	fprintf(fp_log, "\n=== nodal element (mixed, order %d) self test ===\n", TetOrder);
	fprintf(fp_log, "  nodes = %d, tetrahedra = %d, hexahedra = %d (%d-node), prisms = %d (%d-node)\n",
		NNode, NTet, NHex, HexNen, NPrism, PrismNen);
	const double c[6] = {2.0, 3.0, 1.5, 0.4, 0.3, 0.2};
	const double cm[3][3] = {{c[0], c[3], c[5]}, {c[3], c[1], c[4]}, {c[5], c[4], c[2]}};
	const double av[3] = {0.7, -1.3, 0.9};
	double aca = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) aca += av[i] * cm[i][j] * av[j];
	}
	crs_t A;
	crs_alloc(&A);
	crs_zero(&A);
	double vol = 0, vk[3] = {0, 0, 0};
	const int ne = elem3d_count();
	static double ke[ELEM_NMAX][ELEM_NMAX];
	static double kt[TET_NMAX][TET_NMAX], kt2[10][10], kh[HEXP_NMAX][HEXP_NMAX], kp[PRISMP_NMAX][PRISMP_NMAX];
	for (int e = 0; e < ne; e++) {
		const int kind = elem3d_kind(e);
		int32_t nd[ELEM_NMAX];
		const int nen = elem3d_nodes(e, nd);
		double v = 0;
		int bad = 0;
		if (kind == MESHELEM_TET) {
			if (TetOrder >= 3) {
				bad = tet20_element(e, c, kt, &v);
				for (int l = 0; l < nen; l++) for (int m = 0; m < nen; m++) ke[l][m] = kt[l][m];
			}
			else {
				bad = tet10_element(e, c, kt2, &v);
				for (int l = 0; l < nen; l++) for (int m = 0; m < nen; m++) ke[l][m] = kt2[l][m];
			}
		}
		else if (kind == MESHELEM_HEX) {
			bad = hex2_element(e - NTet, c, kh, &v);
			for (int l = 0; l < nen; l++) for (int m = 0; m < nen; m++) ke[l][m] = kh[l][m];
		}
		else {
			bad = prism2_element(e - NTet - NHex, c, kp, &v);
			for (int l = 0; l < nen; l++) for (int m = 0; m < nen; m++) ke[l][m] = kp[l][m];
		}
		if (bad) {
			fprintf(fp_log, "*** element %d is degenerate or inverted\n", e + 1);
			ierr = 1;
			continue;
		}
		vol += v;
		vk[(kind == MESHELEM_TET) ? 0 : (kind == MESHELEM_HEX) ? 1 : 2] += v;
		for (int l = 0; l < nen; l++) {
			for (int m = 0; m < nen; m++) {
				const int64_t p = crs_find(&A, nd[l], nd[m]);
				if (p >= 0) A.val[p] += ke[l][m];
			}
		}
	}
	fprintf(fp_log, "  volume = %.10e (tetrahedra %.4e, hexahedra %.4e, prisms %.4e)\n", vol, vk[0], vk[1], vk[2]);
	double *phi = (double *)malloc((size_t)NNode * sizeof(double));
	for (int i = 0; i < NNode; i++) phi[i] = (av[0] * Xp[i]) + (av[1] * Yp[i]) + (av[2] * Zp[i]);
	double q1 = 0;
	for (int i = 0; i < NNode; i++) q1 += phi[i] * crs_row_dot(&A, i, phi);
	const double want = aca * vol;
	const double err = fabs(q1 - want) / fabs(want);
	fprintf(fp_log, "  linear field : phi^T K phi = %.10e, closed form = %.10e, rel. error = %.2e\n", q1, want, err);
	if (!(err <= 1e-12)) {
		fprintf(fp_log, "*** the linear-field identity failed\n");
		ierr = 1;
	}
	free(phi);
	crs_free(&A);
	return ierr;
}


static int nodal_test_mixed(FILE *fp_log)
{
	int ierr = 0;

	fprintf(fp_log, "\n=== nodal element (mixed) self test ===\n");
	fprintf(fp_log, "  nodes = %d, tetrahedra = %d, hexahedra = %d, prisms = %d, "
		"pyramids = %d\n", NNode, NTet, NHex, NPrism, NPyr);

	const double c[6] = {2.0, 3.0, 1.5, 0.4, 0.3, 0.2};
	const double a[3] = {0.7, -1.3, 0.9};
	double aca = 0;
	{
		const double cm[3][3] = {{c[0], c[3], c[5]}, {c[3], c[1], c[4]}, {c[5], c[4], c[2]}};
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) aca += a[i] * cm[i][j] * a[j];
		}
	}

	crs_t A;
	crs_alloc(&A);
	crs_zero(&A);
	double vol = 0;
	const int ne = elem3d_count();
	for (int e = 0; e < ne; e++) {
		const int kind = elem3d_kind(e);
		int32_t nd[10];
		double ke[10][10], v = 0;
		int nen = 0;
		if (kind == MESHELEM_TET) {
			double g[4][3];
			const int32_t *t = &Tet[e * 4];
			if (tet_grad_pub(t, g, &v)) continue;
			nen = 4;
			for (int l = 0; l < 4; l++) nd[l] = t[l];
			const double cm[3][3] = {{c[0], c[3], c[5]}, {c[3], c[1], c[4]}, {c[5], c[4], c[2]}};
			for (int l = 0; l < 4; l++) {
				for (int m = 0; m < 4; m++) {
					double x = 0;
					for (int i = 0; i < 3; i++) {
						for (int j = 0; j < 3; j++) x += g[l][i] * cm[i][j] * g[m][j];
					}
					ke[l][m] = x * v;
				}
			}
		}
		else if (kind == MESHELEM_HEX) {
			double k8[8][8];
			const int32_t *t = &Hex[(e - NTet) * 8];
			if (hex_element(t, c, k8, &v)) continue;
			nen = 8;
			for (int l = 0; l < 8; l++) nd[l] = t[l];
			for (int l = 0; l < 8; l++) {
				for (int m = 0; m < 8; m++) ke[l][m] = k8[l][m];
			}
		}
		else if (kind == MESHELEM_PRISM) {
			double k6[6][6];
			const int32_t *t = &Prism[(e - NTet - NHex) * 6];
			if (prism_element(t, c, k6, &v)) continue;
			nen = 6;
			for (int l = 0; l < 6; l++) nd[l] = t[l];
			for (int l = 0; l < 6; l++) {
				for (int m = 0; m < 6; m++) ke[l][m] = k6[l][m];
			}
		}
		else {
			double k5[5][5];
			const int32_t *t = &Pyr[(e - NTet - NHex - NPrism) * 5];
			if (pyr_element(t, c, k5, &v)) continue;
			nen = 5;
			for (int l = 0; l < 5; l++) nd[l] = t[l];
			for (int l = 0; l < 5; l++) {
				for (int m = 0; m < 5; m++) ke[l][m] = k5[l][m];
			}
		}
		vol += v;
		for (int l = 0; l < nen; l++) {
			for (int m = 0; m < nen; m++) {
				const int64_t p = crs_find(&A, nd[l], nd[m]);
				if (p >= 0) A.val[p] += ke[l][m];
			}
		}
	}

	double *phi = (double *)malloc((size_t)NNode * sizeof(double));
	for (int i = 0; i < NNode; i++) {
		phi[i] = (a[0] * Xp[i]) + (a[1] * Yp[i]) + (a[2] * Zp[i]);
	}
	double quad = 0;
	for (int i = 0; i < NNode; i++) {
		quad += phi[i] * crs_row_dot(&A, i, phi);
	}
	const double want = aca * vol;
	const double err = ((want != 0) ? (fabs(quad - want) / fabs(want)) : 0);
	fprintf(fp_log, "  volume = %.10e\n", vol);
	fprintf(fp_log, "  linear field : phi^T K phi = %.10e, closed form = %.10e, "
		"rel. error = %.2e\n", quad, want, err);
	if (!(err <= 1e-12)) {		// NaN も不合格にする
		fprintf(fp_log, "*** the linear-field identity failed\n");
		ierr = 1;
	}

	free(phi);
	crs_free(&A);

	return ierr;
}


// 角柱の自己検証 (線形場の恒等式。六面体と同じ考え方)
static int nodal_test_prism(FILE *fp_log)
{
	int ierr = 0;

	fprintf(fp_log, "\n=== nodal element (prism) self test ===\n");
	fprintf(fp_log, "  nodes = %d, prisms = %d, nodes per element = 6\n", NNode, NPrism);

	const double c[6] = {2.0, 3.0, 1.5, 0.4, 0.3, 0.2};
	const double a[3] = {0.7, -1.3, 0.9};
	double aca = 0;
	{
		const double cm[3][3] = {{c[0], c[3], c[5]}, {c[3], c[1], c[4]}, {c[5], c[4], c[2]}};
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) aca += a[i] * cm[i][j] * a[j];
		}
	}

	// ゆがみ : 上下の三角形が平行移動の関係にあると J が ζ に依らなくなる。
	// その形だけで検証すると「J を 1 回だけ評価する」誤りが素通りする
	double warp = 0;
	for (int e = 0; e < NPrism; e++) {
		const int32_t *nd = &Prism[e * 6];
		double h = 0;
		for (int l = 0; l < 6; l++) {
			for (int m = l + 1; m < 6; m++) {
				const double dx = Xp[nd[m]] - Xp[nd[l]];
				const double dy = Yp[nd[m]] - Yp[nd[l]];
				const double dz = Zp[nd[m]] - Zp[nd[l]];
				const double d = sqrt((dx * dx) + (dy * dy) + (dz * dz));
				if (d > h) h = d;
			}
		}
		if (h <= 0) continue;
		for (int cd = 0; cd < 3; cd++) {
			const double *p = ((cd == 0) ? Xp : (cd == 1) ? Yp : Zp);
			// 上下の辺ベクトルの差 (平行移動なら 0)
			for (int l = 0; l < 2; l++) {
				const double d1 = p[nd[l + 1]] - p[nd[l]];
				const double d2 = p[nd[l + 4]] - p[nd[l + 3]];
				const double d = fabs(d2 - d1) / h;
				if (d > warp) warp = d;
			}
		}
	}
	fprintf(fp_log, "  warp = %.3e (0 = the two triangles are a translation of "
		"each other)\n", warp);
	if (warp < 1e-9) {
		fprintf(fp_log, "*** warning : this mesh has no distorted element, so the "
			"isoparametric mapping is not really exercised\n");
	}

	double v6 = 0, vhi = 0;
	for (int e = 0; e < NPrism; e++) {
		const int32_t *nd = &Prism[e * 6];
		for (int q = 0; q < 6; q++) {
			double g[6][3], det;
			if (prism_shape(nd, q, g, &det)) continue;
			v6 += (1.0 / 6) * ((det > 0) ? det : -det);
		}
		vhi += prism_volume_hi(nd);
	}
	const double vdif = ((v6 > 0) ? (fabs(vhi - v6) / v6) : 0);
	fprintf(fp_log, "  volume = %.10e (6-point rule), %.10e (18-point rule), "
		"rel. diff = %.2e\n", v6, vhi, vdif);
	if (vdif > 1e-12) {
		fprintf(fp_log, "*** the two quadrature rules disagree on the volume\n");
		ierr = 1;
	}

	crs_t A;
	crs_alloc(&A);
	crs_zero(&A);
	for (int e = 0; e < NPrism; e++) {
		const int32_t *nd = &Prism[e * 6];
		double ke[6][6], vol;
		if (prism_element(nd, c, ke, &vol)) continue;
		for (int l = 0; l < 6; l++) {
			for (int m = 0; m < 6; m++) {
				const int64_t p = crs_find(&A, nd[l], nd[m]);
				if (p >= 0) A.val[p] += ke[l][m];
			}
		}
	}
	double *phi = (double *)malloc((size_t)NNode * sizeof(double));
	for (int i = 0; i < NNode; i++) {
		phi[i] = (a[0] * Xp[i]) + (a[1] * Yp[i]) + (a[2] * Zp[i]);
	}
	double quad = 0;
	for (int i = 0; i < NNode; i++) {
		quad += phi[i] * crs_row_dot(&A, i, phi);
	}
	const double want = aca * v6;
	const double err = ((want != 0) ? (fabs(quad - want) / fabs(want)) : 0);
	fprintf(fp_log, "  linear field : phi^T K phi = %.10e, closed form = %.10e, "
		"rel. error = %.2e\n", quad, want, err);
	if (!(err <= 1e-12)) {		// NaN も不合格にする
		fprintf(fp_log, "*** the linear-field identity failed\n");
		ierr = 1;
	}

	free(phi);
	crs_free(&A);

	return ierr;
}


// ピラミッドの自己検証 (線形場の恒等式。六面体・角柱と同じ考え方)
static int nodal_test_pyr(FILE *fp_log)
{
	int ierr = 0;

	fprintf(fp_log, "\n=== nodal element (pyramid) self test ===\n");
	fprintf(fp_log, "  nodes = %d, pyramids = %d, nodes per element = 5\n", NNode, NPyr);

	const double c[6] = {2.0, 3.0, 1.5, 0.4, 0.3, 0.2};
	const double a[3] = {0.7, -1.3, 0.9};
	double aca = 0;
	{
		const double cm[3][3] = {{c[0], c[3], c[5]}, {c[3], c[1], c[4]}, {c[5], c[4], c[2]}};
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) aca += a[i] * cm[i][j] * a[j];
		}
	}

	// ゆがみ : 底面が平行四辺形から外れている度合い (符号つき頂点和)。
	// **剛性が厳密なのは底面が平行四辺形のとき** (基準ピラミッドのアフィン像。
	// 頂点の位置は任意でよい) なので、そこから外れた形が無いと有理式の
	// 積分近似が一度も実行されない。ヤコビアン自体はどんな形でも要素内で
	// 変化する ((1-ζ)² の因子) ため、六面体と違い「J を 1 回だけ評価する」
	// 誤りはゆがみが無くても落ちる
	double warp = 0;
	for (int e = 0; e < NPyr; e++) {
		const int32_t *nd = &Pyr[e * 5];
		double h = 0;
		for (int l = 0; l < 5; l++) {
			for (int m = l + 1; m < 5; m++) {
				const double dx = Xp[nd[m]] - Xp[nd[l]];
				const double dy = Yp[nd[m]] - Yp[nd[l]];
				const double dz = Zp[nd[m]] - Zp[nd[l]];
				const double d = sqrt((dx * dx) + (dy * dy) + (dz * dz));
				if (d > h) h = d;
			}
		}
		if (h <= 0) continue;
		for (int cd = 0; cd < 3; cd++) {
			const double *p = ((cd == 0) ? Xp : (cd == 1) ? Yp : Zp);
			const double d = fabs(p[nd[0]] - p[nd[1]] + p[nd[2]] - p[nd[3]]) / h;
			if (d > warp) warp = d;
		}
	}
	fprintf(fp_log, "  warp = %.3e (0 = every base is a parallelogram, where the "
		"stiffness is exact)\n", warp);
	if (warp < 1e-9) {
		fprintf(fp_log, "*** warning : every base is a parallelogram, so the "
			"rational (non-exact) regime of the quadrature is not exercised\n");
	}

	// 体積 : det J は各方向 2 次までなので、8 点則はどんな 5 節点ピラミッド
	// でも厳密。27 点則との差は積分点・重みの誤りだけを検出する
	double v8 = 0, vhi = 0;
	for (int e = 0; e < NPyr; e++) {
		const int32_t *nd = &Pyr[e * 5];
		for (int q = 0; q < 8; q++) {
			double g[5][3], det;
			if (pyr_shape(nd, q, g, &det)) continue;
			v8 += 0.5 * ((det > 0) ? det : -det);
		}
		vhi += pyr_volume_hi(nd);
	}
	const double vdif = ((v8 > 0) ? (fabs(vhi - v8) / v8) : 0);
	fprintf(fp_log, "  volume = %.10e (8-point rule), %.10e (27-point rule), "
		"rel. diff = %.2e\n", v8, vhi, vdif);
	if (vdif > 1e-12) {
		fprintf(fp_log, "*** the two quadrature rules disagree on the volume\n");
		ierr = 1;
	}

	crs_t A;
	crs_alloc(&A);
	crs_zero(&A);
	for (int e = 0; e < NPyr; e++) {
		const int32_t *nd = &Pyr[e * 5];
		double ke[5][5], vol;
		if (pyr_element(nd, c, ke, &vol)) continue;
		for (int l = 0; l < 5; l++) {
			for (int m = 0; m < 5; m++) {
				const int64_t p = crs_find(&A, nd[l], nd[m]);
				if (p >= 0) A.val[p] += ke[l][m];
			}
		}
	}
	double *phi = (double *)malloc((size_t)NNode * sizeof(double));
	for (int i = 0; i < NNode; i++) {
		phi[i] = (a[0] * Xp[i]) + (a[1] * Yp[i]) + (a[2] * Zp[i]);
	}
	double quad = 0;
	for (int i = 0; i < NNode; i++) {
		quad += phi[i] * crs_row_dot(&A, i, phi);
	}
	const double want = aca * v8;
	const double err = ((want != 0) ? (fabs(quad - want) / fabs(want)) : 0);
	fprintf(fp_log, "  linear field : phi^T K phi = %.10e, closed form = %.10e, "
		"rel. error = %.2e\n", quad, want, err);
	if (!(err <= 1e-12)) {		// NaN も不合格にする
		fprintf(fp_log, "*** the linear-field identity failed\n");
		ierr = 1;
	}

	free(phi);
	crs_free(&A);

	return ierr;
}


// ---- 2 次の六面体・角柱の自己検証 (analysis = P) ----
//
// 検査は 3 つ:
//  (1) 体積を 3 点則と 4 点則の両方で積分して一致を見る (積分点・重みの誤り)。
//      **まっすぐな要素 (中間節点が 1 次の写像の位置にある) なら det J は
//      1 次要素と同じ多項式**なので両方とも厳密で、機械精度で一致する。
//  (2) 線形場の恒等式 φᵀKφ = (aᵀCa)V。等パラメトリックなのでどんな形でも成立。
//      **1 次の φ では ∇φ が定数になり、2 次特有の基底 (辺・面・体心) の係数の
//      誤りの多くはここを素通りする** (tet10 と同じ構図)。
//  (3) 2 次場の恒等式 φ = a・r + rᵀBr/2 (**これだけが 2 次の基底と積分則を
//      検査する**)。2 次の φ が空間に入るのは要素がアフィンのときだけなので、
//      軸に平行な直方体 (六面体) / 鉛直な直角柱 (角柱) の要素に限って実行し、
//      期待値は直方体・三角形の閉形式モーメントから組む (組み立てとは独立)。
//      それ以外の要素が混ざる格子では (3) を飛ばして (1)(2) だけ見る。

// 4 点 Gauss-Legendre (3 点則の検算に使う独立な積分)



// 中間節点が「1 次の写像の位置」からどれだけずれているか (最大、辺長比)
static double hex2_curve(void)
{
	double curve = 0;

	for (int e = 0; e < NHex; e++) {
		int32_t nd[ELEM_NMAX];
		const int nen = hex_nodes(e, nd);
		double h = 0;
		for (int l = 0; l < 8; l++) {
			for (int m = l + 1; m < 8; m++) {
				const double dx = Xp[nd[m]] - Xp[nd[l]];
				const double dy = Yp[nd[m]] - Yp[nd[l]];
				const double dz = Zp[nd[m]] - Zp[nd[l]];
				const double d = sqrt((dx * dx) + (dy * dy) + (dz * dz));
				if (d > h) h = d;
			}
		}
		if (h <= 0) continue;
		for (int a = 8; a < nen; a++) {
			// 三重線形の写像での位置 (角の重み付き平均)
			double px = 0, py = 0, pz = 0;
			for (int cn = 0; cn < 8; cn++) {
				const double w = (1 + (HEX_SGN[cn][0] * HEX2_XI[a][0]))
				               * (1 + (HEX_SGN[cn][1] * HEX2_XI[a][1]))
				               * (1 + (HEX_SGN[cn][2] * HEX2_XI[a][2])) / 8.0;
				px += w * Xp[nd[cn]];
				py += w * Yp[nd[cn]];
				pz += w * Zp[nd[cn]];
			}
			const double dx = Xp[nd[a]] - px, dy = Yp[nd[a]] - py, dz = Zp[nd[a]] - pz;
			const double d = sqrt((dx * dx) + (dy * dy) + (dz * dz)) / h;
			if (d > curve) curve = d;
		}
	}

	return curve;
}


// 2 次場 φ = a・r + rᵀBr/2 の勾配 g = a + Br に対する ∫gᵀCg dV を、
// モーメント m1 = ∫r dV/V、m2 = ∫r rᵀ dV/V から組む (要素の形の閉形式)
static double quad_energy(const double cm[3][3], const double av[3],
	const double bm[3][3], const double m1[3], const double m2[3][3], double vol)
{
	double x = 0;

	for (int k = 0; k < 3; k++) {
		for (int l = 0; l < 3; l++) {
			double gg = av[k] * av[l];
			for (int p = 0; p < 3; p++) {
				gg += ((av[k] * bm[l][p]) + (av[l] * bm[k][p])) * m1[p];
				for (int q = 0; q < 3; q++) {
					gg += bm[k][p] * bm[l][q] * m2[p][q];
				}
			}
			x += cm[k][l] * gg * vol;
		}
	}

	return x;
}


static int nodal_test_hex2(FILE *fp_log)
{
	int ierr = 0;

	fprintf(fp_log, "\n=== nodal element (hexahedron, order 2) self test ===\n");
	fprintf(fp_log, "  nodes = %d, hexahedra = %d, nodes per element = %d\n",
		NNode, NHex, HexNen);

	const double c[6] = {2.0, 3.0, 1.5, 0.4, 0.3, 0.2};
	const double cm[3][3] = {{c[0], c[3], c[5]}, {c[3], c[1], c[4]}, {c[5], c[4], c[2]}};
	const double av[3] = {0.7, -1.3, 0.9};
	const double bm[3][3] = {{ 3.0, -1.1,  0.6},
	                         {-1.1,  2.2,  1.7},
	                         { 0.6,  1.7, -0.9}};
	double aca = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) aca += av[i] * cm[i][j] * av[j];
	}

	const double curve = hex2_curve();
	fprintf(fp_log, "  mid-node offset = %.3e (0 = straight elements)\n", curve);

	// ゆがみ (平行六面体からのずれ)。0 だと (2) の線形場が
	// 「J を 1 回だけ評価する」誤りを検出できない (1 次の六面体と同じ)
	double warp = 0;
	for (int e = 0; e < NHex; e++) {
		const int32_t *nd = &Hex[e * 8];
		double h = 0;
		for (int l = 0; l < 8; l++) {
			for (int m = l + 1; m < 8; m++) {
				const double dx = Xp[nd[m]] - Xp[nd[l]];
				const double dy = Yp[nd[m]] - Yp[nd[l]];
				const double dz = Zp[nd[m]] - Zp[nd[l]];
				const double d = sqrt((dx * dx) + (dy * dy) + (dz * dz));
				if (d > h) h = d;
			}
		}
		if (h <= 0) continue;
		for (int cd = 0; cd < 3; cd++) {
			const double *p = ((cd == 0) ? Xp : (cd == 1) ? Yp : Zp);
			double sum = 0;
			for (int l = 0; l < 8; l++) {
				const int sg = HEX_SGN[l][0] * HEX_SGN[l][1] * HEX_SGN[l][2];
				sum += sg * p[nd[l]];
			}
			const double d = fabs(sum) / h;
			if (d > warp) warp = d;
		}
	}
	fprintf(fp_log, "  warp = %.3e (0 = parallelepiped)\n", warp);

	// (1) 体積 : 3 点則と 4 点則の一致
	double v3 = 0, v4 = 0;
	for (int e = 0; e < NHex; e++) {
		int32_t nd[ELEM_NMAX];
		const int nen = hex_nodes(e, nd);
		for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
		for (int k = 0; k < 3; k++) {
			double g[HEXP_NMAX][3], det;
			if (hex2_grad(nd, nen, GX3[i], GX3[j], GX3[k], g, &det)) continue;
			v3 += GW3[i] * GW3[j] * GW3[k] * fabs(det);
		}
		}
		}
		for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
		for (int k = 0; k < 4; k++) {
			double g[HEXP_NMAX][3], det;
			if (hex2_grad(nd, nen, GX4[i], GX4[j], GX4[k], g, &det)) continue;
			v4 += GW4[i] * GW4[j] * GW4[k] * fabs(det);
		}
		}
		}
	}
	const double vdif = ((v3 > 0) ? (fabs(v4 - v3) / v3) : 0);
	fprintf(fp_log, "  volume = %.10e (3-point rule), %.10e (4-point rule), "
		"rel. diff = %.2e\n", v3, v4, vdif);
	if (vdif > 1e-12) {
		fprintf(fp_log, "*** the two quadrature rules disagree on the volume\n");
		ierr = 1;
	}

	// 剛性行列 (この試験用の異方性 c で直接組む)
	crs_t A;
	crs_alloc(&A);
	crs_zero(&A);
	for (int e = 0; e < NHex; e++) {
		int32_t nd[ELEM_NMAX];
		const int nen = hex_nodes(e, nd);
		double ke[HEXP_NMAX][HEXP_NMAX], vol;
		if (hex2_element(e, c, ke, &vol)) continue;
		for (int l = 0; l < nen; l++) {
			for (int m = 0; m < nen; m++) {
				const int64_t p = crs_find(&A, nd[l], nd[m]);
				if (p >= 0) A.val[p] += ke[l][m];
			}
		}
	}

	double *phi = (double *)malloc((size_t)NNode * sizeof(double));

	// (2) 線形場の恒等式 (どんな形でも成立)
	for (int i = 0; i < NNode; i++) {
		phi[i] = (av[0] * Xp[i]) + (av[1] * Yp[i]) + (av[2] * Zp[i]);
	}
	double quad = 0;
	for (int i = 0; i < NNode; i++) {
		quad += phi[i] * crs_row_dot(&A, i, phi);
	}
	{
		const double want = aca * v3;
		const double err = ((want != 0) ? (fabs(quad - want) / fabs(want)) : 0);
		fprintf(fp_log, "  linear field : phi^T K phi = %.10e, closed form = %.10e, "
			"rel. error = %.2e\n", quad, want, err);
		if (!(err <= 1e-12)) {		// NaN も不合格にする
			fprintf(fp_log, "*** the linear-field identity failed\n");
			ierr = 1;
		}
	}

	// (3) 2 次場の恒等式 (軸に平行な直方体の要素に限る)
	{
		int nbox = 0;
		double exact = 0;
		for (int e = 0; e < NHex; e++) {
			const int32_t *nd = &Hex[e * 8];
			double lo[3], hi[3];
			for (int d = 0; d < 3; d++) {
				const double *p = ((d == 0) ? Xp : (d == 1) ? Yp : Zp);
				lo[d] = hi[d] = p[nd[0]];
				for (int l = 1; l < 8; l++) {
					if (p[nd[l]] < lo[d]) lo[d] = p[nd[l]];
					if (p[nd[l]] > hi[d]) hi[d] = p[nd[l]];
				}
			}
			// 8 頂点がすべて min/max の組み合わせにあるか (= 直方体)
			int box = 1;
			double h = 0;
			for (int d = 0; d < 3; d++) {
				if ((hi[d] - lo[d]) > h) h = hi[d] - lo[d];
			}
			for (int l = 0; l < 8; l++) {
				for (int d = 0; d < 3; d++) {
					const double *p = ((d == 0) ? Xp : (d == 1) ? Yp : Zp);
					const double x = p[nd[l]];
					if ((fabs(x - lo[d]) > (1e-9 * h)) && (fabs(x - hi[d]) > (1e-9 * h))) box = 0;
				}
			}
			if (!box || (curve > 1e-9)) continue;
			nbox++;
			double m1[3], m2[3][3], vol = 1;
			for (int d = 0; d < 3; d++) {
				m1[d] = (lo[d] + hi[d]) / 2;
				vol *= (hi[d] - lo[d]);
			}
			for (int p = 0; p < 3; p++) {
				for (int q = 0; q < 3; q++) {
					m2[p][q] = m1[p] * m1[q];
					if (p == q) {
						const double w = hi[p] - lo[p];
						m2[p][p] += (w * w) / 12;
					}
				}
			}
			exact += quad_energy(cm, av, bm, m1, m2, vol);
		}
		if (nbox == NHex) {
			for (int i = 0; i < NNode; i++) {
				const double r[3] = {Xp[i], Yp[i], Zp[i]};
				double f = (av[0] * r[0]) + (av[1] * r[1]) + (av[2] * r[2]);
				for (int k = 0; k < 3; k++) {
					for (int l = 0; l < 3; l++) f += bm[k][l] * r[k] * r[l] / 2;
				}
				phi[i] = f;
			}
			double q2 = 0;
			for (int i = 0; i < NNode; i++) {
				q2 += phi[i] * crs_row_dot(&A, i, phi);
			}
			const double err = fabs(q2 - exact) / ((exact != 0) ? fabs(exact) : 1);
			fprintf(fp_log, "  quadratic field : phi^T K phi = %.10e, exact = %.10e, "
				"rel. error = %.2e\n", q2, exact, err);
			if (!(err <= 1e-10)) {		// NaN も不合格にする
				fprintf(fp_log, "*** the quadratic-field identity failed "
					"(this is the only test that exercises the order-2 basis)\n");
				ierr = 1;
			}
		}
		else {
			fprintf(fp_log, "  quadratic field : skipped (%d of %d elements are "
				"axis-aligned straight boxes; the identity needs all of them)\n",
				nbox, NHex);
		}
	}

	free(phi);
	crs_free(&A);

	return ierr;
}


static int nodal_test_prism2(FILE *fp_log)
{
	int ierr = 0;

	fprintf(fp_log, "\n=== nodal element (prism, order 2) self test ===\n");
	fprintf(fp_log, "  nodes = %d, prisms = %d, nodes per element = %d\n",
		NNode, NPrism, PrismNen);

	const double c[6] = {2.0, 3.0, 1.5, 0.4, 0.3, 0.2};
	const double cm[3][3] = {{c[0], c[3], c[5]}, {c[3], c[1], c[4]}, {c[5], c[4], c[2]}};
	const double av[3] = {0.7, -1.3, 0.9};
	const double bm[3][3] = {{ 3.0, -1.1,  0.6},
	                         {-1.1,  2.2,  1.7},
	                         { 0.6,  1.7, -0.9}};
	double aca = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) aca += av[i] * cm[i][j] * av[j];
	}

	// (1) 体積 : 面内 Duffy 9 点 x ζ 3 点と、面内 Duffy 16 点 x ζ 4 点の一致
	double v3 = 0, v4 = 0;
	for (int e = 0; e < NPrism; e++) {
		int32_t nd[PRISMP_NMAX];
		const int nen = prism_nodes(e, nd);
		double lam[NQTRI][3], wq[NQTRI];
		tri_quad(lam, wq);
		for (int q = 0; q < NQTRI; q++) {
			for (int k = 0; k < 3; k++) {
				double g[PRISMP_NMAX][3], det;
				if (prism2_grad(nd, nen, lam[q][1], lam[q][2], GX3[k], g, &det)) continue;
				v3 += wq[q] * GW3[k] * fabs(det);
			}
		}
		// 独立な高次則 : Duffy 4x4 (面内) x 4 点 (ζ)
		for (int a = 0; a < 4; a++) {
		for (int b = 0; b < 4; b++) {
			const double u = (1 + GX4[a]) / 2, w = (1 + GX4[b]) / 2;
			const double l1 = u, l2 = w * (1 - u);
			const double jw = (GW4[a] / 2) * (GW4[b] / 2) * (1 - u);
			for (int k = 0; k < 4; k++) {
				double g[PRISMP_NMAX][3], det;
				if (prism2_grad(nd, nen, l1, l2, GX4[k], g, &det)) continue;
				v4 += jw * GW4[k] * fabs(det);
			}
		}
		}
	}
	const double vdif = ((v3 > 0) ? (fabs(v4 - v3) / v3) : 0);
	fprintf(fp_log, "  volume = %.10e (9x3-point rule), %.10e (16x4-point rule), "
		"rel. diff = %.2e\n", v3, v4, vdif);
	if (vdif > 1e-12) {
		fprintf(fp_log, "*** the two quadrature rules disagree on the volume\n");
		ierr = 1;
	}

	crs_t A;
	crs_alloc(&A);
	crs_zero(&A);
	for (int e = 0; e < NPrism; e++) {
		int32_t nd[PRISMP_NMAX];
		const int nen = prism_nodes(e, nd);
		double ke[PRISMP_NMAX][PRISMP_NMAX], vol;
		if (prism2_element(e, c, ke, &vol)) continue;
		for (int l = 0; l < nen; l++) {
			for (int m = 0; m < nen; m++) {
				const int64_t p = crs_find(&A, nd[l], nd[m]);
				if (p >= 0) A.val[p] += ke[l][m];
			}
		}
	}

	double *phi = (double *)malloc((size_t)NNode * sizeof(double));

	// (2) 線形場の恒等式
	for (int i = 0; i < NNode; i++) {
		phi[i] = (av[0] * Xp[i]) + (av[1] * Yp[i]) + (av[2] * Zp[i]);
	}
	double quad = 0;
	for (int i = 0; i < NNode; i++) {
		quad += phi[i] * crs_row_dot(&A, i, phi);
	}
	{
		const double want = aca * v3;
		const double err = ((want != 0) ? (fabs(quad - want) / fabs(want)) : 0);
		fprintf(fp_log, "  linear field : phi^T K phi = %.10e, closed form = %.10e, "
			"rel. error = %.2e\n", quad, want, err);
		if (!(err <= 1e-12)) {		// NaN も不合格にする
			fprintf(fp_log, "*** the linear-field identity failed\n");
			ierr = 1;
		}
	}

	// (3) 2 次場の恒等式 (鉛直に押し出したまっすぐな角柱の要素に限る)
	{
		int nright = 0;
		double exact = 0;
		for (int e = 0; e < NPrism; e++) {
			int32_t nd[PRISMP_NMAX];
			const int nen = prism_nodes(e, nd);
			const int32_t *v6 = &Prism[e * 6];
			double h = 0;
			for (int l = 0; l < 6; l++) {
				for (int m = l + 1; m < 6; m++) {
					const double dx = Xp[v6[m]] - Xp[v6[l]];
					const double dy = Yp[v6[m]] - Yp[v6[l]];
					const double dz = Zp[v6[m]] - Zp[v6[l]];
					const double d = sqrt((dx * dx) + (dy * dy) + (dz * dz));
					if (d > h) h = d;
				}
			}
			if (h <= 0) continue;
			// 鉛直な直角柱か : 上面が下面の (0,0,dz) 平行移動で、下面の z が一定
			int right = 1;
			const double dz = Zp[v6[3]] - Zp[v6[0]];
			for (int l = 0; l < 3; l++) {
				if (fabs(Xp[v6[3 + l]] - Xp[v6[l]]) > (1e-9 * h)) right = 0;
				if (fabs(Yp[v6[3 + l]] - Yp[v6[l]]) > (1e-9 * h)) right = 0;
				if (fabs((Zp[v6[3 + l]] - Zp[v6[l]]) - dz) > (1e-9 * h)) right = 0;
				if (fabs(Zp[v6[l]] - Zp[v6[0]]) > (1e-9 * h)) right = 0;
			}
			// 中間節点がまっすぐか (対応する頂点の平均にあるか)
			for (int a = 6; a < nen; a++) {
				// PRISM2_T / PRISM2_Z から「どの頂点の平均か」を作る
				const int t = PRISM2_T[a], z = PRISM2_Z[a];
				int vs[4], nv = 0;
				if (t < 3) { vs[nv++] = t; }
				else { vs[nv++] = TRI6_ED[t - 3][0]; vs[nv++] = TRI6_ED[t - 3][1]; }
				const int base = nv;
				if (z == 0) {
					for (int q = 0; q < base; q++) vs[nv++] = vs[q] + 3;
				}
				else if (z > 0) {
					for (int q = 0; q < base; q++) vs[q] += 3;
				}
				double px = 0, py = 0, pz = 0;
				for (int q = 0; q < nv; q++) {
					px += Xp[v6[vs[q]]] / nv;
					py += Yp[v6[vs[q]]] / nv;
					pz += Zp[v6[vs[q]]] / nv;
				}
				const double ddx = Xp[nd[a]] - px, ddy = Yp[nd[a]] - py, ddz = Zp[nd[a]] - pz;
				if ((sqrt((ddx * ddx) + (ddy * ddy) + (ddz * ddz)) / h) > 1e-9) right = 0;
			}
			if (!right) continue;
			nright++;
			// モーメント : 三角形 (面内) x 一様 (z)。面積は下面の外積から
			const double area = fabs(((Xp[v6[1]] - Xp[v6[0]]) * (Yp[v6[2]] - Yp[v6[0]]))
			                       - ((Xp[v6[2]] - Xp[v6[0]]) * (Yp[v6[1]] - Yp[v6[0]]))) / 2;
			double s[2] = {0, 0}, sq[2][2];
			for (int k = 0; k < 2; k++) {
				const double *p = ((k == 0) ? Xp : Yp);
				for (int a = 0; a < 3; a++) s[k] += p[v6[a]];
			}
			for (int k = 0; k < 2; k++) {
				for (int l = 0; l < 2; l++) {
					const double *pk = ((k == 0) ? Xp : Yp);
					const double *pl = ((l == 0) ? Xp : Yp);
					double t2 = 0;
					for (int a = 0; a < 3; a++) t2 += pk[v6[a]] * pl[v6[a]];
					sq[k][l] = ((s[k] * s[l]) + t2) / 12;
				}
			}
			const double zc = Zp[v6[0]] + (dz / 2);
			double m1[3], m2[3][3];
			m1[0] = s[0] / 3;
			m1[1] = s[1] / 3;
			m1[2] = zc;
			for (int k = 0; k < 2; k++) {
				for (int l = 0; l < 2; l++) m2[k][l] = sq[k][l];
			}
			m2[0][2] = m2[2][0] = m1[0] * zc;
			m2[1][2] = m2[2][1] = m1[1] * zc;
			m2[2][2] = (zc * zc) + ((dz * dz) / 12);
			exact += quad_energy(cm, av, bm, m1, m2, area * fabs(dz));
		}
		if (nright == NPrism) {
			for (int i = 0; i < NNode; i++) {
				const double r[3] = {Xp[i], Yp[i], Zp[i]};
				double f = (av[0] * r[0]) + (av[1] * r[1]) + (av[2] * r[2]);
				for (int k = 0; k < 3; k++) {
					for (int l = 0; l < 3; l++) f += bm[k][l] * r[k] * r[l] / 2;
				}
				phi[i] = f;
			}
			double q2 = 0;
			for (int i = 0; i < NNode; i++) {
				q2 += phi[i] * crs_row_dot(&A, i, phi);
			}
			const double err = fabs(q2 - exact) / ((exact != 0) ? fabs(exact) : 1);
			fprintf(fp_log, "  quadratic field : phi^T K phi = %.10e, exact = %.10e, "
				"rel. error = %.2e\n", q2, exact, err);
			if (!(err <= 1e-10)) {		// NaN も不合格にする
				fprintf(fp_log, "*** the quadratic-field identity failed "
					"(this is the only test that exercises the order-2 basis)\n");
				ierr = 1;
			}
		}
		else {
			fprintf(fp_log, "  quadratic field : skipped (%d of %d elements are "
				"straight vertical prisms; the identity needs all of them)\n",
				nright, NPrism);
		}
	}

	free(phi);
	crs_free(&A);

	return ierr;
}


/*
3 次元の非構造格子の組み立て (種別の混在を含む)。

**要素の走査順は elem3d_* の連番 (四面体 -> 六面体 -> 角柱)。** 純粋な
四面体格子では従来の assemble_tet と同じ順・同じ演算になるので、
浮動小数の丸めまで含めて答えが変わらない (既存の回帰で確認する)。
*/
void assemble_elem3d(crs_t *A, int mode)
{
	crs_zero(A);

	const int ne = elem3d_count();
	const int p2 = (TetOrder >= 2);

	for (int e = 0; e < ne; e++) {
		double c[6];
		material_coef_pub(elem3d_mat(e), mode, c);
		if ((c[0] <= 0) && (c[1] <= 0) && (c[2] <= 0)) continue;

		const int kind = elem3d_kind(e);

		if (kind == MESHELEM_TET) {
			if (TetOrder >= 3) {
				int32_t nd[TET_NMAX];
				double ke[TET_NMAX][TET_NMAX], vol;
				const int nen = tet_nodes(e, nd);
				if (tet20_element(e, c, ke, &vol)) continue;
				for (int l = 0; l < nen; l++) {
					for (int m = 0; m < nen; m++) {
						const int64_t p = crs_find(A, nd[l], nd[m]);
						if (p >= 0) A->val[p] += ke[l][m];
					}
				}
				continue;
			}
			if (p2) {
				int32_t nd[TET_NMAX];
				double ke[10][10], vol;
				tet_nodes(e, nd);
				if (tet10_element(e, c, ke, &vol)) continue;
				for (int l = 0; l < 10; l++) {
					for (int m = 0; m < 10; m++) {
						const int64_t p = crs_find(A, nd[l], nd[m]);
						if (p >= 0) A->val[p] += ke[l][m];
					}
				}
				continue;
			}
			const int32_t *nd = &Tet[e * 4];
			double g[4][3], vol;
			if (tet_grad_pub(nd, g, &vol)) continue;
			for (int l = 0; l < 4; l++) {
				for (int m = 0; m < 4; m++) {
					const double v = (c[0] * g[l][0] * g[m][0])
					               + (c[1] * g[l][1] * g[m][1])
					               + (c[2] * g[l][2] * g[m][2])
					               + (c[3] * ((g[l][0] * g[m][1]) + (g[l][1] * g[m][0])))
					               + (c[4] * ((g[l][1] * g[m][2]) + (g[l][2] * g[m][1])))
					               + (c[5] * ((g[l][2] * g[m][0]) + (g[l][0] * g[m][2])));
					const int64_t p = crs_find(A, nd[l], nd[m]);
					if (p >= 0) A->val[p] += v * vol;
				}
			}
			continue;
		}

		if (kind == MESHELEM_HEX) {
			if (HexNen > 8) {
				int32_t nd[ELEM_NMAX];
				const int nen = hex_nodes(e - NTet, nd);
				double ke[HEXP_NMAX][HEXP_NMAX], vol;
				if (hex2_element(e - NTet, c, ke, &vol)) continue;
				for (int l = 0; l < nen; l++) {
					for (int m = 0; m < nen; m++) {
						const int64_t p = crs_find(A, nd[l], nd[m]);
						if (p >= 0) A->val[p] += ke[l][m];
					}
				}
				continue;
			}
			const int32_t *nd = &Hex[(e - NTet) * 8];
			double ke[8][8], vol;
			if (hex_element(nd, c, ke, &vol)) continue;
			for (int l = 0; l < 8; l++) {
				for (int m = 0; m < 8; m++) {
					const int64_t p = crs_find(A, nd[l], nd[m]);
					if (p >= 0) A->val[p] += ke[l][m];
				}
			}
			continue;
		}

		if (kind == MESHELEM_PRISM) {
			if (PrismNen > 6) {
				int32_t nd[PRISMP_NMAX];
				const int nen = prism_nodes(e - NTet - NHex, nd);
				double ke[PRISMP_NMAX][PRISMP_NMAX], vol;
				if (prism2_element(e - NTet - NHex, c, ke, &vol)) continue;
				for (int l = 0; l < nen; l++) {
					for (int m = 0; m < nen; m++) {
						const int64_t p = crs_find(A, nd[l], nd[m]);
						if (p >= 0) A->val[p] += ke[l][m];
					}
				}
				continue;
			}
			const int32_t *nd = &Prism[(e - NTet - NHex) * 6];
			double ke[6][6], vol;
			if (prism_element(nd, c, ke, &vol)) continue;
			for (int l = 0; l < 6; l++) {
				for (int m = 0; m < 6; m++) {
					const int64_t p = crs_find(A, nd[l], nd[m]);
					if (p >= 0) A->val[p] += ke[l][m];
				}
			}
			continue;
		}

		{
			const int32_t *nd = &Pyr[(e - NTet - NHex - NPrism) * 5];
			double ke[5][5], vol;
			if (pyr_element(nd, c, ke, &vol)) continue;
			for (int l = 0; l < 5; l++) {
				for (int m = 0; m < 5; m++) {
					const int64_t p = crs_find(A, nd[l], nd[m]);
					if (p >= 0) A->val[p] += ke[l][m];
				}
			}
		}
	}
}


void assemble_hex(crs_t *A, int mode)
{
	crs_zero(A);

	for (int e = 0; e < NHex; e++) {
		double c[6];
		material_coef_pub(HexMat[e], mode, c);
		if ((c[0] <= 0) && (c[1] <= 0) && (c[2] <= 0)) continue;

		const int32_t *nd = &Hex[e * 8];
		double ke[8][8], vol;
		if (hex_element(nd, c, ke, &vol)) continue;

		for (int l = 0; l < 8; l++) {
			for (int m = 0; m < 8; m++) {
				const int64_t p = crs_find(A, nd[l], nd[m]);
				if (p >= 0) A->val[p] += ke[l][m];
			}
		}
	}
}


void assemble_tet(crs_t *A, int mode)
{
	crs_zero(A);

	const int p2 = (TetOrder >= 2);

	for (int e = 0; e < NTet; e++) {
		double c[6];
		material_coef_pub(TetMat[e], mode, c);
		if ((c[0] <= 0) && (c[1] <= 0) && (c[2] <= 0)) continue;

		if (TetOrder >= 3) {
			int32_t nd[TET_NMAX];
			double ke[TET_NMAX][TET_NMAX], vol;
			const int nen = tet_nodes(e, nd);
			if (tet20_element(e, c, ke, &vol)) continue;
			for (int l = 0; l < nen; l++) {
				for (int m = 0; m < nen; m++) {
					const int64_t p = crs_find(A, nd[l], nd[m]);
					if (p >= 0) A->val[p] += ke[l][m];
				}
			}
			continue;
		}
		if (p2) {
			int32_t nd[TET_NMAX];
			double ke[10][10], vol;
			tet_nodes(e, nd);
			if (tet10_element(e, c, ke, &vol)) continue;
			for (int l = 0; l < 10; l++) {
				for (int m = 0; m < 10; m++) {
					const int64_t p = crs_find(A, nd[l], nd[m]);
					if (p >= 0) A->val[p] += ke[l][m];
				}
			}
			continue;
		}

		const int32_t *nd = &Tet[e * 4];
		double g[4][3], vol;
		if (tet_grad_pub(nd, g, &vol)) continue;

		for (int l = 0; l < 4; l++) {
			for (int m = 0; m < 4; m++) {
				const double v = (c[0] * g[l][0] * g[m][0])
				               + (c[1] * g[l][1] * g[m][1])
				               + (c[2] * g[l][2] * g[m][2])
				               + (c[3] * ((g[l][0] * g[m][1]) + (g[l][1] * g[m][0])))
				               + (c[4] * ((g[l][1] * g[m][2]) + (g[l][2] * g[m][1])))
				               + (c[5] * ((g[l][2] * g[m][0]) + (g[l][0] * g[m][2])));
				const int64_t p = crs_find(A, nd[l], nd[m]);
				if (p >= 0) A->val[p] += v * vol;
			}
		}
	}
}


// ---- 節点要素 (P1 / P2) の自己検証 (analysis = P) ----
//
// 剛性行列を組んで、多項式の再現性を厳密な恒等式で検査する。
//
//   φ(r) = a・r + (1/2) r^T B r   (B は対称)
//   ∇φ = a + B r  (r について 1 次)
//   φ^T K φ = ∫ (∇φ)^T C (∇φ) dV
//
// 右辺は要素毎に閉形式で書ける。r が 1 次なので被積分関数は 2 次で、
//   ∫ r_k dV = V (Σ_a p_a,k)/4
//   ∫ r_k r_l dV = (V/20)[(Σ_a p_a,k)(Σ_b p_b,l) + Σ_a p_a,k p_a,l]
// (∫λ_aλ_b dV = V(1+δ_ab)/20 と r = Σ λ_a p_a から)。組み立てとは独立に
// 計算するので、形状関数・数値積分・組み立てのどれが壊れても落ちる。
//
// **検査の効き方**
//  ・1 次要素は φ が 1 次のときだけ補間が厳密。2 次の φ は落ちて当然なので
//    次数に応じて実行する検査を変える。
//  ・1 次の φ では ∇φ が要素内で一定になり、どんな数値積分でも厳密になる。
//    積分則の誤り (点数不足・重みの誤り) を捕まえるのは 2 次の φ だけ。
//  ・材料は異方性にすること。等方性だけだと C の非対角成分が死ぬ。
//  ・中間節点が辺の中点に無い (曲がった) 格子では 2 次の φ の補間が厳密で
//    なくなるので、その検査は飛ばして体積だけ見る。曲面の等パラメトリック
//    写像は「積分した体積が解析値と合うか」で別に検証する。
/*
六面体の自己検証 (analysis = P)。

**線形場の恒等式**を機械精度で見る:

	φ(r) = a・r  は等パラメトリック三重線形要素で**厳密に補間できる**
	(φ_h = Σ N_a (a・x_a) = a・Σ N_a x_a = a・x(ξ))。したがって ∇φ_h = a が
	要素内で一定になり、
		φᵀKφ = ∫ (∇φ)ᵀ C (∇φ) dV = (aᵀ C a) V
	が厳密に成り立つ。

**ゆがんだ六面体でないと意味がない。** 平行六面体ではヤコビアンが要素内で
一定になるため、「J を要素中心で 1 回だけ評価する」という誤り (等パラメトリック
要素で最も踏みやすい手抜き) が厳密に正しくなってしまう。実測: 剛体回転した
直方体でも同軸 (回転対称) でも、この誤りは答えを 1 桁も動かさなかった
(同軸は φ が θ に依らないので、要素行列が 3% 違っても радиальный 成分だけが
効いて相殺する)。ゆがんだ格子でだけ 1.6% ずれる。

体積は **2 点則と 3 点則の両方**で積分して一致を見る。detJ は各方向 2 次までの
多項式なので 2 点 Gauss で厳密であり、点数を落とす誤りはここで落ちる
(線形場の恒等式は左右が同じ体積を使うので、それだけでは検出できない)。
*/
static int nodal_test_hex(FILE *fp_log)
{
	int ierr = 0;

	fprintf(fp_log, "\n=== nodal element (hexahedron) self test ===\n");
	fprintf(fp_log, "  nodes = %d, hexahedra = %d, nodes per element = 8\n", NNode, NHex);

	// 異方性テンソル (等方性だと C の非対角成分が一度も実行されない)
	const double c[6] = {2.0, 3.0, 1.5, 0.4, 0.3, 0.2};
	// 軸に平行でない試験場 (軸平行だと C の 1 成分しか見ない)
	const double a[3] = {0.7, -1.3, 0.9};

	double aca = 0;
	{
		const double cm[3][3] = {{c[0], c[3], c[5]}, {c[3], c[1], c[4]}, {c[5], c[4], c[2]}};
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) aca += a[i] * cm[i][j] * a[j];
		}
	}

	// ゆがみ (要素が平行六面体からどれだけ離れているか)。
	// 平行六面体では対角の和 x0 - x1 + x2 - x3 (下面) 等が 0 になる
	double warp = 0, hmax = 0;
	for (int e = 0; e < NHex; e++) {
		const int32_t *nd = &Hex[e * 8];
		double h = 0;
		for (int l = 0; l < 8; l++) {
			for (int m = l + 1; m < 8; m++) {
				const double dx = Xp[nd[m]] - Xp[nd[l]];
				const double dy = Yp[nd[m]] - Yp[nd[l]];
				const double dz = Zp[nd[m]] - Zp[nd[l]];
				const double d = sqrt((dx * dx) + (dy * dy) + (dz * dz));
				if (d > h) h = d;
			}
		}
		if (h > hmax) hmax = h;
		if (h <= 0) continue;
		// 三重線形写像の非アフィン成分 (符号つきの頂点和)
		for (int cdim = 0; cdim < 3; cdim++) {
			const double *p = ((cdim == 0) ? Xp : (cdim == 1) ? Yp : Zp);
			double sum = 0;
			for (int l = 0; l < 8; l++) {
				const int sg = HEX_SGN[l][0] * HEX_SGN[l][1] * HEX_SGN[l][2];
				sum += sg * p[nd[l]];
			}
			const double d = fabs(sum) / h;
			if (d > warp) warp = d;
		}
	}
	fprintf(fp_log, "  warp = %.3e (0 = parallelepiped; the test needs a distorted mesh)\n", warp);
	if (warp < 1e-9) {
		fprintf(fp_log, "*** warning : this mesh has no distorted element, so the "
			"isoparametric mapping is not really exercised\n");
	}

	// 体積 : 2 点則と 3 点則で一致すること (detJ は 2 点で厳密)
	double v2 = 0, v3 = 0;
	for (int e = 0; e < NHex; e++) {
		const int32_t *nd = &Hex[e * 8];
		for (int q = 0; q < 8; q++) {
			double g[8][3], det;
			if (hex_shape(nd, q, g, &det)) continue;
			v2 += ((det > 0) ? det : -det);
		}
		v3 += hex_volume3(nd);
	}
	const double vdif = ((v2 > 0) ? (fabs(v3 - v2) / v2) : 0);
	fprintf(fp_log, "  volume = %.10e (2-point rule), %.10e (3-point rule), rel. diff = %.2e\n",
		v2, v3, vdif);
	if (vdif > 1e-12) {
		fprintf(fp_log, "*** the two quadrature rules disagree on the volume\n");
		ierr = 1;
	}

	// 線形場の恒等式
	crs_t A;
	crs_alloc(&A);
	{
		// material_coef_pub を通さずに直接この c で組む
		crs_zero(&A);
		for (int e = 0; e < NHex; e++) {
			const int32_t *nd = &Hex[e * 8];
			double ke[8][8], vol;
			if (hex_element(nd, c, ke, &vol)) continue;
			for (int l = 0; l < 8; l++) {
				for (int m = 0; m < 8; m++) {
					const int64_t p = crs_find(&A, nd[l], nd[m]);
					if (p >= 0) A.val[p] += ke[l][m];
				}
			}
		}
	}
	double *phi = (double *)malloc((size_t)NNode * sizeof(double));
	for (int i = 0; i < NNode; i++) {
		phi[i] = (a[0] * Xp[i]) + (a[1] * Yp[i]) + (a[2] * Zp[i]);
	}
	double quad = 0;
	for (int i = 0; i < NNode; i++) {
		quad += phi[i] * crs_row_dot(&A, i, phi);
	}
	const double want = aca * v2;
	const double err = ((want != 0) ? (fabs(quad - want) / fabs(want)) : 0);
	fprintf(fp_log, "  linear field : phi^T K phi = %.10e, closed form = %.10e, "
		"rel. error = %.2e\n", quad, want, err);
	if (!(err <= 1e-12)) {		// NaN も不合格にする
		fprintf(fp_log, "*** the linear-field identity failed\n");
		ierr = 1;
	}

	free(phi);
	crs_free(&A);

	return ierr;
}


/*
3 次の六面体 (64 節点) / 角柱 (40 節点) の自己検証 (analysis = P)。
  (1) 線形場 φ = a・r : どんな形でも φᵀKφ = aᵀCa V (V は要素の積分で得た体積)
  (2) 3 次の場 φ = a・r + rᵀBr/2 + γ1 xyz + γ2 x³ + γ3 y²z + γ4 xz² + γ5 z³ (全次数 3) :
      アフィンな要素では Q3 / P3×P3 の空間に厳密に入るので φᵀKφ = Σ_e ∫ ∇φᵀC∇φ。右辺は解析的な
      勾配を**角の 1 次の写像の上で 5 点則**で積分する (形状関数を使わない独立な計算)。
      **3 次特有の基底 (辺の 2 節点・面の 4 節点・体心) を検査できるのはこれだけ** (線形場と平行平板は
      2 次以下に退行しても通る)
*/
static void cubic_grad(const double r[3], const double av[3], const double bm[3][3], const double gm[9], double g[3])
{
	const double x = r[0], y = r[1], z = r[2];
	for (int k = 0; k < 3; k++) {
		g[k] = av[k];
		for (int l = 0; l < 3; l++) g[k] += bm[k][l] * r[l];
	}
	g[0] += (gm[0] * y * z) + (3 * gm[1] * x * x) + (gm[3] * z * z);
	g[1] += (gm[0] * x * z) + (2 * gm[2] * y * z);
	g[2] += (gm[0] * x * y) + (gm[2] * y * y) + (2 * gm[3] * x * z) + (3 * gm[4] * z * z);
	// 4 次の項 (4 次の要素の検査だけで使う。3 次では係数 0)
	g[0] += (4 * gm[5] * x * x * x) + (2 * gm[6] * x * y * y) + (gm[8] * y * z * z);
	g[1] += (2 * gm[6] * x * x * y) + (gm[7] * z * z * z) + (gm[8] * x * z * z);
	g[2] += (3 * gm[7] * y * z * z) + (2 * gm[8] * x * y * z);
}


// 5 次の六面体・角柱の自己検証の 5 次の項。各軸をその軸の格子の大きさ ax[] で割った X, Y, Z で
//   L (0.9 X²Y²Z − 0.5 Y⁵ + 0.7 X Z⁴ + 3 X Z⁵ (+ 3 X⁵ Y : 六面体だけ。角柱の面内 P5 に入らない))
// (L は全体の大きさ、勾配が線形の項と同じ桁になる)。値を返し、勾配を g に足す。**X Z⁵ と X⁵ Y が要素の行列の積分点の
// 数を検査する** : (∂φ/∂x)² が z の 10 次、(∂φ/∂y)² が x の 10 次なので各方向 6 点で厳密、5 点に落とすと外れる。
// 軸ごとに割るのは、薄い箱 (z が 0.2mm) で全体の大きさで割ると Z⁵ の項が 6e-6 しか効かず、5 点に落とす誤りが
// 素通りしたから (実測)
static double quintic_term(double l, const double ax[3], int prism, const double r[3], double g[3])
{
	const double x = r[0] / ax[0], y = r[1] / ax[1], z = r[2] / ax[2];
	const double a = 0.9 * l, b = -0.5 * l, c = 0.7 * l, d = 3.0 * l, e = (prism ? 0 : 3.0 * l);
	const double z4 = z * z * z * z, x4 = x * x * x * x;
	g[0] += ((2 * a * x * y * y * z) + (c * z4) + (d * z4 * z) + (5 * e * x4 * y)) / ax[0];
	g[1] += ((2 * a * x * x * y * z) + (5 * b * y * y * y * y) + (e * x4 * x)) / ax[1];
	g[2] += ((a * x * x * y * y) + (4 * c * x * z * z * z) + (5 * d * x * z4)) / ax[2];
	return (a * x * x * y * y * z) + (b * y * y * y * y * y) + (c * x * z4) + (d * x * z4 * z) + (e * x4 * x * y);
}

static double cubic_val(const double r[3], const double av[3], const double bm[3][3], const double gm[9])
{
	const double x = r[0], y = r[1], z = r[2];
	double f = (av[0] * x) + (av[1] * y) + (av[2] * z);
	for (int k = 0; k < 3; k++) {
		for (int l = 0; l < 3; l++) f += bm[k][l] * r[k] * r[l] / 2;
	}
	return f + (gm[0] * x * y * z) + (gm[1] * x * x * x) + (gm[2] * y * y * z) + (gm[3] * x * z * z) + (gm[4] * z * z * z)
	         + (gm[5] * x * x * x * x) + (gm[6] * x * x * y * y) + (gm[7] * y * z * z * z) + (gm[8] * x * y * z * z);
}


static int nodal_test_ho3(FILE *fp_log, int prism)
{
	int ierr = 0;
	const int ne = (prism ? NPrism : NHex);
	const int po = (prism ? prismp_order(PrismNen) : hexp_order(HexNen));
	fprintf(fp_log, "\n=== nodal element (%s, order %d) self test ===\n", (prism ? "prism" : "hexahedron"), po);
	fprintf(fp_log, "  nodes = %d, elements = %d, nodes per element = %d\n", NNode, ne, (prism ? PrismNen : HexNen));

	const double c[6] = {2.0, 3.0, 1.5, 0.4, 0.3, 0.2};
	const double cm[3][3] = {{c[0], c[3], c[5]}, {c[3], c[1], c[4]}, {c[5], c[4], c[2]}};
	const double av[3] = {0.7, -1.3, 0.9};
	const double bm[3][3] = {{3.0, -1.1, 0.6}, {-1.1, 2.2, 1.7}, {0.6, 1.7, -0.9}};
	// 3 次の項は座標の大きさ (格子の寸法 L) に合わせて振る (L³ の項が線形の項と同じ桁になるように)
	double lmax = 0;
	for (int i = 0; i < NNode; i++) {
		const double m = fabs(Xp[i]) + fabs(Yp[i]) + fabs(Zp[i]);
		if (m > lmax) lmax = m;
	}
	if (lmax <= 0) lmax = 1;
	double axl[3] = {0, 0, 0};		// 軸ごとの大きさ (5 次の項 quintic_term)
	for (int i = 0; i < NNode; i++) {
		if (fabs(Xp[i]) > axl[0]) axl[0] = fabs(Xp[i]);
		if (fabs(Yp[i]) > axl[1]) axl[1] = fabs(Yp[i]);
		if (fabs(Zp[i]) > axl[2]) axl[2] = fabs(Zp[i]);
	}
	for (int k = 0; k < 3; k++) {
		if (axl[k] <= 0) axl[k] = lmax;
	}
	double gm[9] = {1.3, -0.8, 1.1, 0.6, -0.7, 0.9, -1.2, 0.5, 1.4};
	double bs[3][3];
	for (int k = 0; k < 5; k++) gm[k] /= (lmax * lmax);
	for (int k = 5; k < 9; k++) gm[k] = ((po >= 4) ? (gm[k] / (lmax * lmax * lmax)) : 0);
	for (int k = 0; k < 3; k++) {
		for (int l = 0; l < 3; l++) bs[k][l] = bm[k][l] / lmax;
	}

	crs_t A;
	crs_alloc(&A);
	crs_zero(&A);
	double vol = 0;
	for (int e = 0; e < ne; e++) {
		int32_t nd[ELEM_NMAX];
		const int nen = (prism ? prism_nodes(e, nd) : hex_nodes(e, nd));
		static double kh[HEXP_NMAX][HEXP_NMAX], kp[PRISMP_NMAX][PRISMP_NMAX];
		double v;
		if (prism ? prism2_element(e, c, kp, &v) : hex2_element(e, c, kh, &v)) {
			fprintf(fp_log, "*** element %d is degenerate or inverted\n", e + 1);
			ierr = 1;
			continue;
		}
		vol += v;
		for (int l = 0; l < nen; l++) {
			for (int m = 0; m < nen; m++) {
				const int64_t p = crs_find(&A, nd[l], nd[m]);
				if (p >= 0) A.val[p] += (prism ? kp[l][m] : kh[l][m]);
			}
		}
	}
	double *phi = (double *)malloc((size_t)NNode * sizeof(double));

	// (1) 線形場
	double aca = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) aca += av[i] * cm[i][j] * av[j];
	}
	for (int i = 0; i < NNode; i++) phi[i] = (av[0] * Xp[i]) + (av[1] * Yp[i]) + (av[2] * Zp[i]);
	double q1 = 0;
	for (int i = 0; i < NNode; i++) q1 += phi[i] * crs_row_dot(&A, i, phi);
	{
		const double want = aca * vol;
		const double err = fabs(q1 - want) / fabs(want);
		fprintf(fp_log, "  linear field : phi^T K phi = %.10e, closed form = %.10e, rel. error = %.2e\n", q1, want, err);
		if (!(err <= 1e-12)) {
			fprintf(fp_log, "*** the linear-field identity failed\n");
			ierr = 1;
		}
	}

	// (2) 3 次の場 (アフィンな要素だけ。1 つでも外れたら飛ばす)
	int naff = 0;
	double exact = 0;
	// 右辺の積分は各方向 p + 2 点 (要素の行列より 1 点多い。p 次の場の勾配の 2 乗は各方向 2p - 2 次)
	const int nr = po + 2;
	double G5[GL_NMAX], W5[GL_NMAX];
	gauss01(nr, G5, W5);
	for (int k = 0; k < nr; k++) { G5[k] = (2 * G5[k]) - 1; W5[k] *= 2; }
	for (int e = 0; e < ne; e++) {
		int32_t nd[ELEM_NMAX];
		const int nen = (prism ? prism_nodes(e, nd) : hex_nodes(e, nd));
		const int nc = (prism ? 6 : 8);
		double h = 0;
		for (int l = 0; l < nc; l++) {
			for (int m = l + 1; m < nc; m++) {
				const double dx = Xp[nd[m]] - Xp[nd[l]], dy = Yp[nd[m]] - Yp[nd[l]], dz = Zp[nd[m]] - Zp[nd[l]];
				const double d = sqrt((dx * dx) + (dy * dy) + (dz * dz));
				if (d > h) h = d;
			}
		}
		// アフィン : 角の 1 次の写像が平行移動の分を除いて一定 (六面体は平行六面体、角柱は上下が平行移動)、
		// かつ全節点がその写像の上の位置にある (曲がっていない)
		int aff = 1;
		for (int cd = 0; cd < 3; cd++) {
			const double *p = ((cd == 0) ? Xp : (cd == 1) ? Yp : Zp);
			if (prism) {
				for (int t = 0; t < 3; t++) {
					const double dz0 = p[nd[3]] - p[nd[0]], dzt = p[nd[3 + t]] - p[nd[t]];
					if (fabs(dzt - dz0) > (1e-9 * h)) aff = 0;
				}
			}
			else {
				double sum = 0;
				for (int l = 0; l < 8; l++) sum += HEX_SGN[l][0] * HEX_SGN[l][1] * HEX_SGN[l][2] * p[nd[l]];
				// 双 1 次の交差項 (3 組) もアフィンでは 0
				double s01 = 0, s12 = 0, s02 = 0;
				for (int l = 0; l < 8; l++) {
					s01 += HEX_SGN[l][0] * HEX_SGN[l][1] * p[nd[l]];
					s12 += HEX_SGN[l][1] * HEX_SGN[l][2] * p[nd[l]];
					s02 += HEX_SGN[l][0] * HEX_SGN[l][2] * p[nd[l]];
				}
				if ((fabs(sum) + fabs(s01) + fabs(s12) + fabs(s02)) > (1e-9 * h)) aff = 0;
			}
		}
		for (int a = nc; aff && (a < nen); a++) {
			double r0[3] = {0, 0, 0};
			if (prism) {
				const signed char *al = prismp_a(po, a);
				const double lam[3] = {al[0] / (double)po, al[1] / (double)po, al[2] / (double)po};
				const double t = al[3] / (double)po;
				for (int k = 0; k < 3; k++) {
					const double wb = lam[k] * (1 - t), wt = lam[k] * t;
					r0[0] += (wb * Xp[nd[k]]) + (wt * Xp[nd[3 + k]]);
					r0[1] += (wb * Yp[nd[k]]) + (wt * Yp[nd[3 + k]]);
					r0[2] += (wb * Zp[nd[k]]) + (wt * Zp[nd[3 + k]]);
				}
			}
			else {
				const signed char *ix = hexp_i(po, a);
				const double xi[3] = {-1 + (2.0 * ix[0] / po), -1 + (2.0 * ix[1] / po), -1 + (2.0 * ix[2] / po)};
				for (int l = 0; l < 8; l++) {
					const double w = (1 + (HEX_SGN[l][0] * xi[0])) * (1 + (HEX_SGN[l][1] * xi[1])) * (1 + (HEX_SGN[l][2] * xi[2])) / 8;
					r0[0] += w * Xp[nd[l]];
					r0[1] += w * Yp[nd[l]];
					r0[2] += w * Zp[nd[l]];
				}
			}
			const double dx = Xp[nd[a]] - r0[0], dy = Yp[nd[a]] - r0[1], dz = Zp[nd[a]] - r0[2];
			if (sqrt((dx * dx) + (dy * dy) + (dz * dz)) > (1e-9 * h)) aff = 0;
		}
		if (!aff) continue;
		naff++;
		// 角の 1 次の写像の上で 5 点則 (アフィンなので det J は一定)
		if (prism) {
			double lq[QTRI_MAX][3], wq[QTRI_MAX];
			const int nq = tri_quad_pub(nr, lq, wq);
			double j1[3], j2[3], j3[3];
			for (int k = 0; k < 3; k++) {
				const double *p = ((k == 0) ? Xp : (k == 1) ? Yp : Zp);
				j1[k] = p[nd[1]] - p[nd[0]];
				j2[k] = p[nd[2]] - p[nd[0]];
				j3[k] = (p[nd[3]] - p[nd[0]]) / 2;
			}
			const double det = fabs((j1[0] * ((j2[1] * j3[2]) - (j2[2] * j3[1]))) - (j1[1] * ((j2[0] * j3[2]) - (j2[2] * j3[0])))
			                      + (j1[2] * ((j2[0] * j3[1]) - (j2[1] * j3[0]))));
			for (int q = 0; q < nq; q++) {
				for (int k = 0; k < nr; k++) {
					const double t = (G5[k] + 1) / 2;
					double r[3] = {0, 0, 0};
					for (int cc = 0; cc < 3; cc++) {
						const double *p = ((cc == 0) ? Xp : (cc == 1) ? Yp : Zp);
						for (int m = 0; m < 3; m++) r[cc] += lq[q][m] * (((1 - t) * p[nd[m]]) + (t * p[nd[3 + m]]));
					}
					double g[3];
					cubic_grad(r, av, bs, gm, g);
					if (po >= 5) (void)quintic_term(lmax, axl, prism, r, g);
					double gcg = 0;
					for (int i = 0; i < 3; i++) {
						for (int j = 0; j < 3; j++) gcg += g[i] * cm[i][j] * g[j];
					}
					exact += wq[q] * W5[k] * det * gcg;
				}
			}
		}
		else {
			for (int i = 0; i < nr; i++) {
			for (int j = 0; j < nr; j++) {
			for (int k = 0; k < nr; k++) {
				const double xi[3] = {G5[i], G5[j], G5[k]};
				double r[3] = {0, 0, 0}, jm[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
				for (int l = 0; l < 8; l++) {
					const double f[3] = {1 + (HEX_SGN[l][0] * xi[0]), 1 + (HEX_SGN[l][1] * xi[1]), 1 + (HEX_SGN[l][2] * xi[2])};
					const double w = f[0] * f[1] * f[2] / 8;
					const double dw[3] = {HEX_SGN[l][0] * f[1] * f[2] / 8, HEX_SGN[l][1] * f[0] * f[2] / 8, HEX_SGN[l][2] * f[0] * f[1] / 8};
					const double p[3] = {Xp[nd[l]], Yp[nd[l]], Zp[nd[l]]};
					for (int cc = 0; cc < 3; cc++) {
						r[cc] += w * p[cc];
						for (int d = 0; d < 3; d++) jm[cc][d] += p[cc] * dw[d];
					}
				}
				const double det = fabs((jm[0][0] * ((jm[1][1] * jm[2][2]) - (jm[1][2] * jm[2][1])))
				                      - (jm[0][1] * ((jm[1][0] * jm[2][2]) - (jm[1][2] * jm[2][0])))
				                      + (jm[0][2] * ((jm[1][0] * jm[2][1]) - (jm[1][1] * jm[2][0]))));
				double g[3];
				cubic_grad(r, av, bs, gm, g);
				if (po >= 5) (void)quintic_term(lmax, axl, prism, r, g);
				double gcg = 0;
				for (int a = 0; a < 3; a++) {
					for (int b = 0; b < 3; b++) gcg += g[a] * cm[a][b] * g[b];
				}
				exact += W5[i] * W5[j] * W5[k] * det * gcg;
			}
			}
			}
		}
	}
	if (naff == ne) {
		for (int i = 0; i < NNode; i++) {
			const double r[3] = {Xp[i], Yp[i], Zp[i]};
			double gd[3] = {0, 0, 0};
			phi[i] = cubic_val(r, av, bs, gm) + ((po >= 5) ? quintic_term(lmax, axl, prism, r, gd) : 0);
		}
		double q3 = 0;
		for (int i = 0; i < NNode; i++) q3 += phi[i] * crs_row_dot(&A, i, phi);
		const double err = fabs(q3 - exact) / fabs(exact);
		fprintf(fp_log, "  %s field : phi^T K phi = %.10e, exact = %.10e, rel. error = %.2e\n", ((po >= 5) ? "quintic" : (po >= 4) ? "quartic" : "cubic"), q3, exact, err);
		// 許容 1e-12 (実測 1e-15〜7e-14)。1e-10 では 5 次の積分点を 5 点に落とす誤り (1e-11 しか動かない) が素通りした
		if (!(err <= 1e-12)) {
			fprintf(fp_log, "*** the order-%d field identity failed (this is the only test that exercises the order-%d basis)\n", po, po);
			ierr = 1;
		}
	}
	else {
		fprintf(fp_log, "  %s field : skipped (%d of %d elements are affine and straight; the identity needs all of them)\n",
			((po >= 5) ? "quintic" : (po >= 4) ? "quartic" : "cubic"),
			naff, ne);
	}
	free(phi);
	crs_free(&A);
	return ierr;
}


// 自己検証の 6 次 / 7 次の項 (次数 o の単項式 3 つ、係数は 1/L^(o-1) 倍) : 次数 omax 以下の 6 次以上の項の値を返し、
// 勾配を g に足す。6 次 : 0.7 x⁴yz − 0.4 y⁶ + 1.2 x²z⁴、7 次 : 0.6 x³y²z² − 0.3 z⁷ + 1.1 x y⁵ z
static double ho_extra(int omax, double lref, const double r[3], double g[3])
{
	static const int ex[6][4] = {{6, 4, 1, 1}, {6, 0, 6, 0}, {6, 2, 0, 4}, {7, 3, 2, 2}, {7, 0, 0, 7}, {7, 1, 5, 1}};
	static const double cf[6] = {0.7, -0.4, 1.2, 0.6, -0.3, 1.1};
	double f = 0;
	for (int t = 0; t < 6; t++) {
		if (ex[t][0] > omax) continue;
		double c = cf[t];
		for (int k = 1; k < ex[t][0]; k++) c /= lref;
		double pw[3][8];
		for (int a = 0; a < 3; a++) {
			pw[a][0] = 1;
			for (int k = 1; k < 8; k++) pw[a][k] = pw[a][k - 1] * r[a];
		}
		const int *e = &ex[t][1];
		f += c * pw[0][e[0]] * pw[1][e[1]] * pw[2][e[2]];
		for (int a = 0; a < 3; a++) {
			if (e[a] == 0) continue;
			double d = c * e[a];
			for (int b = 0; b < 3; b++) d *= pw[b][(b == a) ? (e[b] - 1) : e[b]];
			g[a] += d;
		}
	}
	return f;
}


int solve_nodal_test(FILE *fp_log)
{
	int ierr = 0;

	if (MeshElem == MESHELEM_HEX) {
		if (hexp_order(HexNen)) return nodal_test_ho3(fp_log, 0);
		return ((HexNen > 8) ? nodal_test_hex2(fp_log) : nodal_test_hex(fp_log));
	}
	if (MeshElem == MESHELEM_PRISM) {
		if (prismp_order(PrismNen)) return nodal_test_ho3(fp_log, 1);
		return ((PrismNen > 6) ? nodal_test_prism2(fp_log) : nodal_test_prism(fp_log));
	}
	if (MeshElem == MESHELEM_PYR)   return nodal_test_pyr(fp_log);
	if (MeshElem == MESHELEM_MIXED) return ((TetOrder >= 2) ? nodal_test_mixed_ho(fp_log) : nodal_test_mixed(fp_log));

	fprintf(fp_log, "\n=== nodal element (P%d) self test ===\n", TetOrder);
	fprintf(fp_log, "  nodes = %d, tetrahedra = %d, nodes per element = %d\n",
		NNode, NTet, tet_nen());

	// 中間節点が辺の中点からどれだけずれているか (曲がった格子の判定)
	double curve = 0;
	if (TetOrder >= 3) {
		// 3 次 / 4 次 : 頂点以外の節点の、直線要素での位置 Σ (a_v/p) x_v からのずれ
		for (int e = 0; e < NTet; e++) {
			int32_t nd[TET_NMAX];
			const int nen = tet_nodes(e, nd);
			double h = 0;
			for (int a = 0; a < 4; a++) {
				for (int b = a + 1; b < 4; b++) {
					const double dx = Xp[nd[b]] - Xp[nd[a]], dy = Yp[nd[b]] - Yp[nd[a]];
					const double dz = Zp[nd[b]] - Zp[nd[a]];
					const double len = sqrt((dx * dx) + (dy * dy) + (dz * dz));
					if (len > h) h = len;
				}
			}
			if (h <= 0) continue;
			for (int l = 4; l < nen; l++) {
				double r[3] = {0, 0, 0};
				const int *a = tet_ho_a(l);
				for (int v = 0; v < 4; v++) {
					r[0] += a[v] * Xp[nd[v]] / TetOrder;
					r[1] += a[v] * Yp[nd[v]] / TetOrder;
					r[2] += a[v] * Zp[nd[v]] / TetOrder;
				}
				const double dx = Xp[nd[l]] - r[0], dy = Yp[nd[l]] - r[1], dz = Zp[nd[l]] - r[2];
				const double d = sqrt((dx * dx) + (dy * dy) + (dz * dz)) / h;
				if (d > curve) curve = d;
			}
		}
		fprintf(fp_log, "  mid-node offset = %.3e (relative to the edge length)\n", curve);
	}
	else if (TetOrder >= 2) {
		static const int ed[6][2] = {{0, 1}, {1, 2}, {2, 0}, {3, 0}, {3, 2}, {3, 1}};
		for (int e = 0; e < NTet; e++) {
			int32_t nd[TET_NMAX];
			tet_nodes(e, nd);
			double h = 0;
			for (int l = 0; l < 6; l++) {
				const int32_t a = nd[ed[l][0]], b = nd[ed[l][1]];
				const double dx = Xp[b] - Xp[a], dy = Yp[b] - Yp[a], dz = Zp[b] - Zp[a];
				const double len = sqrt((dx * dx) + (dy * dy) + (dz * dz));
				if (len > h) h = len;
			}
			if (h <= 0) continue;
			for (int l = 0; l < 6; l++) {
				const int32_t a = nd[ed[l][0]], b = nd[ed[l][1]], c = nd[4 + l];
				const double dx = Xp[c] - ((Xp[a] + Xp[b]) / 2);
				const double dy = Yp[c] - ((Yp[a] + Yp[b]) / 2);
				const double dz = Zp[c] - ((Zp[a] + Zp[b]) / 2);
				const double d = sqrt((dx * dx) + (dy * dy) + (dz * dz)) / h;
				if (d > curve) curve = d;
			}
		}
		fprintf(fp_log, "  mid-node offset = %.3e (relative to the edge length)\n", curve);
	}
	const int straight = (curve < 1e-9);

	// 積分した体積 (2 次では等パラメトリック写像のヤコビアンの検証になる)
	double vsum = 0;
	double *vole = (double *)malloc((size_t)NTet * sizeof(double));
	{
		const double c1[6] = {1, 1, 1, 0, 0, 0};
		for (int e = 0; e < NTet; e++) {
			double vol = 0;
			if (TetOrder >= 3) {
				double ke[TET_NMAX][TET_NMAX];
				if (tet20_element(e, c1, ke, &vol)) {
					fprintf(fp_log, "*** tetrahedron %d is degenerate or inverted\n", e + 1);
					free(vole);
					return 1;
				}
			}
			else if (TetOrder >= 2) {
				double ke[10][10];
				if (tet10_element(e, c1, ke, &vol)) {
					fprintf(fp_log, "*** tetrahedron %d is degenerate or inverted\n", e + 1);
					free(vole);
					return 1;
				}
			}
			else {
				double g[4][3];
				if (tet_grad_pub(&Tet[e * 4], g, &vol)) vol = 0;
			}
			vole[e] = vol;
			vsum += vol;
		}
	}
	fprintf(fp_log, "  integrated volume = %.10e [m^3]\n", vsum);

	crs_t A;
	crs_alloc_tet(&A);
	assemble_tet(&A, 0);			// 誘電率テンソル (材料に anisotropy を持たせること)
	fprintf(fp_log, "  matrix : %lld nonzeros (%.1f per row)\n",
		(long long)A.nnz, (double)A.nnz / ((NNode > 0) ? NNode : 1));
	fflush(fp_log);

	const int n = NNode;
	double *phi = (double *)malloc((size_t)n * sizeof(double));
	double *y = (double *)malloc((size_t)n * sizeof(double));

	double dmax = 0;
	for (int i = 0; i < n; i++) {
		for (int64_t p = A.rowptr[i]; p < A.rowptr[i + 1]; p++) {
			if ((A.col[p] == i) && (fabs(A.val[p]) > dmax)) dmax = fabs(A.val[p]);
		}
	}

	// (a) 定数の零空間 : K 1 = 0 (どんな格子・次数でも成り立つ)
	{
		for (int i = 0; i < n; i++) phi[i] = 1;
		crs_spmv(&A, phi, y, NULL);
		double amax = 0;
		for (int i = 0; i < n; i++) {
			if (fabs(y[i]) > amax) amax = fabs(y[i]);
		}
		const double rel = amax / ((dmax > 0) ? dmax : 1);
		fprintf(fp_log, "  (a) constant null space : max|K 1| / max|Kii| = %.3e\n", rel);
		if (!(rel <= 1e-10)) {		// NaN も不合格にする
			fprintf(fp_log, "*** the stiffness matrix does not annihilate constants\n");
			ierr = 1;
		}
	}

	// (b) 1 次の φ、(c) 2 次の φ
	// 軸に平行でない向き・非対角成分の入った B を選ぶ (テンソルの成分順序や
	// 係数 2 の誤りを見逃さないため)
	const double av[3] = {0.7, -1.3, 2.1};
	const double bm[3][3] = {{ 3.0, -1.1,  0.6},
	                         {-1.1,  2.2,  1.7},
	                         { 0.6,  1.7, -0.9}};

	// 3 次の項 κxyz + μx³ + ν y²z (係数は領域の大きさ L で O(1) の勾配になるよう 1/L² 倍)
	double lref = 0;
	for (int i = 0; i < n; i++) {
		const double r = fabs(Xp[i]) + fabs(Yp[i]) + fabs(Zp[i]);
		if (r > lref) lref = r;
	}
	if (lref <= 0) lref = 1;
	const double ck = 2.3 / (lref * lref), cmu = 1.1 / (lref * lref), cnu = -1.7 / (lref * lref);
	// 4 次の項 ρ x²y² + τ z⁴ + ω x y z² (1/L³ 倍)
	const double crho = 0.9 / (lref * lref * lref), ctau = -0.6 / (lref * lref * lref);
	const double comg = 1.4 / (lref * lref * lref);
	// 5 次の項 β x³y² + γ z⁵ + η x y z³ (1/L⁴ 倍)
	const double l4 = lref * lref * lref * lref;
	const double cbet = 0.8 / l4, cgam = -0.5 / l4, ceta = 1.3 / l4;

	for (int cs = 0; cs < 7; cs++) {
		const int quad = (cs >= 1);		// 0 : φ は 1 次、1 : 2 次、2 : 3 次、3 : 4 次、4 : 5 次、5 : 6 次、6 : 7 次
		const int cubic = (cs >= 2);
		const int quart = (cs >= 3);
		const int quint = (cs >= 4);
		const int ohi = ((cs >= 5) ? (cs + 1) : 0);		// 6 次以上の項の次数 (ho_extra)
		if (quad && ((TetOrder < 2) || !straight)) continue;
		if (cubic && (TetOrder < 3)) continue;
		if (quart && (TetOrder < 4)) continue;
		if (quint && (TetOrder < 5)) continue;
		if (ohi > TetOrder) continue;

		// 節点値
		for (int i = 0; i < n; i++) {
			const double r[3] = {Xp[i], Yp[i], Zp[i]};
			double f = (av[0] * r[0]) + (av[1] * r[1]) + (av[2] * r[2]);
			if (quad) {
				for (int k = 0; k < 3; k++) {
					for (int l = 0; l < 3; l++) f += bm[k][l] * r[k] * r[l] / 2;
				}
			}
			if (cubic) {
				f += (ck * r[0] * r[1] * r[2]) + (cmu * r[0] * r[0] * r[0])
				   + (cnu * r[1] * r[1] * r[2]);
			}
			if (quart) {
				f += (crho * r[0] * r[0] * r[1] * r[1]) + (ctau * r[2] * r[2] * r[2] * r[2])
				   + (comg * r[0] * r[1] * r[2] * r[2]);
			}
			if (quint) {
				f += (cbet * r[0] * r[0] * r[0] * r[1] * r[1]) + (cgam * r[2] * r[2] * r[2] * r[2] * r[2])
				   + (ceta * r[0] * r[1] * r[2] * r[2] * r[2]);
			}
			if (ohi) {
				double gd[3] = {0, 0, 0};
				f += ho_extra(ohi, lref, r, gd);
			}
			phi[i] = f;
		}

		// 厳密値 Σ_e ∫ (∇φ)^T C (∇φ) dV
		//
		// 体積は「直線要素の閉形式」を使う。こうすると 1 次の φ の検査でも
		// 数値積分の体積が独立に検証される。曲がった格子ではその閉形式が
		// 内接多角形の体積になってしまう (実測 -4.5%) ので、そこだけは
		// 数値積分の体積を使い、検査の意味を「∇N が定数場を再現するか」に
		// 絞る (曲がった写像の体積は解析値との比較で別に見る)
		double exact = 0;
		for (int e = 0; e < NTet; e++) {
			const int32_t *nd = &Tet[e * 4];
			double g[4][3], vol;
			if (tet_grad_pub(nd, g, &vol)) continue;
			if (!straight) vol = vole[e];
			double cf[6];
			material_coef_pub(TetMat[e], 0, cf);
			const double cm[3][3] = {{cf[0], cf[3], cf[5]},
			                         {cf[3], cf[1], cf[4]},
			                         {cf[5], cf[4], cf[2]}};

			// 座標のモーメント
			double s[3] = {0, 0, 0}, sq[3][3];
			for (int k = 0; k < 3; k++) {
				for (int a = 0; a < 4; a++) {
					const double *p = ((k == 0) ? Xp : (k == 1) ? Yp : Zp);
					s[k] += p[nd[a]];
				}
			}
			for (int k = 0; k < 3; k++) {
				for (int l = 0; l < 3; l++) {
					const double *pk = ((k == 0) ? Xp : (k == 1) ? Yp : Zp);
					const double *pl = ((l == 0) ? Xp : (l == 1) ? Yp : Zp);
					double t = 0;
					for (int a = 0; a < 4; a++) t += pk[nd[a]] * pl[nd[a]];
					sq[k][l] = ((s[k] * s[l]) + t) / 20;		// ∫ r_k r_l dV / V
				}
			}
			// ∫ g_k g_l dV / V  (g = a + B r、2 次の φ でないときは B = 0)
			for (int k = 0; k < 3; k++) {
				for (int l = 0; l < 3; l++) {
					double gg = av[k] * av[l];
					if (quad) {
						for (int p = 0; p < 3; p++) {
							gg += ((av[k] * bm[l][p]) + (av[l] * bm[k][p])) * s[p] / 4;
							for (int q = 0; q < 3; q++) {
								gg += bm[k][p] * bm[l][q] * sq[p][q];
							}
						}
					}
					exact += cm[k][l] * gg * vol;
				}
			}
		}

		// 3 次の φ : ∫ (∇φ)^T C (∇φ) dV を直線要素の上で**別の積分則** (Duffy + 5 点、
		// 125 点) で厳密に積分する (被積分関数は λ の 4 次)。要素行列 (4 点) とも
		// 形状関数とも独立に、解析的な ∇φ だけを使う。4 次の φ (被積分関数は 6 次、
		// Duffy で u の 8 次) は 6 点 (216 点) で、要素行列の 5 点とは別の則
		// 5 次の φ (被積分関数は 8 次、Duffy で u の 10 次) は 7 点 (343 点)、要素行列の 6 点とは別の則
		if (cubic) {
			static double lq[NQTETMAX][4], wq5[NQTETMAX];
			// 6 次 / 7 次の φ (被積分関数は 10 / 12 次、Duffy で u の 12 / 14 次) は 8 / 9 点
			const int nq5 = tet_quad_n((ohi ? (ohi + 2) : quint ? 7 : quart ? 6 : 5), lq, wq5);
			exact = 0;
			for (int e = 0; e < NTet; e++) {
				const int32_t *nd = &Tet[e * 4];
				double g[4][3], vol;
				if (tet_grad_pub(nd, g, &vol)) continue;
				double cf[6];
				material_coef_pub(TetMat[e], 0, cf);
				const double cm[3][3] = {{cf[0], cf[3], cf[5]},
				                         {cf[3], cf[1], cf[4]},
				                         {cf[5], cf[4], cf[2]}};
				for (int qq = 0; qq < nq5; qq++) {
					double r[3] = {0, 0, 0};
					for (int v = 0; v < 4; v++) {
						r[0] += lq[qq][v] * Xp[nd[v]];
						r[1] += lq[qq][v] * Yp[nd[v]];
						r[2] += lq[qq][v] * Zp[nd[v]];
					}
					double gr[3];
					for (int k = 0; k < 3; k++) {
						gr[k] = av[k];
						for (int l = 0; l < 3; l++) gr[k] += bm[k][l] * r[l];
					}
					gr[0] += (ck * r[1] * r[2]) + (3 * cmu * r[0] * r[0]);
					gr[1] += (ck * r[0] * r[2]) + (2 * cnu * r[1] * r[2]);
					gr[2] += (ck * r[0] * r[1]) + (cnu * r[1] * r[1]);
					if (quart) {
						gr[0] += (2 * crho * r[0] * r[1] * r[1]) + (comg * r[1] * r[2] * r[2]);
						gr[1] += (2 * crho * r[0] * r[0] * r[1]) + (comg * r[0] * r[2] * r[2]);
						gr[2] += (4 * ctau * r[2] * r[2] * r[2]) + (2 * comg * r[0] * r[1] * r[2]);
					}
					if (quint) {
						gr[0] += (3 * cbet * r[0] * r[0] * r[1] * r[1]) + (ceta * r[1] * r[2] * r[2] * r[2]);
						gr[1] += (2 * cbet * r[0] * r[0] * r[0] * r[1]) + (ceta * r[0] * r[2] * r[2] * r[2]);
						gr[2] += (5 * cgam * r[2] * r[2] * r[2] * r[2]) + (3 * ceta * r[0] * r[1] * r[2] * r[2]);
					}
					if (ohi) (void)ho_extra(ohi, lref, r, gr);
					double gg = 0;
					for (int k = 0; k < 3; k++) {
						for (int l = 0; l < 3; l++) gg += cm[k][l] * gr[k] * gr[l];
					}
					exact += 6 * vol * wq5[qq] * gg;
				}
			}
		}

		crs_spmv(&A, phi, y, NULL);
		double q = 0;
		for (int i = 0; i < n; i++) q += phi[i] * y[i];
		const double rel = fabs(q - exact) / ((exact != 0) ? fabs(exact) : 1);
		fprintf(fp_log, "  (%c) %s field : phi^T K phi = %.10e, exact = %.10e, err = %.3e\n",
			((ohi == 7) ? 'j' : (ohi == 6) ? 'i' : quint ? 'h' : quart ? 'g' : cubic ? 'f' : quad ? 'c' : 'b'),
			((ohi == 7) ? "septic   " : (ohi == 6) ? "sextic   " : quint ? "quintic  " : quart ? "quartic  " : cubic ? "cubic    " : quad ? "quadratic" : "linear   "),
			q, exact, rel);
		if (!(rel <= 1e-10)) {		// NaN も不合格にする
			fprintf(fp_log, "*** the %s stiffness matrix is wrong\n",
				((TetOrder == 7) ? "P7" : (TetOrder == 6) ? "P6" : (TetOrder == 5) ? "P5" : (TetOrder == 4) ? "P4"
				 : (TetOrder == 3) ? "P3" : (TetOrder >= 2) ? "P2" : "P1"));
			ierr = 1;
		}

		// (e) 重心での勾配 (場の出力が使う経路)
		//
		// 場の出力は要素あたり 1 本のベクトルなので重心で評価する。この評価点は
		// **一様な場では検証できない** (要素内のどこで取っても同じ値になる)。
		// 実測 : 平行平板の場の出力だけでは、評価点を頂点にずらす変異が素通りした。
		// 2 次の φ なら ∇φ = a + B r が場所で変わるので、重心以外を選ぶと落ちる。
		// 直線要素では重心の物理座標が頂点の平均に一致する (形状関数の重みが
		// 頂点 -1/8、中間節点 +1/4 で、中間節点が中点にあるとき Σ = 平均になる)
		{
			double gmax = 0, amax2 = 0;
			for (int e = 0; e < NTet; e++) {
				double gn[TET_NMAX][3];
				int nen = 0;
				if (tet_grad_center(e, gn, &nen)) continue;
				int32_t nd[TET_NMAX];
				tet_nodes(e, nd);

				double rc[3] = {0, 0, 0};
				for (int a = 0; a < 4; a++) {
					rc[0] += Xp[nd[a]] / 4;
					rc[1] += Yp[nd[a]] / 4;
					rc[2] += Zp[nd[a]] / 4;
				}
				double hg[3] = {0, 0, 0};
				if (ohi) (void)ho_extra(ohi, lref, rc, hg);
				for (int k = 0; k < 3; k++) {
					double gh = 0;
					for (int a = 0; a < nen; a++) gh += phi[nd[a]] * gn[a][k];
					double ge = av[k] + hg[k];
					if (quad) {
						for (int l = 0; l < 3; l++) ge += bm[k][l] * rc[l];
					}
					if (cubic) {
						ge += ((k == 0) ? ((ck * rc[1] * rc[2]) + (3 * cmu * rc[0] * rc[0]))
						     : (k == 1) ? ((ck * rc[0] * rc[2]) + (2 * cnu * rc[1] * rc[2]))
						     : ((ck * rc[0] * rc[1]) + (cnu * rc[1] * rc[1])));
					}
					if (quart) {
						ge += ((k == 0) ? ((2 * crho * rc[0] * rc[1] * rc[1]) + (comg * rc[1] * rc[2] * rc[2]))
						     : (k == 1) ? ((2 * crho * rc[0] * rc[0] * rc[1]) + (comg * rc[0] * rc[2] * rc[2]))
						     : ((4 * ctau * rc[2] * rc[2] * rc[2]) + (2 * comg * rc[0] * rc[1] * rc[2])));
					}
					if (quint) {
						ge += ((k == 0) ? ((3 * cbet * rc[0] * rc[0] * rc[1] * rc[1]) + (ceta * rc[1] * rc[2] * rc[2] * rc[2]))
						     : (k == 1) ? ((2 * cbet * rc[0] * rc[0] * rc[0] * rc[1]) + (ceta * rc[0] * rc[2] * rc[2] * rc[2]))
						     : ((5 * cgam * rc[2] * rc[2] * rc[2] * rc[2]) + (3 * ceta * rc[0] * rc[1] * rc[2] * rc[2])));
					}
					const double dg = fabs(gh - ge);
					if ((dg > gmax) || (dg != dg)) gmax = dg;		// NaN は残す (以後の比較で上書きされない)
					if (fabs(ge) > amax2) amax2 = fabs(ge);
				}
			}
			const double relg = gmax / ((amax2 > 0) ? amax2 : 1);
			fprintf(fp_log, "  (e) centroid gradient (%s field) : max err = %.3e\n",
				((ohi == 7) ? "septic" : (ohi == 6) ? "sextic" : quint ? "quintic" : quart ? "quartic" : cubic ? "cubic" : quad ? "quadratic" : "linear"), relg);
			if (!(relg <= 1e-10)) {		// NaN も不合格にする
				fprintf(fp_log, "*** the element-centre gradient is wrong "
					"(this is what fieldout writes)\n");
				ierr = 1;
			}
		}
	}

	if ((TetOrder >= 2) && !straight) {
		fprintf(fp_log, "  (c) quadratic field : skipped (the mesh is curved, so the "
			"P2 interpolant of a quadratic is not exact)\n");
	}

	// (d) 対称性
	{
		double amax = 0, adif = 0;
		for (int i = 0; i < n; i++) {
			for (int64_t p = A.rowptr[i]; p < A.rowptr[i + 1]; p++) {
				const int32_t j = A.col[p];
				const int64_t pj = crs_find(&A, j, (int32_t)i);
				if (pj < 0) continue;
				if (fabs(A.val[p]) > amax) amax = fabs(A.val[p]);
				const double d = fabs(A.val[p] - A.val[pj]);
				if (d > adif) adif = d;
			}
		}
		const double rel = adif / ((amax > 0) ? amax : 1);
		fprintf(fp_log, "  (d) symmetry : max|Kij - Kji| / max|Kij| = %.3e\n", rel);
		if (!(rel <= 1e-12)) {		// NaN も不合格にする
			fprintf(fp_log, "*** the stiffness matrix is not symmetric\n");
			ierr = 1;
		}
	}

	free(phi);
	free(y);
	free(vole);
	crs_free(&A);

	fprintf(fp_log, "  result : %s\n", (ierr ? "FAILED" : "passed"));

	return ierr;
}


// ---- 断面 2 次元の三角形要素 (MeshDim == 2) ----
//
// M / F は「伝送線路軸 t に垂直な断面での 2 次元問題」なので、四面体ではなく
// 三角形で切った格子に載る。未知数は Az (軸方向のベクトルポテンシャル成分) の
// 節点値で、面内の 2 軸 (p, q) だけが微分に効く。
//
// 1 次三角形は ∇λ が要素内一定なので、剛性は面積を掛けるだけで厳密:
//   K_ij = ν S (∇λ_i・∇λ_j)
// 質量は ∫λ_iλ_j dS = S(1+δ_ij)/12 で厳密:
//   M_ij = σ S (1 + δ_ij)/12
//
// 単位長あたりの量として扱うので「体積」は面積そのもの (線路長 1 m 相当)。

// 三角形の局所節点をまとめて取り出す。戻り値は節点数 (1 次 3、2 次 6、3 次 10)。
// 並びは Gmsh の tri6 と同じ (頂点 3 個のあと辺 (0,1) (1,2) (2,0) の中間節点)。
// 3 次は gmsh の tri10 (多重指数は TRI3_A、実測)
int tri_nodes(int e, int32_t nd[TRI_NMAX])
{
	for (int l = 0; l < 3; l++) nd[l] = Tri[(e * 3) + l];
	if (TetOrder < 2) return 3;
	if ((TetOrder >= 3) && (MeshDim == 2)) {
		const int nho = (((TetOrder + 1) * (TetOrder + 2)) / 2) - 3;
		for (int l = 0; l < nho; l++) nd[3 + l] = Tri3[(e * TRI_HO_MAX) + l];
		return 3 + nho;
	}
	for (int l = 0; l < 3; l++) nd[3 + l] = Tri2[(e * 3) + l];

	return 6;
}


// 面内の 2 軸 (伝送線路軸 t の次の 2 つ)
void tri_axes(int *p, int *q)
{
	const int t = ((Tline == 'X') ? 0 : (Tline == 'Y') ? 1 : 2);

	*p = (t + 1) % 3;
	*q = (t + 2) % 3;
}


// 三角形の面内勾配 ∇λ (要素内一定) と面積。戻り値 : 0 = 正常、1 = 退化
int tri_grad(const int32_t nd[3], double g[3][2], double *area)
{
	int p, q;
	tri_axes(&p, &q);
	const double *cp = ((p == 0) ? Xp : (p == 1) ? Yp : Zp);
	const double *cq = ((q == 0) ? Xp : (q == 1) ? Yp : Zp);

	const double p0 = cp[nd[0]], q0 = cq[nd[0]];
	const double p1 = cp[nd[1]] - p0, q1 = cq[nd[1]] - q0;
	const double p2 = cp[nd[2]] - p0, q2 = cq[nd[2]] - q0;

	const double det = (p1 * q2) - (p2 * q1);
	if (det == 0) return 1;

	// ∇λ1, ∇λ2 は J^-1 の行、∇λ0 = -(∇λ1 + ∇λ2)
	g[1][0] =  q2 / det;  g[1][1] = -p2 / det;
	g[2][0] = -q1 / det;  g[2][1] =  p1 / det;
	g[0][0] = -(g[1][0] + g[2][0]);
	g[0][1] = -(g[1][1] + g[2][1]);
	*area = fabs(det) / 2;

	return 0;
}


// 行 row の中で列 col の位置 (crs_find の公開版)
int64_t crs_find_tri(const crs_t *A, int32_t row, int32_t col)
{
	return crs_find(A, row, col);
}


// 三角形の連結から節点の隣接関係を作る
void crs_alloc_tri(crs_t *A)
{
	const int n = NNode;

	const int nen = ((TetOrder >= 2) ? (((TetOrder + 1) * (TetOrder + 2)) / 2) : 3);
	int *cnt = (int *)malloc((size_t)n * sizeof(int));
	memset(cnt, 0, (size_t)n * sizeof(int));
	for (int e = 0; e < NTri; e++) {
		int32_t nd[TRI_NMAX];
		tri_nodes(e, nd);
		for (int l = 0; l < nen; l++) cnt[nd[l]]++;
	}
	int64_t *nptr = (int64_t *)malloc(((size_t)n + 1) * sizeof(int64_t));
	nptr[0] = 0;
	for (int i = 0; i < n; i++) nptr[i + 1] = nptr[i] + cnt[i];
	int32_t *nlist = (int32_t *)malloc((size_t)nptr[n] * sizeof(int32_t));
	memset(cnt, 0, (size_t)n * sizeof(int));
	for (int e = 0; e < NTri; e++) {
		int32_t nd[TRI_NMAX];
		tri_nodes(e, nd);
		for (int l = 0; l < nen; l++) {
			const int32_t i = nd[l];
			nlist[nptr[i] + cnt[i]] = (int32_t)e;
			cnt[i]++;
		}
	}

	A->n = n;
	A->rowptr = (int64_t *)malloc(((size_t)n + 1) * sizeof(int64_t));
	int cap = 64;
	int32_t *work = (int32_t *)malloc((size_t)cap * sizeof(int32_t));
	int *rown = (int *)malloc((size_t)n * sizeof(int));

	for (int pass = 0; pass < 2; pass++) {
		if (pass == 1) {
			A->rowptr[0] = 0;
			for (int i = 0; i < n; i++) A->rowptr[i + 1] = A->rowptr[i] + rown[i];
			A->nnz = A->rowptr[n];
			A->col = (int32_t *)malloc((size_t)A->nnz * sizeof(int32_t));
			A->val = (double *)malloc((size_t)A->nnz * sizeof(double));
		}
		for (int i = 0; i < n; i++) {
			const int64_t p0 = nptr[i], p1 = nptr[i + 1];
			const int need = (int)(p1 - p0) * nen;
			if (need > cap) {
				cap = need;
				work = (int32_t *)realloc(work, (size_t)cap * sizeof(int32_t));
			}
			int m = 0;
			for (int64_t p = p0; p < p1; p++) {
				int32_t nd[TRI_NMAX];
				tri_nodes(nlist[p], nd);
				for (int l = 0; l < nen; l++) work[m++] = nd[l];
			}
			qsort(work, (size_t)m, sizeof(int32_t), cmp_int32);
			if (pass == 0) {
				int u = 0;
				for (int q = 0; q < m; q++) {
					if ((q == 0) || (work[q] != work[q - 1])) u++;
				}
				rown[i] = u;
			}
			else {
				int64_t w = A->rowptr[i];
				for (int q = 0; q < m; q++) {
					if ((q == 0) || (work[q] != work[q - 1])) A->col[w++] = work[q];
				}
			}
		}
	}

	free(cnt);
	free(nptr);
	free(nlist);
	free(work);
	free(rown);

	crs_zero(A);
}


// ---- 2 次三角形 (6 節点、等パラメトリック) ----
//
// 積分は四面体と同じ考え方で Duffy 変換 + 各方向 3 点 Gauss-Legendre (9 点):
//   λ1 = u, λ2 = v(1-u),  λ0 = 1 - λ1 - λ2,  dλ1 dλ2 = (1-u) du dv
// 重みは全部正で、剛性 (λ について 2 次) も質量 (4 次) も直線要素なら厳密になる。
#define NQTRI (9)

static void tri_quad(double lam[NQTRI][3], double wq[NQTRI])
{
	const double gp[3] = {0.5 - (0.5 * 0.77459666924148337704),
	                      0.5,
	                      0.5 + (0.5 * 0.77459666924148337704)};
	const double gw[3] = {5.0 / 18, 8.0 / 18, 5.0 / 18};

	int q = 0;
	for (int a = 0; a < 3; a++) {
	for (int b = 0; b < 3; b++) {
		const double u = gp[a], v = gp[b];
		lam[q][1] = u;
		lam[q][2] = v * (1 - u);
		lam[q][0] = 1 - lam[q][1] - lam[q][2];
		wq[q] = gw[a] * gw[b] * (1 - u);
		q++;
	}
	}
}


// ---- 3 次三角形 (10 節点、断面 2 次元、等パラメトリック) ----
//
// 局所節点の並びは **gmsh 4.15.2 の tri10 の出力を重心座標で分類した実測の表**
// (多重指数 3λ。tri10_2d.msh の 14 要素がすべてこの並び)。四面体の TET3_A の面と同じ
// 規約で、形状関数も同じ次数によらない Lagrange の式 (lagr_tet に 4 番目の指数 0 を
// 渡す) で作る。積分は Duffy + 各方向 4 点 (16 点)。∇N・∇N は λ の 4 次で (1-u) の
// 因子を含めて u の 5 次、質量 N N は 6 次 + 1 = 7 次なので、直線要素ならどちらも厳密
static const int TRI3_A[10][3] = {
	{3, 0, 0}, {0, 3, 0}, {0, 0, 3},
	{2, 1, 0}, {1, 2, 0},		// 辺 (0,1)
	{0, 2, 1}, {0, 1, 2},		// 辺 (1,2)
	{1, 0, 2}, {2, 0, 1},		// 辺 (2,0)
	{1, 1, 1}					// 面心
};
// 4 次 (15 節点、gmsh の tri15 = 型 23) : gmsh 4.15.2 の出力 (box_p4.msh の電極面 36 枚が
// すべて同じ並び) を分類した実測の表。頂点、辺 (0,1) (1,2) (2,0) に 3 個ずつ、内部 3 個
static const int TRI4_A[15][3] = {
	{4, 0, 0}, {0, 4, 0}, {0, 0, 4},
	{3, 1, 0}, {2, 2, 0}, {1, 3, 0},
	{0, 3, 1}, {0, 2, 2}, {0, 1, 3},
	{1, 0, 3}, {2, 0, 2}, {3, 0, 1},
	{2, 1, 1}, {1, 2, 1}, {1, 1, 2}
};
// 5 次 (21 節点、gmsh の tri21 = 型 25) : gmsh 4.15.2 の出力 (立方体の 540 枚がすべて同じ並び) を分類した実測の表。
// 頂点、辺 (0,1) (1,2) (2,0) に 4 個ずつ、内部 6 個
static const int TRI5_A[21][3] = {
	{5, 0, 0}, {0, 5, 0}, {0, 0, 5},
	{4, 1, 0}, {3, 2, 0}, {2, 3, 0},
	{1, 4, 0}, {0, 4, 1}, {0, 3, 2},
	{0, 2, 3}, {0, 1, 4}, {1, 0, 4},
	{2, 0, 3}, {3, 0, 2}, {4, 0, 1},
	{3, 1, 1}, {1, 3, 1}, {1, 1, 3},
	{2, 2, 1}, {1, 2, 2}, {2, 1, 2}
};
#define NQTRI3 (QTRI_MAX)

// 三角形の再帰の規則 (tri_rule) が実測の表 TRI3_A / TRI4_A / TRI5_A と全項一致するか (0 = 一致)
static int tri_rule_check(void)
{
	int t3[TRI_NMAX][3];
	for (int q = 3; q <= 5; q++) {
		const int nf = tri_rule(q, t3);
		const int (*mf)[3] = ((q == 5) ? TRI5_A : (q == 4) ? TRI4_A : TRI3_A);
		if (nf != (((q + 1) * (q + 2)) / 2)) return 1;
		for (int i = 0; i < nf; i++) {
			for (int c = 0; c < 3; c++) {
				if (t3[i][c] != mf[i][c]) return 1;
			}
		}
	}
	return 0;
}

// 積分点 (2 次は 9 点、3 次は 16 点)。戻り値は点の数
static int tri_quad_ord(double lam[NQTRI3][3], double wq[NQTRI3])
{
	if (TetOrder < 3) {
		double l9[NQTRI][3], w9[NQTRI];
		tri_quad(l9, w9);
		for (int q = 0; q < NQTRI; q++) {
			for (int a = 0; a < 3; a++) lam[q][a] = l9[q][a];
			wq[q] = w9[q];
		}
		return NQTRI;
	}
	// 3 次 : 4 点 (16 点)。4 次 : ∇N・∇N は 6 次で u の 7 次、質量は 8 次で u の 9 次 -> 5 点 (25 点)
	const int ng = TetOrder + 1;		// 5 次 : 6 点 (36 点)
	double gp[GL_NMAX], gw[GL_NMAX];
	gauss01(ng, gp, gw);
	int q = 0;
	for (int a = 0; a < ng; a++) {
	for (int b = 0; b < ng; b++) {
		const double u = gp[a], v = gp[b];
		lam[q][1] = u;
		lam[q][2] = v * (1 - u);
		lam[q][0] = 1 - lam[q][1] - lam[q][2];
		wq[q] = gw[a] * gw[b] * (1 - u);
		q++;
	}
	}
	return q;
}


// 6 節点三角形の形状関数と λ 微分 dl[i][a] = ∂N_i/∂λ_a
static void tri6_shape(const double lam[3], double n[6], double dl[6][3])
{
	static const int ed[3][2] = {{0, 1}, {1, 2}, {2, 0}};	// Gmsh の tri6

	memset(dl, 0, sizeof(double) * 6 * 3);
	for (int a = 0; a < 3; a++) {
		n[a] = lam[a] * ((2 * lam[a]) - 1);
		dl[a][a] = (4 * lam[a]) - 1;
	}
	for (int l = 0; l < 3; l++) {
		const int a = ed[l][0], b = ed[l][1];
		n[3 + l] = 4 * lam[a] * lam[b];
		dl[3 + l][a] = 4 * lam[b];
		dl[3 + l][b] = 4 * lam[a];
	}
}


// 6 / 10 節点三角形の形状関数と λ 微分 (nen で切り替え)
static void trin_shape(int nen, const double lam[3], double n[TRI_NMAX], double dl[TRI_NMAX][3])
{
	if (nen == 6) {
		double n6[6], d6[6][3];
		tri6_shape(lam, n6, d6);
		for (int i = 0; i < 6; i++) {
			n[i] = n6[i];
			for (int a = 0; a < 3; a++) dl[i][a] = d6[i][a];
		}
		return;
	}
	const double l4[4] = {lam[0], lam[1], lam[2], 0};
	const int p = ((nen == 36) ? 7 : (nen == 28) ? 6 : (nen == 21) ? 5 : (nen == 15) ? 4 : 3);
	for (int i = 0; i < nen; i++) {
		const int *ta = ((p >= 6) ? TRIG_A[i] : (p == 5) ? TRI5_A[i] : (p == 4) ? TRI4_A[i] : TRI3_A[i]);
		const int a4[4] = {ta[0], ta[1], ta[2], 0};
		double dn[4];
		n[i] = lagr_tet(p, a4, l4, dn);
		for (int a = 0; a < 3; a++) dl[i][a] = dn[a];
	}
}


// 積分点 1 点での面内勾配 ∇N_i とヤコビアン行列式 (nen = 6 / 10)。戻り値 : 0 = 正常
static int tri6_grad(const int32_t *nd, int nen, const double lam[3],
	double gn[TRI_NMAX][2], double *det)
{
	int p, q;
	tri_axes(&p, &q);
	const double *cp = ((p == 0) ? Xp : (p == 1) ? Yp : Zp);
	const double *cq = ((q == 0) ? Xp : (q == 1) ? Yp : Zp);

	double n[TRI_NMAX], dl[TRI_NMAX][3], dx[TRI_NMAX][2];
	trin_shape(nen, lam, n, dl);
	// 参照座標 ξ = (λ1, λ2) についての微分 (λ0 = 1 - λ1 - λ2)
	for (int i = 0; i < nen; i++) {
		for (int k = 0; k < 2; k++) dx[i][k] = dl[i][k + 1] - dl[i][0];
	}

	double j[2][2] = {{0, 0}, {0, 0}};
	for (int k = 0; k < 2; k++) {
		for (int i = 0; i < nen; i++) {
			j[0][k] += cp[nd[i]] * dx[i][k];
			j[1][k] += cq[nd[i]] * dx[i][k];
		}
	}
	const double d = (j[0][0] * j[1][1]) - (j[0][1] * j[1][0]);
	if (d == 0) return 1;

	// Jinv[k][r] = ∂ξ_k/∂x_r
	const double ji[2][2] = {{ j[1][1] / d, -j[0][1] / d},
	                         {-j[1][0] / d,  j[0][0] / d}};
	for (int i = 0; i < nen; i++) {
		for (int r = 0; r < 2; r++) {
			gn[i][r] = (dx[i][0] * ji[0][r]) + (dx[i][1] * ji[1][r]);
		}
	}
	*det = d;

	return 0;
}


// 断面 2 次元のヒステリシス (J-A) の積分点。1 次は重心の 1 点 (B が要素内一定なので厳密、重み = 面積)、高次は剛性と
// 同じ則 (tri_quad_ord) の点ごとに状態を持つ。状態の枠は三角形あたり tri_ja_stride() 個 (8 の倍数、ja_alloc の単位)
static double JqLam[NQTRI3][3], JqW[NQTRI3];
static int JqN = 0, JqOrd = -1;

int tri_ja_nip(void)
{
	if (TetOrder < 2) return 1;
	if (JqOrd != TetOrder) {
		JqN = tri_quad_ord(JqLam, JqW);
		JqOrd = TetOrder;
	}
	return JqN;
}

int tri_ja_stride(void)
{
	return 8 * ((tri_ja_nip() + 7) / 8);
}

// 三角形 e の積分点 q : 節点 nd、面内勾配 gn (∇N、p, q 成分)、重み w (|det J| 込み)、面内の位置 pos。戻り値は節点数 (0 = 退化)
int tri_ja_ip(int e, int q, int32_t nd[TRI_NMAX], double gn[TRI_NMAX][2], double *w, double pos[2])
{
	int pa, qa;
	tri_axes(&pa, &qa);
	const double *cp = ((pa == 0) ? Xp : (pa == 1) ? Yp : Zp);
	const double *cq = ((qa == 0) ? Xp : (qa == 1) ? Yp : Zp);
	if (tri_ja_nip() == 1) {
		double g[3][2], area;
		for (int l = 0; l < 3; l++) nd[l] = Tri[(e * 3) + l];
		if (tri_grad(nd, g, &area)) return 0;
		for (int l = 0; l < 3; l++) { gn[l][0] = g[l][0]; gn[l][1] = g[l][1]; }
		*w = area;
		pos[0] = (cp[nd[0]] + cp[nd[1]] + cp[nd[2]]) / 3;
		pos[1] = (cq[nd[0]] + cq[nd[1]] + cq[nd[2]]) / 3;
		return 3;
	}
	const int nen = tri_nodes(e, nd);
	double det, n[TRI_NMAX], dl[TRI_NMAX][3];
	if (tri6_grad(nd, nen, JqLam[q], gn, &det)) return 0;
	trin_shape(nen, JqLam[q], n, dl);
	*w = JqW[q] * fabs(det);
	pos[0] = pos[1] = 0;
	for (int l = 0; l < nen; l++) {
		pos[0] += n[l] * cp[nd[l]];
		pos[1] += n[l] * cq[nd[l]];
	}
	return nen;
}


// 2 次 / 3 次三角形の面積 (等パラメトリック写像で積分する)
double tri6_area(int e)
{
	int32_t nd[TRI_NMAX];
	const int nen = tri_nodes(e, nd);
	if (nen < 6) return 0;

	double lam[NQTRI3][3], wq[NQTRI3], a = 0;
	const int nq = tri_quad_ord(lam, wq);
	for (int q = 0; q < nq; q++) {
		double gn[TRI_NMAX][2], det;
		if (tri6_grad(nd, nen, lam[q], gn, &det)) return 0;
		a += wq[q] * fabs(det);
	}

	return a;
}


// 3 次三角形の節点の重み ∫N_i dA (等パラメトリック写像で積分) と重心での N_i
// (直線要素では頂点 A/30、辺 3A/40、面心 9A/20 と、重心で頂点・辺 0、面心 1)。
// 戻り値は節点数 (3 次でなければ 0 を返し、呼び出し側は従来の閉形式を使う)
int tri_weights(int e, double w[TRI_NMAX], double wc[TRI_NMAX])
{
	int32_t nd[TRI_NMAX];
	const int nen = tri_nodes(e, nd);
	if (nen < 10) return 0;
	double lam[NQTRI3][3], wq[NQTRI3];
	const int nq = tri_quad_ord(lam, wq);
	for (int l = 0; l < nen; l++) w[l] = 0;
	for (int q = 0; q < nq; q++) {
		double gn[TRI_NMAX][2], det, n[TRI_NMAX], dl[TRI_NMAX][3];
		if (tri6_grad(nd, nen, lam[q], gn, &det)) return 0;
		trin_shape(nen, lam[q], n, dl);
		for (int l = 0; l < nen; l++) w[l] += wq[q] * fabs(det) * n[l];
	}
	const double lc[3] = {1.0 / 3, 1.0 / 3, 1.0 / 3};
	double dl[TRI_NMAX][3];
	trin_shape(nen, lc, wc, dl);
	return nen;
}


// 要素の重心での面内 ∇N_i (場の出力に使う)。1 次では要素内一定なので同じ値
int tri_grad_center(int e, double gn[TRI_NMAX][2], int *nen)
{
	int32_t nd[TRI_NMAX];
	const int n = tri_nodes(e, nd);

	*nen = n;
	if (n == 3) {
		double g[3][2], area;
		if (tri_grad(nd, g, &area)) return 1;
		for (int i = 0; i < 3; i++) {
			gn[i][0] = g[i][0];
			gn[i][1] = g[i][1];
		}
		return 0;
	}

	const double lam[3] = {1.0 / 3, 1.0 / 3, 1.0 / 3};
	double det;

	return tri6_grad(nd, n, lam, gn, &det);
}


// 剛性行列 K_ij = ∫ ν ∇λ_i・∇λ_j dS  (面内の異方性テンソルも扱う)
// nucell が非 NULL なら要素毎の ν を使う (等方性、非線形解析用)。
void assemble_nu_tri(crs_t *A, const double *nucell)
{
	int p, q;
	tri_axes(&p, &q);
	crs_zero(A);

	const int p2 = (TetOrder >= 2);
	double lam[NQTRI3][3], wq[NQTRI3];
	const int nq = (p2 ? tri_quad_ord(lam, wq) : 0);

	for (int e = 0; e < NTri; e++) {
		int32_t nd6[TRI_NMAX];
		const int nen = tri_nodes(e, nd6);
		const int32_t *nd = nd6;
		double g[3][2], area;
		if (!p2 && tri_grad(nd, g, &area)) continue;

		// 面内 2x2 の磁気抵抗率 [[cpp, cpq], [cpq, cqq]]
		double cpp, cqq, cpq;
		if (nucell != NULL) {
			cpp = cqq = nucell[e];
			cpq = 0;
		}
		else {
			// mode 4 = 「∇Az の基底での ν」。B = ∇×(Az ê_t) なので
			// 面内 2 成分の入れ替えと非対角の符号反転が要る (assemble.c 参照)
			double c[6];
			material_coef_pub(TriMat[e], 4, c);
			const double cm[3][3] = {{c[0], c[3], c[5]},
			                         {c[3], c[1], c[4]},
			                         {c[5], c[4], c[2]}};
			cpp = cm[p][p];
			cqq = cm[q][q];
			cpq = cm[p][q];
		}

		if (p2) {
			for (int q = 0; q < nq; q++) {
				double gn[TRI_NMAX][2], det;
				if (tri6_grad(nd, nen, lam[q], gn, &det)) break;
				const double dw = wq[q] * fabs(det);
				for (int l = 0; l < nen; l++) {
					for (int m = 0; m < nen; m++) {
						const double v = (cpp * gn[l][0] * gn[m][0])
						               + (cqq * gn[l][1] * gn[m][1])
						               + (cpq * ((gn[l][0] * gn[m][1]) + (gn[l][1] * gn[m][0])));
						const int64_t s = crs_find(A, nd[l], nd[m]);
						if (s >= 0) A->val[s] += v * dw;
					}
				}
			}
			continue;
		}

		for (int l = 0; l < 3; l++) {
			for (int m = 0; m < 3; m++) {
				const double v = (cpp * g[l][0] * g[m][0])
				               + (cqq * g[l][1] * g[m][1])
				               + (cpq * ((g[l][0] * g[m][1]) + (g[l][1] * g[m][0])));
				const int64_t s = crs_find(A, nd[l], nd[m]);
				if (s >= 0) A->val[s] += v * area;
			}
		}
	}
}


/*
非線形静磁場 (M + bh) の Newton 反復の三角形版 : ヤコビアン J と内力 res = K(ν(|∇A|)) A を作る (assemble_newton と同じ式)
  R_i  = Σ_g w ν(B) (∇N_i・∇A) |det J|
  J_ij = Σ_g w [ ν (∇N_i・∇N_j) + 2 (dν/dB²) (∇N_i・∇A)(∇N_j・∇A) ] |det J|
ν は積分点ごと (高次の三角形は assemble_nu_tri と同じ積分則、1 次は ∇A が要素内一定なので 1 点で厳密)。軸ごとの曲線
(bhaniso) は構造格子と同じく B_p = ∂A/∂q, B_q = −∂A/∂p で軸ごとに分離する (rank-1 項は出ない)。線形の材料は
dν/dB² = 0 で J = K(ν)、R = K(ν)A に退化する
*/
void assemble_newton_tri(crs_t *J, const double *az, double *res)
{
	int p, q;
	tri_axes(&p, &q);
	crs_zero(J);
	memset(res, 0, (size_t)NNode * sizeof(double));

	const int p2 = (TetOrder >= 2);
	double lam[NQTRI3][3], wq[NQTRI3];
	const int nq = (p2 ? tri_quad_ord(lam, wq) : 1);
	double *je = (double *)malloc((size_t)TRI_NMAX * TRI_NMAX * sizeof(double));

	for (int e = 0; e < NTri; e++) {
		int32_t nd[TRI_NMAX];
		const int nen = tri_nodes(e, nd);
		const material_t *mt = &Material[TriMat[e]];
		const int aniso = (mt->bhaniso && (mt->nbh[0] > 0));
		double re[TRI_NMAX];
		for (int l = 0; l < nen; l++) {
			re[l] = 0;
			for (int m = 0; m < nen; m++) je[(l * TRI_NMAX) + m] = 0;
		}
		int bad = 0;
		for (int g = 0; g < nq; g++) {
			double gn[TRI_NMAX][2], dw;
			if (p2) {
				double det;
				if (tri6_grad(nd, nen, lam[g], gn, &det)) { bad = 1; break; }
				dw = wq[g] * fabs(det);
			}
			else {
				double g1[3][2], area;
				if (tri_grad(nd, g1, &area)) { bad = 1; break; }
				for (int l = 0; l < 3; l++) { gn[l][0] = g1[l][0]; gn[l][1] = g1[l][1]; }
				dw = area;
			}
			double gr[2] = {0, 0};
			for (int l = 0; l < nen; l++) {
				gr[0] += az[nd[l]] * gn[l][0];
				gr[1] += az[nd[l]] * gn[l][1];
			}
			if (aniso) {
				//   R_i  = ∫ [ H_p(B_p) ∂N_i/∂q - H_q(B_q) ∂N_i/∂p ] dS
				//   J_ij = ∫ [ (dH_p/dB) ∂N_i/∂q ∂N_j/∂q + (dH_q/dB) ∂N_i/∂p ∂N_j/∂p ] dS
				const double bp = gr[1], bq = -gr[0];
				double hp, dhp, hq, dhq;
				bh_h_dh_pub(mt, p, fabs(bp), &hp, &dhp);
				bh_h_dh_pub(mt, q, fabs(bq), &hq, &dhq);
				if (bp < 0) hp = -hp;
				if (bq < 0) hq = -hq;
				for (int l = 0; l < nen; l++) {
					re[l] += ((hp * gn[l][1]) - (hq * gn[l][0])) * dw;
					for (int m = 0; m < nen; m++) {
						je[(l * TRI_NMAX) + m] += ((dhp * gn[l][1] * gn[m][1]) + (dhq * gn[l][0] * gn[m][0])) * dw;
					}
				}
				continue;
			}
			double nu, dnudb2;
			bh_nu_dnu_pub(mt, sqrt((gr[0] * gr[0]) + (gr[1] * gr[1])), &nu, &dnudb2);
			double u[TRI_NMAX];
			for (int l = 0; l < nen; l++) {
				u[l] = (gn[l][0] * gr[0]) + (gn[l][1] * gr[1]);
				re[l] += nu * u[l] * dw;
			}
			for (int l = 0; l < nen; l++) {
				for (int m = 0; m < nen; m++) {
					const double kk = (gn[l][0] * gn[m][0]) + (gn[l][1] * gn[m][1]);
					je[(l * TRI_NMAX) + m] += ((nu * kk) + (2 * dnudb2 * u[l] * u[m])) * dw;
				}
			}
		}
		if (bad) continue;		// 退化した三角形は setup で弾いてある (assemble_nu_tri と同じ扱い)
		for (int l = 0; l < nen; l++) {
			res[nd[l]] += re[l];
			for (int m = 0; m < nen; m++) {
				const int64_t s = crs_find(J, nd[l], nd[m]);	// mpisolve = 2 : 他の rank の行は -1
				if (s >= 0) J->val[s] += je[(l * TRI_NMAX) + m];
			}
		}
	}
	free(je);
}


// ベクトル J-A の固定点の行列 (1 次の三角形)。ten6 は要素ごとの ∇A の基底のテンソル (xx, yy, zz, xy, yz, zx、assemble_ja_tri が書く)
void assemble_nu_tri6(crs_t *A, const double *ten6)
{
	int p, q;
	tri_axes(&p, &q);
	const int off = (((p + q) == 1) ? 3 : ((p + q) == 3) ? 4 : 5);		// (0,1) -> xy、(1,2) -> yz、(0,2) -> zx
	crs_zero(A);
	for (int e = 0; e < NTri; e++) {
		const int32_t *nd = &Tri[e * 3];
		double g[3][2], area;
		if (tri_grad(nd, g, &area)) continue;
		const double *t = &ten6[(int64_t)e * 6];
		const double cpp = t[p], cqq = t[q], cpq = t[off];
		for (int l = 0; l < 3; l++) {
			for (int m = 0; m < 3; m++) {
				const double v = (cpp * g[l][0] * g[m][0]) + (cqq * g[l][1] * g[m][1])
				               + (cpq * ((g[l][0] * g[m][1]) + (g[l][1] * g[m][0])));
				const int64_t s2 = crs_find(A, nd[l], nd[m]);
				if (s2 >= 0) A->val[s2] += v * area;
			}
		}
	}
}


// 質量行列 M_ij = ∫ σ λ_i λ_j dS = σ S (1 + δ_ij)/12 (導体要素のみ)
void assemble_mass_tri(crs_t *A)
{
	assemble_mass_tri_sig(A, NULL);
}


// 要素ごとの σ (sig[NTri]、NULL なら導体ごとの CondSigma) で組む。熱連成 (thermal + F) で
// σ(T) を要素ごとに渡すのに使う。NULL のときの演算は従来と同じ
static void mass_tri_add(crs_t *A, const double *sig, int all);

void assemble_mass_tri_sig(crs_t *A, const double *sig)
{
	mass_tri_add(A, sig, 0);
}


// 静磁場の過渡解析 (M + transient) の渦電流の質量 : 要素ごとの σ (sig[NTri]) を導体でない要素にも掛ける
// (構造格子の assemble_mass_cell と同じ扱い。導体の要素は呼び出し側で 0 にする)
void assemble_mass_tri_elem(crs_t *A, const double *sig)
{
	mass_tri_add(A, sig, 1);
}


static void mass_tri_add(crs_t *A, const double *sig, int all)
{
	crs_zero(A);

	for (int e = 0; e < NTri; e++) {
		const int id = TriCond[e];
		if ((id < 0) && !all) continue;
		const double sg = ((sig || all) ? sig[e] : CondSigma[id]);
		if (sg <= 0) continue;

		int32_t nd[TRI_NMAX];
		const int nen = tri_nodes(e, nd);
		if (nen >= 6) {
			// 2 次 / 3 次三角形 : ∫N_i N_j dA を Gauss 積分する
			double lam[NQTRI3][3], wq[NQTRI3];
			const int nq = tri_quad_ord(lam, wq);
			for (int q = 0; q < nq; q++) {
				double gn[TRI_NMAX][2], det, n[TRI_NMAX], dl[TRI_NMAX][3];
				if (tri6_grad(nd, nen, lam[q], gn, &det)) break;
				trin_shape(nen, lam[q], n, dl);
				const double dw = sg * wq[q] * fabs(det);
				for (int l = 0; l < nen; l++) {
					for (int m = 0; m < nen; m++) {
						const int64_t s = crs_find(A, nd[l], nd[m]);
						if (s >= 0) A->val[s] += dw * n[l] * n[m];
					}
				}
			}
			continue;
		}
		const double w = sg * TriArea[e] / 12;
		for (int l = 0; l < 3; l++) {
			for (int m = 0; m < 3; m++) {
				const int64_t s = crs_find(A, nd[l], nd[m]);
				if (s >= 0) A->val[s] += w * ((l == m) ? 2 : 1);
			}
		}
	}
}
