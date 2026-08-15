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

#define ARRAY_INC (100000)


// ---- Gmsh ASCII 2.2 の読み込み ----

static int read_nodes(FILE *fp, int32_t **idmap, int32_t *maxid)
{
	char line[BUFSIZ];

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

static void elem_reset(void)
{
	NTet = 0;
	NTri = 0;
	Tet = NULL;
	TetTag = NULL;
	Tri = NULL;
	TriTag = NULL;
	Tet2 = NULL;
	Tri2 = NULL;
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


// 要素 1 個を格納する。type は Gmsh の要素型
// (4/11 = 四面体、2/9 = 三角形、5 = 六面体、6 = 角柱、7 = ピラミッド、3 = 四角形)
static int elem_store(int type, int tag, const int32_t *nd)
{
	const int order = ((type == 11) || (type == 9)) ? 2 : 1;

	if ((type == 5) || (type == 12) || (type == 17)) {
		// 六面体 (8 / 27 / 20 節点)。局所の並びは Gmsh のまま使う。
		// 頂点 8 個を Hex に、2 次の残り (辺 12 + 面 6 + 体心 1) を Hex2 に
		// 入れる。次数は最初の 1 個で決め、以後は混在を許さない
		const int hn = ((type == 5) ? 8 : (type == 12) ? 27 : 20);
		if (HexNen == 0) HexNen = hn;
		if (HexNen != hn) {
			printf("*** mesh : mixed hexahedron orders (%d-node after %d-node)\n",
				hn, HexNen);
			return 1;
		}
		if (NHex % ARRAY_INC == 0) {
			Hex = (int32_t *)realloc(Hex, (size_t)(NHex + ARRAY_INC) * 8 * sizeof(int32_t));
			HexTag = (int *)realloc(HexTag, (size_t)(NHex + ARRAY_INC) * sizeof(int));
			if (hn > 8) {
				Hex2 = (int32_t *)realloc(Hex2, (size_t)(NHex + ARRAY_INC) * 19 * sizeof(int32_t));
			}
		}
		for (int l = 0; l < 8; l++) Hex[(NHex * 8) + l] = nd[l];
		for (int l = 8; l < hn; l++) Hex2[(NHex * 19) + (l - 8)] = nd[l];
		HexTag[NHex] = tag;
		NHex++;

		return 0;
	}
	if ((type == 6) || (type == 13) || (type == 18)) {
		// 角柱 (6 / 18 / 15 節点)。局所の並びは Gmsh のまま使う。
		// 頂点 6 個を Prism に、2 次の残り (辺 9 + 四角形面 3) を Prism2 に入れる
		const int pn = ((type == 6) ? 6 : (type == 13) ? 18 : 15);
		if (PrismNen == 0) PrismNen = pn;
		if (PrismNen != pn) {
			printf("*** mesh : mixed prism orders (%d-node after %d-node)\n",
				pn, PrismNen);
			return 1;
		}
		if (NPrism % ARRAY_INC == 0) {
			Prism = (int32_t *)realloc(Prism, (size_t)(NPrism + ARRAY_INC) * 6 * sizeof(int32_t));
			PrismTag = (int *)realloc(PrismTag, (size_t)(NPrism + ARRAY_INC) * sizeof(int));
			if (pn > 6) {
				Prism2 = (int32_t *)realloc(Prism2, (size_t)(NPrism + ARRAY_INC) * 12 * sizeof(int32_t));
			}
		}
		for (int l = 0; l < 6; l++) Prism[(NPrism * 6) + l] = nd[l];
		for (int l = 6; l < pn; l++) Prism2[(NPrism * 12) + (l - 6)] = nd[l];
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
	if ((type == 3) || (type == 16) || (type == 10)) {
		// 四角形 (境界面、4 / 8 / 9 節点)。六面体格子の電極面はこれになる。
		// 2 次の中間節点 (辺 4 + 中心) は Quad2 に入れる。無い位置は -1 で
		// 埋める (realloc した領域は未初期化なので、必ず全 5 要素を書くこと)
		const int qn = ((type == 3) ? 4 : (type == 16) ? 8 : 9);
		if (NQuad % ARRAY_INC == 0) {
			Quad = (int32_t *)realloc(Quad, (size_t)(NQuad + ARRAY_INC) * 4 * sizeof(int32_t));
			QuadTag = (int *)realloc(QuadTag, (size_t)(NQuad + ARRAY_INC) * sizeof(int));
			Quad2 = (int32_t *)realloc(Quad2, (size_t)(NQuad + ARRAY_INC) * 5 * sizeof(int32_t));
		}
		for (int l = 0; l < 4; l++) Quad[(NQuad * 4) + l] = nd[l];
		for (int l = 0; l < 5; l++) {
			Quad2[(NQuad * 5) + l] = (((4 + l) < qn) ? nd[4 + l] : -1);
		}
		QuadTag[NQuad] = tag;
		NQuad++;

		return 0;
	}

	if ((type == 4) || (type == 11)) {
		// 四面体。次数は最初の 1 個で決め、以後は混在を許さない
		if (TetOrder == 0) TetOrder = order;
		if (TetOrder != order) {
			printf("*** mesh : mixed element orders (order %d after order %d)\n",
				order, TetOrder);
			return 1;
		}
		if (NTet % ARRAY_INC == 0) {
			Tet = (int32_t *)realloc(Tet, (size_t)(NTet + ARRAY_INC) * 4 * sizeof(int32_t));
			TetTag = (int *)realloc(TetTag, (size_t)(NTet + ARRAY_INC) * sizeof(int));
			if (order == 2) {
				Tet2 = (int32_t *)realloc(Tet2, (size_t)(NTet + ARRAY_INC) * 6 * sizeof(int32_t));
			}
		}
		for (int l = 0; l < 4; l++) Tet[(NTet * 4) + l] = nd[l];
		if (order == 2) {
			for (int l = 0; l < 6; l++) Tet2[(NTet * 6) + l] = nd[4 + l];
		}
		TetTag[NTet] = tag;
		NTet++;
	}
	else {
		if (NTri % ARRAY_INC == 0) {
			Tri = (int32_t *)realloc(Tri, (size_t)(NTri + ARRAY_INC) * 3 * sizeof(int32_t));
			TriTag = (int *)realloc(TriTag, (size_t)(NTri + ARRAY_INC) * sizeof(int));
			Tri2 = (int32_t *)realloc(Tri2, (size_t)(NTri + ARRAY_INC) * 3 * sizeof(int32_t));
		}
		for (int l = 0; l < 3; l++) Tri[(NTri * 3) + l] = nd[l];
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


// 全要素を読んだあとの整合性検査 (格子の次元と次数)
static int elem_finish(void)
{
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
			*/
			if ((NTet > 0) && (NHex > 0) && (NPyr < 1)) {
				printf("%s\n", "*** mesh : hexahedra and tetrahedra cannot share a "
					"face conformingly (the quadrilateral face carries a bilinear "
					"trace, the triangular face a linear one); put prisms or "
					"pyramids between them, or use one element type");
				return 1;
			}
			if ((TetOrder >= 2) || (HexNen > 8) || (PrismNen > 6)) {
				printf("%s\n", "*** mesh : mixed element types must be first order "
					"(10-node tetrahedra and second-order hexahedra / prisms "
					"cannot be mixed with other kinds)");
				return 1;
			}
			MeshElem = MESHELEM_MIXED;
			MeshDim = 3;
			TetOrder = 1;
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
		TetOrder = ((PrismNen > 6) ? 2 : 1);
		if ((NTri < 1) && (NQuad < 1)) {
			printf("%s\n", "*** mesh : a prism mesh needs triangular or "
				"quadrilateral boundary faces for the electrodes");
			return 1;
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
				if (Quad2[(t * 5)] < 0) {
					printf("%s\n", "*** mesh : the prisms are order 2 but a "
						"quadrilateral is order 1 (regenerate the mesh with -order 2)");
					return 1;
				}
				// 18 節点の角柱の側面は面心つきの 9 節点四角形 (六面体と同じ穴)
				if ((PrismNen == 18) && (Quad2[(t * 5) + 4] < 0)) {
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
		TetOrder = ((HexNen > 8) ? 2 : 1);
		if (NQuad < 1) {
			printf("%s\n", "*** mesh : a hexahedral mesh needs quadrilateral "
				"boundary faces for the electrodes");
			return 1;
		}
		// 2 次の六面体に 1 次の四角形が混ざると電極の中間節点が固定されない。
		// **27 節点の六面体は面心節点まで要る** (8 節点四角形では辺の中間節点が
		// 埋まるので、先頭だけ見る検査を素通りして C が -41% になる。実測)
		if (TetOrder >= 2) {
			for (int t = 0; t < NQuad; t++) {
				if (Quad2[(t * 5)] < 0) {
					printf("%s\n", "*** mesh : the hexahedra are order 2 but a "
						"quadrilateral is order 1 (regenerate the mesh with -order 2)");
					return 1;
				}
				if ((HexNen == 27) && (Quad2[(t * 5) + 4] < 0)) {
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
		// 三角形の次数は「中間節点が頂点と違うか」で決まる (1 次では頂点で埋めてある)
		TetOrder = ((Tri2[0] != Tri[0]) ? 2 : 1);
		for (int t = 0; t < NTri; t++) {
			const int o = ((Tri2[(t * 3)] != Tri[(t * 3)]) ? 2 : 1);
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

	return 0;
}


static int read_elements(FILE *fp, const int32_t *idmap, int32_t maxid)
{
	char line[BUFSIZ];

	if (fgets(line, sizeof(line), fp) == NULL) return 1;
	const int ne = atoi(line);

	elem_reset();

	for (int e = 0; e < ne; e++) {
		if (fgets(line, sizeof(line), fp) == NULL) return 1;

		// id type ntags tag1 ... tagn node1 ...
		int nv = 0;
		long v[64];
		char *p = line;
		while ((nv < 64) && (*p != '\0')) {
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
		const int nn = ((type == 4) ? 4 : (type == 11) ? 10
		              : (type == 2) ? 3 : (type == 9) ? 6
		              : (type == 5) ? 8 : (type == 3) ? 4 : (type == 6) ? 6
		              : (type == 7) ? 5
		              : (type == 12) ? 27 : (type == 17) ? 20
		              : (type == 13) ? 18 : (type == 18) ? 15
		              : (type == 10) ? 9 : (type == 16) ? 8 : 0);
		if (nn == 0) continue;			// 点・線分など、使わない要素型
		if (nv < off + nn) continue;

		// 節点番号を先に解決する (未解決なら打ち切り。ここで continue すると
		// 配列の該当要素が未初期化のまま確定し、setup_unstruct() が
		// それを添字に使って領域外書き込みになる)
		int32_t nd[27];
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

typedef struct {
	int dim;
	long tag;
	int phys;			// 物理タグ (無ければ 0)
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
			int phys = 0;
			for (long i = 0; i < np; i++) {
				long pt = 0;
				if (fscanf(fp, "%ld", &pt) != 1) return 1;
				if (i == 0) phys = (int)pt;		// 複数あれば最初のものを使う
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
			if (NEnt < MAXENT) {
				Ent[NEnt].dim = d;
				Ent[NEnt].tag = tag;
				Ent[NEnt].phys = phys;
				NEnt++;
			}
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

		const int nn = ((type == 4) ? 4 : (type == 11) ? 10
		              : (type == 2) ? 3 : (type == 9) ? 6
		              : (type == 5) ? 8 : (type == 3) ? 4 : (type == 6) ? 6
		              : (type == 7) ? 5
		              : (type == 12) ? 27 : (type == 17) ? 20
		              : (type == 13) ? 18 : (type == 18) ? 15
		              : (type == 10) ? 9 : (type == 16) ? 8 : 0);
		// 物理タグはエンティティ側にある (2.2 と違い要素の行には無い)
		const int phys = ent_phys((int)dim, tag);

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
			int32_t nd[27];
			for (int l = 0; l < nn; l++) {
				long g = 0;
				if (fscanf(fp, "%ld", &g) != 1) return 1;
				if ((g < 0) || (g > maxid) || (idmap[g] < 0)) {
					printf("*** mesh : element %ld refers to an unknown node %ld\n", etag, g);
					return 1;
				}
				nd[l] = idmap[g];
			}
			if (elem_store((int)type, phys, nd)) return 1;
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
	char line[BUFSIZ];

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
	char line[BUFSIZ];

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
		              || (type == 5) || (type == 3) || (type == 6) || (type == 7)
		              || (type == 12) || (type == 17) || (type == 13) || (type == 18)
		              || (type == 10) || (type == 16));

		for (long e = 0; e < cnt; e++) {
			int32_t etag = 0;
			if (rd_i32(fp, &etag)) return 1;
			int tag = 0;
			for (long t = 0; t < ntag; t++) {
				int32_t v = 0;
				if (rd_i32(fp, &v)) return 1;
				if (t == 0) tag = (int)v;		// 1 つ目が物理タグ
			}
			int32_t nd[27];
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
			int phys = 0;
			for (int64_t i = 0; i < np; i++) {
				int32_t pt = 0;
				if (rd_i32(fp, &pt)) return 1;
				if (i == 0) phys = (int)pt;		// 複数あれば最初のものを使う
			}
			if (d > 0) {
				int64_t nbd = 0;
				if (rd_i64(fp, &nbd)) return 1;
				for (int64_t i = 0; i < nbd; i++) {
					int32_t dummy = 0;
					if (rd_i32(fp, &dummy)) return 1;
				}
			}
			if (NEnt < MAXENT) {
				Ent[NEnt].dim = d;
				Ent[NEnt].tag = tag;
				Ent[NEnt].phys = phys;
				NEnt++;
			}
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
		const int phys = ent_phys((int)dim, tag);
		const int use = ((type == 4) || (type == 11) || (type == 2) || (type == 9)
		              || (type == 5) || (type == 3) || (type == 6) || (type == 7)
		              || (type == 12) || (type == 17) || (type == 13) || (type == 18)
		              || (type == 10) || (type == 16));

		for (int64_t e = 0; e < cnt; e++) {
			int64_t etag = 0;
			if (rd_i64(fp, &etag)) return 1;
			int32_t nd[27];
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
			if (use && elem_store((int)type, phys, nd)) return 1;
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
static int mesh_read_binary(const char *fname, int v41)
{
	char line[BUFSIZ];
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
			ierr = (v41 ? read_entities_bin41(fp) : 0);
			if (ierr) {
				printf("%s\n", "*** mesh : cannot parse the binary $Entities");
				break;
			}
		}
		else if (!strncmp(line, "$Nodes", 6)) {
			ierr = (v41 ? read_nodes_bin41(fp, &idmap, &maxid)
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
			ierr = (v41 ? read_elements_bin41(fp, idmap, maxid)
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
	char line[BUFSIZ];
	int32_t *idmap = NULL;
	int32_t maxid = 0;
	int ierr = 0;
	int have_nodes = 0, have_elements = 0;
	int v41 = 0;

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
			else if (line[0] == '4') v41 = 1;
			else {
				printf("*** mesh : only Gmsh ASCII format 2.x / 4.x is supported (got %s)", line);
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
						return mesh_read_binary(fname, v41);
					}
				}
			}
		}
		else if (!strncmp(line, "$Entities", 9)) {
			ierr = read_entities_v41(fp);
			if (ierr) {
				printf("%s\n", "*** mesh : cannot parse $Entities");
				break;
			}
		}
		else if (!strncmp(line, "$Nodes", 6)) {
			ierr = (v41 ? read_nodes_v41(fp, &idmap, &maxid)
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
			ierr = (v41 ? read_elements_v41(fp, idmap, maxid)
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


// 要素 e の局所節点をまとめて取り出す。戻り値は節点数 (1 次 4、2 次 10)。
// 並びは Gmsh の tet10 と同じ (頂点 4 個のあと辺上の中間節点 6 個)
int tet_nodes(int e, int32_t nd[10])
{
	for (int l = 0; l < 4; l++) nd[l] = Tet[(e * 4) + l];
	if (TetOrder < 2) return 4;
	for (int l = 0; l < 6; l++) nd[4 + l] = Tet2[(e * 6) + l];

	return 10;
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
		int32_t nd[27];
		const int nen = get(e, nd);
		for (int l = 0; l < nen; l++) cnt[nd[l]]++;
	}
	int64_t *nptr = (int64_t *)malloc(((size_t)n + 1) * sizeof(int64_t));
	nptr[0] = 0;
	for (int i = 0; i < n; i++) nptr[i + 1] = nptr[i] + cnt[i];
	int32_t *nlist = (int32_t *)malloc((size_t)nptr[n] * sizeof(int32_t));
	memset(cnt, 0, (size_t)n * sizeof(int));
	for (int e = 0; e < nelem; e++) {
		int32_t nd[27];
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
			int32_t nd[27];
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
			int32_t nd[27];
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
int hex_nodes(int e, int32_t nd[27])
{
	for (int l = 0; l < 8; l++) nd[l] = Hex[(e * 8) + l];
	for (int l = 8; l < HexNen; l++) nd[l] = Hex2[(e * 19) + (l - 8)];

	return HexNen;
}


// 角柱 e の局所節点をまとめて取り出す。戻り値は節点数 (6 / 15 / 18)
int prism_nodes(int e, int32_t nd[18])
{
	for (int l = 0; l < 6; l++) nd[l] = Prism[(e * 6) + l];
	for (int l = 6; l < PrismNen; l++) nd[l] = Prism2[(e * 12) + (l - 6)];

	return PrismNen;
}


// 要素 e の節点。戻り値は節点数
// (四面体 4 or 10、六面体 8/20/27、角柱 6/15/18、ピラミッド 5)
int elem3d_nodes(int e, int32_t nd[27])
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
	crs_alloc_conn(A, NTet, ((TetOrder >= 2) ? 10 : 4), get_tet);
}


// 3 次元の非構造格子 (種別の混在を含む)
void crs_alloc_elem3d(crs_t *A)
{
	// work バッファの見積りに使う要素あたりの最大節点数 (大きい分には無害)
	int nenmax = ((TetOrder >= 2) ? 10 : 8);
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
	int32_t nd[10];
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


// 要素の重心での ∇N_i (場の出力に使う。1 次では要素内一定なので同じ値)
int tet_grad_center(int e, double gn[10][3], int *nen)
{
	int32_t nd[10];
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

	return tet10_grad(nd, lam, gn, &det);
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


// 局所座標での微分 dn[a][k] = ∂N_a/∂(ξ,η,ζ)
static void hex2_dshape(int nen, double x, double y, double z, double dn[27][3])
{
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
	double g[27][3], double *det)
{
	double dn[27][3];
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


// 3 点 Gauss-Legendre (各方向) の位置と重み
static const double GX3[3] = {-0.7745966692414834, 0.0, 0.7745966692414834};
static const double GW3[3] = {5.0 / 9, 8.0 / 9, 5.0 / 9};


// 2 次六面体の要素行列 (3x3x3 Gauss)。vol には体積 Σ w |det J| を返す
static int hex2_element(int e, const double c[6], double ke[27][27], double *vol)
{
	int32_t nd[27];
	const int nen = hex_nodes(e, nd);

	memset(ke, 0, sizeof(double) * 27 * 27);
	double v = 0, sgn = 0;
	for (int i = 0; i < 3; i++) {
	for (int j = 0; j < 3; j++) {
	for (int k = 0; k < 3; k++) {
		double g[27][3], det;
		if (hex2_grad(nd, nen, GX3[i], GX3[j], GX3[k], g, &det)) return 1;
		if (sgn == 0) sgn = ((det > 0) ? 1 : -1);
		else if (det * sgn <= 0) return 1;		// 要素内で符号が変わる = 裏返り
		const double dw = GW3[i] * GW3[j] * GW3[k] * fabs(det);
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
int hex2_grad_center(int e, double g[27][3], int *nen)
{
	int32_t nd[27];
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
static void prism2_dshape(int nen, double u, double v, double ze, double dn[18][3])
{
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


static int prism2_grad(const int32_t *nd, int nen, double u, double v, double ze,
	double g[18][3], double *det)
{
	double dn[18][3];
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
static int prism2_element(int e, const double c[6], double ke[18][18], double *vol)
{
	int32_t nd[18];
	const int nen = prism_nodes(e, nd);

	double lam[NQTRI][3], wq[NQTRI];
	tri_quad(lam, wq);

	memset(ke, 0, sizeof(double) * 18 * 18);
	double v = 0, sgn = 0;
	for (int q = 0; q < NQTRI; q++) {
	for (int k = 0; k < 3; k++) {
		double g[18][3], det;
		if (prism2_grad(nd, nen, lam[q][1], lam[q][2], GX3[k], g, &det)) return 1;
		if (sgn == 0) sgn = ((det > 0) ? 1 : -1);
		else if (det * sgn <= 0) return 1;
		const double dw = wq[q] * GW3[k] * fabs(det);
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
int prism2_grad_center(int e, double g[18][3], int *nen)
{
	int32_t nd[18];
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
	if (err > 1e-12) {
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
	if (err > 1e-12) {
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
	if (err > 1e-12) {
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
static const double GX4[4] = {-0.8611363115940526, -0.3399810435848563,
                               0.3399810435848563,  0.8611363115940526};
static const double GW4[4] = {0.3478548451374538, 0.6521451548625461,
                              0.6521451548625461, 0.3478548451374538};


// 中間節点が「1 次の写像の位置」からどれだけずれているか (最大、辺長比)
static double hex2_curve(void)
{
	double curve = 0;

	for (int e = 0; e < NHex; e++) {
		int32_t nd[27];
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
		int32_t nd[27];
		const int nen = hex_nodes(e, nd);
		for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
		for (int k = 0; k < 3; k++) {
			double g[27][3], det;
			if (hex2_grad(nd, nen, GX3[i], GX3[j], GX3[k], g, &det)) continue;
			v3 += GW3[i] * GW3[j] * GW3[k] * fabs(det);
		}
		}
		}
		for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
		for (int k = 0; k < 4; k++) {
			double g[27][3], det;
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
		int32_t nd[27];
		const int nen = hex_nodes(e, nd);
		double ke[27][27], vol;
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
		if (err > 1e-12) {
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
			if (err > 1e-10) {
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
		int32_t nd[18];
		const int nen = prism_nodes(e, nd);
		double lam[NQTRI][3], wq[NQTRI];
		tri_quad(lam, wq);
		for (int q = 0; q < NQTRI; q++) {
			for (int k = 0; k < 3; k++) {
				double g[18][3], det;
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
				double g[18][3], det;
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
		int32_t nd[18];
		const int nen = prism_nodes(e, nd);
		double ke[18][18], vol;
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
		if (err > 1e-12) {
			fprintf(fp_log, "*** the linear-field identity failed\n");
			ierr = 1;
		}
	}

	// (3) 2 次場の恒等式 (鉛直に押し出したまっすぐな角柱の要素に限る)
	{
		int nright = 0;
		double exact = 0;
		for (int e = 0; e < NPrism; e++) {
			int32_t nd[18];
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
			if (err > 1e-10) {
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
			if (p2) {
				int32_t nd[10];
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
				int32_t nd[27];
				const int nen = hex_nodes(e - NTet, nd);
				double ke[27][27], vol;
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
				int32_t nd[18];
				const int nen = prism_nodes(e - NTet - NHex, nd);
				double ke[18][18], vol;
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

		if (p2) {
			int32_t nd[10];
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
	if (err > 1e-12) {
		fprintf(fp_log, "*** the linear-field identity failed\n");
		ierr = 1;
	}

	free(phi);
	crs_free(&A);

	return ierr;
}


int solve_nodal_test(FILE *fp_log)
{
	int ierr = 0;

	if (MeshElem == MESHELEM_HEX) {
		return ((HexNen > 8) ? nodal_test_hex2(fp_log) : nodal_test_hex(fp_log));
	}
	if (MeshElem == MESHELEM_PRISM) {
		return ((PrismNen > 6) ? nodal_test_prism2(fp_log) : nodal_test_prism(fp_log));
	}
	if (MeshElem == MESHELEM_PYR)   return nodal_test_pyr(fp_log);
	if (MeshElem == MESHELEM_MIXED) return nodal_test_mixed(fp_log);

	fprintf(fp_log, "\n=== nodal element (P%d) self test ===\n", TetOrder);
	fprintf(fp_log, "  nodes = %d, tetrahedra = %d, nodes per element = %d\n",
		NNode, NTet, ((TetOrder >= 2) ? 10 : 4));

	// 中間節点が辺の中点からどれだけずれているか (曲がった格子の判定)
	double curve = 0;
	if (TetOrder >= 2) {
		static const int ed[6][2] = {{0, 1}, {1, 2}, {2, 0}, {3, 0}, {3, 2}, {3, 1}};
		for (int e = 0; e < NTet; e++) {
			int32_t nd[10];
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
			if (TetOrder >= 2) {
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
		if (rel > 1e-10) {
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

	for (int cs = 0; cs < 2; cs++) {
		const int quad = cs;			// 0 : φ は 1 次、1 : φ は 2 次
		if (quad && ((TetOrder < 2) || !straight)) continue;

		// 節点値
		for (int i = 0; i < n; i++) {
			const double r[3] = {Xp[i], Yp[i], Zp[i]};
			double f = (av[0] * r[0]) + (av[1] * r[1]) + (av[2] * r[2]);
			if (quad) {
				for (int k = 0; k < 3; k++) {
					for (int l = 0; l < 3; l++) f += bm[k][l] * r[k] * r[l] / 2;
				}
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

		crs_spmv(&A, phi, y, NULL);
		double q = 0;
		for (int i = 0; i < n; i++) q += phi[i] * y[i];
		const double rel = fabs(q - exact) / ((exact != 0) ? fabs(exact) : 1);
		fprintf(fp_log, "  (%c) %s field : phi^T K phi = %.10e, exact = %.10e, err = %.3e\n",
			(quad ? 'c' : 'b'), (quad ? "quadratic" : "linear   "), q, exact, rel);
		if (rel > 1e-10) {
			fprintf(fp_log, "*** the %s stiffness matrix is wrong\n",
				((TetOrder >= 2) ? "P2" : "P1"));
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
				double gn[10][3];
				int nen = 0;
				if (tet_grad_center(e, gn, &nen)) continue;
				int32_t nd[10];
				tet_nodes(e, nd);

				double rc[3] = {0, 0, 0};
				for (int a = 0; a < 4; a++) {
					rc[0] += Xp[nd[a]] / 4;
					rc[1] += Yp[nd[a]] / 4;
					rc[2] += Zp[nd[a]] / 4;
				}
				for (int k = 0; k < 3; k++) {
					double gh = 0;
					for (int a = 0; a < nen; a++) gh += phi[nd[a]] * gn[a][k];
					double ge = av[k];
					if (quad) {
						for (int l = 0; l < 3; l++) ge += bm[k][l] * rc[l];
					}
					if (fabs(gh - ge) > gmax) gmax = fabs(gh - ge);
					if (fabs(ge) > amax2) amax2 = fabs(ge);
				}
			}
			const double relg = gmax / ((amax2 > 0) ? amax2 : 1);
			fprintf(fp_log, "  (e) centroid gradient (%s field) : max err = %.3e\n",
				(quad ? "quadratic" : "linear"), relg);
			if (relg > 1e-10) {
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
		if (rel > 1e-12) {
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

// 三角形の局所節点をまとめて取り出す。戻り値は節点数 (1 次 3、2 次 6)。
// 並びは Gmsh の tri6 と同じ (頂点 3 個のあと辺 (0,1) (1,2) (2,0) の中間節点)
int tri_nodes(int e, int32_t nd[6])
{
	for (int l = 0; l < 3; l++) nd[l] = Tri[(e * 3) + l];
	if (TetOrder < 2) return 3;
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

	const int nen = ((TetOrder >= 2) ? 6 : 3);
	int *cnt = (int *)malloc((size_t)n * sizeof(int));
	memset(cnt, 0, (size_t)n * sizeof(int));
	for (int e = 0; e < NTri; e++) {
		int32_t nd[6];
		tri_nodes(e, nd);
		for (int l = 0; l < nen; l++) cnt[nd[l]]++;
	}
	int64_t *nptr = (int64_t *)malloc(((size_t)n + 1) * sizeof(int64_t));
	nptr[0] = 0;
	for (int i = 0; i < n; i++) nptr[i + 1] = nptr[i] + cnt[i];
	int32_t *nlist = (int32_t *)malloc((size_t)nptr[n] * sizeof(int32_t));
	memset(cnt, 0, (size_t)n * sizeof(int));
	for (int e = 0; e < NTri; e++) {
		int32_t nd[6];
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
				int32_t nd[6];
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


// 積分点 1 点での面内勾配 ∇N_i とヤコビアン行列式。戻り値 : 0 = 正常
static int tri6_grad(const int32_t nd[6], const double lam[3],
	double gn[6][2], double *det)
{
	int p, q;
	tri_axes(&p, &q);
	const double *cp = ((p == 0) ? Xp : (p == 1) ? Yp : Zp);
	const double *cq = ((q == 0) ? Xp : (q == 1) ? Yp : Zp);

	double n[6], dl[6][3], dx[6][2];
	tri6_shape(lam, n, dl);
	// 参照座標 ξ = (λ1, λ2) についての微分 (λ0 = 1 - λ1 - λ2)
	for (int i = 0; i < 6; i++) {
		for (int k = 0; k < 2; k++) dx[i][k] = dl[i][k + 1] - dl[i][0];
	}

	double j[2][2] = {{0, 0}, {0, 0}};
	for (int k = 0; k < 2; k++) {
		for (int i = 0; i < 6; i++) {
			j[0][k] += cp[nd[i]] * dx[i][k];
			j[1][k] += cq[nd[i]] * dx[i][k];
		}
	}
	const double d = (j[0][0] * j[1][1]) - (j[0][1] * j[1][0]);
	if (d == 0) return 1;

	// Jinv[k][r] = ∂ξ_k/∂x_r
	const double ji[2][2] = {{ j[1][1] / d, -j[0][1] / d},
	                         {-j[1][0] / d,  j[0][0] / d}};
	for (int i = 0; i < 6; i++) {
		for (int r = 0; r < 2; r++) {
			gn[i][r] = (dx[i][0] * ji[0][r]) + (dx[i][1] * ji[1][r]);
		}
	}
	*det = d;

	return 0;
}


// 2 次三角形の面積 (等パラメトリック写像で積分する)
double tri6_area(int e)
{
	int32_t nd[6];
	if (tri_nodes(e, nd) != 6) return 0;

	double lam[NQTRI][3], wq[NQTRI], a = 0;
	tri_quad(lam, wq);
	for (int q = 0; q < NQTRI; q++) {
		double gn[6][2], det;
		if (tri6_grad(nd, lam[q], gn, &det)) return 0;
		a += wq[q] * fabs(det);
	}

	return a;
}


// 要素の重心での面内 ∇N_i (場の出力に使う)。1 次では要素内一定なので同じ値
int tri_grad_center(int e, double gn[6][2], int *nen)
{
	int32_t nd[6];
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

	return tri6_grad(nd, lam, gn, &det);
}


// 剛性行列 K_ij = ∫ ν ∇λ_i・∇λ_j dS  (面内の異方性テンソルも扱う)
// nucell が非 NULL なら要素毎の ν を使う (等方性、非線形解析用)。
void assemble_nu_tri(crs_t *A, const double *nucell)
{
	int p, q;
	tri_axes(&p, &q);
	crs_zero(A);

	const int p2 = (TetOrder >= 2);
	double lam[NQTRI][3], wq[NQTRI];
	if (p2) tri_quad(lam, wq);

	for (int e = 0; e < NTri; e++) {
		int32_t nd6[6];
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
			for (int q = 0; q < NQTRI; q++) {
				double gn[6][2], det;
				if (tri6_grad(nd, lam[q], gn, &det)) break;
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


// 質量行列 M_ij = ∫ σ λ_i λ_j dS = σ S (1 + δ_ij)/12 (導体要素のみ)
void assemble_mass_tri(crs_t *A)
{
	crs_zero(A);

	for (int e = 0; e < NTri; e++) {
		const int id = TriCond[e];
		if (id < 0) continue;
		const double sg = CondSigma[id];
		if (sg <= 0) continue;

		int32_t nd[6];
		const int nen = tri_nodes(e, nd);
		if (nen == 6) {
			// 2 次三角形 : ∫N_i N_j dA を Gauss 積分する
			double lam[NQTRI][3], wq[NQTRI];
			tri_quad(lam, wq);
			for (int q = 0; q < NQTRI; q++) {
				double gn[6][2], det, n[6], dl[6][3];
				if (tri6_grad(nd, lam[q], gn, &det)) break;
				tri6_shape(lam[q], n, dl);
				const double dw = sg * wq[q] * fabs(det);
				for (int l = 0; l < 6; l++) {
					for (int m = 0; m < 6; m++) {
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
