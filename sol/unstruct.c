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


static int read_elements(FILE *fp, const int32_t *idmap, int32_t maxid)
{
	char line[BUFSIZ];

	if (fgets(line, sizeof(line), fp) == NULL) return 1;
	const int ne = atoi(line);

	NTet = 0;
	NTri = 0;
	Tet = NULL;
	TetTag = NULL;
	Tri = NULL;
	TriTag = NULL;
	Tet2 = NULL;
	Tri2 = NULL;
	TetOrder = 0;			// 最初に出た四面体で決まる

	for (int e = 0; e < ne; e++) {
		if (fgets(line, sizeof(line), fp) == NULL) return 1;

		// id type ntags tag1 ... tagn node1 ...
		int nv = 0;
		long v[32];
		char *p = line;
		while ((nv < 32) && (*p != '\0')) {
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
		const int order = ((type == 11) || (type == 9)) ? 2 : 1;
		const int nn = ((type == 4) ? 4 : (type == 11) ? 10
		              : (type == 2) ? 3 : (type == 9) ? 6 : 0);
		if (nn == 0) continue;			// 点・線分など、使わない要素型
		if (nv < off + nn) continue;

		// 節点番号を先に解決する (未解決なら打ち切り。ここで continue すると
		// 配列の該当要素が未初期化のまま確定し、setup_unstruct() が
		// それを添字に使って領域外書き込みになる)
		int32_t nd[10];
		for (int l = 0; l < nn; l++) {
			const long g = v[off + l];
			if ((g < 0) || (g > maxid) || (idmap[g] < 0)) {
				printf("*** mesh : element %d refers to an unknown node %ld\n", e + 1, g);
				return 1;
			}
			nd[l] = idmap[g];
		}

		if ((type == 4) || (type == 11)) {
			// 四面体。次数は最初の 1 個で決め、以後は混在を許さない
			if (TetOrder == 0) TetOrder = order;
			if (TetOrder != order) {
				printf("*** mesh : mixed element orders (element %d is order %d, "
					"the mesh started as order %d)\n", e + 1, order, TetOrder);
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
			// 三角形 (電極面の指定に使う)。次数は四面体に合わせる
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
	if ((TetOrder == 2) && (NTri > 0)) {
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


int mesh_read(const char *fname)
{
	char line[BUFSIZ];
	int32_t *idmap = NULL;
	int32_t maxid = 0;
	int ierr = 0;
	int have_nodes = 0, have_elements = 0;

	FILE *fp = fopen(fname, "r");
	if (fp == NULL) {
		printf("*** mesh file %s open error.\n", fname);
		return 1;
	}

	while (fgets(line, sizeof(line), fp) != NULL) {
		if      (!strncmp(line, "$MeshFormat", 11)) {
			if (fgets(line, sizeof(line), fp) == NULL) break;
			if (line[0] != '2') {
				printf("*** mesh : only Gmsh ASCII format 2.x is supported (got %s)", line);
				ierr = 1;
				break;
			}
		}
		else if (!strncmp(line, "$Nodes", 6)) {
			ierr = read_nodes(fp, &idmap, &maxid);
			if (ierr) break;
			have_nodes = 1;
		}
		else if (!strncmp(line, "$Elements", 9)) {
			if (!have_nodes) {
				printf("%s\n", "*** mesh : $Elements before $Nodes");
				ierr = 1;
				break;
			}
			ierr = read_elements(fp, idmap, maxid);
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


// 四面体の連結から節点の隣接関係を作る (対角成分を含む)
void crs_alloc_tet(crs_t *A)
{
	const int n = NNode;
	const int nen = ((TetOrder >= 2) ? 10 : 4);		// 要素あたりの節点数

	// 節点毎の要素数を数える
	int *cnt = (int *)malloc((size_t)n * sizeof(int));
	memset(cnt, 0, (size_t)n * sizeof(int));
	for (int e = 0; e < NTet; e++) {
		int32_t nd[10];
		tet_nodes(e, nd);
		for (int l = 0; l < nen; l++) cnt[nd[l]]++;
	}
	int64_t *nptr = (int64_t *)malloc(((size_t)n + 1) * sizeof(int64_t));
	nptr[0] = 0;
	for (int i = 0; i < n; i++) nptr[i + 1] = nptr[i] + cnt[i];
	int32_t *nlist = (int32_t *)malloc((size_t)nptr[n] * sizeof(int32_t));
	memset(cnt, 0, (size_t)n * sizeof(int));
	for (int e = 0; e < NTet; e++) {
		int32_t nd[10];
		tet_nodes(e, nd);
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
		const int need = (int)(p1 - p0) * nen;
		if (need > cap) {
			cap = need;
			work = (int32_t *)realloc(work, (size_t)cap * sizeof(int32_t));
		}
		int m = 0;
		for (int64_t p = p0; p < p1; p++) {
			int32_t nd[10];
			tet_nodes(nlist[p], nd);
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
			int32_t nd[10];
			tet_nodes(nlist[p], nd);
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
int solve_nodal_test(FILE *fp_log)
{
	int ierr = 0;

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
