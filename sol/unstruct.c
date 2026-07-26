/*
unstruct.c

非構造格子 (4 節点四面体) の読み込み・CRS 構築・要素行列。

格子は Gmsh ASCII 2.2 形式 (.msh) で与える。物理タグで領域と電極を指定し、
`region` / `electrode` キーで材料番号・導体番号に対応づける。

構造格子との違いは「節点の並びと隣接関係」だけなので、Dirichlet の扱い・
反作用からの電荷抽出・反復解法は構造格子版をそのまま使える。
現時点では静電系 (C / L / R) のみ対応する (M / F は断面 2 次元の定式化なので
構造格子専用)。
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

		if (type == 4) {
			// 四面体
			if (nv < off + 4) continue;
			if (NTet % ARRAY_INC == 0) {
				Tet = (int32_t *)realloc(Tet, (size_t)(NTet + ARRAY_INC) * 4 * sizeof(int32_t));
				TetTag = (int *)realloc(TetTag, (size_t)(NTet + ARRAY_INC) * sizeof(int));
			}
			for (int l = 0; l < 4; l++) {
				const long g = v[off + l];
				if ((g < 0) || (g > maxid) || (idmap[g] < 0)) {
					printf("*** mesh : element %d refers to an unknown node %ld\n", e + 1, g);
					return 1;
				}
				Tet[(NTet * 4) + l] = idmap[g];
			}
			TetTag[NTet] = tag;
			NTet++;
		}
		else if (type == 2) {
			// 三角形 (電極面の指定に使う)
			if (nv < off + 3) continue;
			if (NTri % ARRAY_INC == 0) {
				Tri = (int32_t *)realloc(Tri, (size_t)(NTri + ARRAY_INC) * 3 * sizeof(int32_t));
				TriTag = (int *)realloc(TriTag, (size_t)(NTri + ARRAY_INC) * sizeof(int));
			}
			// 四面体と同じく未解決の節点番号は打ち切る。ここで continue すると
			// Tri[] の該当要素が未初期化のまま確定し、setup_unstruct() が
			// それを添字に使って領域外書き込みになる
			for (int l = 0; l < 3; l++) {
				const long g = v[off + l];
				if ((g < 0) || (g > maxid) || (idmap[g] < 0)) {
					printf("*** mesh : element %d refers to an unknown node %ld\n", e + 1, g);
					return 1;
				}
				Tri[(NTri * 3) + l] = idmap[g];
			}
			TriTag[NTri] = tag;
			NTri++;
		}
	}

	if (NTet < 1) {
		printf("%s\n", "*** mesh : no tetrahedron found");
		return 1;
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


// 四面体の連結から節点の隣接関係を作る (対角成分を含む)
void crs_alloc_tet(crs_t *A)
{
	const int n = NNode;

	// 節点毎の要素数を数える
	int *cnt = (int *)malloc((size_t)n * sizeof(int));
	memset(cnt, 0, (size_t)n * sizeof(int));
	for (int e = 0; e < NTet; e++) {
		for (int l = 0; l < 4; l++) cnt[Tet[(e * 4) + l]]++;
	}
	int64_t *nptr = (int64_t *)malloc(((size_t)n + 1) * sizeof(int64_t));
	nptr[0] = 0;
	for (int i = 0; i < n; i++) nptr[i + 1] = nptr[i] + cnt[i];
	int32_t *nlist = (int32_t *)malloc((size_t)nptr[n] * sizeof(int32_t));
	memset(cnt, 0, (size_t)n * sizeof(int));
	for (int e = 0; e < NTet; e++) {
		for (int l = 0; l < 4; l++) {
			const int32_t nd = Tet[(e * 4) + l];
			nlist[nptr[nd] + cnt[nd]] = (int32_t)e;
			cnt[nd]++;
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
		const int need = (int)(p1 - p0) * 4;
		if (need > cap) {
			cap = need;
			work = (int32_t *)realloc(work, (size_t)cap * sizeof(int32_t));
		}
		int m = 0;
		for (int64_t p = p0; p < p1; p++) {
			const int32_t e = nlist[p];
			for (int l = 0; l < 4; l++) work[m++] = Tet[(e * 4) + l];
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
			const int32_t e = nlist[p];
			for (int l = 0; l < 4; l++) work[m++] = Tet[(e * 4) + l];
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


// 全体行列の作成 (非構造格子)
//   K_ij = V (∇N_i)^T C (∇N_j)   1 次四面体では被積分関数が一定なので厳密
void assemble_tet(crs_t *A, int mode)
{
	crs_zero(A);

	for (int e = 0; e < NTet; e++) {
		const int32_t *nd = &Tet[e * 4];
		double g[4][3], vol;
		if (tet_grad_pub(nd, g, &vol)) continue;

		double c[6];
		material_coef_pub(TetMat[e], mode, c);
		if ((c[0] <= 0) && (c[1] <= 0) && (c[2] <= 0)) continue;

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
