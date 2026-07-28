/*
fieldout.c

解いた場を可視化用に書き出す (VTK legacy ASCII、`fieldout = 1`)。

これまで出力は集中定数 (R/L/C 行列、rlc.csv、SPICE) だけで、解いた場そのものを
見る手段が無かった。外部ライブラリに依存せずに ParaView / VisIt で開ける形式が
要るので、レガシー VTK の ASCII を選ぶ (書式が単純で、読み手が最も多い)。

  構造格子   -> DATASET RECTILINEAR_GRID
  非構造格子 -> DATASET UNSTRUCTURED_GRID (VTK_TETRA = 10)

**節点の並びに注意。** VTK の構造格子は x が最内 (x, y, z の順に遅くなる) だが、
本コードの node_index() は z が最内なので、書き出しでは並べ替える。
セルも同じ (VTK は i が最内、CellMaterial は k が最内)。

解は solve() の途中でしか手元に無いので、その場で field_add_* に渡して溜め、
最後に field_write() で一度に書く。fieldout = 0 (既定) なら何も溜めず何も書かない。
*/

#include "fem.h"
#include "fem_prototype.h"


int64_t num_cell(void)
{
	if (MeshMode) return ((MeshDim == 2) ? (int64_t)NTri : (int64_t)NTet);

	return ((int64_t)Nx * Ny * Nz);
}


// 構造格子のセル番号 (k が最内)
static int64_t cell_index(int i, int j, int k)
{
	return (((int64_t)i * Ny) + j) * Nz + k;
}


// ---- 蓄積 ----

static void name_copy(char *dst, const char *src)
{
	// VTK の配列名に空白は使えないので '_' に落とす
	int n = 0;
	for (; (src[n] != '\0') && (n < 62); n++) {
		dst[n] = ((src[n] == ' ') ? '_' : src[n]);
	}
	dst[n] = '\0';
}


void field_add_node(const char *name, const double *u)
{
	if (!FieldOut || (NFieldN >= MAXFIELD)) return;

	const int64_t n = num_node();
	FieldN[NFieldN] = (double *)malloc((size_t)n * sizeof(double));
	memcpy(FieldN[NFieldN], u, (size_t)n * sizeof(double));
	name_copy(FieldNName[NFieldN], name);
	NFieldN++;
}


void field_add_cellvec(const char *name, const double *v)
{
	if (!FieldOut || (NFieldC >= MAXFIELD)) return;

	const int64_t nc = num_cell();
	FieldC[NFieldC] = (double *)malloc((size_t)nc * 3 * sizeof(double));
	memcpy(FieldC[NFieldC], v, (size_t)nc * 3 * sizeof(double));
	name_copy(FieldCName[NFieldC], name);
	NFieldC++;
}


// 節点スカラー u からセル中心の -∇u (kind 0) または
// 2 次元断面の回転 B = ∇×(u ẑ) (kind 1、伝送線路軸 Tline に垂直な 2 軸) を作る。
// 構造格子では中心での 3 重線形補間の微分 = 各方向 4 本の辺差分の平均 (厳密)、
// 四面体では ∇λ が要素内一定なので厳密。
void field_add_grad(const char *name, const double *u, int kind)
{
	if (!FieldOut) return;

	const int64_t nc = num_cell();
	double *v = (double *)malloc((size_t)nc * 3 * sizeof(double));

	if (MeshMode && (MeshDim == 2)) {
		// 断面 2 次元 : 面内 2 軸だけが微分に効く
		int p, q;
		tri_axes(&p, &q);
		for (int e = 0; e < NTri; e++) {
			const int32_t *nd = &Tri[e * 3];
			double g[3][2], area, d[3] = {0, 0, 0};
			if (!tri_grad(nd, g, &area)) {
				for (int a = 0; a < 3; a++) {
					d[p] += u[nd[a]] * g[a][0];
					d[q] += u[nd[a]] * g[a][1];
				}
			}
			for (int c = 0; c < 3; c++) v[(e * 3) + c] = -d[c];
		}
	}
	else if (MeshMode) {
		for (int e = 0; e < NTet; e++) {
			double gn[10][3], d[3] = {0, 0, 0};
			int nen = 0;
			if (!tet_grad_center(e, gn, &nen)) {
				int32_t nd[10];
				tet_nodes(e, nd);
				for (int a = 0; a < nen; a++) {
					const double ua = u[nd[a]];
					for (int c = 0; c < 3; c++) d[c] += ua * gn[a][c];
				}
			}
			for (int c = 0; c < 3; c++) v[(e * 3) + c] = -d[c];
		}
	}
	else {
		int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
		for (int k = 0; k < Nz; k++) {
			const double dx = Xn[i + 1] - Xn[i];
			const double dy = Yn[j + 1] - Yn[j];
			const double dz = Zn[k + 1] - Zn[k];
			double p[8];
			for (int l = 0; l < 8; l++) {
				p[l] = u[node_index(i + ((l >> 2) & 1), j + ((l >> 1) & 1), k + (l & 1))];
			}
			// 中心での偏微分 (各方向 4 本の辺差分の平均)
			const double gx = ((p[4] - p[0]) + (p[5] - p[1]) + (p[6] - p[2]) + (p[7] - p[3])) / (4 * dx);
			const double gy = ((p[2] - p[0]) + (p[3] - p[1]) + (p[6] - p[4]) + (p[7] - p[5])) / (4 * dy);
			const double gz = ((p[1] - p[0]) + (p[3] - p[2]) + (p[5] - p[4]) + (p[7] - p[6])) / (4 * dz);
			const int64_t c = cell_index(i, j, k);
			v[(c * 3) + 0] = -gx;
			v[(c * 3) + 1] = -gy;
			v[(c * 3) + 2] = -gz;
		}
		}
		}
	}

	if (kind == 1) {
		// B = ∇×(u ẑ_t)。伝送線路軸 t に対し B_p = ∂u/∂q, B_q = -∂u/∂p
		// (v には -∇u が入っているので符号に注意)
		const int t = ((Tline == 'X') ? 0 : (Tline == 'Y') ? 1 : 2);
		const int p = (t + 1) % 3, q = (t + 2) % 3;
		for (int64_t c = 0; c < nc; c++) {
			const double gp = -v[(c * 3) + p];		// ∂u/∂p
			const double gq = -v[(c * 3) + q];		// ∂u/∂q
			v[(c * 3) + t] = 0;
			v[(c * 3) + p] = gq;
			v[(c * 3) + q] = -gp;
		}
	}

	field_add_cellvec(name, v);
	free(v);
}


void field_free(void)
{
	for (int i = 0; i < NFieldN; i++) free(FieldN[i]);
	for (int i = 0; i < NFieldC; i++) free(FieldC[i]);
	NFieldN = 0;
	NFieldC = 0;
}


// ---- 書き出し ----

static void write_grid(FILE *fp)
{
	if (MeshMode && (MeshDim == 2)) {
		fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
		fprintf(fp, "POINTS %d double\n", NNode);
		for (int i = 0; i < NNode; i++) {
			fprintf(fp, "%.9e %.9e %.9e\n", Xp[i], Yp[i], Zp[i]);
		}
		fprintf(fp, "\nCELLS %d %d\n", NTri, 4 * NTri);
		for (int e = 0; e < NTri; e++) {
			fprintf(fp, "3 %d %d %d\n", Tri[e * 3], Tri[(e * 3) + 1], Tri[(e * 3) + 2]);
		}
		fprintf(fp, "\nCELL_TYPES %d\n", NTri);
		for (int e = 0; e < NTri; e++) fprintf(fp, "5\n");	// VTK_TRIANGLE
	}
	else if (MeshMode) {
		fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
		fprintf(fp, "POINTS %d double\n", NNode);
		for (int i = 0; i < NNode; i++) {
			fprintf(fp, "%.9e %.9e %.9e\n", Xp[i], Yp[i], Zp[i]);
		}
		// VTK_QUADRATIC_TETRA (24) の中間節点の並びは Gmsh の tet10 と違う。
		//   VTK  : (0,1) (1,2) (0,2) (0,3) (1,3) (2,3)
		//   Gmsh : (0,1) (1,2) (2,0) (3,0) (3,2) (3,1)
		// 頂点からの並べ替えは最後の 2 つの入れ替えになる
		static const int g2v[6] = {0, 1, 2, 3, 5, 4};
		const int nen = ((TetOrder >= 2) ? 10 : 4);
		fprintf(fp, "\nCELLS %d %d\n", NTet, (nen + 1) * NTet);
		for (int e = 0; e < NTet; e++) {
			int32_t nd[10];
			tet_nodes(e, nd);
			fprintf(fp, "%d %d %d %d %d", nen, nd[0], nd[1], nd[2], nd[3]);
			for (int l = 4; l < nen; l++) fprintf(fp, " %d", nd[4 + g2v[l - 4]]);
			fprintf(fp, "\n");
		}
		fprintf(fp, "\nCELL_TYPES %d\n", NTet);
		// VTK_TETRA = 10, VTK_QUADRATIC_TETRA = 24
		for (int e = 0; e < NTet; e++) fprintf(fp, "%d\n", ((nen == 10) ? 24 : 10));
	}
	else {
		fprintf(fp, "DATASET RECTILINEAR_GRID\n");
		fprintf(fp, "DIMENSIONS %d %d %d\n", Nx + 1, Ny + 1, Nz + 1);
		fprintf(fp, "X_COORDINATES %d double\n", Nx + 1);
		for (int i = 0; i <= Nx; i++) fprintf(fp, "%.9e\n", Xn[i]);
		fprintf(fp, "Y_COORDINATES %d double\n", Ny + 1);
		for (int j = 0; j <= Ny; j++) fprintf(fp, "%.9e\n", Yn[j]);
		fprintf(fp, "Z_COORDINATES %d double\n", Nz + 1);
		for (int k = 0; k <= Nz; k++) fprintf(fp, "%.9e\n", Zn[k]);
	}
}


// 節点データ 1 本 (構造格子は VTK の並び x 最内に直す)
static void write_node_array(FILE *fp, const char *name, const double *u)
{
	fprintf(fp, "SCALARS %s double 1\nLOOKUP_TABLE default\n", name);
	if (MeshMode) {
		for (int i = 0; i < NNode; i++) fprintf(fp, "%.9e\n", u[i]);
	}
	else {
		for (int k = 0; k <= Nz; k++) {
		for (int j = 0; j <= Ny; j++) {
		for (int i = 0; i <= Nx; i++) {
			fprintf(fp, "%.9e\n", u[node_index(i, j, k)]);
		}
		}
		}
	}
}


static void write_cell_vec(FILE *fp, const char *name, const double *v)
{
	fprintf(fp, "VECTORS %s double\n", name);
	if (MeshMode) {
		const int64_t ne = num_cell();
		for (int64_t e = 0; e < ne; e++) {
			fprintf(fp, "%.9e %.9e %.9e\n", v[e * 3], v[(e * 3) + 1], v[(e * 3) + 2]);
		}
	}
	else {
		for (int k = 0; k < Nz; k++) {
		for (int j = 0; j < Ny; j++) {
		for (int i = 0; i < Nx; i++) {
			const int64_t c = cell_index(i, j, k);
			fprintf(fp, "%.9e %.9e %.9e\n", v[c * 3], v[(c * 3) + 1], v[(c * 3) + 2]);
		}
		}
		}
	}
}


static void write_cell_int(FILE *fp, const char *name, int structured_from_cell)
{
	fprintf(fp, "SCALARS %s int 1\nLOOKUP_TABLE default\n", name);
	if (MeshMode && (MeshDim == 2)) {
		// 2 次元格子では導体番号も要素データとして持っている
		for (int e = 0; e < NTri; e++) {
			fprintf(fp, "%d\n", (structured_from_cell ? (int)TriCond[e] : (int)TriMat[e]));
		}
	}
	else if (MeshMode) {
		for (int e = 0; e < NTet; e++) fprintf(fp, "%d\n", (int)TetMat[e]);
	}
	else {
		for (int k = 0; k < Nz; k++) {
		for (int j = 0; j < Ny; j++) {
		for (int i = 0; i < Nx; i++) {
			const int64_t c = cell_index(i, j, k);
			fprintf(fp, "%d\n", (structured_from_cell
				? (int)CellConductor[c] : (int)CellMaterial[c]));
		}
		}
		}
	}
}


int field_write(FILE *fp_log)
{
	if (!FieldOut) return 0;
	if ((NFieldN == 0) && (NFieldC == 0)) {
		fprintf(fp_log, "*** warning : fieldout = 1 but the selected analysis "
			"produced no field to write\n");
		return 0;
	}

	FILE *fp = fopen(FN_field, "w");
	if (fp == NULL) {
		fprintf(fp_log, "*** %s open error\n", FN_field);
		return 1;
	}

	fprintf(fp, "# vtk DataFile Version 3.0\n");
	fprintf(fp, "%s (%s)\n", ((Title[0] != '\0') ? Title : PROGRAM), PROGRAM);
	fprintf(fp, "ASCII\n");
	write_grid(fp);

	if (NFieldN > 0) {
		fprintf(fp, "\nPOINT_DATA %lld\n", (long long)num_node());
		for (int i = 0; i < NFieldN; i++) {
			write_node_array(fp, FieldNName[i], FieldN[i]);
		}
	}

	fprintf(fp, "\nCELL_DATA %lld\n", (long long)num_cell());
	write_cell_int(fp, "material", 0);
	if (MeshMode ? (MeshDim == 2) : (CellConductor != NULL)) {
		write_cell_int(fp, "conductor", 1);
	}
	for (int i = 0; i < NFieldC; i++) {
		write_cell_vec(fp, FieldCName[i], FieldC[i]);
	}

	fclose(fp);
	fprintf(fp_log, "  field output : %s (%lld points, %lld cells, "
		"%d node arrays, %d cell vectors)\n",
		FN_field, (long long)num_node(), (long long)num_cell(), NFieldN, NFieldC);

	return 0;
}
