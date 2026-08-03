/*
hdf5out.c

掃引・履歴の**系列**を HDF5 (`ofe_series.h5`) に書き出す (`hdf5 = 1`)。

これまで系列の出口は 2 つに分かれていた: 周波数掃引は `ofe_sweep.csv` に
1 行 1 点、電流掃引 (J-A ヒステリシス) の履歴は **`ofe.log` のテキスト表だけ**で
機械可読な形が無く、場は掃引の最後の 1 点しか残らなかった。表示側で
「周波数を動かしながら R(f) と場を見る」「B-H ループを描く」ことができない。
1 つの自己記述的なファイルにまとめて、その用途に使えるようにする。

**HDF5 は任意依存**。`-DWITH_HDF5=ON` のときだけ実体がコンパイルされ、
既定 (OFF) では下の `#else` のスタブになる。スタブは `hdf5 = 1` を
**黙って無視せず入力エラーにする** (書いたつもりでファイルが無い、が一番困る)。

書き方は**追記 (ストリーミング)**。掃引の各点で溜め込まずにその場で追記して
解放するので、メモリの使い方は従来と変わらない。第 1 次元を無制限にした
チャンク化データセットなので、途中で解を諦めてもファイルの長さは
「実際に計算できた点数」を正しく表す (0 で埋めた嘘が残らない)。

配置 (すべて属性で単位と軸を自己記述する):

  /                    format / version / title / program
  /sweep               frequency [nf]、C L G R Rf Lf [nf][np][np]
  /hysteresis          step [ns]、current H B L [ns]、iterations [ns]
  /mesh                構造格子 : x y z / 非構造格子 : points cells
                       material conductor [ncell]
  /field               axis [np]、節点スカラー [np][nnode]、
                       要素ベクトル [np][ncell][3]

節点・要素の並びは **VTK と同じ (構造格子は x が最内)** にそろえる。
`/mesh` の属性に `node_order` / `cell_order` と `dims` を書いてあるので、
表示側は reshape だけで格子に戻せる。ここを VTK と変えると、同じ解を
2 つの出口で見たときに食い違う。
*/

#include "fem.h"
#include "fem_prototype.h"


#ifdef OFE_HDF5

#include <hdf5.h>

#define SERIES_VERSION (1)

static hid_t Hfile = -1;
static hid_t Gsweep = -1, Ghyst = -1, Gfield = -1;
static int MeshDone = 0;
static int NptSweep = 0, NptHyst = 0, NptField = 0;
static double *Pack = NULL;			// 並べ替え用の作業配列
static int64_t NPack = 0;


int h5_enabled(void)
{
	return 1;
}


int h5_nfield(void)
{
	return NptField;
}


// ---- 小道具 ----

static void attr_str(hid_t obj, const char *name, const char *val)
{
	const hid_t sp = H5Screate(H5S_SCALAR);
	const hid_t tp = H5Tcopy(H5T_C_S1);
	H5Tset_size(tp, strlen(val) + 1);
	const hid_t at = H5Acreate2(obj, name, tp, sp, H5P_DEFAULT, H5P_DEFAULT);
	if (at >= 0) {
		H5Awrite(at, tp, val);
		H5Aclose(at);
	}
	H5Tclose(tp);
	H5Sclose(sp);
}


static void attr_int(hid_t obj, const char *name, const int *val, int n)
{
	const hsize_t dim = (hsize_t)n;
	const hid_t sp = H5Screate_simple(1, &dim, NULL);
	const hid_t at = H5Acreate2(obj, name, H5T_NATIVE_INT, sp, H5P_DEFAULT, H5P_DEFAULT);
	if (at >= 0) {
		H5Awrite(at, H5T_NATIVE_INT, val);
		H5Aclose(at);
	}
	H5Sclose(sp);
}


/*
系列の 1 点を追記する。データセットが無ければ作る。

第 1 次元 (点) を無制限にしたチャンク化データセットで、1 チャンク = 1 点。
ndim は 1 点あたりの次元数 (0 = スカラー、1 = ベクトル、2 = 行列)。
*/
static int append_row(hid_t grp, const char *name, hid_t type, const void *data,
	int ndim, const hsize_t *dims, const char *units)
{
	const int rank = ndim + 1;
	hsize_t cur[3], start[3], count[3];

	count[0] = 1;
	for (int i = 0; i < ndim; i++) {
		if (dims[i] < 1) return 0;			// 空の点は書かない
		count[i + 1] = dims[i];
	}

	hid_t ds = -1;
	if (H5Lexists(grp, name, H5P_DEFAULT) > 0) {
		ds = H5Dopen2(grp, name, H5P_DEFAULT);
		if (ds < 0) return 1;
		const hid_t sp = H5Dget_space(ds);
		H5Sget_simple_extent_dims(sp, cur, NULL);
		H5Sclose(sp);
		cur[0] += 1;
		if (H5Dset_extent(ds, cur) < 0) {
			H5Dclose(ds);
			return 1;
		}
	}
	else {
		hsize_t maxd[3], chunk[3];
		cur[0] = 1;
		maxd[0] = H5S_UNLIMITED;
		chunk[0] = 1;
		for (int i = 0; i < ndim; i++) {
			cur[i + 1] = maxd[i + 1] = chunk[i + 1] = dims[i];
		}
		const hid_t sp = H5Screate_simple(rank, cur, maxd);
		const hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
		H5Pset_chunk(dcpl, rank, chunk);
		ds = H5Dcreate2(grp, name, type, sp, H5P_DEFAULT, dcpl, H5P_DEFAULT);
		H5Pclose(dcpl);
		H5Sclose(sp);
		if (ds < 0) return 1;
		if (units != NULL) attr_str(ds, "units", units);
	}

	start[0] = cur[0] - 1;
	for (int i = 0; i < ndim; i++) start[i + 1] = 0;

	const hid_t fsp = H5Dget_space(ds);
	H5Sselect_hyperslab(fsp, H5S_SELECT_SET, start, NULL, count, NULL);
	const hid_t msp = H5Screate_simple(rank, count, NULL);
	const herr_t st = H5Dwrite(ds, type, msp, fsp, H5P_DEFAULT, data);
	H5Sclose(msp);
	H5Sclose(fsp);
	H5Dclose(ds);

	return ((st < 0) ? 1 : 0);
}


// 掃引しない量 (格子) を 1 回だけ書く
static int write_once(hid_t grp, const char *name, hid_t type, const void *data,
	int ndim, const hsize_t *dims, const char *units)
{
	const hid_t sp = H5Screate_simple(ndim, dims, NULL);
	const hid_t ds = H5Dcreate2(grp, name, type, sp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Sclose(sp);
	if (ds < 0) return 1;
	if (units != NULL) attr_str(ds, "units", units);
	const herr_t st = H5Dwrite(ds, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
	H5Dclose(ds);

	return ((st < 0) ? 1 : 0);
}


// ---- 格子 ----

// 節点スカラーを VTK と同じ並び (構造格子は x が最内) に詰め直す
static const double *pack_node(const double *u)
{
	if (MeshMode) return u;

	for (int k = 0; k <= Nz; k++) {
	for (int j = 0; j <= Ny; j++) {
	for (int i = 0; i <= Nx; i++) {
		Pack[(((int64_t)k * (Ny + 1)) + j) * (Nx + 1) + i] = u[node_index(i, j, k)];
	}
	}
	}

	return Pack;
}


// 要素ベクトルを VTK と同じ並び (構造格子は i が最内) に詰め直す
static const double *pack_cell(const double *v, int nc3)
{
	if (MeshMode) return v;

	int64_t p = 0;
	for (int k = 0; k < Nz; k++) {
	for (int j = 0; j < Ny; j++) {
	for (int i = 0; i < Nx; i++) {
		const int64_t c = (((int64_t)i * Ny) + j) * Nz + k;
		for (int d = 0; d < nc3; d++) Pack[p++] = v[(c * nc3) + d];
	}
	}
	}

	return Pack;
}


static int write_mesh(void)
{
	if (MeshDone) return 0;

	const hid_t g = H5Gcreate2(Hfile, "/mesh", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (g < 0) return 1;

	const int64_t nn = num_node();
	const int64_t nc = num_cell();
	int ierr = 0;
	hsize_t d[2];

	if (MeshMode) {
		attr_str(g, "type", "unstructured");
		attr_str(g, "node_order", "as read from the mesh file");
		attr_str(g, "cell_order", "as read from the mesh file");
		double *pt = (double *)malloc((size_t)NNode * 3 * sizeof(double));
		for (int i = 0; i < NNode; i++) {
			pt[(i * 3) + 0] = Xp[i];
			pt[(i * 3) + 1] = Yp[i];
			pt[(i * 3) + 2] = Zp[i];
		}
		d[0] = (hsize_t)NNode;
		d[1] = 3;
		ierr |= write_once(g, "points", H5T_NATIVE_DOUBLE, pt, 2, d, "m");
		free(pt);

		// 要素の節点表。並びは Gmsh のまま (VTK 用の入れ替えはしない) なので、
		// 属性に「Gmsh の並び」と明記しておく
		const int p2 = (TetOrder >= 2);
		const int nen = ((MeshDim == 2) ? (p2 ? 6 : 3) : (p2 ? 10 : 4));
		int32_t *cl = (int32_t *)malloc((size_t)nc * nen * sizeof(int32_t));
		for (int64_t e = 0; e < nc; e++) {
			int32_t nd[10];
			if (MeshDim == 2) tri_nodes((int)e, nd);
			else              tet_nodes((int)e, nd);
			for (int l = 0; l < nen; l++) cl[(e * nen) + l] = nd[l];
		}
		d[0] = (hsize_t)nc;
		d[1] = (hsize_t)nen;
		ierr |= write_once(g, "cells", H5T_NATIVE_INT32, cl, 2, d, NULL);
		free(cl);
		attr_str(g, "cell_node_order", "Gmsh");
		{
			const int v = nen;
			attr_int(g, "nodes_per_cell", &v, 1);
		}
	}
	else {
		attr_str(g, "type", "rectilinear");
		attr_str(g, "node_order", "x fastest : reshape to [nz+1][ny+1][nx+1]");
		attr_str(g, "cell_order", "x fastest : reshape to [nz][ny][nx]");
		d[0] = (hsize_t)(Nx + 1);  ierr |= write_once(g, "x", H5T_NATIVE_DOUBLE, Xn, 1, d, "m");
		d[0] = (hsize_t)(Ny + 1);  ierr |= write_once(g, "y", H5T_NATIVE_DOUBLE, Yn, 1, d, "m");
		d[0] = (hsize_t)(Nz + 1);  ierr |= write_once(g, "z", H5T_NATIVE_DOUBLE, Zn, 1, d, "m");
		const int dims[3] = {Nx + 1, Ny + 1, Nz + 1};
		attr_int(g, "dims", dims, 3);
	}

	// 材料・導体番号 (VTK の CELL_DATA と同じ内容・同じ並び)
	int32_t *m = (int32_t *)malloc((size_t)nc * sizeof(int32_t));
	int32_t *cd = (int32_t *)malloc((size_t)nc * sizeof(int32_t));
	int havecond = 0;
	if (MeshMode && (MeshDim == 2)) {
		for (int64_t e = 0; e < nc; e++) {
			m[e] = (int32_t)TriMat[e];
			cd[e] = (int32_t)TriCond[e];
		}
		havecond = 1;
	}
	else if (MeshMode) {
		for (int64_t e = 0; e < nc; e++) m[e] = (int32_t)TetMat[e];
	}
	else {
		int64_t p = 0;
		for (int k = 0; k < Nz; k++) {
		for (int j = 0; j < Ny; j++) {
		for (int i = 0; i < Nx; i++) {
			const int64_t c = (((int64_t)i * Ny) + j) * Nz + k;
			m[p] = (int32_t)CellMaterial[c];
			if (CellConductor != NULL) cd[p] = (int32_t)CellConductor[c];
			p++;
		}
		}
		}
		havecond = (CellConductor != NULL);
	}
	d[0] = (hsize_t)nc;
	ierr |= write_once(g, "material", H5T_NATIVE_INT32, m, 1, d, NULL);
	if (havecond) ierr |= write_once(g, "conductor", H5T_NATIVE_INT32, cd, 1, d, NULL);
	free(m);
	free(cd);

	H5Gclose(g);

	// 並べ替えの作業配列 (節点スカラーと要素ベクトルの大きい方)
	NPack = ((nn > (nc * 3)) ? nn : (nc * 3));
	Pack = (double *)malloc((size_t)NPack * sizeof(double));
	if (Pack == NULL) ierr = 1;

	MeshDone = 1;

	return ierr;
}


// ---- 開閉 ----

int h5_open(FILE *fp_log)
{
	if (!Hdf5Out) return 0;

	Hfile = H5Fcreate(FN_series, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if (Hfile < 0) {
		fprintf(fp_log, "*** %s create error\n", FN_series);
		printf("*** %s create error\n", FN_series);
		return 1;
	}

	char an[16];
	int na = 0;
	if (Analysis & ANALYSIS_C) an[na++] = 'C';
	if (Analysis & ANALYSIS_L) an[na++] = 'L';
	if (Analysis & ANALYSIS_R) an[na++] = 'R';
	if (Analysis & ANALYSIS_M) an[na++] = 'M';
	if (Analysis & ANALYSIS_F) an[na++] = 'F';
	if (Analysis & ANALYSIS_E) an[na++] = 'E';
	if (Analysis & ANALYSIS_A) an[na++] = 'A';
	if (Analysis & ANALYSIS_P) an[na++] = 'P';
	an[na] = '\0';

	const int ver = SERIES_VERSION;
	attr_str(Hfile, "format", "OpenFEM series");
	attr_int(Hfile, "version", &ver, 1);
	attr_str(Hfile, "program", PROGRAM);
	attr_str(Hfile, "title", ((Title[0] != '\0') ? Title : ""));
	attr_str(Hfile, "analysis", an);
	attr_int(Hfile, "nport", &NPort, 1);

	fprintf(fp_log, "  series output : %s (HDF5)\n", FN_series);

	return 0;
}


void h5_close(FILE *fp_log)
{
	if (Hfile < 0) return;

	if (Gsweep >= 0) H5Gclose(Gsweep);
	if (Ghyst >= 0)  H5Gclose(Ghyst);
	if (Gfield >= 0) H5Gclose(Gfield);
	H5Fclose(Hfile);
	Hfile = Gsweep = Ghyst = Gfield = -1;
	free(Pack);
	Pack = NULL;

	if (fp_log != NULL) {
		fprintf(fp_log, "  series : %s (%d sweep points, %d hysteresis steps, "
			"%d field snapshots)\n", FN_series, NptSweep, NptHyst, NptField);
	}
	MeshDone = 0;
	NptSweep = NptHyst = NptField = 0;
}


// ---- 系列の追記 ----

// 行列 1 本 (無い量は書かない。あるかどうかは点ごとに変わり得る)
static int add_mat(hid_t g, const char *name, int have, const double *mat, const char *units)
{
	if (!have) return 0;

	hsize_t d[2];
	d[0] = d[1] = (hsize_t)NPort;

	return append_row(g, name, H5T_NATIVE_DOUBLE, mat, 2, d, units);
}


int h5_add_sweep(double freq)
{
	if (Hfile < 0) return 0;

	if (Gsweep < 0) {
		Gsweep = H5Gcreate2(Hfile, "/sweep", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (Gsweep < 0) return 1;
		attr_str(Gsweep, "axis", "frequency");
		attr_str(Gsweep, "layout", "[point][port][port]");
	}

	// 単位長あたりか否か (tline) で単位が変わる。rlc.csv と同じ規約
	const int pul = (Tline != 0);
	int ierr = 0;
	hsize_t d1 = 1;
	ierr |= append_row(Gsweep, "frequency", H5T_NATIVE_DOUBLE, &freq, 0, &d1, "Hz");
	ierr |= add_mat(Gsweep, "C",  HaveC, Cmat,  (pul ? "F/m" : "F"));
	ierr |= add_mat(Gsweep, "L",  HaveL, Lmat,  (pul ? "H/m" : "H"));
	ierr |= add_mat(Gsweep, "G",  HaveR, Gmat,  (pul ? "S/m" : "S"));
	ierr |= add_mat(Gsweep, "R",  HaveR, Rmat,  (pul ? "ohm.m" : "ohm"));
	ierr |= add_mat(Gsweep, "M",  HaveM, Mmat,  (pul ? "H/m" : "H"));
	ierr |= add_mat(Gsweep, "Rf", HaveF, Rfmat, (pul ? "ohm/m" : "ohm"));
	ierr |= add_mat(Gsweep, "Lf", HaveF, Lfmat, (pul ? "H/m" : "H"));
	NptSweep++;

	return ierr;
}


int h5_add_hysteresis(int step, double cur, double h, double b, double l, int iter)
{
	if (Hfile < 0) return 0;

	if (Ghyst < 0) {
		Ghyst = H5Gcreate2(Hfile, "/hysteresis", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (Ghyst < 0) return 1;
		attr_str(Ghyst, "axis", "step");
		attr_str(Ghyst, "note", "history dependent : the order of the points is the "
			"magnetisation history and must not be sorted");
	}

	int ierr = 0;
	hsize_t d1 = 1;
	ierr |= append_row(Ghyst, "step", H5T_NATIVE_INT, &step, 0, &d1, NULL);
	ierr |= append_row(Ghyst, "current", H5T_NATIVE_DOUBLE, &cur, 0, &d1, "A");
	ierr |= append_row(Ghyst, "H", H5T_NATIVE_DOUBLE, &h, 0, &d1, "A/m");
	ierr |= append_row(Ghyst, "B", H5T_NATIVE_DOUBLE, &b, 0, &d1, "T");
	ierr |= append_row(Ghyst, "L", H5T_NATIVE_DOUBLE, &l, 0, &d1, (Tline ? "H/m" : "H"));
	ierr |= append_row(Ghyst, "iterations", H5T_NATIVE_INT, &iter, 0, &d1, NULL);
	NptHyst++;

	return ierr;
}


int h5_add_field(double axis)
{
	if (Hfile < 0) return 0;
	if ((NFieldN == 0) && (NFieldC == 0)) return 0;

	if (write_mesh()) return 1;

	if (Gfield < 0) {
		Gfield = H5Gcreate2(Hfile, "/field", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (Gfield < 0) return 1;
		attr_str(Gfield, "layout", "[point][node] / [point][cell][3]");
		attr_str(Gfield, "note", "the node and cell order is the one described in /mesh");
	}

	int ierr = 0;
	hsize_t d1 = 1;
	// 軸の値 (周波数掃引なら Hz、電流掃引なら A)。系列の並びは追記した順
	ierr |= append_row(Gfield, "axis", H5T_NATIVE_DOUBLE, &axis, 0, &d1, NULL);

	hsize_t dn = (hsize_t)num_node();
	for (int i = 0; i < NFieldN; i++) {
		ierr |= append_row(Gfield, FieldNName[i], H5T_NATIVE_DOUBLE,
			pack_node(FieldN[i]), 1, &dn, NULL);
	}
	hsize_t dc[2];
	dc[0] = (hsize_t)num_cell();
	dc[1] = 3;
	for (int i = 0; i < NFieldC; i++) {
		ierr |= append_row(Gfield, FieldCName[i], H5T_NATIVE_DOUBLE,
			pack_cell(FieldC[i], 3), 2, dc, NULL);
	}
	NptField++;

	return ierr;
}


#else		// OFE_HDF5 が無いビルド


int h5_enabled(void)
{
	return 0;
}


int h5_nfield(void)
{
	return 0;
}


int h5_open(FILE *fp_log)
{
	(void)fp_log;

	return 0;
}


void h5_close(FILE *fp_log)
{
	(void)fp_log;
}


int h5_add_sweep(double freq)
{
	(void)freq;

	return 0;
}


int h5_add_hysteresis(int step, double cur, double h, double b, double l, int iter)
{
	(void)step; (void)cur; (void)h; (void)b; (void)l; (void)iter;

	return 0;
}


int h5_add_field(double axis)
{
	(void)axis;

	return 0;
}


#endif		// OFE_HDF5
