// fem_prototype.h

#ifndef _FEM_PROTOTYPE_H_
#define _FEM_PROTOTYPE_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdint.h>
#include "fem.h"

// input_data.c
int input_data(FILE *);
void input_warn(const char *, ...);
int material_freq(void);

// setup.c
int setup(void);
void memfree(void);
int64_t node_index(int, int, int);
int64_t num_node(void);

// ingeometry.c
int ingeometry(double, double, double, int, const double *, double);

// unstruct.c
int mesh_read(const char *);
int tet_nodes(int, int32_t [TET_NMAX]);
int tet_nen(void);
int tet20_element(int, const double [6], double [TET_NMAX][TET_NMAX], double *);
int hex_nodes(int, int32_t [ELEM_NMAX]);
int prism_nodes(int, int32_t [PRISMP_NMAX]);
void crs_alloc_tet(crs_t *);
void crs_alloc_hex(crs_t *);
void crs_alloc_elem3d(crs_t *);
void assemble_elem3d(crs_t *, int);
int elem3d_count(void);
int elem3d_kind(int);
int elem3d_nodes(int, int32_t [ELEM_NMAX]);
int elem3d_mat(int);
int elem3d_tag(int);
int elem3d_nip(int);
int tet_ho_point(int, const double [4], double [ELEM_NMAX], double [ELEM_NMAX][3], double *);
int elem3d_ip(int, int, double [ELEM_NMAX], double [ELEM_NMAX][3], double *);
int tet_face_shape(int, const int [3], const double [3], double [TET_NMAX], double [3], double [3]);
int elem3d_face_point(int, const int32_t *, int, double, double, double [ELEM_NMAX], double [3], double [3], double [3]);
int face_quad_pub(int, int, double [][2], double *);
int tri_quad_pub(int, double [][3], double *);
int elem3d_face_nip(int);
int elem3d_face_ip(int, const int32_t *, int, int, double [ELEM_NMAX], double *);
void assemble_hex(crs_t *, int);
int hex_grad_center(int, double [8][3]);
int hex2_grad_center(int, double [HEXP_NMAX][3], int *);
void crs_alloc_prism(crs_t *);
void assemble_prism(crs_t *, int);
int prism_grad_center(int, double [6][3]);
int prism2_grad_center(int, double [PRISMP_NMAX][3], int *);
int pyr_grad_center(int, double [5][3]);
void assemble_tet(crs_t *, int);
int tet10_element(int, const double [6], double [10][10], double *);
int tet_grad_center(int, double [TET_NMAX][3], int *);
int solve_nodal_test(FILE *);
int tri_nodes(int, int32_t [TRI_NMAX]);
int tri_weights(int, double [TRI_NMAX], double [TRI_NMAX]);
void tri_axes(int *, int *);
int tri_grad(const int32_t [3], double [3][2], double *);
int64_t crs_find_tri(const crs_t *, int32_t, int32_t);
void crs_alloc_tri(crs_t *);
void assemble_nu_tri(crs_t *, const double *);
void assemble_mass_tri(crs_t *);
void assemble_mass_tri_sig(crs_t *, const double *);
int eddy2d_heat_solve(FILE *, const double *, double, double *, double *, double *, int);
double tri6_area(int);
int tri_grad_center(int, double [TRI_NMAX][2], int *);

// edge.c
void edge_build(void);
int edge_elem_count(void);
int edge_elem_kind(int);
int edge_elem_nedge(int);
int edge_elem_nen(int);
const int32_t *edge_elem_nodes(int);
const int (*edge_elem_table(int))[2];
int edge_elem_mat(int);
void edge_elem_matrices(int, const double [6], double, double [12][12], double [12][12]);
void edge_elem_center(int, double [12][3], double [12][3]);
int64_t edge_id(int32_t, int32_t);
void edge_free(void);
void crs_alloc_edge(crs_t *);
void edge_element(int, const double [6], double, double [6][6], double [6][6]);
void assemble_edge(crs_t *, crs_t *);
int solve_edge_test(FILE *);
int edge_tree(unsigned char *);
void edge_grad(const double *, double *);
void edge_gradT(const double *, double *);
void edge_nodal_aux(const crs_t *, crs_t *);
int solver_cg_edge(const crs_t *, const crs_t *, const double *, double *,
	const unsigned char *, int, int, int, double, FILE *, const char *);
int tet_grad_pub(const int32_t [4], double [4][3], double *);

// eddy3d.c
void crs_alloc_edge_node(crs_t *);
void assemble_eddy3d(crs_t *, crs_t *, double, const double *);
int solve_eddy3d(FILE *);
int solve_eddy3d_transient(FILE *);
int solve_eddy3d_p2(FILE *);
int solve_edge_test_p2(FILE *);
int edge_tree_gauge(const unsigned char *, unsigned char *, int *);
int solve_eddy3d_p2_transient(FILE *);
int eddy3d_p2_heat_open(FILE *);
int eddy3d_p2_heat_solve(FILE *, const double *, double, double *, double *, double *, int);
void eddy3d_p2_heat_close(void);
int32_t *p2_faces(int *, int32_t **);
int p2_face_id(const int32_t *, int, const int32_t [3]);

// eddy3d_pr.c (3 次以上の Nedelec)
int solve_eddy3d_pr(FILE *);
int solve_eddy3d_pr_transient(FILE *);
int solve_edge_test_pr(FILE *);
int eddy3d_pr_heat_open(FILE *);
int eddy3d_pr_heat_solve(FILE *, const double *, double, double *, double *, double *, int);
void eddy3d_pr_heat_close(void);
const int *tet_ho_index(int);

// assemble.c (公開した材料係数)
void material_coef_pub(int, int, double [6]);

// fieldout.c
int64_t num_cell(void);
void field_add_node(const char *, const double *);
void field_add_cellvec(const char *, const double *);
void field_add_grad(const char *, const double *, int);
void field_cell_grad(const double *, int, double *);
void field_free(void);
int field_write(FILE *);

// hdf5out.c
int  h5_enabled(void);
int  h5_open(FILE *);
void h5_close(FILE *);
int  h5_add_sweep(double);
int  h5_add_hysteresis(int, double, double, double, double, int);
int  h5_add_field(double);
int  h5_nfield(void);

// crs.c
void crs_alloc(crs_t *);
void crs_free(crs_t *);
void crs_alloc_dist(crs_t *);
int crs_distribute(crs_t *);
int crs_distribute_pair(FILE *, crs_t *, crs_t *);
int crs_owns(const crs_t *, int64_t);
int crs_dist_info(const crs_t *, int64_t *, int64_t *);
void crs_zero(crs_t *);
int64_t crs_offset(int64_t, int, int, int, int, int, int);
void crs_spmv(const crs_t *, const double *, double *, const unsigned char *);
void crs_spmv_c(const crs_t *, const crs_t *, double,
	const double *, const double *, double *, double *, const unsigned char *,
	double *, double *, double *, double *);
void crs_diag(const crs_t *, double *);
double crs_row_dot(const crs_t *, int64_t, const double *);

// assemble.c
void element_matrix(double, double, double, const double [6], double [8][8]);
int tensor6_inverse(const double [6], double [6]);
void assemble(crs_t *, int);
void assemble_nu(crs_t *, const double *);
void assemble_cell(crs_t *, const double *);
void assemble_mass_cell(crs_t *, const double *);
double temp_factor(const material_t *, int);
double temp_factor_at(const material_t *, int, double);
double table_interp(int, const double *, const double *, double);
void joule_load(const double *, const double *, double *);
void assemble_mass(crs_t *);
double bh_nu(const material_t *, int, double);
void bh_nu_dnu_pub(const material_t *, double, double *, double *);
double bh_b_of_h(const material_t *, int, double);
void assemble_newton(crs_t *, const double *, double *);
void assemble_ja(const double *, double *, double *);
void ja_eval(const ja_t *, double, double, double, double, double *, double *, double *);
void ja_eval_vec(const ja_t *, const double [2], const double [2], const double [2], const double [2], const double [2], double [2], double [2], double [2][2]);
void ja_eval_vec3(const ja_t *, const double [3], const double [3], const double [3],
	const double [3], const double [3], double [3], double [3], double [6]);
double ja_vec3_selftest(const ja_t *);
void ja_commit(int64_t);
void assemble_nu6(crs_t *, const double *);
double ja_linear_energy(const double *);

// solver_cg.c
int solver_cg(const crs_t *, const double *, double *, const unsigned char *,
	int, int, double, FILE *, const char *);

// solver_direct.c
int solver_direct(const crs_t *, const double *, double *, const unsigned char *,
	FILE *, const char *);
int solver_direct_c(const crs_t *, const crs_t *, double,
	const double *, const double *, double *, double *, const unsigned char *,
	FILE *, const char *);

// solve.c
int solve(FILE *);
int solve_port(const crs_t *, const unsigned char *, int, double *, double *,
	double *, double *, FILE *, const char *);

// thermal.c
int solve_thermal(FILE *, double *, double *);
int eddy3d_heat_open(FILE *);
int eddy3d_heat_solve(FILE *, const double *, double, double *, double *, double *, int);
void eddy3d_heat_close(void);

// mpiwrap.c (WITH_MPI=OFF ではスタブ。mpi.h はこのモジュールの外に出さない)
int mpi_init(int *, char ***);
int mpi_rank(void);
int mpi_size(void);
void mpi_finalize(void);
void mpi_send_dbl(const double *, int, int, int);
void mpi_wait_sends(void);
void mpi_recv_dbl(double *, int, int, int);

// solver_cocg.c
int solver_cocg(const crs_t *, const crs_t *, double,
	const double *, const double *, double *, double *, const unsigned char *,
	int, int, double, FILE *, const char *);

// matutil.c
int mat_inverse(const double *, double *, int);
int mat_inverse_c(const double *, const double *, double *, double *, int);

// outputRLC.c
void outputRLC(FILE *);
void writeout(FILE *);

// utils.c
double cputime(void);
void monitor1(FILE *, const char *);
void monitor2(FILE *, int);
void comline(int, char *[], int *, int *, char *);

#ifdef __cplusplus
}
#endif

// matutil.c : Anderson 加速 (非線形渦電流の実効磁気抵抗率の反復)
#define ANDERSON_MAX (8)
void anderson_step(int64_t, double *, const double *, double *, double *, double *, double *, int, int *, double);

// mpiwrap.c (行ブロックの行列ベクトル積・出力)
void mpi_rows_enable(int);
int mpi_enabled(void);
void mpi_row_range(int, int *, int *);
void mpi_allgather_rows(double *, int);
int mpi_row_owner(int, int);
void mpi_bcast_dbl(double *, int, int);
FILE *mpi_fopen_out(const char *, const char *);

// expr.c
int expr_compile(const char *, expr_t *);
double expr_eval(const expr_t *, double, double, double);
int expr_uses_t(const expr_t *);
double expr_eval_t(const expr_t *, double, double, double, double, double *);

#endif		// _FEM_PROTOTYPE_H_
