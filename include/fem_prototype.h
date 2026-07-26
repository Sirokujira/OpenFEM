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

// setup.c
int setup(void);
void memfree(void);
int64_t node_index(int, int, int);
int64_t num_node(void);

// ingeometry.c
int ingeometry(double, double, double, int, const double *, double);

// unstruct.c
int mesh_read(const char *);
void crs_alloc_tet(crs_t *);
void assemble_tet(crs_t *, int);

// edge.c
void edge_build(void);
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

// assemble.c (公開した材料係数)
void material_coef_pub(int, int, double [6]);

// crs.c
void crs_alloc(crs_t *);
void crs_free(crs_t *);
void crs_zero(crs_t *);
int64_t crs_offset(int64_t, int, int, int, int, int, int);
void crs_spmv(const crs_t *, const double *, double *, const unsigned char *);
void crs_diag(const crs_t *, double *);
double crs_row_dot(const crs_t *, int64_t, const double *);

// assemble.c
void element_matrix(double, double, double, const double [6], double [8][8]);
int tensor6_inverse(const double [6], double [6]);
void assemble(crs_t *, int);
void assemble_nu(crs_t *, const double *);
void assemble_mass(crs_t *);
double bh_nu(const material_t *, int, double);
void assemble_newton(crs_t *, const double *, double *);
void assemble_ja(const double *, double *, double *);
void ja_eval(const ja_t *, double, double, double, double, double *, double *, double *);

// solver_cg.c
int solver_cg(const crs_t *, const double *, double *, const unsigned char *,
	int, int, double, FILE *, const char *);

// solve.c
int solve(FILE *);

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

#endif		// _FEM_PROTOTYPE_H_
