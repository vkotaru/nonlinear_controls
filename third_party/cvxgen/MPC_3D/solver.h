/* Produced by CVXGEN, 2020-07-28 18:32:12 -0400.  */
/* CVXGEN is Copyright (C) 2006-2017 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2017 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: solver.h. */
/* Description: Header file with relevant definitions. */
#ifndef SOLVER_H_MPC3D
#define SOLVER_H_MPC3D
#include <stdio.h>
#include <math.h>
#define pm(A, m, n) printmatrix(#A, A, m, n, 1)
typedef struct Params_t {
  double x_0[6];
  double Q[36];
  double R[9];
  double Q_final[36];
  double A[36];
  double B[18];
  double u_max[1];
  double S[1];
  double *x[1];
} Params;
typedef struct Vars_t {
  double *u_0; /* 3 rows. */
  double *x_1; /* 6 rows. */
  double *u_1; /* 3 rows. */
  double *x_2; /* 6 rows. */
  double *u_2; /* 3 rows. */
  double *x_3; /* 6 rows. */
  double *u_3; /* 3 rows. */
  double *x_4; /* 6 rows. */
  double *u_4; /* 3 rows. */
  double *x_5; /* 6 rows. */
  double *u_5; /* 3 rows. */
  double *x_6; /* 6 rows. */
  double *u_6; /* 3 rows. */
  double *x_7; /* 6 rows. */
  double *u_7; /* 3 rows. */
  double *x_8; /* 6 rows. */
  double *u_8; /* 3 rows. */
  double *x_9; /* 6 rows. */
  double *u_9; /* 3 rows. */
  double *x_10; /* 6 rows. */
  double *u_10; /* 3 rows. */
  double *x_11; /* 6 rows. */
  double *t_01; /* 3 rows. */
  double *t_02; /* 3 rows. */
  double *t_03; /* 3 rows. */
  double *t_04; /* 3 rows. */
  double *t_05; /* 3 rows. */
  double *t_06; /* 3 rows. */
  double *t_07; /* 3 rows. */
  double *t_08; /* 3 rows. */
  double *t_09; /* 3 rows. */
  double *t_10; /* 3 rows. */
  double *t_11; /* 3 rows. */
  double *t_12; /* 1 rows. */
  double *t_13; /* 1 rows. */
  double *t_14; /* 1 rows. */
  double *t_15; /* 1 rows. */
  double *t_16; /* 1 rows. */
  double *t_17; /* 1 rows. */
  double *t_18; /* 1 rows. */
  double *t_19; /* 1 rows. */
  double *t_20; /* 1 rows. */
  double *t_21; /* 1 rows. */
  double *u[11];
  double *x[12];
} Vars;
typedef struct Workspace_t {
  double h[169];
  double s_inv[169];
  double s_inv_z[169];
  double b[66];
  double q[142];
  double rhs[546];
  double x[546];
  double *s;
  double *z;
  double *y;
  double lhs_aff[546];
  double lhs_cc[546];
  double buffer[546];
  double buffer2[546];
  double KKT[1783];
  double L[2351];
  double d[546];
  double v[546];
  double d_inv[546];
  double gap;
  double optval;
  double ineq_resid_squared;
  double eq_resid_squared;
  double block_33[1];
  /* Pre-op symbols. */
  double quad_645199597568[1];
  int converged;
} Workspace;
typedef struct Settings_t {
  double resid_tol;
  double eps;
  int max_iters;
  int refine_steps;
  int better_start;
  /* Better start obviates the need for s_init and z_init. */
  double s_init;
  double z_init;
  int verbose;
  /* Show extra details of the iterative refinement steps. */
  int verbose_refinement;
  int debug;
  /* For regularization. Minimum value of abs(D_ii) in the kkt D factor. */
  double kkt_reg;
} Settings;
extern Vars vars;
extern Params params;
extern Workspace work;
extern Settings settings;
/* Function definitions in ldl.c: */
void ldl_solve(double *target, double *var);
void ldl_factor(void);
double check_factorization(void);
void matrix_multiply(double *result, double *source);
double check_residual(double *target, double *multiplicand);
void fill_KKT(void);

/* Function definitions in matrix_support.c: */
void multbymA(double *lhs, double *rhs);
void multbymAT(double *lhs, double *rhs);
void multbymG(double *lhs, double *rhs);
void multbymGT(double *lhs, double *rhs);
void multbyP(double *lhs, double *rhs);
void fillq(void);
void fillh(void);
void fillb(void);
void pre_ops(void);

/* Function definitions in solver.c: */
double eval_gap(void);
void set_defaults(void);
void setup_pointers(void);
void setup_indexed_params(void);
void setup_indexed_optvars(void);
void setup_indexing(void);
void set_start(void);
double eval_objv(void);
void fillrhs_aff(void);
void fillrhs_cc(void);
void refine(double *target, double *var);
double calc_ineq_resid_squared(void);
double calc_eq_resid_squared(void);
void better_start(void);
void fillrhs_start(void);
long solve(void);

/* Function definitions in testsolver.c: */
int main(int argc, char **argv);
void load_default_data(void);

/* Function definitions in util.c: */
void tic(void);
float toc(void);
float tocq(void);
void printmatrix(char *name, double *A, int m, int n, int sparse);
double unif(double lower, double upper);
float ran1(long*idum, int reset);
float randn_internal(long *idum, int reset);
double randn(void);
void reset_rand(void);

#endif
