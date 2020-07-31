/* Produced by CVXGEN, 2020-07-30 04:42:42 -0400.  */
/* CVXGEN is Copyright (C) 2006-2017 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2017 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: matrix_support.c. */
/* Description: Support functions for matrix multiplication and vector filling. */
#include "solver.h"
void multbymA(double *lhs, double *rhs) {
}
void multbymAT(double *lhs, double *rhs) {
  lhs[0] = 0;
  lhs[1] = 0;
  lhs[2] = 0;
}
void multbymG(double *lhs, double *rhs) {
  lhs[0] = -rhs[0]*(params.A[0])-rhs[1]*(params.A[1])-rhs[2]*(params.A[2]);
  lhs[1] = -rhs[0]*(-1);
  lhs[2] = -rhs[1]*(-1);
  lhs[3] = -rhs[2]*(-1);
  lhs[4] = -rhs[0]*(1);
  lhs[5] = -rhs[1]*(1);
  lhs[6] = -rhs[2]*(1);
}
void multbymGT(double *lhs, double *rhs) {
  lhs[0] = -rhs[0]*(params.A[0])-rhs[1]*(-1)-rhs[4]*(1);
  lhs[1] = -rhs[0]*(params.A[1])-rhs[2]*(-1)-rhs[5]*(1);
  lhs[2] = -rhs[0]*(params.A[2])-rhs[3]*(-1)-rhs[6]*(1);
}
void multbyP(double *lhs, double *rhs) {
  /* TODO use the fact that P is symmetric? */
  /* TODO check doubling / half factor etc. */
  lhs[0] = rhs[0]*(2*params.Q[0])+rhs[1]*(2*params.Q[3])+rhs[2]*(2*params.Q[6]);
  lhs[1] = rhs[0]*(2*params.Q[1])+rhs[1]*(2*params.Q[4])+rhs[2]*(2*params.Q[7]);
  lhs[2] = rhs[0]*(2*params.Q[2])+rhs[1]*(2*params.Q[5])+rhs[2]*(2*params.Q[8]);
}
void fillq(void) {
  work.q[0] = params.c[0];
  work.q[1] = params.c[1];
  work.q[2] = params.c[2];
}
void fillh(void) {
  work.h[0] = params.b[0];
  work.h[1] = -params.xlb[0];
  work.h[2] = -params.xlb[1];
  work.h[3] = -params.xlb[2];
  work.h[4] = params.xub[0];
  work.h[5] = params.xub[1];
  work.h[6] = params.xub[2];
}
void fillb(void) {
}
void pre_ops(void) {
}
