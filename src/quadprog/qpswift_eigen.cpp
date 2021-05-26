//
// Created by kotaru on 5/26/21.
//
#include "quadprog/qpswift_eigen.h"

namespace nonlinear_controls {

QPSwiftEigen::QPSwiftEigen(const int &n, const int &m, const int &p)
    : QuadProg(n, m, p) {}

int QPSwiftEigen::solve() {
  QP *myQP;

  // Setup Function //
  myQP = QP_SETUP_dense(this->nv, this->ni, this->ne, this->H.data(),
                        this->Aeq.data(), this->A.data(), this->f.data(),
                        this->b.data(), this->beq.data(), NULL);

  /****************************************
   *	After this, you can change the solver settings like this
   *	myQP->options->maxit  = 30   (to change the maximum number of
   iterations *to 30; default is 100)
   *
   * myQP->options->reltol = 1e-3 (to change the Relative tolerance to 1e-3;
   *default is 1e-6)
   *
   * myQP->options->abstol  = 1e-3 (to change the Absolute
   *tolerance to 1e-3; default is 1e-6)
   *
   * myQP->options->SIGMA  = 50 (to change the SIGMA to 50; default is 100;
   *recommended not to change this)
   *
   *myQP->options->VERBOSE  = 0 (displays no output when set to 0; default
   is *1 which corresponds to complete verbose mode)
   ******************************************/

  /* The Solution can be found as real pointer in myQP->x;It is an array of
   * Dimension n */

  qp_int ExitCode = QP_SOLVE(myQP);

  std::cout << "Solution" << std::endl;

  for (int i = 0; i < 3; ++i) {
    this->xOpt(i) = myQP->x[i];
  }
  std::cout << "xOpt: " << this->xOpt.transpose() << std::endl;

  QP_CLEANUP_dense(myQP);
  return 1;
}

} // namespace nonlinear_controls
