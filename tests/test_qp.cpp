#include "common/qpoases_eigen.hpp"
#include "controls/clf_qp.h"
#include "quadprog/qpswift_eigen.h"
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <qpOASES.hpp>

namespace nlc = nonlinear_controls;

#define nlc_real double

int main() {
  USING_NAMESPACE_QPOASES

  /* Setup data of first QP. */
  real_t H[2 * 2] = {1.0, 0.0, 0.0, 0.5};
  real_t g[2] = {1.5, 1.0};
  real_t lb[2] = {0.5, -2.0};
  real_t ub[2] = {5.0, 2.0};

  /* Setup data of second QP. */
  real_t H_new[2 * 2] = {1.0, 0.5, 0.5, 0.5};
  real_t A_new[2 * 2] = {1.0, 5.0, 0.0, 0.0};
  real_t g_new[2] = {1.0, 1.5};
  real_t lb_new[2] = {-0.5, -1.0};
  real_t ub_new[2] = {5.0, -0.5};
  real_t lbA_new[2] = {-INFINITY, -INFINITY};
  real_t ubA_new[2] = {1.0, INFINITY};

  /* Setting up SQProblem object. */
  SQProblem *example;
  // example = new SQProblem( 2, 1);
  // real_t A[1*2] = { 1.0, 1.0};
  // real_t lbA[1] ={ -INFINITY};
  // real_t ubA[1] = { 1.0 };

  example = new SQProblem(2, 2);
  real_t A[2 * 2] = {1.0, 1.0, 0.0, 0.0};
  real_t lbA[2] = {-INFINITY, -INFINITY};
  real_t ubA[2] = {1.0, -INFINITY};

  /* Solve first QP. */
  int_t nWSR = 10;
  example->init(H, g, A, lb, ub, lbA, ubA, nWSR, 0);
  real_t xOpt[2];
  example->getPrimalSolution(xOpt);
  printf("\nxOpt = [ %e, %e ];  objVal = %e\n\n", xOpt[0], xOpt[1],
         example->getObjVal());

  /* Solve second QP. */
  nWSR = 10;
  example = new SQProblem(2, 2);
  real_t A_new2[2 * 2] = {1.0, 5.0, 1.0, -1.0};
  real_t lbA_new2[2] = {-INFINITY, -INFINITY};
  real_t ubA_new2[2] = {1.0, 0.75};
  example->init(H_new, g_new, A_new2, lb_new, ub_new, lbA_new2, ubA_new2, nWSR,
                0);

  /* Get and print solution of second QP. */
  example->getPrimalSolution(xOpt);
  printf("\nxOpt = [ %e, %e ];  objVal = %e\n\n", xOpt[0], xOpt[1],
         example->getObjVal());

  // nlc::ClfQP testqp;

  // testqp.problem_.H = 2 * Eigen::Matrix<nlc_real, 3, 3>::Identity();
  // // testqp.problem_.H(1,0) = -1;
  // // testqp.problem_.H(2,2) = 4;
  // testqp.problem_.f = Eigen::Matrix<nlc_real, 3, 1>::Zero();
  // testqp.problem_.A = Eigen::Matrix<nlc_real, 1, 3>::Ones();
  // testqp.problem_.b = -3 * Eigen::Matrix<nlc_real, 1, 1>::Ones();
  // testqp.problem_.xlb = -100 * Eigen::Matrix<nlc_real, 3, 1>::Ones();
  // testqp.problem_.xub = 100 * Eigen::Matrix<nlc_real, 3, 1>::Ones();

  // testqp.setup();
  // testqp.solve();
  // std::cout << testqp.getOptimizer() << std::endl;

  // USING_NAMESPACE_QPOASES
  // Eigen::Matrix<nlc_real, 1, 1> lbAmat, ubAmat; lbAmat << -100; ubAmat << -3;

  // /* Setup data of first QP. */
  // // real_t H[3*3], A[3], g[3], lb[3], ub[3],
  // // real_t lbA[1], ubA[1];
  // real_t *H, *A, *g, *lb, *ub, *lbA, *ubA;
  // H = testqp.problem_.H.transpose().data();
  // A = testqp.problem_.A.transpose().data();
  // g = testqp.problem_.f.transpose().data();
  // lb = testqp.problem_.xlb.transpose().data();
  // ub = testqp.problem_.xub.transpose().data();
  // lbA = lbAmat.transpose().data();
  // ubA = ubAmat.transpose().data();
  // // Eigen::Map<Eigen::Matrix<double, 3, 3>>(H, 3, 3) =
  // testqp.problem_.H.transpose(); for (int i = 0; i <9 ; i++) printf("H[%d]=
  // %f\n", i, H[i]);
  // // Eigen::Map<Eigen::Matrix<double, 3, 1>>(A, 3, 1) =
  // testqp.problem_.A.transpose();
  // // Eigen::Map<Eigen::Matrix<double, 1, 3>>(g, 1, 3) =
  // testqp.problem_.f.transpose();
  // // Eigen::Map<Eigen::Matrix<double, 1, 3>>(lb, 1, 3) =
  // testqp.problem_.xlb.transpose();
  // // Eigen::Map<Eigen::Matrix<double, 1, 3>>(ub, 1, 3) =
  // testqp.problem_.xub.transpose();
  // // lbA[0] = -100;
  // // ubA[0] = -3;
  // /* Setting up QProblem object. */
  // QProblem example(3, 1);

  // Options options;
  // example.setOptions(options);

  // /* Solve first QP. */
  // int_t nWSR = 10;
  // example.init(H, g, A, lb, ub, lbA, ubA, nWSR);

  // /* Get and print solution of first QP. */
  // real_t xOpt[3];
  // real_t yOpt[3 + 1];
  // example.getPrimalSolution(xOpt);
  // example.getDualSolution(yOpt);
  // printf("\nxOpt = [ %e, %e, %e ];  yOpt = [ %e, %e, %e ];  objVal = %e\n\n",
  //        xOpt[0], xOpt[1], xOpt[2], yOpt[0], yOpt[1], yOpt[2],
  //        example.getObjVal());
  // /*************************************/

  nlc::QPOasesEigen example2(3, 1);
  example2.data_.H = 2 * Eigen::Matrix<real_t, 3, 3>::Identity();
  example2.data_.g = Eigen::Matrix<real_t, 3, 1>::Zero();
  example2.data_.A = Eigen::Matrix<real_t, 1, 3>::Ones();
  example2.data_.lbA = -100 * Eigen::Matrix<real_t, 1, 1>::Ones();
  example2.data_.ubA = -3 * Eigen::Matrix<real_t, 1, 1>::Ones();
  example2.data_.lb = -100 * Eigen::Matrix<real_t, 3, 1>::Ones();
  example2.data_.ub = 100 * Eigen::Matrix<real_t, 3, 1>::Ones();
  example2.setup();
  example2.solve();
  std::cout << "nlc::QPOasesEigen" << std::endl;
  std::cout << example2.getOptimizer().transpose() << std::endl;

  nlc::QPSwiftEigen ex_swft(3, 8, 0);
  ex_swft.H = example2.data_.H;

  //     std::cout << "H: \n"
  //               << H << "\nf: \n"
  //               << f.transpose() << "\n A: \n"
  //               << A << "\n b: " << b.transpose() << "\n Aeq: \n"
  //               << Aeq << "\n beq: " << beq.transpose() << "\n lb: \n"
  //               << lb.transpose() << "\n ub: " << ub.transpose() <<
  //               std::endl;

  //     // Eigen::Map<Eigen::Matrix<nlc_real, Eigen::Dynamic,
  //     // Eigen::Dynamic>>(&nlc::clf3D::params.Q, H.rows(), H.cols()) = H;
  //     // nlc::clf3D::params.Q = (float) H.data();
  //     //   memcpy(&nlc::clf3D::params.Q, H.data(), sizeof(H.data()));

  //     // double Q[9];
  //     double *p1 = &nlc::clf3D::params.Q[0];
  //     Eigen::Map<Eigen::Matrix<double, 3, 3> >(&nlc::clf3D::params.Q[0],3,3)
  //     = H.cast<double>(); Eigen::Map<Eigen::Matrix<double, 3, 1>
  //     >(&nlc::clf3D::params.c[0],3,1) = f.cast<double>();
  //     Eigen::Map<Eigen::Matrix<double, 1, 3> >(&nlc::clf3D::params.A[0],1,3)
  //     = A.cast<double>(); H(2,2) = 4;

  //     // Eigen::Map<Eigen::Matrix<double, 3, 3>
  //     >(&nlc::clf3D::params.Q[0],3,3) = b.cast<double>();

  //     for (int i = 0; i < 3; i++) {
  //         for (int j = 0; j < 3; j++) {
  //             // nlc::clf3D::params.Q[i + j * 3] = H(i, j);  // Q[i+j*3] =
  //             H(i,j);
  //         }
  //         // nlc::clf3D::params.c[i] = 0;
  //         // nlc::clf3D::params.A[i] = A(i);
  //         nlc::clf3D::params.xlb[i] = -100;
  //         nlc::clf3D::params.xub[i] = 100;
  //     }
  //     nlc::clf3D::params.b[0] = -b(0);

  //     for (int i = 0; i < 3; i++) {
  //         for (int j = 0; j < 3; j++) {
  //             printf("Q[%d, %d] = %f ", i, j, nlc::clf3D::params.Q[i + j *
  //             3]);
  //         }
  //         printf("c[%d]=%f ", i, nlc::clf3D::params.c[i]);
  //         printf("A[%d]=%f ", i, nlc::clf3D::params.A[i]);
  //         printf("xlb[%d]=%f ", i, nlc::clf3D::params.xlb[i]);
  //         printf("xub[%d]=%f\n", i, nlc::clf3D::params.xub[i]);
  //     }
  //     printf("b=%f\n", nlc::clf3D::params.b[0]);

  //   nlc::clf3D::set_defaults();
  //   nlc::clf3D::setup_indexing();
  // //   nlc::clf3D::load_default_data();

  //     nlc::clf3D::settings.verbose = 1;
  //     nlc::clf3D::vars.x[0] = 0.0;
  //     nlc::clf3D::vars.x[1] = 0.0;
  //     nlc::clf3D::vars.x[2] = 0.0;

  //     long num_iters = nlc::clf3D::solve();

  //     for (int i = 0; i < 3; i++) {
  //         std::cout << nlc::clf3D::vars.x[i] << " x[" << i << "]" <<
  //         std::endl;
  //     }

  //     double Q[9];
  //     double *p = &Q[0];
  //     Eigen::Map<Eigen::Matrix<double, 3, 3> >(p,3,3) = H.cast<double>();

  //     std::cout << "+++++++++++++++++++++++++++" << std::endl;
  //     std::cout << "Q" << std::endl;
  //     for (int i = 0; i < 9; ++i) {
  //         std::cout << "params.Q " << nlc::clf3D::params.Q[i] << " <---> Q "
  //         << *(Q+i) << std::endl;
  //     }
  //     std::cout << "\n";
}
