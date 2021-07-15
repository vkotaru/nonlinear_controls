#include "common/qpoases_eigen.hpp"
#include "controls/clf_qp.h"
#include "quadprog/qpswift_eigen.h"
#include <chrono>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <qpOASES.hpp>

namespace nlc = nonlinear_controls;

#define nlc_real double

int main() {
  USING_NAMESPACE_QPOASES
  //////////////////////////////////////////////////////////////////
  nlc::QPOasesEigen example2(3, 1);
  example2.data_.H = 2 * Eigen::Matrix<real_t, 3, 3>::Identity();
  example2.data_.g = Eigen::Matrix<real_t, 3, 1>::Zero();
  example2.data_.A << 1, 3, 4;
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
  ex_swft.f = example2.data_.g;
  ex_swft.c = 0;
  ex_swft.A << Eigen::Matrix3d::Identity(), -Eigen::Matrix3d::Identity(),
      example2.data_.A, -example2.data_.A;
  ex_swft.b << example2.data_.ub, -example2.data_.lb, example2.data_.ubA,
      -example2.data_.lbA;
  std::cout << "Ax<=b\nA:\n" << ex_swft.A << "\nb:\n" << ex_swft.b << "\n";
  ex_swft.solve();
  //////////////////////////////////////////////////////////////////
  printf("...............................\n"
         "test case 3\n");
  printf("QPOases: \n");
  /* Setup data of first QP. */
  real_t H[2 * 2] = {1.0, 0.0, 0.0, 0.5};
  real_t A[1 * 2] = {1.0, 1.0};
  real_t g[2] = {1.5, 1.0};
  real_t lb[2] = {0.5, -2.0};
  real_t ub[2] = {5.0, 2.0};
  real_t lbA[1] = {-1.0};
  real_t ubA[1] = {2.0};
  /* Setting up QProblem object. */
  QProblem example(2, 1);
  Options options;
  options.printLevel = PL_LOW;
  example.setOptions(options);
  /* Solve first QP. */
  int_t nWSR = 10;
  example.init(H, g, A, lb, ub, lbA, ubA, nWSR);
  /* Get and print solution of first QP. */
  real_t xOpt[2];
  real_t yOpt[2 + 1];
  example.getPrimalSolution(xOpt);
  example.getDualSolution(yOpt);
  printf("\nxOpt = [ %e, %e ];\n  yOpt = [ %e, %e, %e ];  objVal = %e\n\n",
         xOpt[0], xOpt[1], yOpt[0], yOpt[1], yOpt[2], example.getObjVal());

  nlc::QPOasesEigen ex3a(2, 1);
  ex3a.data_.H << 1.0, 0.0, 0.0, 0.5;
  ex3a.data_.g << 1.5, 1.0;
  ex3a.data_.A << 1, 1;
  ex3a.data_.lbA << -1.0;
  ex3a.data_.ubA << 2.0;
  ex3a.data_.lb << 0.5, -2.0;
  ex3a.data_.ub << 5.0, 2.0;

  ex3a.options.setToFast();
  ex3a.options.printLevel = qpOASES::PL_LOW;
  ex3a.setup();
  auto start = std::chrono::high_resolution_clock::now();
  ex3a.solve();
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  std::cout << "Time taken by qpoases: " << duration.count() << " microseconds"
            << std::endl;

  std::cout << "nlc::QPOasesEigen: xOpt: ";
  std::cout << ex3a.getOptimizer().transpose() << std::endl;

  nlc::QPSwiftEigen ex3s(2, 6, 0);
  ex3s.H = ex3a.data_.H;
  ex3s.f = ex3a.data_.g;
  ex3s.c = 0;
  ex3s.A << Eigen::Matrix2d::Identity(), -Eigen::Matrix2d::Identity(),
      ex3a.data_.A, -ex3a.data_.A;
  ex3s.b << ex3a.data_.ub, -ex3a.data_.lb, ex3a.data_.ubA, -ex3a.data_.lbA;
  std::cout << "Ax<=b\nA:\n" << ex_swft.A << "\nb:\n" << ex_swft.b << "\n";

  start = std::chrono::high_resolution_clock::now();
  ex3s.solve();
  stop = std::chrono::high_resolution_clock::now();
  duration =
      std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  std::cout << "Time taken by qpswift: " << duration.count() << " microseconds"
            << std::endl;
}
