#ifndef QPOASES_EIGEN_H
#define QPOASES_EIGEN_H

#include <eigen3/Eigen/Dense>
#include <iostream>
#include <qpOASES.hpp>
#include <common/qpoases_eigen.hpp>

namespace nonlinear_controls {
typedef Eigen::Matrix<qpOASES::real_t, Eigen::Dynamic, Eigen::Dynamic> QPMatrix;

QPMatrix qpOases_solve(const int n, const int m, QPOasesData &data) {
    qpOASES::QProblem problem(n, m);
    qpOASES::Options options;
    qpOASES::real_t *H, *A, *g, *lb, *ub, *lbA, *ubA;
    qpOASES::int_t nWSR; 

    // data transfer
    H       = data_.H.transpose().data();
    A       = data_.A.transpose().data();
    g       = data_.g.transpose().data();
    lb      = data_.lb.transpose().data();
    ub      = data_.ub.transpose().data();
    lbA     = data_.lbA.transpose().data();
    ubA     = data_.ubA.transpose().data();
    nWSR    = data_.nWSR;

    options.printLevel = qpOASES::PL_NONE;
    problem.setOptions(options);
    problem.init(H, g, A, lb, ub, lbA, ubA, nWSR);

}

// class QPOasesEigen2 {
// protected:
//   qpOASES::QProblemB *problem;
//   int nvars, ncons;
//   bool has_constraints = false;
//   int qp_type = 1;

//   qpOASES::real_t *H, *A, *g, *lb, *ub, *lbA, *ubA;
//   qpOASES::int_t nWSR;

//   void eigen2array() {
//     H = data_.H.transpose().data();
//     g = data_.g.transpose().data();
//     lb = data_.lb.transpose().data();
//     ub = data_.ub.transpose().data();

//     if (has_constraints) {
//       A = data_.A.transpose().data();
//       lbA = data_.lbA.transpose().data();
//       ubA = data_.ubA.transpose().data();
//     }
//     nWSR = data_.nWSR;
//   }

// public:
//   QPOasesEigen2(const int n) { QPOasesEigen(n, 0); }

//   QPOasesEigen2(const int n, const int m) : nvars(n), ncons(m) {
//     data_.H.resize(n, n);
//     data_.A.resize(m, n);
//     data_.g.resize(n, 1);
//     data_.lb.resize(n, 1);
//     data_.ub.resize(n, 1);
//     data_.lbA.resize(m, 1);
//     data_.ubA.resize(m, 1);
//     data_.nWSR = 10;
//     reset_data();

//     // instantiate the qp
//     if (m == 0) {
//       has_constraints = false;
//       problem = new qpOASES::QProblemB(n);
//     } else {
//       has_constraints = true;
//       problem = new qpOASES::SQProblem(n, m);
//     }
//   }
//   ~QPOasesEigen2() { delete problem; }
//   void setup() { problem->setOptions(options); }
//   void solve() {
//     eigen2array();
//     if (has_constraints) {
//       problem->init(H, g, A, lb, ub, lbA, ubA, nWSR, 0);
//     } else {
//       problem->init(H, g, lb, ub, nWSR, 0);
//     }
//   }
//   Eigen::Matrix<double, Eigen::Dynamic, 1> getOptimizer() {
//     qpOASES::real_t xOpt[nvars];
//     problem->getPrimalSolution(xOpt);
//     Eigen::Matrix<double, Eigen::Dynamic, 1> xOptVec;
//     xOptVec = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>>(xOpt, nvars, 1);
//     return xOptVec;
//   }

//   QPOasesData data_;
//   qpOASES::Options options;

//   void print() {
//     std::cout << "--------------------------------------------------" << std::endl;
//     std::cout << "*           QP setup       *" << std::endl;
//     std::cout << "--------------------------------------------------" << std::endl;
//     std::cout << "H: \n"
//               << data_.H << "\nf: \n"
//               << data_.g.transpose() << std::endl;
//     std::cout << "\nlbA: \n"
//               << data_.lbA.transpose() << "\nA: \n"
//               << data_.A << "\nubA: \n"
//               << data_.ubA.transpose() << std::endl;
//     std::cout << "\nlb: " << data_.lb.transpose()
//               << "\nub: " << data_.ub.transpose() << std::endl;
//     std::cout << "-------------------*****---------------------------" << std::endl;
//   }
//   void reset_data() {
//     data_.H.setIdentity();
//     data_.A.setZero();
//     data_.g.setZero();
//     data_.lb.setZero();
//     data_.ub.setZero();
//     data_.lbA.setZero();
//     data_.ubA.setZero();
//   }
// };

} // namespace nonlinear_control

#endif // QPOASES_EIGEN_H