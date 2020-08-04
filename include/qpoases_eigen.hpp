#ifndef QPOASES_EIGEN_H
#define QPOASES_EIGEN_H

#include <eigen3/Eigen/Dense>
#include <iostream>
#include <qpOASES.hpp>

namespace nonlinear_control {

struct QPOasesData {
    Eigen::Matrix<qpOASES::real_t, Eigen::Dynamic, Eigen::Dynamic> H, A;
    Eigen::Matrix<qpOASES::real_t, Eigen::Dynamic, Eigen::Dynamic> g, lb, ub, lbA, ubA;
    qpOASES::int_t nWSR;
};

class QPOasesEigen {
   protected:
    qpOASES::QProblem problem;
    qpOASES::Options options;
    int nvars, ncons;

    qpOASES::real_t *H, *A, *g, *lb, *ub, *lbA, *ubA;
    qpOASES::int_t nWSR; 

    void eigen2array() {
        H       = data_.H.transpose().data();
        A       = data_.A.transpose().data();
        g       = data_.g.transpose().data();
        lb      = data_.lb.transpose().data();
        ub      = data_.ub.transpose().data();
        lbA     = data_.lbA.transpose().data();
        ubA     = data_.ubA.transpose().data();
        nWSR    = data_.nWSR;
    }

   public:
    QPOasesEigen(const int n, const int m) : problem(n, m), nvars(n), ncons(m) {
        data_.H.resize(n,n);
        data_.A.resize(m,n);
        data_.g.resize(n,1);
        data_.lb.resize(n,1);
        data_.ub.resize(n,1);
        data_.lbA.resize(n,1);
        data_.ubA.resize(n,1);
        data_.nWSR = 10;
    }
    ~QPOasesEigen() {
    }
    void setup() {
        problem.setOptions(options);
    }
    void solve() {
        eigen2array();
        problem.init(H, g, A, lb, ub, lbA, ubA, nWSR);
    }
    Eigen::Matrix<double, Eigen::Dynamic, 1> getOptimizer() {
        qpOASES::real_t xOpt[nvars];
        problem.getPrimalSolution(xOpt);
        Eigen::Matrix<double, Eigen::Dynamic, 1> xOptVec;
        xOptVec = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>>(xOpt, nvars, 1);
        return xOptVec;
    }

    QPOasesData data_;
};

}  // namespace nonlinear_control

#endif  // QPOASES_EIGEN_H