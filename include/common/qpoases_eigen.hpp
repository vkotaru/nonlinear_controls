#ifndef NONLINEAR_CONTROLS_COMMON_QPOASES_EIGEN_HPP
#define NONLINEAR_CONTROLS_COMMON_QPOASES_EIGEN_HPP

#include <eigen3/Eigen/Dense>
#include <iostream>
#include <qpOASES.hpp>
#include <sstream>

#include "common/log.hpp"

namespace nonlinear_controls {

struct QPOasesData {
  Eigen::Matrix<qpOASES::real_t, Eigen::Dynamic, Eigen::Dynamic> H, A;
  Eigen::Matrix<qpOASES::real_t, Eigen::Dynamic, Eigen::Dynamic> g, lb, ub, lbA, ubA;
  qpOASES::int_t nWSR{};
};

class QPOasesEigen {
protected:
  qpOASES::QProblem problem;
  int nvars, ncons;

  qpOASES::real_t *H{}, *A{}, *g{}, *lb{}, *ub{}, *lbA{}, *ubA{};
  qpOASES::int_t nWSR{};

  void eigen_to_array() {
    H = data_.H.transpose().data();
    A = data_.A.transpose().data();
    g = data_.g.transpose().data();
    lb = data_.lb.transpose().data();
    ub = data_.ub.transpose().data();
    lbA = data_.lbA.transpose().data();
    ubA = data_.ubA.transpose().data();
    nWSR = data_.nWSR;
  }

public:
  QPOasesEigen(const int n, const int m) : problem(n, m), nvars(n), ncons(m) {
    data_.H.resize(n, n);
    data_.g.resize(n, 1);
    data_.lb.resize(n, 1);
    data_.ub.resize(n, 1);

    data_.A.resize(m, n);
    data_.lbA.resize(m, 1);
    data_.ubA.resize(m, 1);

    data_.nWSR = 10;
    reset_data();

    this->options.setToFast();
    this->options.printLevel = qpOASES::PL_LOW;
  }
  ~QPOasesEigen() = default;

  void setup() { problem.setOptions(options); }

  qpOASES::returnValue solve() {
    eigen_to_array();
    qpOASES::returnValue sol = problem.init(H, g, A, lb, ub, lbA, ubA, nWSR);
    if (sol == qpOASES::SUCCESSFUL_RETURN) {
      Logger::SUCCESS("QP solved: SUCCESSFUL_RETURN");
    } else if (sol == qpOASES::RET_MAX_NWSR_REACHED) {
      Logger::WARN("QP solver: RET_MAX_NWSR_REACHED");
    } else if (sol == qpOASES::RET_INIT_FAILED) {
      Logger::ERROR("QP solver: RET_INIT_FAILED");
    } else {
      Logger::ERROR("QP solver: Unknown return value");
    }
    return sol;
  }
  void hotstart() {
    eigen_to_array();
    problem.hotstart(g, lb, ub, lbA, ubA, nWSR);
  }
  Eigen::Matrix<double, Eigen::Dynamic, 1> get_optimizer() {
    qpOASES::real_t x_opt[nvars];
    problem.getPrimalSolution(x_opt);
    Eigen::Matrix<double, Eigen::Dynamic, 1> x_opt_vec;
    x_opt_vec = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>>(x_opt, nvars, 1);
    return x_opt_vec;
  }

  QPOasesData data_;
  qpOASES::Options options;

  void print() {
    std::ostringstream oss;
    oss << "--------------------------------------------------\n"
        << "*           QP setup       *\n"
        << "--------------------------------------------------\n"
        << "H:\n"
        << data_.H << "\nf:\n"
        << data_.g.transpose() << "\n\nlbA:\n"
        << data_.lbA.transpose() << "\nA:\n"
        << data_.A << "\nubA:\n"
        << data_.ubA.transpose() << "\n\nlb: " << data_.lb.transpose()
        << "\nub: " << data_.ub.transpose()
        << "\n-------------------*****---------------------------";
    Logger::INFO(oss.str());
  }

  void reset_data() {
    data_.H.setIdentity();
    data_.A.setZero();
    data_.g.setZero();
    data_.lb.setZero();
    data_.ub.setZero();
    data_.lbA.setZero();
    data_.ubA.setZero();
  }
};

}  // namespace nonlinear_controls

#endif  // NONLINEAR_CONTROLS_COMMON_QPOASES_EIGEN_HPP
