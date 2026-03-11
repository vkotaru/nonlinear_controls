#include "controls/SO3_clf.h"

#include <sstream>

namespace nonlinear_controls {

SO3Clf::SO3Clf(const Eigen::Matrix3d& J) : SO3Controller(), inertia_(J) {
  Eigen::EigenSolver<Eigen::Matrix3d> s(this->inertia_);
  min_eigval_inertia_ = s.eigenvalues().real().minCoeff();
  inertia_scaled_ = this->inertia_ * (1 / min_eigval_inertia_);
  // 0.00489992 is min eigenvalue of inertia matrix

  /// QP setup
  //    solver = new qpOASES::SQProblem(4, 1);
  options.setToMPC();
  options.printLevel = qpOASES::PL_LOW;
  solver.setOptions(options);

  for (int i = 0; i < 16; i++) {
    H[i] = 0.0;
    H_new[i] = 0.0;
  }
  for (int j = 0; j < 4; j++) {
    H[4 * j + j] = 1.0;
    A[j] = 0.0;
    g[j] = 0.0;
    lb[j] = -1000;
    ub[j] = 1000;

    H_new[4 * j + j] = 1.0;
    A_new[j] = 0.0;
    g_new[j] = 0.0;
    lb_new[j] = -1000;
    ub_new[j] = 1000;
  }
  H[15] = 4e2;
  A[3] = -1;
  lbA[0] = -1. * INFINITY;
  ubA[0] = INFINITY;
  nWSR = 10;
  H_new[15] = 4e2;
  A_new[3] = -1;
  lbA_new[0] = -1. * INFINITY;
  ubA_new[0] = INFINITY;
  nWSR_new = 10;
}

SO3Clf::~SO3Clf() = default;

void SO3Clf::print_qp_setup() {
  std::ostringstream oss;
  oss << "--------------------------------------------------\n"
      << "*           QP setup       *\n"
      << "--------------------------------------------------\n"
      << "Cost function:\n";
  for (int i = 0; i < 4; i++) {
    oss << "[";
    for (int j = 0; j < 4; j++) {
      oss << H[4 * i + j] << "\t";
    }
    oss << "]\t[" << g[i] << "]\n";
  }
  oss << "Bounds:\n";
  for (int i = 0; i < 4; i++) {
    oss << "[" << lb[i] << "]\t<x" << (i + 1) << "<\t[" << ub[i] << "]\n";
  }
  oss << "Constraints:\n[" << lbA[0] << "]<[";
  for (int i = 0; i < 4; i++) {
    oss << "\t" << A[i];
  }
  oss << "]x<[" << ubA[0] << "]\n"
      << "-------------------*****---------------------------";
  Logger::INFO(oss.str());
}

void SO3Clf::print_qp2_setup() {
  std::ostringstream oss;
  oss << "--------------------------------------------------\n"
      << "*           QP_new setup       *\n"
      << "--------------------------------------------------\n"
      << "Cost function:\n";
  for (int i = 0; i < 4; i++) {
    oss << "[";
    for (int j = 0; j < 4; j++) {
      oss << H_new[4 * i + j] << "\t";
    }
    oss << "]\t[" << g_new[i] << "]\n";
  }
  oss << "Bounds:\n";
  for (int i = 0; i < 4; i++) {
    oss << "[" << lb_new[i] << "]\t<x" << (i + 1) << "<\t[" << ub_new[i] << "]\n";
  }
  oss << "Constraints:\n[" << lbA_new[0] << "]<[";
  for (int i = 0; i < 4; i++) {
    oss << "\t" << A_new[i];
  }
  oss << "]x<[" << ubA_new[0] << "]\n"
      << "-------------------*****---------------------------";
  Logger::INFO(oss.str());
}

void SO3Clf::init() {
  std::ostringstream oss;
  oss << "/////////////////////////////////////////\n"
      << "///////////////// init //////////////////\n"
      << "inertia\n"
      << this->inertia_;
  Logger::INFO(oss.str());
  solver.init(H, g, A, lb, ub, lbA, ubA, nWSR);
  print_qp_setup();
  Logger::INFO("/////////////////////////////////////////");
}

void SO3Clf::run(double dt, TSO3 x, TSO3 xd, Eigen::Vector3d& u) {
  /// computing errors
  Eigen::Matrix<double, 6, 1> error = x - xd;
  eR = error.block(0, 0, 3, 1);
  eOmega = error.block(3, 0, 3, 1);

  dR = x.R * utils::hat(x.Omega);
  dRc = xd.R * utils::hat(xd.Omega);
  Eigen::Matrix3d m1, m2;
  m1 = (dRc.transpose() * x.R - x.R.transpose() * dRc);
  m2 = (xd.R.transpose() * dR - dR.transpose() * xd.R);
  deR = 0.5 * (utils::vee(m1) + utils::vee(m2));

  /// setting up the CLF-QP
  V2 = (eOmega.transpose() * inertia_scaled_ * eOmega * 0.5 + epsilon2 * (eR.transpose() * eOmega) +
        c2 * (eR.transpose() * eR) * 0.5)
           .value();
  LgV2 = (eOmega.transpose() * inertia_scaled_ + epsilon2 * eR.transpose());
  LfV2 = ((epsilon2 * eOmega.transpose() + c2 * eR.transpose()) * deR -
          LgV2 * (dR.transpose() * xd.R * xd.Omega + x.R.transpose() * xd.R * xd.d_omega))
             .value();

  for (int i = 0; i < 3; i++) {
    A_new[i] = LgV2(0, i);
  }
  ubA_new[0] = -LfV2 - eta2 * V2;

  print_qp2_setup();
  /// solving the QP
  qpOASES::int_t nwsr = 1000;
  qpOASES::real_t cpu_time = 1 / 500;  // << modify this >>
  qpOASES::returnValue sol_info =
      solver.hotstart(H_new, g_new, A_new, lb_new, ub_new, lbA_new, ubA_new, nwsr);
  if (sol_info == qpOASES::SUCCESSFUL_RETURN) {
    this->pause = false;
    Logger::SUCCESS("Optimal solution found");
  } else {
    this->pause = true;
    Logger::ERROR("Optimal solution NOT found");
  }
  solver.getPrimalSolution(x_opt);
  Eigen::Vector3d d_omega;
  d_omega << x_opt[0], x_opt[1], x_opt[2];
  std::ostringstream oss;
  oss << "x_opt = [ " << x_opt[0] << ", " << x_opt[1] << ", " << x_opt[2]
      << " ];  objVal = " << solver.getObjVal();
  Logger::INFO(oss.str());

  /// computing the input
  u = this->inertia_ * d_omega + utils::hat(x.Omega) * this->inertia_ * x.Omega;
}

}  // namespace nonlinear_controls
