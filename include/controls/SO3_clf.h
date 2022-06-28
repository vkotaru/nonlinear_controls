#ifndef NONLINEAR_CONTROLS_SO3_CLF_H
#define NONLINEAR_CONTROLS_SO3_CLF_H

#include <eigen3/Eigen/Dense>
#include <iostream>
#include <qpOASES.hpp>

#include "data_types/data_types.hpp"
#include "geometric_control.h"
#include "SO3_control.h"
#include "linear_mpc.h"
#include "common/qpoases_eigen.hpp"
#include "common/utils.hpp"

namespace nonlinear_controls {

class SO3Clf : public SO3Controller {
protected:
  Eigen::Matrix3d inertia_, inertia_scaled_;
  double min_eigval_inertia_;
  Eigen::Matrix3d dR, dRc;
  Eigen::Vector3d eR, deR, eOmega;

  VectorXd err_state, uOpt;
  TSO3 state_des;

  /// CLF
  double eta2 = 100.0;
  double epsilon2 = 4.0;
  double c2 = 20.0;
  double V2, LfV2;
  Eigen::Matrix<double, 1, 3> LgV2;

  /// Quadratic Programming
  qpOASES::real_t H[16], H_new[16];
  qpOASES::real_t g[4], A[4], g_new[4], A_new[4];
  qpOASES::real_t lb[4], ub[4], lb_new[4], ub_new[4];
  qpOASES::real_t lbA[1], ubA[1], lbA_new[1], ubA_new[1];
  qpOASES::real_t xOpt[4];

  qpOASES::SQProblem solver{4, 1};
  qpOASES::int_t nWSR = 10, nWSR_new;
  qpOASES::Options options;
  qpOASES::real_t cpu_time_limit;

  void print_qp_setup();
  void print_qp2_setup();
  bool pause = false;

public:
  //    SO3Clf();
  ~SO3Clf();
  SO3Clf(const Eigen::Matrix3d &J);

  virtual void init() override;
  virtual void run(double dt) {
    std::cout << "This is not the correct way to call this function" << std::endl;
  }
  void run(double dt, TSO3 x, TSO3 xd, Eigen::Vector3d &u);
};

}

#endif // NONLINEAR_CONTROLS_SO3_CLF_H
