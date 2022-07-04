#ifndef NONLINEAR_CONTROLS_SO3_VBLMPC_H
#define NONLINEAR_CONTROLS_SO3_VBLMPC_H

#include <eigen3/Eigen/Dense>
#include <iostream>
#include <qpOASES.hpp>

#include "common/qpoases_eigen.hpp"
#include "common/utils.hpp"
#include "data_types/data_types.hpp"
#include "geometric_control.h"
#include "mpc_qpoases.hpp"
#include "SO3_control.h"

namespace nonlinear_controls {

class SO3VblMPC : public SO3Controller {
protected:
  bool IS_LTI = true; // is linear time invariant
  double dt;

  const double g{G_SI};
  int N{5}, nx{6}, nu{3};
  Eigen::Matrix3d inertia_, inertia_inv_;
  TSO3 state_;

  LinearMPC *mpcSolver;
  VectorXd err_state, uOpt;
  TSO3 state_des;

  MatrixXd A, B;                                   // Dynamics
  MatrixXd Q, P, R;                                // gains
  VectorXd state_lb, state_ub, input_lb, input_ub; // bounds

  bool verbose = true;

public:
  SO3VblMPC(bool islti, int _N, double _dt);
  SO3VblMPC(bool islti, int _N, double _dt, Eigen::Matrix3d _J);
  ~SO3VblMPC();

  void generate_dynamics(Eigen::Vector3d Omd, MatrixXd &A);
  void init_dynamics(Eigen::Vector3d Omd);
  void init_dynamics(TSO3 attd);
  void construct() { mpcSolver->construct(); }

  void init() override {}
  void run(double dt) override {}
  VectorXd run(const VectorXd &_err_state) {
    uOpt = mpcSolver->run(_err_state);
    return uOpt.block(0, 0, nu, 1);
  }
  VectorXd run(TSO3 state) {
    uOpt = mpcSolver->run(state - state_des);
    return uOpt.block(0, 0, nu, 1);
  }
  void run(double dt, TSO3 x, TSO3 xd, Eigen::Vector3d &u);
  void set_gains(const MatrixXd &_Q, const MatrixXd &_P,
                 const MatrixXd &_R);
  void set_input_bounds(VectorXd lb, VectorXd ub);
  void set_state_bounds(VectorXd lb, VectorXd ub);

  /// function to integrate dynamics
  void updateState(Eigen::Matrix<double, 6, 1> &state, const Eigen::Vector3d &input) {
    std::cout << "updateState: not implemented" << std::endl;
  }
};

} // namespace nonlinear_controls
#endif // NONLINEAR_CONTROLS_SO3_VBLMPC_H
