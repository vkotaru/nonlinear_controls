#ifndef NONLINEAR_CONTROLS_DOUBLE_INT_MPC_HPP
#define NONLINEAR_CONTROLS_DOUBLE_INT_MPC_HPP

#include "linear_mpc.h"
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <qpOASES.hpp>

namespace nonlinear_controls {

template<unsigned int D> class DoubleIntMPC {
protected:
  const double g{G_SI};
  int N, nx, nu;
  double dt;

  LinearMPC *mpcSolver;
  VectorXd err_state, uOpt;

  MatrixXd A, B;                                   // Dynamics
  MatrixXd Q, P, R;                                // gains
  VectorXd state_lb, state_ub, input_lb, input_ub; // bounds

public:
  DoubleIntMPC(const int N, const double _dt) : N(N), dt(_dt) {
    nx = 2 * D;
    nu = D;
    mpcSolver = new LinearMPC(N, nx, nu);
    /////////////////////////////////////////////////
    /// setting up discrete-translational dynamics
    /// x_{k+1} = x_{k} + dt*v_{k} + 0.5*dt*dt*a_{k}
    /// v_{k+1} = v_{k} + dt*a_{k}
    /////////////////////////////////////////////////
    A.resize(nx, nx);
    B.resize(nx, nu);
    A.setIdentity();
    A.topRightCorner(3, 3) += dt * Eigen::Matrix<double, D, D>::Identity();
    B << 0.5 * dt * dt * Eigen::Matrix<double, D, D>::Identity(),
        dt * Eigen::Matrix<double, D, D>::Identity();
    mpcSolver->init_dynamics(A, B);

    /// gains
    Q.resize(nx, nx);
    P.resize(nx, nx);
    R.resize(nu, nu);
    Q << 1000 * Eigen::Matrix<double, D, D>::Identity(),
        Eigen::Matrix<double, D, D>::Zero(), Eigen::Matrix<double, D, D>::Zero(),
        100 * Eigen::Matrix<double, D, D>::Identity();
    P << 8.1314, 0.0000, -0.0000, 0.6327, -0.0000, -0.0000, 0.0000, 8.1314,
        0.0000, 0.0000, 0.6327, 0.0000, -0.0000, 0.0000, 8.1314, -0.0000,
        0.0000, 0.6327, 0.6327, 0.0000, -0.0000, 0.2606, -0.0000, -0.0000,
        -0.0000, 0.6327, 0.0000, -0.0000, 0.2606, 0.0000, -0.0000, 0.0000,
        0.6327, -0.0000, 0.0000, 0.2606;
    P = 1e4 * P;
    R = 1 * Eigen::Matrix<double, D, D>::Identity();
    mpcSolver->set_mpc_gains(Q, P, R);

    /// bounds
    input_lb.resize(nu, 1);
    input_ub.resize(nu, 1);
    state_lb.resize(nx, 1);
    state_ub.resize(nx, 1);
    input_lb << -g, -g, -g;
    input_ub << g, g, g;
    state_lb.setOnes();
    state_lb = -100 * state_lb;
    state_ub = -state_lb;
    mpcSolver->set_input_bounds(input_lb, input_ub);
    mpcSolver->set_state_bounds(state_lb, state_ub);

    uOpt.resize(nu * N, 1);

    ///
    mpcSolver->construct();
  }

  ~DoubleIntMPC() = default;

  void set_gains(const MatrixXd &_Q, const MatrixXd &_P,
                 const MatrixXd &_R) {
    this->Q = _Q;
    this->P = _P;
    this->R = _R;

    std::cout << "SO3VblMPC: setting gains\n Q:\n" << std::endl;
    std::cout << _Q << std::endl;
    std::cout << "P:\n";
    std::cout << _P << std::endl;
    std::cout << "R:\n" << std::endl;
    std::cout << _R << std::endl;

    mpcSolver->set_mpc_gains(_Q, _P, _R);
  }

  void set_input_bounds(VectorXd lb, VectorXd ub) {
    this->input_lb = lb;
    this->input_ub = ub;
    mpcSolver->set_input_bounds(lb, ub);
  }

  void set_state_bounds(VectorXd lb, VectorXd ub) {
    this->state_lb = lb;
    this->state_ub = ub;
    mpcSolver->set_state_bounds(lb, ub);
  }

  void construct() {
    /// call this function after any update in the MPC parameters
    mpcSolver->construct();
  }

  VectorXd run(VectorXd _err_state) {
    uOpt = mpcSolver->run(_err_state);
    return uOpt.block(0, 0, nu, 1);
  }

  VectorXd updateState(VectorXd state, VectorXd input) {
    return this->A * state + this->B * input;
  }
};

} // namespace nonlinear_controls
#endif // NONLINEAR_CONTROLS_DOUBLE_INT_MPC_HPP
