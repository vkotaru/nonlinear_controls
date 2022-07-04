//
// Created by kotaru on 5/26/21.
//

#ifndef NONLINEAR_CONTROLS_TRACKING_MPC_H
#define NONLINEAR_CONTROLS_TRACKING_MPC_H

#include "controls/mpc_qpoases.hpp"
#include "deque"
#include "quadprog/qpswift_eigen.h"
#include "quadprog/quadprog.h"
#include "vector"

namespace nonlinear_controls {

class TrackingMPC {
protected:
  /// \brief Time Horizon
  int N;
  /// \brief state-space size
  int nx;
  /// \brief input-space size
  int nu;

  /// \brief QP # variables, # equality constraints, # inequality constraints;
  int nv, ne, ni;

  /// \brief flags
  bool _VERBOSE = false;
  bool _TIME_INVARIANT = true;

  /**
   * Using the UC Berkeley, ME231A MPC with substitution notation
   */
  /// \brief MPC gains
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Q, P, R;
  /// \brief System Dynamics
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> A, B;
  /// \brief MPC formulation
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Sx, Su, Qbar, Rbar;
  /// \brief state & input bounds
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Ulb, Uub, Xlb, Xub;

  /// \brief QP formulation
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> H, F;

public:
  TrackingMPC(const int &_N, const int &_nx, const int &_nu);
  ~TrackingMPC();

  /**
   * Function to set MPC cost gains
   * @param _Q
   * @param _P
   * @param _R
   */
  void set_mpc_gains(const MatrixXd &_Q, const MatrixXd &_P,
                     const MatrixXd &_R);
  /**
   * Function to set MPC input bounds
   * @param _ulb
   * @param _uub
   */
  void set_input_bounds(const VectorXd &_ulb, const VectorXd &_uub);
  /**
   * Function o set MPC state bounds
   * @param _xlb
   * @param _xub
   */
  void set_state_bounds(const VectorXd &_xlb, const VectorXd &_xub);
  /**
   * Function to set the linear dynamics matrices
   * @param _A
   * @param _B
   */
  void set_lti_dynamics(const MatrixXd &_A, const MatrixXd &_B);

  void construct();
  void run(const VectorXd &x0);
  void run(const VectorXd &x0, const VectorXd &Xd, const VectorXd &Ud);
};

} // namespace nonlinear_controls

#endif // NONLINEAR_CONTROLS_TRACKING_MPC_H
