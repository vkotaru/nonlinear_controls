//
// Created by kotaru on 5/26/21.
//
#include "controls/tracking_mpc.h"

namespace nonlinear_controls {

TrackingMPC::TrackingMPC(const int &_N, const int &_nx, const int &_nu)
    : N(_N), nx(_nx), nu(_nu) {

  nv = N * nu;
  ni = N * nx * 2 /* state-bounds */;
  ne = 0;

  A.resize(nx, nx);
  B.resize(nx, nu);

  Sx.resize(nx * N, nx);
  Sx.setZero();
  Su.resize(nx * N, N * nu);
  Su.setZero();
  Qbar.resize(nx * N, nx * N);
  Qbar.setIdentity();
  Rbar.resize(nu * N, nu * N);
  Rbar.setIdentity();

  Q.resize(nx, nx);
  Q.setIdentity();
  R.resize(nu, nu);
  R.setIdentity();
  P.resize(nx, nx);
  P.setIdentity();

  H.resize(nv, nv);
  H.setIdentity();
  F.resize(nx, nv);
  F.setZero();

  Ulb.resize(nu * N, 1);
  Ulb.setOnes();
  Ulb = -1e6 * Ulb;
  Uub.resize(nu * N, 1);
  Uub.setOnes();
  Uub = 1e6 * Uub;
  Xlb.resize(nx * N, 1);
  Xlb.setOnes();
  Xlb = -1e6 * Xlb;
  Xub.resize(nx * N, 1);
  Xub.setOnes();
  Xub = 1e6 * Xub;
}

TrackingMPC::~TrackingMPC() = default;

void TrackingMPC::set_mpc_gains(const MatrixXd &_Q, const MatrixXd &_P,
                                const MatrixXd &_R) {
  assert(_Q.rows() == nx && _Q.cols() == nx);
  assert(_P.rows() == nx && _P.cols() == nx);
  assert(_R.rows() == nu && _R.cols() == nu);
  this->Q = _Q;
  this->P = _P;
  this->R = _R;
  for (int i = 0; i < N; ++i) {
    this->Qbar.block(nx * i, nx * i, nx, nx) = Q;
    if (i == N - 1) {
      this->Qbar.block(nx * i, nx * i, nx, nx) = P;
    }
    this->Rbar.block(nu * i, nu * i, nu, nu) = R;
  }
}

void TrackingMPC::set_input_bounds(const VectorXd &_ulb, const VectorXd &_uub) {
  assert(_ulb.rows() == nu && _ulb.cols() == 1);
  assert(_uub.rows() == nu && _uub.cols() == 1);
  Ulb = _ulb.replicate(N, 1);
  Uub = _uub.replicate(N, 1);
}

void TrackingMPC::set_state_bounds(const VectorXd &_xlb, const VectorXd &_xub) {
  assert(_xlb.rows() == nx && _xlb.cols() == 1);
  assert(_xub.rows() == nx && _xub.cols() == 1);
  Xlb = _xlb.replicate(N, 1);
  Xub = _xub.replicate(N, 1);
}

void TrackingMPC::set_lti_dynamics(const MatrixXd &_A, const MatrixXd &_B) {
  assert(_A.rows() == nx && _A.cols() == nx);
  assert(_B.rows() == nx && _B.cols() == nu);
  this->A = _A;
  this->B = _B;
}

void TrackingMPC::construct() {
  ////////////////////////////////////////
  /// Note: remember to call this function
  /// if dynamics, input bounds or
  /// MPC gains are updated
  ////////////////////////////////////////
  for (int i = 0; i < N; ++i) {
    if (i == 0) {
      Sx.block(nx * i, 0, nx, nx) = A;
      Su.block(nx * i, 0, nx, nu) = B;
    } else {
      Sx.block(nx * i, 0, nx, nx) = A * Sx.block(nx * (i - 1), 0, nx, nx);

      Su.block(nx * i, 0, nx, nu) = A * Su.block(nx * (i-1), 0, nx, nu);
      Su.block(nx * i, nu, nx, nu * (N - 1)) = Su.block(nx * (i - 1), 0, nx, nu * (N - 1));

/*      Su.block(nx * i, 0, nx, nu) = A * Su.block(nx * (i - 1), 0, nx, nu);
      Su.block(nx * i, nu, nx, nu * (N - 1)) =
          Su.block(nx * (i - 1), 0, nx, nu * (N - 1)); */

    }
  }
  // Cost function
  H = Su.transpose() * Qbar * Su + Rbar;
  F = Sx.transpose() * Qbar * Su;
}

void TrackingMPC::run(const VectorXd &x0, const VectorXd &Xd,
                      const VectorXd &Ud) {
  assert(x0.rows() == nx && x0.cols() == 1);
  assert(Xd.rows() == nx * N && Xd.cols() == 1);
  assert(Ud.rows() == nu * N && Ud.cols() == 1);
}

} // namespace nonlinear_controls
