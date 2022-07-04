#ifndef _NONLINEAR_CONTROLS_LINEAR_MPC_HPP_
#define _NONLINEAR_CONTROLS_LINEAR_MPC_HPP_

#include "common/qpoases_eigen.hpp"
#include "common/constants.h"
#include <deque>
#include <vector>
#include <stdexcept>
#include <optional>

namespace nonlinear_controls {
//////////////////////////////////////////
/// Eigen Dynamic Matrices and Vectors
//////////////////////////////////////////
using MatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using VectorXd = Eigen::Matrix<double, Eigen::Dynamic, 1>;

#define INFY 1.e6
//////////////////////////////////////////
/// Linear Time Invariant Dynamics MPC
//////////////////////////////////////////
class LinearMPCBase {
public:
  struct Bounds {
    explicit Bounds(const long N) {
      lb.resize(N);
      ub.resize(N);

      lb.setOnes();
      ub.setOnes();

      lb = -INFY * lb;
      ub = INFY * ub;
    }
    VectorXd lb;
    VectorXd ub;
  };
protected:
  long N, nx, nu;
  long nVars{}, nCons{};

  bool qp_initialized{false};
  bool solver_initialized{false};
  double cpu_time_limit{0};
  long MaxIter{0};

  /// stage state cost
  MatrixXd Q;
  /// terminal state cost
  MatrixXd P;
  /// stage input cost
  MatrixXd R;
  /// stage delta input cost
  MatrixXd R2;
  /// dynamics state matrix
  MatrixXd A;
  /// dynamics input mapping
  MatrixXd B;

  Bounds state_bnds_;
  Bounds input_bnds_;

  bool verbose = false;

  /// input state
  VectorXd x0;
  VectorXd xref;

public:
  LinearMPCBase(const long &N, const long &nx, const long &nu)
      : N(N), nx(nx), nu(nu), state_bnds_(nx), input_bnds_(nu) {


    /// Sizing all Eigen Dynamic Matrices
    A.resize(nx, nx);
    B.resize(nx, nu);
    x0.resize(nx);
    xref.resize(nx);

    Q.resize(nx, nx);
    Q.setIdentity();
    R.resize(nu, nu);
    R.setIdentity();
    P.resize(nx, nx);
    P.setIdentity();
  }
  ~LinearMPCBase() = default;

  // setters
  virtual void set_gains(const MatrixXd &q, const MatrixXd &p, const MatrixXd &r) {
    assert(q.rows() == nx && q.cols() == nx);
    assert(p.rows() == nx && p.cols() == nx);
    assert(r.rows() == nu && r.cols() == nu);
    this->Q = q;
    this->P = p;
    this->R = r;
  }
  virtual void set_input_bounds(const VectorXd &lb, const VectorXd &ub) {
    assert(lb.rows() == nu && lb.cols() == 1);
    assert(ub.rows() == nu && ub.cols() == 1);
    /// set only once
    input_bnds_.lb = lb;
    input_bnds_.ub = ub;
  }
  virtual void set_state_bounds(const VectorXd &lb, const VectorXd &ub) {
    assert(lb.rows() == nx && lb.cols() == 1);
    assert(ub.rows() == nx && ub.cols() == 1);
    /// set only once
    state_bnds_.lb = lb;
    state_bnds_.ub = ub;
  }
  virtual void init_dynamics(const MatrixXd &a, const MatrixXd &b) {
    update_dynamics(a, b);
  }
  virtual void update_dynamics(const MatrixXd &a, const MatrixXd &b) {
    this->A = a;
    this->B = b;
  }

  // main functions
  virtual void construct() {
    throw std::invalid_argument("LinearMPCBase::construct not implemented");
  }
  virtual void init(const VectorXd &_x0) {
    throw std::invalid_argument("LinearMPCBase::init not implemented");
  }

  virtual std::optional<MatrixXd> run(const VectorXd &_x0) {
    throw std::invalid_argument("LinearMPCBase::run not implemented");
  }
  virtual std::optional<MatrixXd> run(const VectorXd &_x0, const VectorXd &xd_) {
    throw std::invalid_argument("LinearMPCBase::run not implemented");
  }
  virtual std::optional<MatrixXd> run(const VectorXd &_x0, const MatrixXd &a,
                                      const MatrixXd &b) {
    throw std::invalid_argument("LinearMPCBase::updateCArrays not implemented");
  }

  virtual void updateCArrays(const VectorXd &_x0) {
    throw std::invalid_argument("LinearMPCBase::updateCArrays not implemented");
  }
  virtual void print() {
    throw std::invalid_argument("LinearMPCBase::print not implemented");
  }

  // setters
  void set_max_iter(const int max_iter) { this->MaxIter = max_iter; }
  void set_cpu_time_limit(const qpOASES::real_t t) {
    cpu_time_limit = t;
  }
  virtual Eigen::MatrixXd X() {
    throw std::invalid_argument("LinearMPCBase::print not implemented");
  }
  virtual void reset() {
    throw std::invalid_argument("LinearMPCBase::reset not implemented");
  }
  bool valid_state(const VectorXd &_x0) {
    if ((state_bnds_.lb - _x0).maxCoeff() < 0 && (_x0 - state_bnds_.ub).maxCoeff() < 0)
      return true;
    else
      return false;
  }
  bool valid_input(const VectorXd &_u0) {
    if ((input_bnds_.lb - _u0).maxCoeff() < 0 && (_u0 - input_bnds_.ub).maxCoeff() < 0)
      return true;
    else
      return false;
  }
};

class LinearMPC : public LinearMPCBase {
protected:

  // Using the UC Berkeley, ME231A MPC with substitution notation
  MatrixXd Sx, Su, Qbar, Rbar, H, F;
  VectorXd Xlb, Xub, Ulb, Uub;
  VectorXd Uarray;

public:
  // MPC with state substitution in the dynamics
  // only N inputs are the optimization variables
  LinearMPC(const long &N, const long &nx, const long &nu) : LinearMPCBase(N, nx, nu) {
    nVars = N * nu;
    nCons = (N + 1) * nx;

    Sx.resize(nx * (N + 1), nx);
    Sx.setZero();
    Su.resize(nx * (N + 1), N * nu);
    Su.setZero();

    Qbar.resize(nx * (N + 1), nx * (N + 1));
    Qbar.setIdentity();

    Rbar.resize(nu * N, nu * N);
    Rbar.setIdentity();

    H.resize(nVars, nVars);
    H.setIdentity();
    F.resize(nx, nVars);
    F.setZero();

    for (int i = 0; i < N; ++i) {
      this->Qbar.block(nx * i, nx * i, nx, nx) = Q;
      this->Rbar.block(nu * i, nu * i, nu, nu) = R;
    }
    this->Qbar.block(nx * N, nx * N, nx, nx) = P;

    Ulb.resize(nu * N);
    Ulb.setOnes();
    Ulb = -INFY * Ulb;
    Uub.resize(nu * N);
    Uub.setOnes();
    Uub = INFY * Uub;
    Xlb.resize(nx * (N + 1));
    Xlb.setOnes();
    Xlb = -INFY * Xlb;
    Xub.resize(nx * (N + 1));
    Xub.setOnes();
    Xub = INFY * Xub;

    Uarray.resize(nu * N);
    Uarray.setZero();
  }

  void set_input_bounds(const VectorXd &lb, const VectorXd &ub) override {
    LinearMPCBase::set_input_bounds(lb, ub);
    Ulb = lb.replicate(N, 1);
    Uub = ub.replicate(N, 1);
  }
  void set_state_bounds(const VectorXd &lb, const VectorXd &ub) override {
    LinearMPCBase::set_state_bounds(lb, ub);
    Xlb = lb.replicate(N + 1, 1);
    Xub = ub.replicate(N + 1, 1);
  }
  void set_gains(const MatrixXd &q, const MatrixXd &p, const MatrixXd &r) override {
    LinearMPCBase::set_gains(q, p, r);
    for (int i = 0; i < N; ++i) {
      this->Qbar.block(nx * i, nx * i, nx, nx) = Q;
      this->Rbar.block(nu * i, nu * i, nu, nu) = R;
    }
    this->Qbar.block(nx * N, nx * N, nx, nx) = P;
  }

  void construct() override {
    ////////////////////////////////////////
    /// Note: remember to call this function
    /// if dynamics, input bounds or
    /// MPC gains are updated
    ////////////////////////////////////////
    Sx.block(0, 0, nx, nx).setIdentity();
    Su.block(0, 0, nx, N * nu).setZero();

    for (int i = 0; i < N; ++i) {
      Sx.block(nx * (i + 1), 0, nx, nx) =
          Sx.block(nx * i, 0, nx, nx) * A;
      if (i == 0)
        Su.block(nx * (i + 1), 0, nx, nu) = B;
      else {
        Su.block(nx * (i + 1), 0, nx, nu) = A * Su.block(nx * i, 0, nx, nu);
        Su.block(nx * (i + 1), nu, nx, nu * (N - 1)) = Su.block(nx * i, 0, nx, nu * (N - 1));
      }
    }
    // Cost function
    H = Su.transpose() * Qbar * Su + Rbar;
    F = Sx.transpose() * Qbar * Su;
  }

  void print() override {
    LinearMPCBase::print();
  }

  Eigen::MatrixXd X() override {
    // use this function only after the solve function is implemented
    Eigen::MatrixXd x_traj;
    x_traj.resize(nx, N);
    x_traj.setZero();
    Eigen::VectorXd x = x0;
    for (int i = 0; i < N; ++i) {
      x = A * x + B * Uarray.block(nu * i, 0, nu, 1);
      x_traj.col(i) = x;
    }
    return x_traj;
  }

};

//////////////////////////////////////////
/// Linear Time Varying Dynamics MPC
//////////////////////////////////////////
//class LinearMPCt : public LinearMPC {
//protected:
//  std::deque<MatrixXd> Astorage, Bstorage;
//
//public:
//  LinearMPCt(const int &_N, const int &_nx, const int &_nu);
//  ~LinearMPCt();
//
//  // main functions
//  virtual void construct();
//  virtual VectorXd run(const VectorXd &x0, const MatrixXd &A,
//                       const MatrixXd &B);
//
//  // dynamics
//  virtual void init_dynamics(MatrixXd A, MatrixXd B);
//  virtual void update_dynamics(MatrixXd A, MatrixXd B);
//};

} // namespace nonlinear_controls
#endif //_NONLINEAR_CONTROLS_LINEAR_MPC_HPP_
