#ifndef NONLINEAR_CONTROLS_LINEAR_MPC_H
#define NONLINEAR_CONTROLS_LINEAR_MPC_H

#include "common/qpoases_eigen.hpp"
#include "deque"
#include "vector"

namespace nonlinear_controls {
//////////////////////////////////////////
/// Eigen Dynamic Matrices and Vectors
//////////////////////////////////////////
using MatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using VectorXd = Eigen::Matrix<double, Eigen::Dynamic, 1>;

//////////////////////////////////////////
/// Linear Time Invariant Dynamics MPC
//////////////////////////////////////////
class LinearMPC {
protected:

  qpOASES::SQProblem *solver;
  qpOASES::int_t nWSR = 10;
  qpOASES::Options options;
  qpOASES::real_t cpu_time_limit{};
  QPOasesData data_;
  qpOASES::real_t *Hc{}, *Ac{}, *gc{}, *lbc{}, *ubc{}, *lbAc{}, *ubAc{};
  bool solver_initialized{false};

  int N, nx, nu;
  int nVars, nCons;

  bool QP_INITIALIZED = false;

  // Using the UC Berkeley, ME231A MPC with substitution notation
  MatrixXd Q, P, R; // mpc gains
  MatrixXd A, B;    // dynamics
  MatrixXd Sx, Su, Qbar, Rbar, H, F;
  MatrixXd Ulb, Uub, Xlb, Xub; // bounds

  bool verbose = false;

public:
  LinearMPC(const int &_N, const int &_nx, const int &_nu);
  ~LinearMPC();

  // main functions
  virtual void construct(); // TODO: define an init function
  virtual VectorXd init(const VectorXd &x0);
  virtual VectorXd run(const VectorXd &x0);
  virtual VectorXd run(const VectorXd &x0, const MatrixXd &A,
                       const MatrixXd &B);
  // dynamics
  virtual void init_dynamics(MatrixXd A, MatrixXd B);
  virtual void update_dynamics(MatrixXd A, MatrixXd B);

  void updateCArrays(const VectorXd &x0);

  void print();
  // setters
  void set_max_iter(const int MAX_ITER) { nWSR = MAX_ITER; }
  void set_mpc_gains(MatrixXd Q, MatrixXd P, MatrixXd R);
  void set_input_bounds(VectorXd lb, VectorXd ub);
  void set_state_bounds(VectorXd lb, VectorXd ub);
  void set_cpu_time_limit(const qpOASES::real_t t) {
    cpu_time_limit = t;
  }
};

//////////////////////////////////////////
/// Linear Time Varying Dynamics MPC
//////////////////////////////////////////
class LinearMPCt : public LinearMPC {
protected:
  std::deque<MatrixXd> Astorage, Bstorage;

public:
  LinearMPCt(const int &_N, const int &_nx, const int &_nu);
  ~LinearMPCt();

  // main functions
  virtual void construct();
  virtual VectorXd run(const VectorXd &x0, const MatrixXd &A,
                       const MatrixXd &B);

  // dynamics
  virtual void init_dynamics(MatrixXd A, MatrixXd B);
  virtual void update_dynamics(MatrixXd A, MatrixXd B);
};

} // namespace nonlinear_controls
#endif // NONLINEAR_CONTROLS_LINEAR_MPC_H
