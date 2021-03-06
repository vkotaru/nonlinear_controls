#ifndef NONLINEAR_CONTROLS_LINEAR_MPC_H
#define NONLINEAR_CONTROLS_LINEAR_MPC_H
#include "common/qpoases_eigen.hpp"
#include "deque"
#include "vector"
namespace nonlinear_controls {
//////////////////////////////////////////
/// Eigen Dynamic Matrices and Vectors
//////////////////////////////////////////
template<typename T>
using MatrixX = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
template<typename T> using VectorX = Eigen::Matrix<T, Eigen::Dynamic, 1>;

typedef MatrixX<double> MatrixXd;
typedef MatrixX<float> MatrixXf;
typedef VectorX<double> VectorXd;
typedef VectorX<float> VectorXf;

//////////////////////////////////////////
/// Linear Time Invariant Dynamics MPC
//////////////////////////////////////////
template<typename T> class LinearMPC {
protected:

  qpOASES::SQProblem *solver;
  qpOASES::int_t nWSR = 10;
  qpOASES::Options options;
  qpOASES::real_t cpu_time_limit{};
  QPOasesData data_;
  qpOASES::real_t *Hc{}, *Ac{}, *gc{}, *lbc{}, *ubc{}, *lbAc{}, *ubAc{};

  int N, nx, nu;
  int nVars, nCons;

  bool QP_INITIALIZED = false;

  // Using the UC Berkeley, ME231A MPC with substitution notation
  MatrixX<T> Q, P, R; // mpc gains
  MatrixX<T> A, B;    // dynamics
  MatrixX<T> Sx, Su, Qbar, Rbar, H, F;
  MatrixX<T> Ulb, Uub, Xlb, Xub; // bounds

  bool verbose = false;

public:
  LinearMPC(const int &_N, const int &_nx, const int &_nu);
  ~LinearMPC();

  // main functions
  virtual void construct(); // TODO: define an init function
  virtual VectorX<T> init(const VectorX<T> &x0);
  virtual VectorX<T> run(const VectorX<T> &x0);
  virtual VectorX<T> run(const VectorX<T> &x0, const MatrixX<T> &A,
                         const MatrixX<T> &B);
  // dynamics
  virtual void init_dynamics(MatrixX<T> A, MatrixX<T> B);
  virtual void update_dynamics(MatrixX<T> A, MatrixX<T> B);

  void updateCArrays(const VectorX<T> &x0);

  void print();
  // setters
  void set_max_iter(const int MAX_ITER) { nWSR = MAX_ITER; }
  void set_mpc_gains(MatrixX<T> Q, MatrixX<T> P, MatrixX<T> R);
  void set_input_bounds(VectorX<T> lb, VectorX<T> ub);
  void set_state_bounds(VectorX<T> lb, VectorX<T> ub);
};

//////////////////////////////////////////
/// Linear Time Varying Dynamics MPC
//////////////////////////////////////////
template<typename T> class LinearMPCt : public LinearMPC<T> {
protected:
  std::deque<MatrixX<T>> Astorage, Bstorage;

public:
  LinearMPCt(const int &_N, const int &_nx, const int &_nu);
  ~LinearMPCt();

  // main functions
  virtual void construct() override;
  virtual VectorX<T> run(const VectorX<T> &x0, const MatrixX<T> &A,
                         const MatrixX<T> &B) override;

  // dynamics
  virtual void init_dynamics(MatrixX<T> A, MatrixX<T> B) override;
  virtual void update_dynamics(MatrixX<T> A, MatrixX<T> B) override;
};

} // namespace nonlinear_controls
#endif // NONLINEAR_CONTROLS_LINEAR_MPC_H
