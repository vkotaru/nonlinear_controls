#ifndef NONLINEAR_CONTROLS_MPC_BASE_H
#define NONLINEAR_CONTROLS_MPC_BASE_H
#include "qpoases_eigen.hpp"
#include "vector"
#include "deque"
namespace nonlinear_control {

template <typename T> struct MPCGains {
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Q, R, P;
};

template <typename T> struct MPCBounds {
  Eigen::Matrix<T, Eigen::Dynamic, 1> lb, ub;
};

template <typename T> struct MPCProblem {
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A, B;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> lb, ub;
  int N, nx, nu;
};

template <typename T> class MPCBase {
private:
  QPOasesEigen *qp;

  MPCProblem<T> problem;
  MPCBounds<T> state, input;
  MPCGains<T> gains;

  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Au, Ax, Af, bx, bu, bf;
  std::deque<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> Astorage;
  std::deque<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> Bstorage;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Inx, Inu, Onx, Onu;


  bool isLTI;
  const int nx, nu, N;
  int nVars, nCons;

public:
  MPCBase(const bool _isLTI, const int _N, const int _nx, const int _nu);
  ~MPCBase();

  MPCProblem<T> problem_;

  void setup();
};

} // namespace nonlinear_control

#endif // NONLINEAR_CONTROLS_MPC_BASE_H