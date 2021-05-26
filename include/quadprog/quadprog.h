//
// Created by kotaru on 5/26/21.
//

#ifndef NONLINEAR_CONTROLS_QUADPROG_H
#define NONLINEAR_CONTROLS_QUADPROG_H

#include <eigen3/Eigen/Dense>

namespace nonlinear_controls {

class QuadProg {
protected:
  /// \brief Equality tolerance value
  double equality_tol{1e-6};

  /// \brief Number of decision variables
  int nv;
  /// \brief Number of inequality constraints
  int ni;
  /// \brief Number of equality constraints
  int ne;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  QuadProg(const int &n, const int &m, const int &p);
  ~QuadProg();

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> H;
  Eigen::Matrix<double, Eigen::Dynamic, 1> f;
  double c{};
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> A;
  Eigen::Matrix<double, Eigen::Dynamic, 1> b;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Aeq;
  Eigen::Matrix<double, Eigen::Dynamic, 1> beq;
  Eigen::Matrix<double, Eigen::Dynamic, 1> lb;
  Eigen::Matrix<double, Eigen::Dynamic, 1> ub;
  Eigen::Matrix<double, Eigen::Dynamic, 1> x0;
  Eigen::Matrix<double, Eigen::Dynamic, 1> xOpt;

  virtual int solve();
};

} // namespace nonlinear_controls

#endif // NONLINEAR_CONTROLS_QUADPROG_H
