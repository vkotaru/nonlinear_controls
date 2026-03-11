#ifndef NLC_QP_INTERFACE_H
#define NLC_QP_INTERFACE_H

#include "eigen3/Eigen/Dense"

namespace nonlinear_controls {

/**
 * @brief Data structure for quadratic programming problems.
 *
 * Standard QP form: min 0.5*x'Hx + f'x
 *                   s.t. A*x <= b, Aeq*x = beq, xlb <= x <= xub
 */
struct QuadProgData {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> H;    ///< Hessian matrix
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> A;    ///< Inequality constraint matrix
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Aeq;  ///< Equality constraint matrix
  Eigen::Matrix<double, Eigen::Dynamic, 1> f;                 ///< Linear cost vector
  Eigen::Matrix<double, Eigen::Dynamic, 1> b;                 ///< Inequality constraint bounds
  Eigen::Matrix<double, Eigen::Dynamic, 1> beq;               ///< Equality constraint values
  Eigen::Matrix<double, Eigen::Dynamic, 1> xlb;               ///< Lower bounds on x
  Eigen::Matrix<double, Eigen::Dynamic, 1> xub;               ///< Upper bounds on x
};

/**
 * @brief Abstract interface for quadratic programming solvers.
 *
 * Base class that defines the interface for QP solver implementations.
 * Derived classes should implement setup() and solve() for specific solvers.
 */
class QPInterface {
public:
  QPInterface();
  virtual ~QPInterface();

  QuadProgData problem_;  ///< QP problem data

  /// @brief Configure the solver with problem data
  virtual void setup();
  /// @brief Solve the QP and store result
  virtual void solve();
};

}  // namespace nonlinear_controls

#endif  // NLC_QP_INTERFACE_H