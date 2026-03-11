#ifndef NLC_BASE_CONTROLLER_H
#define NLC_BASE_CONTROLLER_H

#include <eigen3/Eigen/Dense>

namespace nonlinear_controls {

/**
 * @brief Abstract base class for all controllers.
 *
 * Provides a common interface for controller initialization and execution.
 * Derived classes should override init() and run() to implement specific
 * control algorithms.
 */
class BaseController {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  BaseController() = default;
  virtual ~BaseController() = default;

  /// @brief Initialize controller state and parameters
  virtual void init() {}
  /// @brief Execute one control step
  virtual void run() {}
};

}  // namespace nonlinear_controls
#endif  // NLC_BASE_CONTROLLER_H
