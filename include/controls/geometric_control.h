#ifndef NLC_GEOMETRIC_CONTROL_H
#define NLC_GEOMETRIC_CONTROL_H

#include "base_control.h"
#include "data_types/data_types.hpp"

namespace nonlinear_controls {

/**
 * @brief Base class for geometric controllers on Lie groups.
 *
 * Provides foundation for controllers that operate on manifolds such as
 * SO(3) and SE(3). Geometric controllers exploit the structure of these
 * manifolds to achieve globally valid control laws.
 *
 * @see SO3Controller, SE3Controller
 */
class GeometricController : public BaseController {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  GeometricController() = default;
  ~GeometricController() override = default;

  void init() override {}
  void run() override {}
};

}  // namespace nonlinear_controls
#endif  // NLC_GEOMETRIC_CONTROL_H