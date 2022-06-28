#ifndef NLC_GEOMETRIC_CONTROL_H
#define NLC_GEOMETRIC_CONTROL_H

#include "base_control.h"
#include "data_types/data_types.hpp"

namespace nonlinear_controls {

class GeometricController : public BaseController {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  GeometricController(/* args */) = default;
  ~GeometricController() = default;

  void init() override {}
  void run() override {}
};

}  // namespace nonlinear_controls
#endif  // NLC_GEOMETRIC_CONTROL_H