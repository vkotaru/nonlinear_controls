#ifndef NLC_SE3_CONTROL_H
#define NLC_SE3_CONTROL_H

#include "SO3_control.h"
#include "data_types/data_types.hpp"
#include "linear_mpc.h"
#include <iostream>

namespace nonlinear_controls {
class SE3Controller : public SO3Controller {

protected:
  double mass_{};

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  SE3Controller(/* args */);
  ~SE3Controller();

  void init(const Eigen::Matrix3d &J, const double m) {
    this->inertia_ = J;
    mass_ = m;
  }

  Gains pgains_;
  void run(double dt, TSE3 x, TSE3 xd, Wrench &u);
};

} // namespace nonlinear_controls
#endif // NLC_SE3_CONTROL_H
