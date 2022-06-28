#ifndef NLC_SO3_CONTROL_H
#define NLC_SO3_CONTROL_H

#include "data_types/data_types.hpp"
#include "geometric_control.h"
#include <iostream>

namespace nonlinear_controls {

class SO3Controller : public GeometricController {
protected:
  TSO3 state_;
  Eigen::Matrix3d inertia_;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  SO3Controller();
  explicit SO3Controller(const Eigen::Matrix3d &J) {
    inertia_ = J;
    gains_.set_kp((Eigen::Vector3d() << 1000.0, 1000, 1000).finished());
    gains_.set_kd((Eigen::Vector3d() << 100.0, 100, 100).finished());
  }
  ~SO3Controller();

  inline const TSO3 &state() const {
    return state_;
  }
  virtual void init();
  void init(const Eigen::Matrix3d &J) {
    inertia_ = J;
  }

  virtual void run(double dt);
  Gains gains_;
  void run(double dt, TSO3 xd, Eigen::Vector3d &u);
  void run(double dt, TSO3 x, TSO3 xd, Eigen::Vector3d &u);
};

}  // namespace nonlinear_controls

#endif  // NLC_SO3_CONTROL_H
