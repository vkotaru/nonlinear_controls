#ifndef NONLINEAR_CONTROL_SE3_VBLMPC_H
#define NONLINEAR_CONTROL_SE3_VBLMPC_H

#include "SO3_vblmpc.h"
#include "data_types/data_types.hpp"
#include "double_int_mpc.hpp"
#include "geometric_control.h"
#include <iostream>
namespace nonlinear_controls {

class SE3VblMPC : public GeometricController {
protected:
  double mass_;
  Eigen::Matrix3d inertia_;
  bool IS_LTI;
  int N;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  SE3VblMPC(bool islti, const int N, const double dt, const double m,
            const Eigen::Matrix3d &J);
  ~SE3VblMPC();

  DoubleIntMPC<3> pos_mpc_;
  SO3VblMPC att_mpc_;

  void init();
  void run(double dt, TSE3 x, TSE3 xd, Wrench &u);
};

} // namespace nonlinear_controls
#endif // NONLINEAR_CONTROL_SE3_VBLMPC_H
