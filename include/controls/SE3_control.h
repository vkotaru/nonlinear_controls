#ifndef NLC_SE3_CONTROL_H
#define NLC_SE3_CONTROL_H

#include "SO3_control.h"
#include "data_types/data_types.hpp"
#include "linear_mpc.h"
#include <iostream>

namespace nonlinear_controls {
template <typename T> class SE3Controller : public SO3Controller<T> {

protected:
  T mass_;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  SE3Controller(/* args */);
  ~SE3Controller();

  void init(const Eigen::Matrix<T, 3, 3> &J, const T &m) {
    this->inertia_ = J;
    mass_ = m;
  }

  Gains<T> pgains_;
  void run(T dt, TSE3<T> x, TSE3<T> xd, Wrench<T> &u);
};

typedef SE3Controller<double> SE3Controllerd;
typedef SE3Controller<float> SE3Controllerf;

} // namespace nonlinear_controls

#endif // NLC_SE3_CONTROL_H
