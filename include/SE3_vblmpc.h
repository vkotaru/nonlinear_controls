#ifndef NONLINEAR_CONTROL_SE3_VBLMPC_H
#define NONLINEAR_CONTROL_SE3_VBLMPC_H

#include "data_types.hpp"
#include "double_int_mpc.hpp"
#include "geometric_control.h"
#include "SO3_vblmpc.h"
#include<iostream>
namespace nonlinear_control {

template <typename T>
class SE3VblMPC: public GeometricController<T> {

  protected:
    T mass_;
    Eigen::Matrix<T, 3, 3> inertia_;
    bool IS_LTI;
    int N;
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    SE3VblMPC(bool islti, int N, T dt, const T m, const Eigen::Matrix<T, 3, 3>J);
    ~SE3VblMPC();

    DoubleIntMPC<T, 3> pos_mpc_;
    SO3VblMPC<T> att_mpc_;

    void run(T dt, TSE3<T> x, TSE3<T> xd, Wrench<T>& u);
};
} // namespace nonlinear_control
#endif // NONLINEAR_CONTROL_SE3_VBLMPC_H
