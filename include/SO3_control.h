#ifndef NLC_SO3_CONTROL_H
#define NLC_SO3_CONTROL_H

#include "data_types.hpp"
#include "geometric_control.h"
#include <iostream>

namespace nonlinear_control {

template <typename T>
class SO3Controller : public GeometricController<T> {
  protected:
    TSO3<T> state_;
    Eigen::Matrix<T, 3, 3> inertia_;

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    SO3Controller(/* args */);
    SO3Controller(const Eigen::Matrix<T, 3, 3>& J) {
        inertia_ = J;
      gains_.set_kp((Eigen::Matrix<T, 3, 1>()<< 1000.0, 1000, 1000).finished());
      gains_.set_kd((Eigen::Matrix<T, 3, 1>()<< 100.0, 100, 100).finished());
    }
    ~SO3Controller();

    inline const TSO3<T>& state() const {
        return state_;
    }
    virtual void init();
    void init(const Eigen::Matrix<T, 3, 3>& J) {
        inertia_ = J;
    }

    virtual void run(T dt);
    Gains<T> gains_;
    void run(T dt, TSO3<T> xd, Eigen::Matrix<T, 3, 1>& u);
    void run(T dt, TSO3<T> x, TSO3<T> xd, Eigen::Matrix<T, 3, 1>& u);
};

}  // namespace nonlinear_control

#endif  // NLC_SO3_CONTROL_H
