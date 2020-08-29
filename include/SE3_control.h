#ifndef NLC_SE3_CONTROL_H
#define NLC_SE3_CONTROL_H

#include "data_types.hpp"
#include "SO3_control.h"
#include "linear_mpc.h"
#include<iostream>

namespace nonlinear_control {
template <typename T>
class SE3Controller : public SO3Controller<T> {

  protected:
    T mass_;

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    SE3Controller(/* args */);
    ~SE3Controller();

    void init(const Eigen::Matrix<T, 3, 3>& J, const T& m) {
        this->inertia_ = J;
        mass_ = m;
    }

    Gains<T> pgains_;
    void run(T dt, TSE3<T> x, TSE3<T> xd, Wrench<T>& u);
};


}

#endif // SE3_CONTROL_H
