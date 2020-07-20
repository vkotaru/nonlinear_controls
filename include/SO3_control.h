#ifndef NLC_SO3_CONTROL_H
#define NLC_SO3_CONTROL_H

#include "data_types.hpp"
#include "geometric_control.h"

namespace nonlinear_control {

template <typename T>
class SO3Controller : public GeometricController<T> {
   private:
    TSO3<T> state_;
    Gains<T> gains_;

    Eigen::Matrix<T, 3, 3> inertia_;

   public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    SO3Controller(/* args */);
    ~SO3Controller();

    inline const TSO3<T>& state() const { return state_; }
    virtual void init();

    virtual void run(T dt);
    void run(T dt, TSO3<T> xd, Eigen::Matrix<T, 3, 1>& u);
    void run(T dt, TSO3<T> x, TSO3<T> xd, Eigen::Matrix<T, 3, 1>& u);
};

}  // namespace nonlinear_control

#endif  // NLC_SO3_CONTROL_H