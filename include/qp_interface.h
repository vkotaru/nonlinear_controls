#ifndef NLC_QP_INTERFACE_H
#define NLC_QP_INTERFACE_H

#include "eigen3/Eigen/Dense"

namespace nonlinear_control {

template <typename T>
struct QuadProg {
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> H, A, Aeq;
    Eigen::Matrix<T, Eigen::Dynamic, 1> f, b, beq;
    Eigen::Matrix<T, Eigen::Dynamic, 1> lb, ub;
};

template <typename T>
class QPInterface {
   protected:
    QuadProg<T> problem_; 

   public:
    QPInterface(/* args */);
    ~QPInterface();

    virtual void setup();
    virtual void solve();
};

}  // namespace nonlinear_control

#endif  // NLC_QP_INTERFACE_H