#ifndef NLC_BASE_CONTROLLER_H
#define NLC_BASE_CONTROLLER_H

#include <eigen3/Eigen/Dense>

namespace nonlinear_control {

template <typename T>
class BaseController {
   public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    BaseController();
    ~BaseController();

    virtual void init();
    virtual void run();

   private:
};

}  // namespace nonlinear_control
#endif  // NLC_BASE_CONTROLLER_H
