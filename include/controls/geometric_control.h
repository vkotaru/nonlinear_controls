#ifndef NLC_GEOMETRIC_CONTROL_H
#define NLC_GEOMETRIC_CONTROL_H

#include "base_control.h"
#include "data_types/data_types.hpp"

namespace nonlinear_controls {

template <typename T>
class GeometricController : public BaseController<T> {
   public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    GeometricController(/* args */);
    ~GeometricController();

    // TODO: look into static_cast to access the state from the base class itself
    // inline const TSO3<T>& state() const { return state_; }
    virtual void init();
    virtual void run();
};

}  // namespace nonlinear_controls
#endif  // NLC_GEOMETRIC_CONTROL_H