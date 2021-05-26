//
// Created by kotaru on 5/26/21.
//

#ifndef NONLINEAR_CONTROLS_TRACKING_LMPC_H
#define NONLINEAR_CONTROLS_TRACKING_LMPC_H

#include "deque"
#include "quadprog/qpswift_eigen.h"
#include "vector"

namespace nonlinear_controls {

class TrackingLinearMPC {
protected:
public:
  TrackingLinearMPC();
  ~TrackingLinearMPC();
};

} // namespace nonlinear_controls

#endif // NONLINEAR_CONTROLS_TRACKING_LMPC_H
