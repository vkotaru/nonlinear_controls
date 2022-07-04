#ifndef NLC_QP_INTERFACE_H
#define NLC_QP_INTERFACE_H

#include "eigen3/Eigen/Dense"

namespace nonlinear_controls {

struct QuadProgData {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> H, A, Aeq;
  Eigen::Matrix<double, Eigen::Dynamic, 1> f, b, beq;
  Eigen::Matrix<double, Eigen::Dynamic, 1> xlb, xub;
};

class QPInterface {
protected:

public:
  QPInterface(/* args */);
  ~QPInterface();
  QuadProgData problem_;

  virtual void setup();
  virtual void solve();
};

}  // namespace nonlinear_controls

#endif  // NLC_QP_INTERFACE_H