#include "controls/SO3_control.h"

namespace nonlinear_controls {

SO3Controller::SO3Controller() {
  gains_.set_kp((Eigen::Vector3d() << 1000.0, 1000, 1000).finished());
  gains_.set_kd((Eigen::Vector3d() << 100.0, 100, 100).finished());
}

SO3Controller::~SO3Controller() = default;

void SO3Controller::init() {
}

void SO3Controller::run(double dt) {
}

void SO3Controller::run(double dt, TSO3 xd, Eigen::Vector3d &u) {
  this->run(dt, state_, xd, u);
}

void SO3Controller::run(double dt, TSO3 x, TSO3 xd, Eigen::Vector3d &u) {
  Eigen::Matrix<double, 6, 1> errors = x - xd;
  Eigen::Vector3d eR = errors.head(3);
  Eigen::Vector3d eOm = errors.tail(3);  // TODO: remove the temporary variables

  u = -gains_.kp().cwiseProduct(eR) - gains_.kd().cwiseProduct(eOm);
  u += x.Omega.cross(inertia_ * x.Omega);
  u += -inertia_ * (x.Omega.cross(x.R.transpose() * xd.R * xd.Omega) - x.R.transpose() * xd.R * xd.dOmega);
}

}  // namespace nonlinear_controls
