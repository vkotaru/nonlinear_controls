#include "controls/SE3_control.h"

namespace nonlinear_controls {
SE3Controller::SE3Controller() : SO3Controller() {
}

SE3Controller::~SE3Controller() = default;

void SE3Controller::run(double dt, TSE3 x, TSE3 xd, Wrench &u) {
  auto errors = x - xd;
  Eigen::Vector3d ex = errors.head(3);  // TODO: remove temporary variables
  Eigen::Vector3d ev = errors.segment(3, 3);
  Eigen::Vector3d eR = errors.segment(6, 3);
  Eigen::Vector3d eOm = errors.tail(3);
  u.force = -pgains_.kp().cwiseProduct(ex) - pgains_.kd().cwiseProduct(ev);
  u.force += mass_ * (Eigen::Vector3d() << 0.0, 0.0, 9.81).finished();

  u.torque = -this->gains_.kp().cwiseProduct(eR) - this->gains_.kd().cwiseProduct(eOm);
  u.torque += x.Omega.cross(this->inertia_ * x.Omega);
  u.torque += -this->inertia_ * (x.Omega.cross(x.R.transpose() * xd.R * xd.Omega) - x.R.transpose() * xd.R * xd.dOmega);
}

}  // namespace nonlinear_controls
