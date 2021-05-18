#include "controls/SE3_control.h"

namespace nonlinear_controls {
template <typename T>
SE3Controller<T>::SE3Controller() : SO3Controller<T>() {
}

template <typename T>
SE3Controller<T>::~SE3Controller() = default;

template <typename T>
void SE3Controller<T>::run(T dt, TSE3<T> x, TSE3<T> xd, Wrench<T>& u) {
    auto errors = x - xd;
    auto ex = errors.head(3);  // TODO: remove temporary variables
    auto ev = errors.segment(3, 3);
    auto eR = errors.segment(6, 3);
    auto eOm = errors.tail(3);
    u.force = -pgains_.kp().cwiseProduct(ex) - pgains_.kd().cwiseProduct(ev);
    u.force += mass_ * (Eigen::Matrix<T, 3, 1>() << 0.0, 0.0, 9.81).finished();

    u.torque = -this->gains_.kp().cwiseProduct(eR) - this->gains_.kd().cwiseProduct(eOm);
    u.torque += x.Omega.cross(this->inertia_ * x.Omega);
    u.torque += -this->inertia_ * (x.Omega.cross(x.R.transpose() * xd.R * xd.Omega) - x.R.transpose() * xd.R * xd.dOmega);
}

template class SE3Controller<float>;
template class SE3Controller<double>;

}  // namespace nonlinear_controls
