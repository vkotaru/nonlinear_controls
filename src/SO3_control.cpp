#include "SO3_control.h"

namespace nonlinear_control {

template <typename T>
SO3Controller<T>::SO3Controller() {
}

template <typename T>
SO3Controller<T>::~SO3Controller() {
}

template <typename T>
void SO3Controller<T>::init() {
}

template <typename T>
void SO3Controller<T>::run(T dt) {
}

template <typename T>
void SO3Controller<T>::run(T dt, TSO3<T> xd, Eigen::Matrix<T, 3, 1>& u) {
    this->run(dt, state_, xd, u);
}

template <typename T>
void SO3Controller<T>::run(T dt, TSO3<T> x, TSO3<T> xd, Eigen::Matrix<T, 3, 1>& u) {
    auto errors = x - xd;
    auto eR = errors.head(3);
    auto eOm = errors.tail(3);  // TODO: remove the temporary variables

    u = -gains_.kp().cwiseProduct(eR) - gains_.kd().cwiseProduct(eOm);
    u += x.Omega.cross(inertia_ * x.Omega);
    u += -inertia_ * (x.Omega.cross(x.R.transpose() * xd.R * xd.Omega) - x.R.transpose() * xd.R * xd.dOmega);
}

template class SO3Controller<float>;
template class SO3Controller<double>;

}  // namespace nonlinear_control
