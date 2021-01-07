#include "controls/SO3_control.h"

namespace nonlinear_controls {

template <typename T>
SO3Controller<T>::SO3Controller() {
    gains_.set_kp((Eigen::Matrix<T, 3, 1>()<< 1000.0, 1000, 1000).finished());
    gains_.set_kd((Eigen::Matrix<T, 3, 1>()<< 100.0, 100, 100).finished());
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
    Eigen::Matrix<T, 6, 1> errors = x - xd;
    Eigen::Matrix<T, 3, 1> eR = errors.head(3);
    Eigen::Matrix<T, 3, 1> eOm = errors.tail(3);  // TODO: remove the temporary variables

    u = -gains_.kp().cwiseProduct(eR) - gains_.kd().cwiseProduct(eOm);
    u += x.Omega.cross(inertia_ * x.Omega);
    u += -inertia_ * (x.Omega.cross(x.R.transpose() * xd.R * xd.Omega) - x.R.transpose() * xd.R * xd.dOmega);
}

template class SO3Controller<float>;
template class SO3Controller<double>;

}  // namespace nonlinear_controls
