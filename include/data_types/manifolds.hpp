#ifndef NLC_MANIFOLDS_H
#define NLC_MANIFOLDS_H

#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <iostream>

namespace nonlinear_controls {
namespace manifolds {

template <typename T, size_t _Dim>
class SpecialOrthogonal : public Eigen::Matrix<T, _Dim, _Dim> {
public:
  SpecialOrthogonal(void) : Eigen::Matrix<T, _Dim, _Dim>() {}

  template <typename OtherDerived>
  SpecialOrthogonal(const Eigen::MatrixBase<OtherDerived> &other)
      : Eigen::Matrix<T, _Dim, _Dim>(other) {}

  template <typename OtherDerived>
  SpecialOrthogonal &operator=(const Eigen::MatrixBase<OtherDerived> &other) {
    this->Eigen::Matrix<T, _Dim, _Dim>::operator=(other);
    return *this;
  }
  SpecialOrthogonal inverse() { return this->transpose(); }
};

template <typename T, size_t _Dim> class Sn : public Eigen::Matrix<T, _Dim, 1> {
public:
  Sn(void) : Eigen::Matrix<T, _Dim, 1>() {}

  template <typename OtherDerived>
  Sn(const Eigen::MatrixBase<OtherDerived> &other)
      : Eigen::Matrix<T, _Dim, 1>(other) {}

  template <typename OtherDerived>
  Sn &operator=(const Eigen::MatrixBase<OtherDerived> &other) {
    this->Eigen::Matrix<T, _Dim, 1>::operator=(other);
    return *this;
  }

  Sn inverse() = delete;
};

} // namespace manifolds
} // namespace nonlinear_controls
#endif // NLC_MANIFOLDS_H
