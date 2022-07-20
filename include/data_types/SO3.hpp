#ifndef NLC_SO3_HPP
#define NLC_SO3_HPP

#include "manifolds.hpp"

namespace nonlinear_controls {
using namespace manifolds;

class SO3 : public SpecialOrthogonal<double, 3> {
public:
  SO3() : SpecialOrthogonal<double, 3>() {}

  template<typename OtherDerived>
  explicit SO3(const Eigen::MatrixBase<OtherDerived> &other) : SpecialOrthogonal<double, 3>(other) {}

  template<typename OtherDerived>
  SO3 &operator=(const Eigen::MatrixBase<OtherDerived> &other) {
    this->SpecialOrthogonal<double, 3>::operator=(other);
    return *this;
  }

  template<typename T1, typename T2>
  static Eigen::Vector3d error(const Eigen::MatrixBase<T1> &R, const Eigen::MatrixBase<T2> &Rd) {
    // eR = 0.5*vee(Rd'*R-R'*Rd);
    auto eR = 0.5 * (Rd.transpose() * R - R.transpose() * Rd);
    return (Eigen::Vector3d() << eR(2, 1), eR(0, 2), eR(1, 0))
        .finished();
  }
  template<typename T1, typename T2>
  double config_error(const Eigen::MatrixBase<T1> &R, const Eigen::MatrixBase<T2> &Rd) {
    // Psi =  0.5*trace(eye(3)-Rd'*R);
    return 0.5 * (Eigen::Matrix3d::Identity() - Rd.transpose() * R).trace();
  }

  Eigen::Vector3d error(const SO3 &other) {
    return SO3::error(*this, other);
  }

  double config_error(const SO3 &other) {
    return SO3::config_error(*this, other);
  }

  static SO3 Identity() {
    // TODO this should be part of the inherited behaviour
    return SO3(Eigen::Matrix3d::Identity());
  }
};

class TSO3 {
public:
  TSO3(/* args */) {
    R.setIdentity();
    Omega.setZero();
    dOmega.setZero();
  }
  ~TSO3() = default;
  SO3 R;
  Eigen::Vector3d Omega;
  Eigen::Vector3d dOmega; // feed-forward usage

  TSO3 &operator=(const TSO3 &other) {
    this->R = other.R;
    this->Omega = other.Omega;
    this->dOmega = other.dOmega;
    return *this;
  }

  Eigen::Matrix<double, 6, 1> error(const TSO3 &other) {
    Eigen::Matrix<double, 6, 1> err_;
    err_ << this->R.error(other.R),
        this->Omega - this->R.transpose() * other.R * other.Omega;
    return err_;
  }

  Eigen::Matrix<double, 6, 1> operator-(const TSO3 &other) {
    return this->error(other);
  }

  virtual void print() const {
    // TODO overload operator<<
    std::cout << "rotation: " << this->R << std::endl;
    std::cout << "angular velocity: " << this->Omega.transpose() << std::endl;
  }
};

} // namespace nonlinear_controls

#endif // NLC_SO3_HPP