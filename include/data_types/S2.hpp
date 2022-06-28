#ifndef NLC_S2_HPP
#define NLC_S2_HPP

#include "manifolds.hpp"

namespace nonlinear_controls {
using namespace manifolds;

class S2 : public Sn<double, 3> {
public:
  S2() : Sn<double, 3>() {}

  template<typename OtherDerived>
  explicit S2(const Eigen::MatrixBase<OtherDerived> &other) : Sn<double, 3>(other) {}

  template<typename OtherDerived>
  S2 &operator=(const Eigen::MatrixBase<OtherDerived> &other) {
    this->Sn<double, 3>::operator=(other);
    return *this;
  }

  template<typename T1, typename T2>
  static Eigen::Vector3d error(const Eigen::MatrixBase<T1> &q,
                               const Eigen::MatrixBase<T2> &qd) {
    // eR = qd x q
    return (Eigen::Vector3d() << q(1) * qd(2) - q(2) * qd(1),
        q(2) * qd(0) - q(0) * qd(2), q(0) * qd(1) - q(1) * qd(0))
        .finished();
  }

  Eigen::Vector3d error(const S2 &other) {
    return S2::error(*this, other);
  }

  template<typename T1, typename T2>
  double config_error(const Eigen::MatrixBase<T1> &q,
                      const Eigen::MatrixBase<T2> &qd) {
    // Psi =  1 - q.qd
    return (1.0 - q.dot(qd));
  }
};

class TS2 {
public:
  TS2() {
    q << 0., 0., 1;
  }

  ~TS2() = default;

  S2 q{};
  Eigen::Vector3d dq{};
  Eigen::Vector3d omega{};
  Eigen::Vector3d domega{};

  TS2 &operator=(const TS2 &other) {
    this->q = other.q;
    this->dq = other.dq;
    this->omega = other.omega;
    this->domega = other.domega;
  }

  Eigen::Matrix<double, 6, 1> error(const TS2 &other) {
    Eigen::Matrix3d transport_operator;
    transport_operator << q(1) * q(1) + q(2) * q(2), -q(0) * q(1), -q(0) * q(2),
        -q(0) * q(1), q(0) * q(0) + q(2) * q(2), -q(1) * q(2), -q(0) * q(2),
        -q(1) * q(2), q(0) * q(0) + q(1) * q(1);
    Eigen::Matrix<double, 6, 1> err_;
    err_ << this->q.error(other.q),
        this->omega - transport_operator * other.omega;
    return err_;
  }

  Eigen::Matrix<double, 6, 1> operator-(const TS2 &other) {
    return this->error(other);
  }
};

} // namespace nonlinear_controls

#endif // NLC_S2_HPP
