#ifndef NLC_S2_HPP
#define NLC_S2_HPP

#include "manifolds.hpp"

namespace nonlinear_controls {
using namespace manifolds;

template <typename T> class S2 : public Sn<T, 3> {
public:
  S2(void) : Sn<T, 3>() {}

  template <typename OtherDerived>
  S2(const Eigen::MatrixBase<OtherDerived> &other) : Sn<T, 3>(other) {}

  template <typename OtherDerived>
  S2 &operator=(const Eigen::MatrixBase<OtherDerived> &other) {
    this->Sn<T, 3>::operator=(other);
    return *this;
  }

  template <typename T1, typename T2>
  static Eigen::Matrix<T, 3, 1> error(const Eigen::MatrixBase<T1> &q,
                                      const Eigen::MatrixBase<T2> &qd) {
    // eR = qd x q
    return (Eigen::Matrix<T, 3, 1>() << q(1) * qd(2) - q(2) * qd(1),
            q(2) * qd(0) - q(0) * qd(2), q(0) * qd(1) - q(1) * qd(0))
        .finished();
  }

  template <typename T1, typename T2>
  static T config_error(const Eigen::MatrixBase<T1> &q,
                        const Eigen::MatrixBase<T2> &qd) {
    // Psi =  1 - q.qd
    T Psi = static_cast<T>(1.0) - q.dot(qd);
    return Psi;
  }
};

template <typename T> class TS2 {

public:
  TS2() {
    q << static_cast<T>(0.), static_cast<T>(0.), static_cast<T>(1.);
    dq.setZero ();
    omega.setZero();
    domega.setZero();
  }

  ~TS2() = default;

  S2<T> q;
  Eigen::Matrix<T, 3, 1> dq;
  Eigen::Matrix<T, 3, 1> omega;
  Eigen::Matrix<T, 3, 1> domega;

  template <typename OtherDerived>
  TS2 &operator=(const TS2<OtherDerived> &other) {
    this->q = other.q;
    this->dq = other.dq;
    this->omega = other.omega;
    this->domega = other.domega;
  }

  template <typename OtherDerived>
  Eigen::Matrix<T, 6, 1> error(const TS2<OtherDerived> &other) {
    Eigen::Matrix<T, 3, 3> transport_operator;
    transport_operator << q(1) * q(1) + q(2) * q(2), -q(0) * q(1), -q(0) * q(2),
        -q(0) * q(1), q(0) * q(0) + q(2) * q(2), -q(1) * q(2), -q(0) * q(2),
        -q(1) * q(2), q(0) * q(0) + q(1) * q(1);
    Eigen::Matrix<T, 6, 1> err_;
    err_ << this->q.error(other.q),
        this->omega - transport_operator * other.omega;
    return err_;
  }

  template <typename OtherDerived>
  Eigen::Matrix<T, 6, 1> operator-(const TS2<OtherDerived> &other) {
    return this->error(other);
  }
};

typedef S2<float> S2f;
typedef S2<double> S2d;

typedef TS2<float> TS2f;
typedef TS2<double> TS2d;

} // namespace nonlinear_controls

#endif // NLC_S2_HPP
