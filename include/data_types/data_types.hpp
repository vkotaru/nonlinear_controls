#ifndef NLC_DATA_TYPES_HPP
#define NLC_DATA_TYPES_HPP

#include "manifolds.hpp"
#include "SO3.hpp"
#include "SE3.hpp"
#include "S2.hpp"

namespace nonlinear_controls {
using namespace manifolds;

template <typename T> class Gains {
private:
  Eigen::Matrix<T, 3, 1> kp_;
  Eigen::Matrix<T, 3, 1> kd_;
  Eigen::Matrix<T, 3, 1> ki_;

public:
  Gains() = default;
  ~Gains() = default;

  inline Eigen::Matrix<T, 3, 1> kp() const { return kp_; }
  inline Eigen::Matrix<T, 3, 1> ki() const { return ki_; }
  inline Eigen::Matrix<T, 3, 1> kd() const { return kd_; }

  template <typename OtherDerived>
  void set_kp(Eigen::Matrix<OtherDerived, 3, 1> _kp) {
    kp_ = _kp;
  }

  template <typename OtherDerived>
  void set_ki(Eigen::Matrix<OtherDerived, 3, 1> _ki) {
    ki_ = _ki;
  }

  template <typename OtherDerived>
  void set_kd(Eigen::Matrix<OtherDerived, 3, 1> _kd) {
    kd_ = _kd;
  }
};

template <typename T> class Wrench {
public:
  Wrench(/* args */) = default;
  ~Wrench() = default;

  Eigen::Matrix<T, 3, 1> force;
  Eigen::Matrix<T, 3, 1> torque;

  void reset() {
    force.setZero();
    torque.setZero();
  }

  Eigen::Matrix<T, 6, 1> operator()() const {
    return (Eigen::Matrix<T, 6, 1>() << force, torque).finished();
  }
};

typedef Wrench<double> Wrenchd;
typedef Wrench<float> Wrenchf;

} // namespace nonlinear_controls
#endif // NLC_DATA_TYPES_HPP
