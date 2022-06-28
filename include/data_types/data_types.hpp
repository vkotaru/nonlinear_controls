#ifndef NLC_DATA_TYPES_HPP
#define NLC_DATA_TYPES_HPP

#include <utility>

#include "manifolds.hpp"
#include "SO3.hpp"
#include "SE3.hpp"
#include "S2.hpp"

namespace nonlinear_controls {
using namespace manifolds;

class Gains {
private:
  Eigen::Vector3d kp_;
  Eigen::Vector3d kd_;
  Eigen::Vector3d ki_;

public:
  Gains() = default;
  ~Gains() = default;

  inline Eigen::Vector3d kp() const { return kp_; }
  inline Eigen::Vector3d ki() const { return ki_; }
  inline Eigen::Vector3d kd() const { return kd_; }

  void set_kp(Eigen::Vector3d _kp) {
    kp_ = std::move(_kp);
  }

  void set_ki(Eigen::Vector3d _ki) {
    ki_ = std::move(_ki);
  }

  void set_kd(Eigen::Vector3d _kd) {
    kd_ = std::move(_kd);
  }
};

class Wrench {
public:
  Wrench() = default;
  ~Wrench() = default;

  Eigen::Vector3d force{};
  Eigen::Vector3d torque{};

  void reset() {
    force.setZero();
    torque.setZero();
  }

  Eigen::Matrix<double, 6, 1> operator()() const {
    return (Eigen::Matrix<double, 6, 1>() << force, torque).finished();
  }
};

} // namespace nonlinear_controls
#endif // NLC_DATA_TYPES_HPP
