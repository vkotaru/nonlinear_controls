#ifndef NLC_DATA_TYPES_HPP
#define NLC_DATA_TYPES_HPP

#include <utility>

#include "S2.hpp"
#include "SE3.hpp"
#include "SO3.hpp"
#include "manifolds.hpp"

namespace nonlinear_controls {
using namespace manifolds;

/**
 * @brief PID controller gains container.
 *
 * Stores proportional (kp), derivative (kd), and integral (ki) gains
 * as 3D vectors for use in geometric controllers.
 */
class Gains {
private:
  Eigen::Vector3d kp_;  ///< Proportional gains
  Eigen::Vector3d kd_;  ///< Derivative gains
  Eigen::Vector3d ki_;  ///< Integral gains

public:
  Gains() = default;
  ~Gains() = default;

  /// @brief Get proportional gains
  Eigen::Vector3d kp() const { return kp_; }
  /// @brief Get integral gains
  Eigen::Vector3d ki() const { return ki_; }
  /// @brief Get derivative gains
  Eigen::Vector3d kd() const { return kd_; }

  /// @brief Set proportional gains
  void set_kp(Eigen::Vector3d _kp) { kp_ = std::move(_kp); }
  /// @brief Set integral gains
  void set_ki(Eigen::Vector3d _ki) { ki_ = std::move(_ki); }
  /// @brief Set derivative gains
  void set_kd(Eigen::Vector3d _kd) { kd_ = std::move(_kd); }
};

/**
 * @brief 6-DOF wrench (force + torque) container.
 *
 * Represents a spatial force consisting of a 3D force vector
 * and a 3D torque vector, commonly used in rigid body dynamics.
 */
class Wrench {
public:
  Wrench() = default;
  ~Wrench() = default;

  Eigen::Vector3d force{};   ///< Linear force vector [N]
  Eigen::Vector3d torque{};  ///< Torque/moment vector [N·m]

  /// @brief Reset force and torque to zero
  void reset() {
    force.setZero();
    torque.setZero();
  }

  /// @brief Get wrench as 6x1 vector [force; torque]
  Eigen::Matrix<double, 6, 1> operator()() const {
    return (Eigen::Matrix<double, 6, 1>() << force, torque).finished();
  }
};

}  // namespace nonlinear_controls
#endif  // NLC_DATA_TYPES_HPP
