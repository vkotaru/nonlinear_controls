#ifndef NLC_SE3_HPP
#define NLC_SE3_HPP

#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/MatrixFunctions>

#include "manifolds.hpp"
#include "SO3.hpp"

namespace nonlinear_controls {
using namespace manifolds;

class TSE3 : public TSO3 {
public:
  TSE3() : TSO3() {}
  ~TSE3() = default;

  Eigen::Vector3d position{};
  Eigen::Vector3d velocity{};
  Eigen::Vector3d acceleration{}; // feed-forward

  TSE3 &operator=(const TSE3 &other) {
    this->position = other.position;
    this->velocity = other.velocity;
    this->acceleration = other.acceleration;
    this->R = other.R;
    this->Omega = other.Omega;
    this->dOmega = other.dOmega;
    return *this;
  }

  void print() const override {
    std::cout << "position: " << this->position.transpose() << std::endl;
    std::cout << "velocity: " << this->velocity.transpose() << std::endl;
    TSO3::print();
  }

  Eigen::Matrix<double, 12, 1> error(const TSE3 &other) {
//    std::cout << "this " << std::endl;
//    this->print();
//    std::cout << "other" << std::endl;
//    other.print();

    Eigen::Vector3d pos_err = this->position - other.position;
    Eigen::Vector3d vel_err = this->velocity - other.velocity;
    Eigen::Vector3d rot_err = this->R.error(other.R);
    Eigen::Vector3d ang_vel_err =
        this->Omega - this->R.transpose() * other.R * other.Omega;

    return (Eigen::Matrix<double, 12, 1>() << pos_err, vel_err, rot_err, ang_vel_err).finished();
  }

  Eigen::Matrix<double, 12, 1> operator-(const TSE3 &other) {
    return this->error(other);
  }

  TSO3 extractTSO3() {
    TSO3 att;
    att.R = this->R;
    att.Omega = this->Omega;
    return att;
  }
};

} // namespace nonlinear_controls

#endif // NLC_SE3_HPP