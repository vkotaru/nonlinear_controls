#ifndef NLC_SE3_HPP
#define NLC_SE3_HPP


#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/MatrixFunctions>

#include "manifolds.hpp"
#include "SO3.hpp"

namespace nonlinear_controls {
using namespace manifolds;

template <typename T> class TSE3 : public TSO3<T> {
public:
  TSE3(/* args */) : TSO3<T>() {
    position.setZero();
    velocity.setZero();
    acceleration.setZero();
  }
  ~TSE3() = default;

  Eigen::Matrix<T, 3, 1> position;
  Eigen::Matrix<T, 3, 1> velocity;
  Eigen::Matrix<T, 3, 1> acceleration; // feed-forward

  template <typename OtherDerived>
  TSE3 &operator=(const TSE3<OtherDerived> &other) {
    this->position = other.position;
    this->velocity = other.velocity;
    this->acceleration = other.acceleration;
    this->R = other.R;
    this->Omega = other.Omega;
    this->dOmega = other.dOmega;
  }

  void print() const override {
    std::cout << "position: " << this->position.transpose() << std::endl;
    std::cout << "velocity: " << this->velocity.transpose() << std::endl;
    TSO3<T>::print();
  }

  template <typename OtherDerived>
  Eigen::Matrix<T, 12, 1> error(const TSE3<OtherDerived> &other) {
    std::cout << "this " << std::endl;
    this->print();

    std::cout << "other" << std::endl;
    other.print();

    auto pos_err = this->position - other.position;
    auto vel_err = this->velocity - other.velocity;
    auto rot_err = this->R.error(other.R);
    auto ang_vel_err =
        this->Omega - this->R.transpose() * other.R * other.Omega;

    Eigen::Matrix<T, 12, 1> err_;
    err_ << pos_err, vel_err, rot_err, ang_vel_err;
    return err_;
  }

  template <typename OtherDerived>
  Eigen::Matrix<T, 12, 1> operator-(const TSE3<OtherDerived> &other) {
    return (Eigen::Matrix<T, 12, 1>() << position - other.position,
            velocity - other.velocity, this->R.error(other.R),
            this->Omega - this->R.transpose() * other.R * other.Omega)
        .finished();
  }

  TSO3<T> extractTSO3() {
    TSO3<T> att;
    att.R = this->R;
    att.Omega = this->Omega;
    return att;
  }
};


typedef TSE3<float> TSE3f;
typedef TSE3<double> TSE3d;

} // namespace nonlinear_controls

#endif // NLC_SE3_HPP