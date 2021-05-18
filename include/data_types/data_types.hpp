#ifndef NLC_DATA_TYPES_H
#define NLC_DATA_TYPES_H

#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/MatrixFunctions>

#include "manifolds.hpp"

namespace nonlinear_controls {
using namespace manifolds;

template <typename T> class TSO3 {
public:
  TSO3(/* args */) {
    R.setIdentity();
    Omega.setZero();
    dOmega.setZero();
  }
  ~TSO3() = default;
  SO3<T> R;
  Eigen::Matrix<T, 3, 1> Omega;
  Eigen::Matrix<T, 3, 1> dOmega; // feed-forward usage

  template <typename OtherDerived>
  TSO3 &operator=(const TSO3<OtherDerived> &other) {
    this->R = other.R;
    this->Omega = other.Omega;
    this->dOmega = other.dOmega;
  }

  template <typename OtherDerived>
  Eigen::Matrix<T, 6, 1> error(const TSO3<OtherDerived> &other) {
    Eigen::Matrix<T, 6, 1> err_;
    err_ << this->R.error(other.R),
        this->Omega - this->R.transpose() * other.R * other.Omega;
    return err_;
  }

  template <typename OtherDerived>
  Eigen::Matrix<T, 6, 1> operator-(const TSO3<OtherDerived> &other) {
    Eigen::Matrix<T, 6, 1> err_;
    err_ << this->R.error(other.R),
        this->Omega - this->R.transpose() * other.R * other.Omega;
    return err_;
  }

  virtual void print() const {
    std::cout << "rotation: " << this->R << std::endl;
    std::cout << "angular velocity: " << this->Omega.transpose() << std::endl;
  }
};

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

template class Wrench<float>;
template class TSO3<float>;
template class TSE3<float>;

typedef TSO3<float> TSO3f;
typedef TSO3<double> TSO3d;

typedef TSE3<float> TSE3f;
typedef TSE3<double> TSE3d;

typedef Wrench<double> Wrenchd;
typedef Wrench<float> Wrenchf;

} // namespace nonlinear_controls
#endif // NLC_DATA_TYPES_H
