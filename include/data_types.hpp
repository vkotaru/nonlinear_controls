#ifndef NLC_DATA_TYPES_H
#define NLC_DATA_TYPES_H

#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/MatrixFunctions>

#include "manifolds.hpp"

namespace nonlinear_control {
using namespace manifolds;

template <typename T>
class TSO3 {
   private:
   public:
    TSO3(/* args */) {
        R.setIdentity();
        Omega.setZero();
        dOmega.setZero();
    }
    ~TSO3() {}

    template <typename OtherDerived>
    TSO3 &operator=(const TSO3<OtherDerived> &other) {
        this->R = other.R;
        this->Omega = other.Omega;
        this->dOmega = other.dOmega;
    }

    template <typename OtherDerived>
    Eigen::Matrix<T, 6, 1> &operator-(const TSO3<OtherDerived> &other) {
        return (Eigen::Matrix<T, 6, 1>() << R.error(other.R), Omega - R.transpose() * other.R * other.Omega).finished();
    }

    SO3<T> R;
    Eigen::Matrix<T, 3, 1> Omega;
    Eigen::Matrix<T, 3, 1> dOmega;
};

template <typename T>
class Gains {
   private:
    Eigen::Matrix<T, 3, 1> kp_;
    Eigen::Matrix<T, 3, 1> kd_;
    Eigen::Matrix<T, 3, 1> ki_;

   public:
    Gains(void) {}
    ~Gains(void) {}

    inline const Eigen::Matrix<T, 3, 1> kp() const { return kp_; }
    inline const Eigen::Matrix<T, 3, 1> ki() const { return ki_; }
    inline const Eigen::Matrix<T, 3, 1> kd() const { return kd_; }

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

template <typename T>
class Wrench {
   public:
    Wrench(/* args */);
    ~Wrench();

    Eigen::Matrix<T, 3, 1> force;
    Eigen::Matrix<T, 4, 1> torque;

    void reset() {
        force.setZero();
        torque.setZero();
    }

    const Eigen::Matrix<T, 6, 1> operator()() {
        return (Eigen::Matrix<T, 6, 1>() << force, torque).finished();
    }
};

}  // namespace nonlinear_control
#endif  // NLC_DATA_TYPES_H