#ifndef NLC_SO3_HPP
#define NLC_SO3_HPP

#include "manifolds.hpp"

namespace nonlinear_controls {
using namespace manifolds;

template <typename T>
class SO3 : public SpecialOrthogonal<T, 3> {
   public:
    SO3(void) : SpecialOrthogonal<T, 3>() {}

    template <typename OtherDerived>
    SO3(const Eigen::MatrixBase<OtherDerived>& other) : SpecialOrthogonal<T, 3>(other) {}

    template <typename OtherDerived>
    SO3& operator=(const Eigen::MatrixBase<OtherDerived>& other) {
        this->SpecialOrthogonal<T, 3>::operator=(other);
        return *this;
    }

    template <typename T1, typename T2>
    static Eigen::Matrix<T, 3, 1> error(const Eigen::MatrixBase<T1>& R, const Eigen::MatrixBase<T2>& Rd) {
        // eR = 0.5*vee(Rd'*R-R'*Rd);
        auto eR = 0.5 * (Rd.transpose() * R - R.transpose() * Rd);
        return (Eigen::Matrix<T, 3, 1>() << eR(2, 1), eR(0, 2), eR(1, 0))
            .finished();
    }
    template <typename T1, typename T2>
    static T config_error(const Eigen::MatrixBase<T1>& R, const Eigen::MatrixBase<T2>& Rd) {
        // Psi =  0.5*trace(eye(3)-Rd'*R);
        T Psi = 0.5 * (Eigen::Matrix<T, 3, 3>::Identity() - Rd.transpose() * R).trace();
        return Psi;
    }

    template <typename OtherDerived>
    Eigen::Matrix<T, 3, 1> error(const SO3<OtherDerived>& other) {
        return SO3<T>::error(*this, other);
    }
    template <typename OtherDerived>
    T config_error(const SO3<OtherDerived>& other) {
        return SO3<T>::config_error(*this, other);
    }
};


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

typedef SO3<float> SO3f;
typedef SO3<double> SO3d;

typedef TSO3<float> TSO3f;
typedef TSO3<double> TSO3d;

} // namespace nonlinear_controls

#endif // NLC_SO3_HPP