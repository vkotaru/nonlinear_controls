#ifndef NLC_MANIFOLDS_H
#define NLC_MANIFOLDS_H

#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <iostream>

namespace nonlinear_controls {
namespace manifolds {

template <typename T, int _Dim>
class SpecialOrthogonal : public Eigen::Matrix<T, _Dim, _Dim> {
   public:
    SpecialOrthogonal(void) : Eigen::Matrix<T, _Dim, _Dim>() {}

    template <typename OtherDerived>
    SpecialOrthogonal(const Eigen::MatrixBase<OtherDerived>& other) : Eigen::Matrix<T, _Dim, _Dim>(other) {}

    template <typename OtherDerived>
    SpecialOrthogonal& operator=(const Eigen::MatrixBase<OtherDerived>& other) {
        this->Eigen::Matrix<T, _Dim, _Dim>::operator=(other);
        return *this;
    }
    SpecialOrthogonal inverse() {
        return this->transpose();
    }
};

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

typedef SO3<float> SO3f;
typedef SO3<double> SO3d;

}  // namespace manifolds
}  // namespace nonlinear_controls
#endif  // NLC_MANIFOLDS_H
