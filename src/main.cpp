#include <eigen3/Eigen/Dense>
#include <iostream>

class MyVectorType : public Eigen::VectorXd {
  public:
    MyVectorType(void) : Eigen::VectorXd() {}
    template <typename OtherDerived>
    MyVectorType(const Eigen::MatrixBase<OtherDerived>& other) : Eigen::VectorXd(other) {
    }
    template <typename OtherDerived>
    MyVectorType& operator=(const Eigen::MatrixBase<OtherDerived>& other) {
        this->Eigen::VectorXd::operator=(other);
        return *this;
    }
};

template <typename T, int _Dim>
class SpecialOrthogonal : public Eigen::Matrix<T, _Dim, _Dim> {
  public:
    SpecialOrthogonal(void) : Eigen::Matrix<T, _Dim, _Dim>() {}
    template <typename OtherDerived>
    SpecialOrthogonal(const Eigen::MatrixBase<OtherDerived>& other) : Eigen::Matrix<T, _Dim, _Dim>(other) {
    }
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
    SO3(void) : SpecialOrthogonal<T, 3>() {
    }
    template <typename OtherDerived>
    SO3(const Eigen::MatrixBase<OtherDerived>& other) : SpecialOrthogonal<T, 3>(other) {
    }
    template <typename OtherDerived>
    SO3& operator=(const Eigen::MatrixBase<OtherDerived>& other) {
        this->SpecialOrthogonal<T, 3>::operator=(other);
        return *this;
    }

    template <typename T1, typename T2>
    static Eigen::Matrix<T, 3, 1> error(const Eigen::MatrixBase<T1>& R, const Eigen::MatrixBase<T2>& Rd) {
        // eR = 0.5*vee(Rd'*R-R'*Rd);
        auto eR = 0.5 * (Rd.transpose() * R - R.transpose() * Rd);
        return (Eigen::Matrix<T, 3, 1>() << eR(2, 1), eR(0, 2), eR(1, 0)).finished();
    }
    template <typename T1, typename T2>
    static T config_error(const Eigen::MatrixBase<T1>& R, const Eigen::MatrixBase<T2>& Rd) {
        // eR =  0.5*trace(eye(3)-Rd'*R);
        T psi = 0.5 * (Eigen::Matrix<T, 3, 3>::Identity() - Rd.transpose() * R).trace();
        return psi;
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

int main() {
    MyVectorType v = MyVectorType::Ones(4);
    v(2) += 10;
    v = 2 * v;
    std::cout << v.transpose() << std::endl;

    SpecialOrthogonal<float, 4> R = SpecialOrthogonal<float, 4>::Identity();
    std::cout << R << "\n"
              << R.inverse() << std::endl;

    Eigen::Matrix3d r2;
    r2 << 3, 5, 7, 2, 5, 7, 6, 5, 4;
    Eigen::Matrix3d b = r2.inverse();
    std::cout << r2 << "\n"
              << r2.transpose() << "\n"
              << b << std::endl;

    Eigen::VectorXf v2;

    SO3<float> Rtest = SO3<float>::Identity();
    std::cout << Rtest << "\n"
              << Rtest.inverse() << std::endl;

    std::cout << "***********************" << std::endl;
    SO3d R1, R2;
    R1 << 0.8073, -0.5176, 0.2835, 0.5229, 0.8501, 0.0630, -0.2736, 0.0974, 0.9569;
    R2 << 0.4569, -0.0903, 0.8849,
    0.7532, 0.5684, -0.3309,
    -0.4731, 0.8178, 0.3278;

    std::cout << R1.error(R2) << std::endl;
    std::cout << SO3d::error(R1, R2) << std::endl;
    std::cout << R1.config_error(R2) << std::endl;
    std::cout << SO3d::config_error(R1, R2) << std::endl;
}
