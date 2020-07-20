#include <eigen3/Eigen/Dense>
#include <iostream>

#include "SO3_control.h"
using namespace nonlinear_control;
using namespace manifolds;

int main() {
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
    R1 << 0.8073, -0.5176, 0.2835, 0.5229, 0.8501, 0.0630, -0.2736, 0.0974,
        0.9569;
    R2 << 0.4569, -0.0903, 0.8849, 0.7532, 0.5684, -0.3309, -0.4731, 0.8178,
        0.3278;

    std::cout << R1.error(R2) << std::endl;
    std::cout << SO3d::error(R1, R2) << std::endl;
    std::cout << R1.config_error(R2) << std::endl;
    std::cout << SO3d::config_error(R1, R2) << std::endl;

    std::cout << "***********************" << std::endl;
    SO3Controller<float> geo_;

    std::cout << geo_.state().R << "\n"
              << geo_.state().Omega << std::endl;
}
