#include <eigen3/Eigen/Dense>
#include <iostream>

#include "cvxgen_interface.h"
// #include "qp_interface.h"

namespace nlc = nonlinear_control;

#define nlc_real float

int main() {
    /*************************************/
    Eigen::Matrix<nlc_real, Eigen::Dynamic, Eigen::Dynamic> H, A, Aeq;
    Eigen::Matrix<nlc_real, Eigen::Dynamic, 1> f, b, beq, lb, ub;

    H = Eigen::Matrix<nlc_real, 3, 3>::Identity();
    f = Eigen::Matrix<nlc_real, 3, 1>::Zero();

    Aeq = Eigen::Matrix<nlc_real, 1, 3>::Ones();
    beq.resize(1, 1);
    beq << 3.0;

    // ub = 100 * Eigen::Matrix<nlc_real, 3, 1>::Ones();
    // lb = -100 * Eigen::Matrix<nlc_real, 3, 1>::Ones();
    /*************************************/

    std::cout << "H: \n"
              << H << "\nf: \n"
              << f.transpose() << "\n A: \n"
              << A << "\n b: " << b.transpose() << "\n Aeq: \n"
              << Aeq << "\n beq: " << beq.transpose() << "\n lb: \n"
              << lb.transpose() << "\n ub: " << ub.transpose() << std::endl;

    // Eigen::Map<Eigen::Matrix<nlc_real, Eigen::Dynamic,
    // Eigen::Dynamic>>(&nlc::clf3D::params.Q, H.rows(), H.cols()) = H;
    // nlc::clf3D::params.Q = (float) H.data();
    //   memcpy(&nlc::clf3D::params.Q, H.data(), sizeof(H.data()));

    double Q[9];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            nlc::clf3D::params.Q[i+j*3] = H(i,j);
            // Q[i+j*3] = H(i,j);
        }
    }
    Q = reinterpret_cast<double *>( H.data());

    std::cout << "+++++++++++++++++++++++++++" << std::endl;
    std::cout << "Q" << std::endl;
    for (int i = 0; i < 9; ++i) {
        std::cout << "params.Q " << nlc::clf3D::params.Q[i] << " <---> Q " << Q[i] << std::endl;
    }
    std::cout << "\n";
}
