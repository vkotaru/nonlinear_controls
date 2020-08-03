#include <eigen3/Eigen/Dense>
#include <iostream>

#include "clf_qp.h"

namespace nlc = nonlinear_control;

#define nlc_real double

int main() {
    /*************************************/
    // Eigen::Matrix<nlc_real, Eigen::Dynamic, Eigen::Dynamic> H, A, Aeq;
    // Eigen::Matrix<nlc_real, Eigen::Dynamic, 1> f, b, beq, lb, ub;

    // H =2*Eigen::Matrix<nlc_real, 3, 3>::Identity();
    // f = Eigen::Matrix<nlc_real, 3, 1>::Zero();
    // H(1,0) = -1;

    // A = Eigen::Matrix<nlc_real, 1, 3>::Ones();
    // b.resize(1, 1);
    // b << 3;

    // Aeq = Eigen::Matrix<nlc_real, 1, 3>::Ones();
    // beq.resize(1, 1);
    // beq << 3.0;

    nlc::ClfQP testqp;

    testqp.problem_.H = 2*Eigen::Matrix<nlc_real, 3, 3>::Identity();
    // testqp.problem_.H(1,0) = -1;
    // testqp.problem_.H(2,2) = 4;
    testqp.problem_.f = Eigen::Matrix<nlc_real, 3, 1>::Zero();
    testqp.problem_.A = Eigen::Matrix<nlc_real, 1, 3>::Ones();
    testqp.problem_.b = -3*Eigen::Matrix<nlc_real, 1, 1>::Ones();
    testqp.problem_.xlb = -100*Eigen::Matrix<nlc_real, 3, 1>::Ones();
    testqp.problem_.xub = 100*Eigen::Matrix<nlc_real, 3, 1>::Ones();

    testqp.setup();
    testqp.solve();
    std::cout << testqp.getOptimizer() << std::endl;

    // ub = 100 * Eigen::Matrix<nlc_real, 3, 1>::Ones();
    // lb = -100 * Eigen::Matrix<nlc_real, 3, 1>::Ones();
    /*************************************/

//     std::cout << "H: \n"
//               << H << "\nf: \n"
//               << f.transpose() << "\n A: \n"
//               << A << "\n b: " << b.transpose() << "\n Aeq: \n"
//               << Aeq << "\n beq: " << beq.transpose() << "\n lb: \n"
//               << lb.transpose() << "\n ub: " << ub.transpose() << std::endl;

//     // Eigen::Map<Eigen::Matrix<nlc_real, Eigen::Dynamic,
//     // Eigen::Dynamic>>(&nlc::clf3D::params.Q, H.rows(), H.cols()) = H;
//     // nlc::clf3D::params.Q = (float) H.data();
//     //   memcpy(&nlc::clf3D::params.Q, H.data(), sizeof(H.data()));

//     // double Q[9];
//     double *p1 = &nlc::clf3D::params.Q[0];
//     Eigen::Map<Eigen::Matrix<double, 3, 3> >(&nlc::clf3D::params.Q[0],3,3) = H.cast<double>();
//     Eigen::Map<Eigen::Matrix<double, 3, 1> >(&nlc::clf3D::params.c[0],3,1) = f.cast<double>();
//     Eigen::Map<Eigen::Matrix<double, 1, 3> >(&nlc::clf3D::params.A[0],1,3) = A.cast<double>();
//     H(2,2) = 4;

//     // Eigen::Map<Eigen::Matrix<double, 3, 3> >(&nlc::clf3D::params.Q[0],3,3) = b.cast<double>();

//     for (int i = 0; i < 3; i++) {
//         for (int j = 0; j < 3; j++) {
//             // nlc::clf3D::params.Q[i + j * 3] = H(i, j);  // Q[i+j*3] = H(i,j);
//         }
//         // nlc::clf3D::params.c[i] = 0;
//         // nlc::clf3D::params.A[i] = A(i);
//         nlc::clf3D::params.xlb[i] = -100;
//         nlc::clf3D::params.xub[i] = 100;
//     }
//     nlc::clf3D::params.b[0] = -b(0);

//     for (int i = 0; i < 3; i++) {
//         for (int j = 0; j < 3; j++) {
//             printf("Q[%d, %d] = %f ", i, j, nlc::clf3D::params.Q[i + j * 3]);
//         }
//         printf("c[%d]=%f ", i, nlc::clf3D::params.c[i]);
//         printf("A[%d]=%f ", i, nlc::clf3D::params.A[i]);
//         printf("xlb[%d]=%f ", i, nlc::clf3D::params.xlb[i]);
//         printf("xub[%d]=%f\n", i, nlc::clf3D::params.xub[i]);
//     }
//     printf("b=%f\n", nlc::clf3D::params.b[0]);


//   nlc::clf3D::set_defaults();
//   nlc::clf3D::setup_indexing();
// //   nlc::clf3D::load_default_data();

//     nlc::clf3D::settings.verbose = 1;
//     nlc::clf3D::vars.x[0] = 0.0;
//     nlc::clf3D::vars.x[1] = 0.0;
//     nlc::clf3D::vars.x[2] = 0.0;

//     long num_iters = nlc::clf3D::solve();

//     for (int i = 0; i < 3; i++) {
//         std::cout << nlc::clf3D::vars.x[i] << " x[" << i << "]" << std::endl;
//     }


//     double Q[9];
//     double *p = &Q[0];
//     Eigen::Map<Eigen::Matrix<double, 3, 3> >(p,3,3) = H.cast<double>();

//     std::cout << "+++++++++++++++++++++++++++" << std::endl;
//     std::cout << "Q" << std::endl;
//     for (int i = 0; i < 9; ++i) {
//         std::cout << "params.Q " << nlc::clf3D::params.Q[i] << " <---> Q " << *(Q+i) << std::endl;
//     }
//     std::cout << "\n";
}
