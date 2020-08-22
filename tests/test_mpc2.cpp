#include "mpc_base2.h"
namespace nlc = nonlinear_control;
int main() {
    const int N = 5;
    const double dt = double(1 / 200.0);  // 200 Hz
    const int nx = 6;
    const int nu = 3;

    nlc::LinearMPCBase2 mpc(true, N, nx, nu);
    Eigen::Matrix<double, 6, 6> Q, P;
    Eigen::Matrix<double, 3, 3> R;
    Q = 1*Eigen::MatrixXd::Identity(6,6);
    P = 2*Eigen::MatrixXd::Identity(6,6);
    R = 3*Eigen::MatrixXd::Identity(3,3);

    Eigen::Matrix<double, 2, 2> R2;
    R2.setIdentity();

    mpc.set_gains(Q, P, R);

    return 0;
}
