#include "mpc_base2.h"
#include "SO3_vblmpc2.hpp"
#include <fstream>
#include <iostream>

void sinusoidalTraj(double& yd, double& dyd, double& d2yd, const float t, const double freq, const double amp, const double ampPercent, double ampOffset, double freqOffset) {
    yd = amp * (1 - ampPercent + ampPercent * pow(cos(freq * t + freqOffset), 2)) + ampOffset;
    dyd = -amp * ampPercent * freq * sin(2 * freq * t + freqOffset);
    d2yd =  -2 * amp * ampPercent * pow(freq, 2) * cos(2 * freq * t + freqOffset);
}
namespace nlc = nonlinear_control;

#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;
int main() {
    plt::plot({1, 3, 2, 4});
    plt::show();
    //    const int N = 1;
    //    const double dt = double(1 / 200.0);  // 200 Hz
    //    const int nx = 6;
    //    const int nu = 3;
    //    const double T = 20;
    //    nlc::LinearMPCBase2<double> mpc(true, N, nx, nu);
    //    Eigen::Matrix<double, 6, 6> Q, P, A;
    //    Eigen::Matrix<double, 6, 3> B;
    //    Eigen::Matrix<double, 3, 3> R;
    //    Q = 1 * Eigen::MatrixXd::Identity(6, 6);
    //    P = 2 * Eigen::MatrixXd::Identity(6, 6);
    //    R = 3 * Eigen::MatrixXd::Identity(3, 3);
    //    mpc.set_mpc_gains(Q, P, R);

    //    A = 2 * Eigen::MatrixXd::Ones(6, 6);
    //    B = 3 * Eigen::MatrixXd::Ones(6, 3);
    //    mpc.init_dynamics(N, A, B);
    //    mpc.construct();

    //    mpc.update_dynamics(10 * Eigen::MatrixXd::Ones(6, 6), 9 * Eigen::MatrixXd::Ones(6, 3));
    //    mpc.construct();

    //    /***********************************************/
    //    Eigen::Matrix3d J, iJ;
    //    J << 112533, 0, 0, 0, 36203, 0, 0, 0, 42673;
    //    J = J * 1e-6;
    //    iJ = J.inverse();
    //    nlc::SO3VblMPC2<double> att_mpc_(N, dt, J);
    //    Eigen::Vector3d uMPC, uGeo;
    //    uMPC.setZero();

    //    Q << 1000 * Eigen::Matrix3d::Identity(), Eigen::Matrix3d::Zero(),
    //    Eigen::Matrix3d::Zero(), 100 * Eigen::Matrix3d::Identity();
    //    P.setIdentity();
    //    R.setIdentity();
    //    P <<     6.6514,   -0.0000,   -0.0000,    0.0894,   -0.0000,   -0.0000,
    //    -0.0000,    6.5123,   -0.0000,   -0.0000,    0.0440,   -0.0000,
    //    -0.0000,   -0.0000,    6.5231,   -0.0000,   -0.0000,    0.0475,
    //    0.0894,   -0.0000,   -0.0000,    0.0293,   -0.0000,   -0.0000,
    //    -0.0000,    0.0440,   -0.0000,   -0.0000,    0.0141,   -0.0000,
    //    -0.0000,   -0.0000,    0.0475,   -0.0000,   -0.0000,    0.0152;
    //    P = P * 1e4;
    //    att_mpc_.updateGains(Q, P, R);

    //    nlc::TSE3<double> att, attd;
    //    Eigen::Vector3d OmegadFuture;
    //    OmegadFuture.setZero();
    //    att.R.setZero();
    //    att.R(0, 2) = 1;
    //    att.R(1, 1) = 1;
    //    att.R(2, 0) = -1;
    //    double ydp;
    //    double dydp;
    //    double d2ydp;
    //    double freqp = 0.5;

    //    std::cout << "Initializing dynamics ... " << std::endl;
    //    for (int j = 0; j < N; ++j) {
    //        //        sinusoidalTraj(ydp, dydp, d2ydp, dt * j, freqp, 0.25, 1, -0.125, 3.141 / 4);
    //        //        attd.R = nlc::utils::rotmZ(0) * nlc::utils::rotmY(ydp);
    //        //        attd.Omega(1) = dydp;
    //        //        attd.dOmega(1) = d2ydp;
    //        att_mpc_.initDynamics(attd.Omega);
    //    }
    //    att_mpc_.reconstrutMPC();
    //    std::cout << "Done" << std::endl;


    //    //    std::ofstream outfile;
    //    //    outfile.open("afilempc.dat");
    //    for (int j = 0; j < int(T / dt); ++j) {
    //        //        sinusoidalTraj(ydp, dydp, d2ydp, dt * j, freqp, 0.25, 1, -0.125, 3.141 / 4);
    //        //        attd.R = nlc::utils::rotmZ(0) * nlc::utils::rotmY(ydp);
    //        //        attd.Omega(1) = dydp;
    //        //        attd.dOmega(1) = d2ydp;
    //        att_mpc_.run(dt, att.extractTSO3(), attd.extractTSO3(), uMPC);
    //        std::cout << "att error: " << (att - attd).transpose() << " umpc:" << uMPC.transpose() << std::endl;
    //        //        outfile << "att error: " << (att - attd).transpose() << " umpc:" << uMPC.transpose() << std::endl;
    //        att_mpc_.updateState(att, uMPC);
    //        //        sinusoidalTraj(ydp, dydp, d2ydp, dt * (j + N), freqp, 0.25, 1, -0.125, 3.141 / 4);
    //        //        OmegadFuture(1) = dydp;
    //        //        att_mpc_.updateDynamics(OmegadFuture);
    //    }
    //    //    outfile.close();
    //    //    outfile.open("afilegeo.dat");
    //    return 0;
}
