#include <eigen3/Eigen/Dense>
#include <iostream>
#include <qpOASES.hpp>

#include "mpc_base.h"
#include "qpoases_eigen.hpp"
#include "position_mpc3D.hpp"
#include "SO3_vblmpc.h"
#include "SO3_control.h"
#include "data_types.hpp"
#include "eigen3/unsupported/Eigen/MatrixFunctions"
#include "utils.hpp"
#include <fstream>
#include <iostream>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
namespace nlc = nonlinear_control;

#define nlc_real double

void sinusoidalTraj(double& yd, double& dyd, double& d2yd, const float t, const double freq, const double amp, const double ampPercent, double ampOffset, double freqOffset) {
    yd = amp * (1 - ampPercent + ampPercent * pow(cos(freq * t + freqOffset), 2)) + ampOffset;
    dyd = -amp * ampPercent * freq * sin(2 * freq * t + freqOffset);
    d2yd =  -2 * amp * ampPercent * pow(freq, 2) * cos(2 * freq * t + freqOffset);
}

int main() {
    const int N = 5;
    const double dt = double(1 / 200.0);  // 200 Hz
    const int nx = 6;
    const int nu = 3;
    std::vector<double> t, x, y, z, ux, uy, uz;

    nlc::LinearMPCBase mpc(true, N, nx, nu);

    /* setting up discrete-translational dynamics
    * x_{k+1} = x_{k} + dt*v_{k} + 0.5*dt*dt*a_{k}
    * v_{k+1} = v_{k} + dt*a_{k}
    */
    mpc.problem_.A.setIdentity();
    mpc.problem_.A.topRightCorner(3, 3) += dt * Eigen::Matrix<double, 3, 3>::Identity();
    mpc.problem_.B << 0.5 * dt* dt* Eigen::Matrix<double, 3, 3>::Identity(), dt* Eigen::Matrix<double, 3, 3>::Identity();

    mpc.gains.Q << 1000 * Eigen::Matrix3d::Identity(), Eigen::Matrix3d::Zero(),
                Eigen::Matrix3d::Zero(), 100 * Eigen::Matrix3d::Identity();
    mpc.gains.P <<     8.1314,    0.0000,   - 0.0000,    0.6327,   - 0.0000,   - 0.0000,
                0.0000,    8.1314,    0.0000,    0.0000,    0.6327,    0.0000,
                - 0.0000,    0.0000,    8.1314,   - 0.0000,    0.0000,    0.6327,
                0.6327,    0.0000,   - 0.0000,    0.2606,   - 0.0000,   - 0.0000,
                - 0.0000,    0.6327,    0.0000,   - 0.0000,    0.2606,    0.0000,
                - 0.0000,    0.0000,    0.6327,   - 0.0000,    0.0000,    0.2606;
    mpc.gains.P = 1e4 * mpc.gains.P;
    mpc.gains.R = 1 * Eigen::Matrix<double, 3, 3>::Identity();


    const double g = 9.81;
    mpc.input.lb << -5 * g, -5 * g, -5 * g;
    mpc.input.ub << 5 * g, 5 * g, 5 * g;
    mpc.state.lb = -100 * Eigen::Matrix<double, 6, 1>::Ones();
    mpc.state.ub = 100 * Eigen::Matrix<double, 6, 1>::Ones();

    mpc.construct();

    /************************** position MPC **************************************/
    nlc::PositionMPC3D pos_mpc_(N, dt);
    pos_mpc_.init();

    Eigen::Matrix<double, 6, 1> Xgoal, X0;
    Xgoal << 1.0, 2.0, 3.0, 0.0, 0.0, 0.0;
    X0.setZero();
    std::cout << "\n X0: " << X0.transpose() << std::endl;
    std::cout << "\n Xgoal: " << Xgoal.transpose() << std::endl;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> uOpt;

    float T = 50;
    for (int j = 0; j < int(T / dt); ++j) {
        uOpt = pos_mpc_.run(X0 - Xgoal);
        //        std::cout << "X0: " << X0.transpose() << " uOpt: " << uOpt.transpose() << std::endl;
        pos_mpc_.updateState(X0, uOpt);

        t.push_back(dt * j);
        x.push_back(X0(0));
        y.push_back(X0(1));
        z.push_back(X0(2));
        ux.push_back(uOpt(0));
        uy.push_back(uOpt(1));
        uz.push_back(uOpt(2));
    }
    // Set the "super title"
    plt::suptitle("Position MPC");
    plt::subplot(1, 2, 1);
    plt::plot(t, x, "r-");
    plt::plot(t, y, "g--");
    plt::plot(t, z, "b-.");
    plt::grid(true);

    plt::subplot(1, 2, 2);
    plt::plot(t, ux, "r-");
    plt::plot(t, uy, "g--");
    plt::plot(t, uz, "b-.");
    plt::grid(true);

    plt::show();

    /************************** attitude MPC **************************************/
    // open a file in write mode.


    Eigen::Matrix3d J, iJ;
    J << 112533, 0, 0, 0, 36203, 0, 0, 0, 42673;
    J = J * 1e-6;
    iJ = J.inverse();
    nlc::SO3VblMPC<double> att_mpc_(N, dt, J);
    nlc::SO3Controller<double> att_ctrl_(J);
    att_ctrl_.init();
    att_ctrl_.gains_.set_kp((Eigen::Matrix<double, 3, 1>() << 25.1721, 16.4628, 17.9819).finished());
    att_ctrl_.gains_.set_kd((Eigen::Matrix<double, 3, 1>() << 8.3715, 5.3605, 5.8649).finished());
    att_mpc_.init();
    std::cout << "MPC Constructed" << std::endl;

    nlc::TSE3<double> att, attd;
    att.R.setZero();
    att.R(0, 2) = 1;
    att.R(1, 1) = 1;
    att.R(2, 0) = -1;
    att.print();
    attd.print();
    std::cout << att - attd << std::endl;
    Eigen::Vector3d uMoment, uGeo;
    uGeo.setZero();


    Eigen::Matrix<double, 6, 6> Q, P;
    Eigen::Matrix<double, 3, 3> R;
    Q << 1000 * Eigen::Matrix3d::Identity(), Eigen::Matrix3d::Zero(),
    Eigen::Matrix3d::Zero(), 100 * Eigen::Matrix3d::Identity();
    P.setIdentity();
    R.setIdentity();
    P <<     6.6514,   -0.0000,   -0.0000,    0.0894,   -0.0000,   -0.0000,
    -0.0000,    6.5123,   -0.0000,   -0.0000,    0.0440,   -0.0000,
    -0.0000,   -0.0000,    6.5231,   -0.0000,   -0.0000,    0.0475,
    0.0894,   -0.0000,   -0.0000,    0.0293,   -0.0000,   -0.0000,
    -0.0000,    0.0440,   -0.0000,   -0.0000,    0.0141,   -0.0000,
    -0.0000,   -0.0000,    0.0475,   -0.0000,   -0.0000,    0.0152;
    P = P * 1e4;
    double ydp;
    double dydp;
    double d2ydp;
    double freqp = 0.5;
    att_mpc_.updateGains(Q, P, R);
    att_mpc_.reconstructMPC();
    for (int j = 0; j < int(T / dt); ++j) {

        sinusoidalTraj(ydp, dydp, d2ydp, dt * j, freqp, 0.25, 1, -0.125, 3.141 / 4);
        // if (j ==0) {
        //     att.R = nlc::utils::rotmZ(0) * nlc::utils::rotmY(ydp);
        //     att.Omega(1) = dydp;
        //     att.dOmega(1) = d2ydp;
        // }
        attd.R = nlc::utils::rotmZ(0) * nlc::utils::rotmY(ydp);
        attd.Omega(1) = dydp;
        attd.dOmega(1) = d2ydp;
        // att_mpc_.run(att-attd, uMoment);
        att_mpc_.run(dt, att.extractTSO3(), attd.extractTSO3(), uMoment);
        att_ctrl_.run(dt, att.extractTSO3(), attd.extractTSO3(), uGeo);
        // std::cout << "att error: " << (att-attd).transpose() << "\n umpc:" << uMoment.transpose()  <<" ugeo:" << uGeo.transpose()  << std::endl;
        // dynamics integration
        Eigen::Matrix<double, 3, 3> hat_Om = nlc::utils::hatd(att.Omega);
        att.R = att.R * (hat_Om * dt).exp();
        Eigen::Vector3d dOm;
        dOm = iJ * (uMoment - att.Omega.cross(J * att.Omega));
        att.Omega += dOm * dt;
    }

    att.R.setZero();
    att.R(0, 2) = 1;
    att.R(1, 1) = 1;
    att.R(2, 0) = -1;
    att.Omega.setZero();
    att.dOmega.setZero();
    for (int j = 0; j < int(T / dt); ++j) {

        sinusoidalTraj(ydp, dydp, d2ydp, dt * j, freqp, 0.25, 1, -0.125, 3.141 / 4);
        // if (j ==0) {
        //     att.R = nlc::utils::rotmZ(0) * nlc::utils::rotmY(ydp);
        //     att.Omega(1) = dydp;
        //     att.dOmega(1) = d2ydp;
        // }
        attd.R = nlc::utils::rotmZ(0) * nlc::utils::rotmY(ydp);
        attd.Omega(1) = dydp;
        attd.dOmega(1) = d2ydp;
        // att_mpc_.run(att-attd, uMoment);
        att_mpc_.run(dt, att.extractTSO3(), attd.extractTSO3(), uMoment);
        att_ctrl_.run(dt, att.extractTSO3(), attd.extractTSO3(), uGeo);
        // std::cout << "att error: " << (att-attd).transpose() << "\n umpc:" << uMoment.transpose()  <<" ugeo:" << uGeo.transpose()  << std::endl;
        // dynamics integration
        Eigen::Matrix<double, 3, 3> hat_Om = nlc::utils::hatd(att.Omega);
        att.R = att.R * (hat_Om * dt).exp();
        Eigen::Vector3d dOm;
        dOm = iJ * (uGeo - att.Omega.cross(J * att.Omega));
        att.Omega += dOm * dt;
    }

    return 0;
}
