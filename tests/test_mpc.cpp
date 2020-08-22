#include <eigen3/Eigen/Dense>
#include <iostream>
#include <qpOASES.hpp>

#include "mpc_base.h"
#include "qpoases_eigen.hpp"
#include "position_mpc3D.hpp"
#include "SO3_vblmpc.h"
#include "data_types.hpp"
#include "eigen3/unsupported/Eigen/MatrixFunctions"
#include "utils.hpp"
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

    nlc::LinearMPCBase mpc(true, N, nx, nu);

    /* setting up discrete-translational dynamics 
    * x_{k+1} = x_{k} + dt*v_{k} + 0.5*dt*dt*a_{k}
    * v_{k+1} = v_{k} + dt*a_{k}
    */
    mpc.problem_.A.setIdentity();
    mpc.problem_.A.topRightCorner(3, 3) += dt * Eigen::Matrix<double, 3, 3>::Identity();
    mpc.problem_.B << 0.5*dt*dt * Eigen::Matrix<double, 3, 3>::Identity(), dt * Eigen::Matrix<double, 3, 3>::Identity();

    mpc.gains.Q = Eigen::Matrix<double, 6, 6>::Identity();
    mpc.gains.P = 100*Eigen::Matrix<double, 6, 6>::Identity();
    mpc.gains.R = 1e-6*Eigen::Matrix<double, 3, 3>::Identity();

    const double g = 9.81;
    mpc.input.lb << -g, -g, -g;
    mpc.input.ub << g, g, g;
    mpc.state.lb = -100*Eigen::Matrix<double, 6, 1>::Ones();
    mpc.state.ub = 100*Eigen::Matrix<double, 6, 1>::Ones();

    mpc.construct();

    /************************** position MPC **************************************/
    nlc::PositionMPC3D pos_mpc_(N, dt);
    pos_mpc_.init();

    Eigen::Matrix<double, 6, 1> Xgoal, X0;
    Xgoal << 1.0, 1.0, 1.0, 0.0, 0.0, 0.0;
    X0.setZero();
    std::cout << "\n X0: " << X0.transpose() << std::endl;
    std::cout << "\n Xgoal: " << Xgoal.transpose() << std::endl;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> uOpt;

    float T = 20; 
    for (int j = 0; j < int(T/dt); ++j) {
        uOpt = pos_mpc_.run(X0-Xgoal);
        std::cout << "X0: " << X0.transpose() << " uOpt: " << uOpt.transpose() << std::endl;
        pos_mpc_.updateState(X0, uOpt);
    }

    /************************** attitude MPC **************************************/
    Eigen::Matrix3d J, iJ;
    J << 112533, 0, 0, 0, 36203, 0, 0, 0, 42673;
    J = J*1e-6;
    iJ = J.inverse();
    nlc::SO3VblMPC<double> att_mpc_(N, dt, J);
    att_mpc_.init();
    std::cout << "MPC Constructed" << std::endl;

    nlc::TSE3<double> att, attd;
    att.R.setZero();
    att.R(0,2) = 1; att.R(1,1) = 1; att.R(2,0) = -1;
    att.print();
    attd.print();
    std::cout << att-attd << std::endl;
    Eigen::Vector3d uMoment;


    Eigen::Matrix<double, 6, 6> Q, P;
    Eigen::Matrix<double, 3, 3> R;
    Q.setIdentity();
    P.setIdentity();
    R.setIdentity();
    P = P*100;
    R = R*1e-3;

    double ydp;
    double dydp;
    double d2ydp;
    double freqp = 0.5;



    att_mpc_.updateGains(Q, P, R);
    att_mpc_.reconstructMPC();
    for (int j = 0; j < int(T/dt); ++j) {

        sinusoidalTraj(ydp, dydp, d2ydp, dt*j, freqp, 0.25, 1, -0.125, 3.141 / 4);
        attd.R = nlc::utils::rotmZ(0) * nlc::utils::rotmY(ydp);
        attd.Omega(1) = dydp;
        attd.dOmega(1) = d2ydp;
        // att_mpc_.run(att-attd, uMoment);
        att_mpc_.run(dt, att.extractTSO3(), attd.extractTSO3(), uMoment);
        std::cout << "att error: " << (att-attd).transpose() << " uOpt: " << uMoment.transpose()  << std::endl;

        // dynamics integration
        Eigen::Matrix<double, 3, 3> hat_Om = nlc::utils::hatd(att.Omega);
        att.R = att.R*(hat_Om*dt).exp();
        Eigen::Vector3d dOm;
        dOm = iJ*(uMoment -att.Omega.cross(J*att.Omega));
        att.Omega += dOm*dt;
    }


    return 0;
}