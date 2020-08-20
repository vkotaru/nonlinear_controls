#include <eigen3/Eigen/Dense>
#include <iostream>
#include <qpOASES.hpp>

#include "mpc_base.h"
#include "qpoases_eigen.hpp"
#include "position_mpc3D.hpp"
#include "SO3_vbl_mpc.hpp"

namespace nlc = nonlinear_control;

#define nlc_real double

int main() {
    const int N = 5;
    const double dt = double(1 / 200.0);  // 200 Hz
    // const int nx = 6;
    // const int nu = 3;

    //     nlc::MPCBase mpc(true, N, nx, nu);

    //     /* setting up discrete-translational dynamics 
    //    * x_{k+1} = x_{k} + dt*v_{k} + 0.5*dt*dt*a_{k}
    //    * v_{k+1} = v_{k} + dt*a_{k}
    //    */
    //     mpc.problem_.A.setIdentity();
    //     mpc.problem_.A.topRightCorner(3, 3) += dt * Eigen::Matrix<double, 3, 3>::Identity();
    //     mpc.problem_.B << 0.5*dt*dt * Eigen::Matrix<double, 3, 3>::Identity(), dt * Eigen::Matrix<double, 3, 3>::Identity();

    //     mpc.gains.Q = Eigen::Matrix<double, 6, 6>::Identity();
    //     mpc.gains.P = 100*Eigen::Matrix<double, 6, 6>::Identity();
    //     mpc.gains.R = 1e-6*Eigen::Matrix<double, 3, 3>::Identity();

    //     const double g = 9.81;
    //     mpc.input.lb << -g, -g, -g;
    //     mpc.input.ub << g, g, g;
    //     mpc.state.lb = -100*Eigen::Matrix<double, 6, 1>::Ones();
    //     mpc.state.ub = 100*Eigen::Matrix<double, 6, 1>::Ones();

    //     mpc.construct();
    //     mpc.setup();

    nlc::SO3VblMPC att_mpc_(N, dt);

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

    // std::cout << "\n Ulb >>>> \n" << mpc.Ulb.transpose() << "\n Uub >>>> \n" << mpc.Uub.transpose() << std::endl;
    // std::cout << "\n Qbar >>>> \n" << mpc.Qbar << "\n Rbar >>>> \n" << mpc.Rbar << std::endl;
    return 0;
}