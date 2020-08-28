#include "matplotlibcpp.h"
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
namespace plt = matplotlibcpp;

int main() {
    const int N = 5;
    const double dt = double(1 / 200.0);  // 200 Hz
    const double T = 50;
    int nx, nu;
    std::vector<double> t, x, y, z, ux, uy, uz;


    //    /* 1D double integrator */
    //    std::cout << "------------------------------------------" << std::endl;
    //    std::cout << "Testing 1D position MPC with LTI Dynamics... " << std::endl;
    //    nx = 2, nu = 1;
    //    nlc::LinearMPCBase2<double> pos1D_mpc(true, N, nx, nu);
    //    Eigen::Matrix<double, 2, 2> Q, P, A;
    //    Eigen::Matrix<double, 2, 1> B;
    //    Eigen::Matrix<double, 1, 1> R;
    //    Q << 1000, 0, 0, 100;
    //    P << 8.1314 * 1e4, 0.6327 * 1e4, 0.6327 * 1e4, 0.2606 * 1e4;
    //    R << 1e-3;
    //    pos1D_mpc.set_mpc_gains(Q, P, R);

    //    A << 1, dt, 0, 1;
    //    B << 0.5 * dt* dt, dt;

    //    std::cout << "Initializing dynamics ... " << std::endl;
    //    std::cout << "A: \n" << A << std::endl;
    //    std::cout << "B: \n" << B << std::endl;
    //    pos1D_mpc.init_dynamics(N, A, B);

    //    nlc::MPCBounds<double> state_bnds, input_bnds;
    //    state_bnds.lb = -100 * Eigen::Matrix<double, 2, 1>::Ones();
    //    state_bnds.ub = -state_bnds.lb;
    //    input_bnds.lb = -100 * Eigen::Matrix<double, 1, 1>::Ones();
    //    input_bnds.ub = -input_bnds.lb;
    //    pos1D_mpc.set_bounds(state_bnds, input_bnds);

    //    Eigen::Matrix<double, 2, 1> goal_state, state;
    //    goal_state << 1.0, 0.0;
    //    state.setZero();
    //    std::cout << "\n X0: " << state.transpose() << std::endl;
    //    std::cout << "\n Xgoal: " << goal_state.transpose() << std::endl;
    //    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> zOpt, stateMPC, inputMPC;
    //    stateMPC.resize(N * nx, 1);
    //    inputMPC.resize(N * nu, 1);

    //    for (int j = 0; j < int(T / dt); ++j) {
    //        zOpt = pos1D_mpc.update((state - goal_state));
    //        // extracting states and inputs from the decision variables
    //        stateMPC = zOpt.block(0, 0, nx * N, 1);
    //        inputMPC = zOpt.block(nx * N, 0, nu * N, 1);
    //        //        std::cout << "state: " << state.transpose() << " uOpt: " << inputMPC.block(0, 0, nu, 1).transpose() << std::endl;
    //        // integrating dynamics using the input at N = 0
    //        state = A * state + B * inputMPC.block(0, 0, nu, 1);

    //        t.push_back(dt * j);
    //        x.push_back(state(0));
    //        ux.push_back(inputMPC(0));
    //    }
    //    // Set the "super title"
    //    plt::suptitle("Position MPC");
    //    plt::subplot(1, 2, 1);
    //    plt::plot(t, x, "r-");
    //    //    plt::plot(t, y, "g--");
    //    //    plt::plot(t, z, "b-.");
    //    plt::grid(true);

    //    plt::subplot(1, 2, 2);
    //    plt::plot(t, ux, "r-");
    //    //    plt::plot(t, uy, "g--");
    //    //    plt::plot(t, uz, "b-.");
    //    plt::grid(true);
    //    plt::show();
    //    std::cout << "------------------------------------------" << std::endl;

    /* position mpc with LTI dynamics */
    std::cout << "------------------------------------------" << std::endl;
    std::cout << "Testing 3D position MPC with LTI Dynamics... " << std::endl;
    nx = 6, nu = 3;
    nlc::LinearMPCBase2<double> pos_mpc(true, N, nx, nu);
    Eigen::Matrix<double, 6, 6> Q3, P3, A3;
    Eigen::Matrix<double, 6, 3> B3;
    Eigen::Matrix<double, 3, 3> R3;
    Q3 << 1000 * Eigen::Matrix3d::Identity(), Eigen::Matrix3d::Zero(),
    Eigen::Matrix3d::Zero(), 100 * Eigen::Matrix3d::Identity();
    P3 <<     8.1314,    0.0000,   - 0.0000,    0.6327,   - 0.0000,   - 0.0000,
    0.0000,    8.1314,    0.0000,    0.0000,    0.6327,    0.0000,
    - 0.0000,    0.0000,    8.1314,   - 0.0000,    0.0000,    0.6327,
    0.6327,    0.0000,   - 0.0000,    0.2606,   - 0.0000,   - 0.0000,
    - 0.0000,    0.6327,    0.0000,   - 0.0000,    0.2606,    0.0000,
    - 0.0000,    0.0000,    0.6327,   - 0.0000,    0.0000,    0.2606;
    P3 = 1e4 * P3;
    R3 = 1e-3 * Eigen::Matrix<double, 3, 3>::Identity();
    pos_mpc.set_mpc_gains(Q3, P3, R3);

    A3.setIdentity();
    A3.topRightCorner(3, 3) += dt * Eigen::Matrix<double, 3, 3>::Identity();
    B3 << 0.5 * dt* dt* Eigen::Matrix<double, 3, 3>::Identity(), dt* Eigen::Matrix<double, 3, 3>::Identity();
    std::cout << "Initializing dynamics ... " << std::endl;
    std::cout << "A: \n" << A3 << std::endl;
    std::cout << "B: \n" << B3 << std::endl;
    pos_mpc.init_dynamics(N, A3, B3);

    nlc::MPCBounds<double> state_bnds3, input_bnds3;
    state_bnds3.lb = -100 * Eigen::Matrix<double, 6, 1>::Ones();
    state_bnds3.ub = -state_bnds3.lb;
    input_bnds3.lb = -100 * Eigen::Matrix<double, 3, 1>::Ones();
    input_bnds3.ub = -input_bnds3.lb;
    pos_mpc.set_bounds(state_bnds3, input_bnds3);

    Eigen::Matrix<double, 6, 1> goal_state3, state3;
    goal_state3 << 1.0, 1.0, 1.0, 0.0, 0.0, 0.0;
    state3.setZero();
    std::cout << "\n X0: " << state3.transpose() << std::endl;
    std::cout << "\n Xgoal: " << goal_state3.transpose() << std::endl;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> K3, zOpt3, stateMPC3, inputMPC3, uOpt3;
    uOpt3.resize(nu, 1);
    stateMPC3.resize(N * nx, 1);
    inputMPC3.resize(N * nu, 1);
    K3.resize(nu, nx);
    K3 << 30.6287 * Eigen::Matrix3d::Identity(),  12.4527 * Eigen::Matrix3d::Identity();

    t.clear();
    x.clear();
    ux.clear();
    for (int j = 0; j < int(2 / dt); ++j) { // only for 2 seconds

        uOpt3 = -K3 * (state3 - goal_state3);


        //        zOpt3 = pos_mpc.update((state3 - goal_state3));
        //        // extracting states and inputs from the decision variables
        //        stateMPC3 = zOpt3.block(0, 0, nx * N, 1);
        //        inputMPC3 = zOpt3.block(nx * N, 0, nu * N, 1);

        //        std::cout << "state: " << state.transpose() << " uOpt: " << inputMPC.block(0, 0, nu, 1).transpose() << std::endl;

        // integrating dynamics using the input at N = 0
        std::cout << "Input:>>> " << uOpt3.transpose() << std::endl;
        state3 = A3 * state3 + B3 * uOpt3;

        t.push_back(dt * j);
        x.push_back(state3(0));
        y.push_back(state3(1));
        z.push_back(state3(2));
        ux.push_back(uOpt3(0));
        uy.push_back(uOpt3(1));
        uz.push_back(uOpt3(2));
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
    std::cout << "------------------------------------------" << std::endl;


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
