#ifndef NONLINEAR_CONTROLS_SO3_VBLMPC2_HPP
#define NONLINEAR_CONTROLS_SO3_VBLMPC2_HPP
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <iostream>
#include <qpOASES.hpp>

#include "data_types.hpp"
#include "geometric_control.h"
#include "mpc_base2.h"
#include "qpoases_eigen.hpp"
#include "utils.hpp"

namespace nonlinear_control {

template <typename T>
class SO3VblMPC2 : public GeometricController<T> {
  protected:
    T dt;
    bool isLTI; // is linear time invariant
    int N, nx, nu;
    TSO3<T> state_;
    Eigen::Matrix<T, 3, 3> inertia_, inertia_inv_;
    Eigen::Matrix<T, 6, 6> A;
    Eigen::Matrix<T, 6, 3> B;

    MPCBounds<T> state, input;
    MPCGains<T> gains;

    LinearMPCBase2<T>* mpcSolver;
    VectorX<T> uOpt;

  public:
    SO3VblMPC2(const int _N, const double _dt) {
        Eigen::Matrix<T, 3, 3> J = Eigen::Matrix<T, 3, 3>::Identity();
        J(0, 0) = 1e-3;
        J(1, 1) = 1e-3;
        J(2, 2) = 9e-3;
        SO3VblMPC(_N, _dt, J);
    }
    SO3VblMPC2(const int _N, const double _dt, const Eigen::Matrix<T, 3, 3> _J) : N(_N), dt(_dt), inertia_(_J) {
        inertia_inv_ = inertia_.inverse();
        nx = 6;
        nu = 3;
        uOpt.resize(3 * N, 1);
        mpcSolver = new LinearMPCBase2<T>(true, N, nx, nu);

        B << Eigen::Matrix<double, 3, 3>::Zero(), inertia_inv_;
        B = B * dt;

        gains.Q = Eigen::Matrix<T, 6, 6>::Identity();
        gains.P = 100 * Eigen::Matrix<T, 6, 6>::Identity();
        gains.R = 1 * Eigen::Matrix<T, 3, 3>::Identity();
        mpcSolver->set_mpc_gains(gains);

        input.lb = -5 * Eigen::Matrix<T, 3, 1>::Ones();
        input.ub = 5 * Eigen::Matrix<T, 3, 1>::Ones();
        state.lb = -500 * Eigen::Matrix<T, 6, 1>::Ones();
        state.ub = 500 * Eigen::Matrix<T, 6, 1>::Ones();
        mpcSolver->set_bounds(state, input);
    }
    ~SO3VblMPC2() {

    }

    MatrixX<T> generateA(Eigen::Matrix<T, 3, 1> Omd) {
        Eigen::Matrix<T, 6, 6> A = Eigen::Matrix<T, 6, 6>::Identity();
        A << -utils::hat<T>(Omd), Eigen::Matrix<T, 3, 3>::Identity(),
        Eigen::Matrix<T, 3, 3>::Zero(), inertia_inv_ * (utils::hat<T>(inertia_ * Omd) - utils::hat<T>(Omd) * inertia_);
        return MatrixX<T>::Identity(6, 6) + A * dt;
    }
    void updateDynamics(Eigen::Matrix<T, 3, 1> Omd) {
        mpcSolver->update_dynamics(generateA(Omd), B);
    }
    void initDynamics(Eigen::Matrix<T, 3, 1> Omd) {
        mpcSolver->init_dynamics(generateA(Omd), B);
    }
    void reconstrutMPC () {
        mpcSolver->construct();
    }
    virtual void init() {}
    virtual void run(T dt) {}
    void run(T dt, TSO3<T> x, TSO3<T> xd, Eigen::Matrix<T, 3, 1>& u) {
        uOpt = mpcSolver->update(x.error(xd));
        u = uOpt.block(0, 0, nu, 1);
        u += x.Omega.cross(inertia_ * x.Omega);
        u += -inertia_ * (x.Omega.cross(x.R.transpose() * xd.R * xd.Omega) - x.R.transpose() * xd.R * xd.dOmega);
    }
    void run(Eigen::Matrix<T, 6, 1> _err_state, Eigen::Matrix<T, 3, 1>& u) {
        // time-invariant dynamics
        uOpt = mpcSolver->update(_err_state);
        u = uOpt.block(0, 0, nu, 1);
    }

    void updateGains(Eigen::Matrix<T, 6, 6> Q, Eigen::Matrix<T, 6, 6> P, Eigen::Matrix<T, 3, 3> R) {
        mpcSolver->set_mpc_gains(Q, P, R);
    }
    void updateState(Eigen::Matrix<T, 6, 1>& state, Eigen::Matrix<T, 3, 1> input) {
        state = mpcSolver->getCurrentA() * state + mpcSolver->getCurrentB() * input;
    }
    void updateState(TSO3<T>& state, Eigen::Matrix<T, 3, 1> input) {
        state.R = state.R * (utils::hat<T>(state.Omega) * dt).exp();
        state.Omega += inertia_inv_ * (input - state.Omega.cross(inertia_ * state.Omega)) * dt;
    }
};

}  // namespace nonlinear_control
#endif  // NONLINEAR_CONTROLS_SO3_VBLMPC2_HPP

