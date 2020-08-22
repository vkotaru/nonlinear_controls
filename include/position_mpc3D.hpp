#include <eigen3/Eigen/Dense>
#include <iostream>
#include <qpOASES.hpp>

#include "mpc_base.h"
#include "qpoases_eigen.hpp"

namespace nonlinear_control {

class PositionMPC3D {
   protected:
    int N, nx, nu;
    double dt;

    LinearMPCBase *mpcSolver;
    const double g = 9.81;
    Eigen::Matrix<double, 6, 1> err_state;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> uOpt;

   public:
    PositionMPC3D(const int _N, const double _dt) : N(_N), dt(_dt) {
        nx = 6;
        nu = 3;
        mpcSolver = new LinearMPCBase(true, N, nx, nu);

        /* setting up discrete-translational dynamics
         * x_{k+1} = x_{k} + dt*v_{k} + 0.5*dt*dt*a_{k}
         * v_{k+1} = v_{k} + dt*a_{k}
         */
        mpcSolver->problem_.A.setIdentity();
        mpcSolver->problem_.A.topRightCorner(3, 3) += dt * Eigen::Matrix<double, 3, 3>::Identity();
        mpcSolver->problem_.B << 0.5 * dt * dt * Eigen::Matrix<double, 3, 3>::Identity(), dt * Eigen::Matrix<double, 3, 3>::Identity();

        mpcSolver->gains.Q = Eigen::Matrix<double, 6, 6>::Identity();
        mpcSolver->gains.P = 100 * Eigen::Matrix<double, 6, 6>::Identity();
        mpcSolver->gains.R = 1e-6 * Eigen::Matrix<double, 3, 3>::Identity();

        mpcSolver->input.lb << -g, -g, -g;
        mpcSolver->input.ub << g, g, g;
        mpcSolver->state.lb = -100*Eigen::Matrix<double, 6, 1>::Ones();
        mpcSolver->state.ub = 100*Eigen::Matrix<double, 6, 1>::Ones();

        uOpt.resize(3*N,1);
    }

    ~PositionMPC3D() {}

    void init() {
        mpcSolver->construct();
    }

    Eigen::Vector3d run(Eigen::Matrix<double, 6, 1> _err_state) {
        uOpt = mpcSolver->update(false, _err_state);
        return uOpt.block(0,0,3,1);
    }

    void updateState(Eigen::Matrix<double, 6, 1> &state, Eigen::Vector3d input) {
        state = mpcSolver->problem_.A*state + mpcSolver->problem_.B*input;
    }
};

}  // namespace nonlinear_control
