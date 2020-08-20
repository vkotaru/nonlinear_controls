#include <eigen3/Eigen/Dense>
#include <iostream>
#include <qpOASES.hpp>

#include "mpc_base.h"
#include "utils.hpp"
#include "qpoases_eigen.hpp"

namespace nonlinear_control {

class SO3VblMPC {
   protected:
    int N, nx, nu;
    double dt;
    Eigen::Matrix3d J, iJ;

    MPCBase *mpcSolver;
    const double g = 9.81;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> uOpt;

   public:
    SO3VblMPC(const int _N, const double _dt) {
        // temporary inertia matrix
        Eigen::Matrix3d J = Eigen::Matrix3d::Identity();
        J(0, 0) = 1e-3;
        J(1, 1) = 1e-3;
        J(2, 2) = 9e-3;  
        SO3VblMPC(_N, _dt, J);
    }

    SO3VblMPC(const int _N, const double _dt, const Eigen::Matrix3d _J) : N(_N), dt(_dt), J(_J) {
        nx = 6;
        nu = 3;
        iJ = J.inverse();
        mpcSolver = new MPCBase(true, N, nx, nu);

        /* setting up discrete-translational dynamics
         * dot eta_{k} = -hat(Omega_d)*eta_{k} + I*(\delta\Omega_{k})
         * dot I*(\delta\Omega_{k}) = (JQ^{-1}(\hat(JQ*Omegad)-hat(\Omegad)*JQ))*(\delta\Omega_{k})
         */
        // TODO
        Eigen::Matrix<double, 3, 3> Rd = Eigen::Matrix<double, 3, 3>::Identity();
        Eigen::Matrix<double, 3, 1> Omd = Eigen::Matrix<double, 3, 1>::Zero();
        Eigen::Matrix<double, 6, 6> A = Eigen::Matrix<double, 6, 6>::Identity();
        Eigen::Matrix<double, 6, 3> B = Eigen::Matrix<double, 6, 3>::Zero();
        A <<  -utils::hatd(Omd), Eigen::Matrix<double, 3, 3>::Identity(),
                Eigen::Matrix<double, 3, 3>::Zero(), iJ*(utils::hatd(J*Omd)-utils::hatd(Omd)*J);
        B << Eigen::Matrix<double, 3, 3>::Zero(), iJ;

        mpcSolver->problem_.A.setIdentity();
        mpcSolver->problem_.A += dt * A;
        mpcSolver->problem_.B = B;

        mpcSolver->gains.Q = Eigen::Matrix<double, 6, 6>::Identity();
        mpcSolver->gains.P = 100 * Eigen::Matrix<double, 6, 6>::Identity();
        mpcSolver->gains.R = 1e-6 * Eigen::Matrix<double, 3, 3>::Identity();

        mpcSolver->input.lb << -5, -5, -5;
        mpcSolver->input.ub << 5, 5, 5;
        mpcSolver->state.lb = -500 * Eigen::Matrix<double, 6, 1>::Ones();
        mpcSolver->state.ub = 500 * Eigen::Matrix<double, 6, 1>::Ones();

        uOpt.resize(3 * N, 1);
    }

    ~SO3VblMPC() {}

    void init() {
        mpcSolver->construct();
        mpcSolver->setup();
    }

    Eigen::Vector3d run(Eigen::Matrix<double, 6, 1> _err_state) {
        uOpt = mpcSolver->update(false, _err_state);
        return uOpt.block(0, 0, 3, 1);
    }

    void updateState(Eigen::Matrix<double, 6, 1> &state, Eigen::Vector3d input) {
        state = mpcSolver->problem_.A * state + mpcSolver->problem_.B * input;
    }
};

}  // namespace nonlinear_control
