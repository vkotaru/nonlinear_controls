#ifndef NONLINEAR_CONTROLS_SO3_VBLMPC_H
#define NONLINEAR_CONTROLS_SO3_VBLMPC_H

#include <eigen3/Eigen/Dense>
#include <iostream>
#include <qpOASES.hpp>

#include "data_types/data_types.hpp"
#include "geometric_control.h"
#include "linear_mpc.h"
#include "common/qpoases_eigen.hpp"
#include "common/utils.hpp"

namespace nonlinear_controls {

template <typename T>
class SO3VblMPC : public GeometricController<T> {
  protected:
    const double g = 9.81;
    int N, nx, nu;
    T dt;
    Eigen::Matrix<T, 3, 3> inertia_, inertia_inv_;
    TSO3<T> state_;
    bool IS_LTI = true; // is linear time invariant

    LinearMPC<T>* mpcSolver;
    VectorX<T> err_state, uOpt;
    TSO3<T> state_des;

    MatrixX<T> A, B; // Dynamics
    MatrixX<T> Q, P, R; // gains
    VectorX<T> state_lb, state_ub, input_lb, input_ub; // bounds

  public:
    SO3VblMPC(bool islti, int _N, double _dt);
    SO3VblMPC(bool islti, int _N, double _dt, Eigen::Matrix<T, 3, 3> _J);
    ~SO3VblMPC();

    void generate_dynamics(Eigen::Matrix<T, 3, 1> Omd, MatrixX<T>& A);
    void init_dynamics(Eigen::Matrix<T, 3, 1> Omd);
    void init_dynamics(TSO3<T> attd);
    void construct() {
        mpcSolver->construct();
    }

    virtual void init() {}
    virtual void run(T dt) {}
    VectorX<T> run(VectorX<T> _err_state) {
        uOpt = mpcSolver->run(_err_state);
        return uOpt.block(0, 0, nu, 1);
    }
    VectorX<T> run(TSO3<T> state) {
        uOpt = mpcSolver->run(state - state_des);
        return uOpt.block(0, 0, nu, 1);
    }
    void run(T dt, TSO3<T> x, TSO3<T> xd, Eigen::Matrix<T, 3, 1>& u);
    void set_gains(const MatrixX<T> Q, const MatrixX<T> P, const MatrixX<T> R);
    void set_input_bounds(VectorX<T> lb, VectorX<T> ub);
    void set_state_bounds(VectorX<T> lb, VectorX<T> ub);

    /// funtion to integrate dynamics
    void updateState(Eigen::Matrix<double, 6, 1>& state, Eigen::Vector3d input);
};

}  // namespace nonlinear_controls
#endif  // NONLINEAR_CONTROLS_SO3_VBLMPC_H
