#ifndef NONLINEAR_CONTROLS_SO3_VBLMPC_H
#define NONLINEAR_CONTROLS_SO3_VBLMPC_H

#include <eigen3/Eigen/Dense>
#include <iostream>
#include <qpOASES.hpp>

#include "data_types.hpp"
#include "geometric_control.h"
#include "mpc_base.h"
#include "qpoases_eigen.hpp"
#include "utils.hpp"

namespace nonlinear_control {

template <typename T>
class SO3VblMPC : public GeometricController<T> {
   protected:
    int N, nx, nu;
    T dt;
    Eigen::Matrix3d inertia_, inertia_inv_;
    TSO3<T> state_;
    bool isLTI; // is linear time invariant

    LinearMPCBase *mpcSolver;
    const double g = 9.81;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> uOpt;

   public:
    SO3VblMPC(const int _N, const double _dt);
    SO3VblMPC(const int _N, const double _dt, const Eigen::Matrix<T, 3, 3> _J);
    ~SO3VblMPC();

    virtual void init();
    virtual void run(T dt) {}
    void run(T dt, TSO3<T> x, TSO3<T> xd, Eigen::Matrix<T, 3, 1>& u);
    void run(Eigen::Matrix<T, 6, 1> _err_state, Eigen::Matrix<T, 3, 1>& u);
    void updateDynamics(Eigen::Matrix<T, 3, 1> Omd);
    void updateGains(Eigen::Matrix<T, 6, 6> Q, Eigen::Matrix<T, 6, 6> P, Eigen::Matrix<T, 3, 3> R);

    void reconstructMPC();

    // funtion to integrate dynamics
    void updateState(Eigen::Matrix<double, 6, 1> &state, Eigen::Vector3d input);
};

}  // namespace nonlinear_control
#endif  // NONLINEAR_CONTROLS_SO3_VBLMPC_H