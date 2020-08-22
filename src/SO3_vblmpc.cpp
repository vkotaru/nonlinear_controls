#include "SO3_vblmpc.h"

namespace nonlinear_control {

template <typename T>
SO3VblMPC<T>::SO3VblMPC(const int _N, const double _dt) {
    // temporary inertia matrix
    Eigen::Matrix<T, 3, 3> J = Eigen::Matrix<T, 3, 3>::Identity();
    J(0, 0) = 1e-3;
    J(1, 1) = 1e-3;
    J(2, 2) = 9e-3;
    SO3VblMPC(_N, _dt, J);
}

template <typename T>
SO3VblMPC<T>::SO3VblMPC(const int _N, const double _dt, const  Eigen::Matrix<T, 3, 3> _J) : N(_N), dt(_dt){
    nx = 6;
    nu = 3;
    inertia_ = _J.template cast<double>();
    inertia_inv_ = inertia_.inverse();
    mpcSolver = new LinearMPCBase(true, N, nx, nu);

    /* setting up discrete-translational dynamics
   * dot eta_{k} = -hat(Omega_d)*eta_{k} + I*(\delta\Omega_{k})
   * dot I*(\delta\Omega_{k}) =
   * (JQ^{-1}(\hat(JQ*Omegad)-hat(\Omegad)*JQ))*(\delta\Omega_{k})
   */
    // TODO

    mpcSolver->gains.Q = Eigen::Matrix<double, 6, 6>::Identity();
    mpcSolver->gains.P = 100 * Eigen::Matrix<double, 6, 6>::Identity();
    mpcSolver->gains.R = 1 * Eigen::Matrix<double, 3, 3>::Identity();

    mpcSolver->input.lb << -5, -5, -5;
    mpcSolver->input.ub << 5, 5, 5;
    mpcSolver->state.lb = -500 * Eigen::Matrix<double, 6, 1>::Ones();
    mpcSolver->state.ub = 500 * Eigen::Matrix<double, 6, 1>::Ones();

    uOpt.resize(3 * N, 1);

    Eigen::Matrix<T, 3, 1> Omd = Eigen::Matrix<T, 3, 1>::Zero();
    updateDynamics(Omd);
}

template <typename T>
SO3VblMPC<T>::~SO3VblMPC() {}

template <typename T>
void SO3VblMPC<T>::updateGains(Eigen::Matrix<T, 6, 6> Q, Eigen::Matrix<T, 6, 6> P, Eigen::Matrix<T, 3, 3> R) {
    mpcSolver->gains.Q = Q.template cast<double>();
    mpcSolver->gains.P = P.template cast<double>();
    mpcSolver->gains.R = R.template cast<double>();
}

template<typename T>
void SO3VblMPC<T>::reconstructMPC() {
    mpcSolver->construct();
}

template <typename T>
void SO3VblMPC<T>::updateDynamics(Eigen::Matrix<T, 3, 1> Omd) {
    Eigen::Matrix<double, 6, 6> A = Eigen::Matrix<double, 6, 6>::Identity();
    Eigen::Matrix<double, 6, 3> B = Eigen::Matrix<double, 6, 3>::Zero();
    Eigen::Matrix<double, 3, 3> A11, A22;

    Eigen::Matrix<double, 3, 1> Omdd = Omd.template cast<double> ();
    A11 = -utils::hatd(Omdd);
    A22 = inertia_inv_ * (utils::hatd(inertia_ * Omdd) - utils::hatd(Omdd) * inertia_);
    A << A11, Eigen::Matrix<double, 3, 3>::Identity(),
        Eigen::Matrix<double, 3, 3>::Zero(), A22;
    B << Eigen::Matrix<double, 3, 3>::Zero(), inertia_inv_;

    mpcSolver->problem_.A.setIdentity();
    mpcSolver->problem_.A += dt * A;
    mpcSolver->problem_.B = B * dt;
    mpcSolver->construct();
}

template <typename T>
void SO3VblMPC<T>::init() {
}

template <typename T>
void SO3VblMPC<T>::run(T dt, TSO3<T> x, TSO3<T> xd, Eigen::Matrix<T, 3, 1>& u) {
    updateDynamics(xd.Omega);
    uOpt = mpcSolver->update(false, (x.error(xd)).template cast<double>());
    u = uOpt.block(0, 0, 3, 1).cast<T>();
}

template <typename T>
void SO3VblMPC<T>::run(Eigen::Matrix<T, 6, 1> _err_state, Eigen::Matrix<T, 3, 1>& u) {
    uOpt = mpcSolver->update(false, _err_state.template cast<double>());
    u = uOpt.block(0, 0, 3, 1).cast<T>();
}

template <typename T>
void SO3VblMPC<T>::updateState(Eigen::Matrix<double, 6, 1> &state, Eigen::Vector3d input) {
    state = mpcSolver->problem_.A * state + mpcSolver->problem_.B * input;
}

template class  SO3VblMPC<float>;
template class  SO3VblMPC<double>;

}  // namespace nonlinear_control
