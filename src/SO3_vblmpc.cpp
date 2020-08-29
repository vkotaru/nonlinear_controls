#include "SO3_vblmpc.h"

namespace nonlinear_control {

template <typename T>
SO3VblMPC<T>::SO3VblMPC(bool islti, int _N, double _dt) {
    // temporary inertia matrix
    Eigen::Matrix<T, 3, 3> J;
    J << 112533, 0, 0, 0, 36203, 0, 0, 0, 42673;
    J = J * 1e-6;
    SO3VblMPC(islti, _N, _dt, J);
}

template <typename T>
SO3VblMPC<T>::SO3VblMPC(bool islti, int _N, double _dt, Eigen::Matrix<T, 3, 3> _J) : IS_LTI(islti), N(_N), dt(_dt) {
    nx = 6;
    nu = 3;
    inertia_ = _J;
    inertia_inv_ = inertia_.inverse();

    if (IS_LTI) {
        mpcSolver = new LinearMPC<T>(N, nx, nu);
    } else {
        mpcSolver = new LinearMPCt<T>(N, nx, nu);
    }
    /////////////////////////////////////////////////
    /// setting up discrete-translational dynamics
    /// dot eta_{k} = -hat(Omega_d)*eta_{k} + I*(\delta\Omega_{k})
    /// dot I*(\delta\Omega_{k})
    ///     = (JQ^{-1}(\hat(JQ*Omegad)-hat(\Omegad)*JQ))*(\delta\Omega_{k})
    /////////////////////////////////////////////////
    Eigen::Matrix<T, 3, 1> Omd;
    Omd.setZero();
    A.resize(nx, nx);
    B.resize(nx, nu);
    A.setIdentity();
    generate_dynamics(Omd, A);
    B << Eigen::Matrix<T, 3, 3>::Zero(), inertia_inv_;
    B = B * dt;
    std::cout << "Initializing dynamics ... " << std::endl;
    std::cout << "A: \n" << A << std::endl;
    std::cout << "B: \n" << B << std::endl;
    if (IS_LTI) {
        mpcSolver->init_dynamics(A, B);
    }

    /// gains
    Q.resize(nx, nx);
    P.resize(nx, nx);
    R.resize(nu, nu);
    Q << 1000 * Eigen::Matrix<T, 3, 3>::Identity(), Eigen::Matrix<T, 3, 3>::Zero(),
    Eigen::Matrix<T, 3, 3>::Zero(), 100 * Eigen::Matrix<T, 3, 3>::Identity();
    P << 6.6514,   -0.0000,   -0.0000,    0.0894,   -0.0000,   -0.0000,
    -0.0000,    6.5123,   -0.0000,   -0.0000,    0.0440,   -0.0000,
    -0.0000,   -0.0000,    6.5231,   -0.0000,   -0.0000,    0.0475,
    0.0894,   -0.0000,   -0.0000,    0.0293,   -0.0000,   -0.0000,
    -0.0000,    0.0440,   -0.0000,   -0.0000,    0.0141,   -0.0000,
    -0.0000,   -0.0000,    0.0475,   -0.0000,   -0.0000,    0.0152;
    P = 1e4 * P;
    R = 1 * Eigen::Matrix<T, 3, 3>::Identity();
    mpcSolver->set_mpc_gains(Q, P, R);

    /// bounds
    input_lb.resize(nu, 1);
    input_ub.resize(nu, 1);
    state_lb.resize(nx, 1);
    state_ub.resize(nx, 1);
    input_lb << -2, -2, -2;
    input_ub << 2, 2, 2;
    state_lb.setOnes();
    state_lb = -100 * state_lb;
    state_ub = -state_lb;
    mpcSolver->set_input_bounds(input_lb, input_ub);
    mpcSolver->set_state_bounds(state_lb, state_ub);

    //    uOpt.resize(nu * N, 1);

    ///
    if (IS_LTI) {
        mpcSolver->construct();
    }
}

template <typename T>
SO3VblMPC<T>::~SO3VblMPC() {}

template <typename T>
void SO3VblMPC<T>::generate_dynamics(Eigen::Matrix<T, 3, 1> Omd, MatrixX<T>& A) {
    A << -utils::hat<T>(Omd), Eigen::Matrix<T, 3, 3>::Identity(),
    Eigen::Matrix<T, 3, 3>::Zero(), inertia_inv_ * (utils::hat<T>(inertia_ * Omd) - utils::hat<T>(Omd) * inertia_);
    A =  Eigen::Matrix<T, 6, 6>::Identity() + A * dt;
}

template <typename T>
void SO3VblMPC<T>::init_dynamics(Eigen::Matrix<T, 3, 1> Omd) {
    generate_dynamics(Omd, A);
    mpcSolver->init_dynamics(A, B);
}
template <typename T>
void SO3VblMPC<T>::init_dynamics(TSO3<T> attd) {
    generate_dynamics(attd.Omega, A);
    mpcSolver->init_dynamics(A, B);
}

template <typename T>
void SO3VblMPC<T>::run(T dt, TSO3<T> x, TSO3<T> xd, Eigen::Matrix<T, 3, 1>& u) {
    // time-varying tracjectory with only operating point used for the dynamics for the next N steps
    if (IS_LTI) {
        uOpt = mpcSolver->run(x - xd);
    } else {
        generate_dynamics(xd.Omega, A);
        uOpt = mpcSolver->run(x - xd, A, B);
    }
    u = uOpt.block(0, 0, 3, 1);
    u += x.Omega.cross(inertia_ * x.Omega);
    u += -inertia_ * (x.Omega.cross(x.R.transpose() * xd.R * xd.Omega) - x.R.transpose() * xd.R * xd.dOmega);
}

template <typename T>
void SO3VblMPC<T>::set_gains(const MatrixX<T> Q, const MatrixX<T> P, const MatrixX<T> R) {
    this->Q = Q;
    this->P = P;
    this->R = R;
    mpcSolver->set_mpc_gains(Q, P, R);
}
template <typename T>
void SO3VblMPC<T>::set_input_bounds(VectorX<T> lb, VectorX<T> ub) {
    this->input_lb = lb;
    this->input_ub = ub;
    mpcSolver->set_input_bounds(lb, ub);
}
template <typename T>
void SO3VblMPC<T>::set_state_bounds(VectorX<T> lb, VectorX<T> ub) {
    this->state_lb = lb;
    this->state_ub = ub;
    mpcSolver->set_state_bounds(lb, ub);

}

//////////////////////////////////
/// templates
//////////////////////////////////
template class  SO3VblMPC<float>;
template class  SO3VblMPC<double>;
}  // namespace nonlinear_control
