#include "mpc_base2.h"
namespace nonlinear_control {

template <typename T>
LinearMPCBase2<T>::LinearMPCBase2(const bool _isLTI, const int _N, const int _nx, const int _nu) : isLTI(_isLTI), nx(_nx), nu(_nu), N(_N) {
    std::cout << "LinearMPCBase2 Constructor beginning" << std::endl;

    /*************************
     *  setting-the QP solver
     * ***********************/
    nVars = N * (nu + nx);
    nCons = N * nx;

    QP = new QPOasesEigen(nVars, nCons);
    QP->options.setToMPC();
    QP->options.printLevel = qpOASES::PL_LOW;
    QP->setup();
    QP->data_.nWSR = 10e5; //

    // initialize matrices
    problem_.N = N;
    problem_.nx = nx;
    problem_.nu = nu;

    G0.resize(nCons, nVars);
    G0.setIdentity();
    E0.resize(nCons, nx);
    E0.setZero();
    tolCons.resize(nCons, 1);
    tolCons = 1e-6 * VectorX<T>::Ones(nCons, 1);

    problem_.A.resize(nx, nx);
    problem_.B.resize(nx, nu);
    state.lb.resize(nx, 1);
    state.ub.resize(nx, 1);
    input.lb.resize(nu, 1);
    input.ub.resize(nu, 1);

    H.resize(nVars, nVars);
    H.setIdentity();
    g.resize(nVars, 1);
    g.setZero();
    set_mpc_cost(H, g);

    gains.Q.resize(nx, nx);
    gains.Q.setIdentity();
    gains.R.resize(nu, nu);
    gains.R.setIdentity();
    gains.P.resize(nx, nx);
    gains.P.setIdentity();
    set_mpc_gains(gains);

    Ulb.resize(nu * N, 1);
    Ulb.setOnes();
    Ulb = -1e6 * Ulb;
    Uub.resize(nu * N, 1);
    Uub.setOnes();
    Uub = 1e6 * Uub;
    Xlb.resize(nx * N, 1);
    Xlb.setOnes();
    Xlb = -1e6 * Xlb;
    Xub.resize(nx * N, 1);
    Xub.setOnes();
    Xub = 1e6 * Xub;
    set_mpc_bounds();
    std::cout << "LinearMPCBase2 Constructor done" << std::endl;

}

template <typename T>
LinearMPCBase2<T>::~LinearMPCBase2() {}

template <typename T>
void LinearMPCBase2<T>::init_dynamics(MatrixX<T> A, MatrixX<T> B) {
    Astorage.push_back(A);
    Bstorage.push_back(B);
}

template <typename T>
void LinearMPCBase2<T>::init_dynamics(const int N, MatrixX<T> A, MatrixX<T> B) {
    for (int i = 0; i < N; ++i) {
        Astorage.push_back(A);
        Bstorage.push_back(B);
    }
}

template <typename T>
void LinearMPCBase2<T>::update_dynamics(MatrixX<T> A, MatrixX<T> B) {
    Astorage.pop_front();
    Astorage.push_back(A);
    Bstorage.pop_front();
    Bstorage.push_back(B);
}

template <typename T>
void LinearMPCBase2<T>::set_mpc_cost(MatrixX<T> _H) {
    //    H = _H;
    QP->data_.H = _H.template cast<double>();
}

template <typename T>
void LinearMPCBase2<T>::set_mpc_cost(MatrixX<T> _H, VectorX<T> _g) {
    //    H = _H;
    //    g = _g; do I need this?
    QP->data_.H = _H.template cast<double>();
    QP->data_.g = _g.template cast<double>();
}

template <typename T>
void LinearMPCBase2<T>::set_mpc_gains(MatrixX<T> Q, MatrixX<T> P, MatrixX<T> R) {
    assert(Q.rows() == nx && Q.cols() == nx);
    assert(P.rows() == nx && P.cols() == nx);
    assert(R.rows() == nu && R.cols() == nu);

    for (int i = 0; i < N; ++i) {
        H.block(i * nx, i * nx, nx, nx) = Q;
        if (i == N - 1) {
            H.block(i * nx, i * nx, nx, nx) = P;
        }
        H.block(N * nx + i * nu, N * nx + i * nu, nu, nu) = R;
    }
    QP->data_.H = H.template cast<double>();
    QP->data_.g = g.template cast<double>();
}

template <typename T>
void LinearMPCBase2<T>::set_mpc_gains(MPCGains<T> gains) {
    assert(gains.Q.rows() == nx && gains.Q.cols() == nx);
    assert(gains.P.rows() == nx && gains.P.cols() == nx);
    assert(gains.R.rows() == nu && gains.R.cols() == nu);

    for (int i = 0; i < N; ++i) {
        H.block(i * nx, i * nx, nx, nx) = gains.Q;
        if (i == N - 1) {
            H.block(i * nx, i * nx, nx, nx) = gains.P;
        }
        H.block(N * nx + i * nu, N * nx + i * nu, nu, nu) = gains.R;
    }
    QP->data_.H = H.template cast<double>();
    QP->data_.g = g.template cast<double>();
}


template <typename T>
void LinearMPCBase2<T>::set_input_bounds(VectorX<T> lb, VectorX<T> ub) {
    assert(lb.rows() == nu && lb.cols() == 1);
    assert(ub.rows() == nu && ub.cols() == 1);
    input.lb = lb;
    input.ub = ub;
}

template <typename T>
void LinearMPCBase2<T>::set_input_bounds(MPCBounds<T> _input) {
    assert(_input.lb.rows() == nu && _input.lb.cols() == 1);
    assert(_input.ub.rows() == nu && _input.ub.cols() == 1);
    input.lb = _input.lb;
    input.ub = _input.ub;
}

template <typename T>
void LinearMPCBase2<T>::set_state_bounds(VectorX<T> lb, VectorX<T> ub) {
    assert(lb.rows() == nu && lb.cols() == 1);
    assert(ub.rows() == nu && ub.cols() == 1);
    state.lb = lb;
    state.ub = ub;
}

template <typename T>
void LinearMPCBase2<T>::set_state_bounds(MPCBounds<T> _state) {
    assert(_state.lb.rows() == nx && _state.lb.cols() == 1);
    assert(_state.ub.rows() == nx && _state.ub.cols() == 1);
    state.lb = _state.lb;
    state.ub = _state.ub;
}

template <typename T>
void LinearMPCBase2<T>::set_bounds(MPCBounds<T> _state, MPCBounds<T> _input) {
    LinearMPCBase2::set_state_bounds(_state);
    LinearMPCBase2::set_input_bounds(_input);
    LinearMPCBase2::set_mpc_bounds();
}

template <typename T>
void LinearMPCBase2<T>::set_mpc_bounds() {
    for (int i = 0; i < N; ++i) {
        Xlb.block(i * nx, 0, nx, 1) = state.lb;
        Xub.block(i * nx, 0, nx, 1) = state.ub;
        Ulb.block(i * nu, 0, nu, 1) = input.lb;
        Uub.block(i * nu, 0, nu, 1) = input.ub;
    }
    QP->data_.lb << Xlb.template cast<double>(), Ulb.template cast<double>();
    QP->data_.ub << Xub.template cast<double>(), Uub.template cast<double>();
}

template <typename T>
void LinearMPCBase2<T>::construct() {
    assert(Astorage.size() == N);
    assert(Bstorage.size() == N);

    // generating the data
    // Using the MPC notation from ME231A
    E0.block(0, 0, nx, nx) = Astorage.at(0);
    for (int i  = 0; i < N; ++i) {
        G0.block(nx * i, nx * N + nu * i, nx, nu) = -Bstorage.at(i);
        if (i > 0) {
            G0.block(nx * i, nx * (i - 1), nx, nx) = -Astorage.at(i);
        }
    }
    //    std::cout << G0 << std::endl;
    //    std::cout << E0 << std::endl;
    // debug_print();
}

template <typename T>
VectorX<T> LinearMPCBase2<T>::update(const VectorX<T> x0, const MatrixX<T> A, const MatrixX<T> B) {
    update_dynamics(A, B);
    construct();
    return update(x0);
}

template <typename T>
VectorX<T> LinearMPCBase2<T>::update(const VectorX<T> x0) {
    construct();
    QP->data_.A = G0.template cast<double>();
    QP->data_.lbA = (E0 * x0 - tolCons).template cast<double>();
    QP->data_.ubA = (E0 * x0 + tolCons).template cast<double>();

    //    std::cout << "****************************" << std::endl;
    //    std::cout << QP->data_.H << std::endl;
    //    std::cout << QP->data_.g.transpose() << std::endl;
    //    std::cout << QP->data_.lb.transpose() << std::endl;
    //    std::cout << QP->data_.ub.transpose() << std::endl;
    //    std::cout << QP->data_.lbA.transpose() << std::endl;
    //    std::cout << QP->data_.ubA.transpose() << std::endl;
    //    std::cout << QP->data_.A << std::endl;
    //    std::cout << "****************************" << std::endl;

    QP->solve();
    return QP->getOptimizer().template cast<T>();
}

// void LinearMPCBase::debug_print() {
//     std::cout << "LTI Dynamics: \nA: \n"
//               << problem_.A << "\nB:\n"
//               << problem_.B << std::endl;

//     std::cout << "nVars: " << nVars << " nCons: " << nCons << std::endl;
//     std::cout << "Sx: " << Sx.rows() << " x " << Sx.cols() << std::endl;
//     std::cout << "Su: " << Su.rows() << " x " << Su.cols() << std::endl;
//     std::cout << "Qbar: " << Qbar.rows() << " x " << Qbar.cols() << std::endl;
//     std::cout << "Rbar: " << Rbar.rows() << " x " << Rbar.cols() << std::endl;
//     std::cout << "H: " << H.rows() << " x " << H.cols() << std::endl;
//     std::cout << "F: " << F.rows() << " x " << F.cols() << std::endl;
// }

// template <typename T>
// void LinearMPCBase<T>::update() {

// }

template class LinearMPCBase2<float>;
template class LinearMPCBase2<double>;

}  // namespace nonlinear_control
