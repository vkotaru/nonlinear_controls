#include "mpc_base.h"
namespace nonlinear_control {

LinearMPCBase::LinearMPCBase(const bool _isLTI, const int _N, const int _nx, const int _nu)
    : isLTI(_isLTI), nx(_nx), nu(_nu), N(_N) {
    // setup the qp solver
    nVars = N * nu;
    nCons = (N + 1) * nx;  // Fix this later
    qp = new QPOasesEigen(nVars, nCons);
    qp->options.setToMPC();
    qp->options.printLevel = qpOASES::PL_LOW;
    qp->setup();
    qp->data_.nWSR = 20;

    // initialize matrices
    problem_.N = N;
    problem_.nx = nx;
    problem_.nu = nu;

    problem_.A.resize(nx, nx);
    problem_.B.resize(nx, nu);
    state.lb.resize(nx, 1);
    state.ub.resize(nu, 1);
    input.lb.resize(nu, 1);
    input.ub.resize(nu, 1);

    Sx.resize(nx * (N + 1), nx);
    Sx.setZero();
    Su.resize(nx * (N + 1), N * nu);
    Su.setZero();
    Qbar.resize(nx * (N + 1), nx * (N + 1));
    Qbar.setIdentity();
    Rbar.resize(nu * N, nu * N);
    Rbar.setIdentity();

    gains.Q.resize(nx, nx);
    gains.Q.setIdentity();
    gains.R.resize(nu, nu);
    gains.R.setIdentity();
    gains.P.resize(nu, nu);
    gains.P.setIdentity();

    H.resize(nVars, nVars);
    H.setIdentity();
    F.resize(nx, nVars);
    F.setZero();

    Ulb.resize(nu * N, 1);
    Ulb.setOnes();
    Ulb = -INFINITY * Ulb;
    Uub.resize(nu * N, 1);
    Uub.setOnes();
    Uub = INFINITY * Uub;
    Xlb.resize(nx * (N + 1), 1);
    Xlb.setOnes();
    Xlb = -INFINITY * Xlb;
    Xub.resize(nx * (N + 1), 1);
    Xub.setOnes();
    Xub = INFINITY * Xub;
}

LinearMPCBase::~LinearMPCBase() {}

void LinearMPCBase::construct() {
    // generating the data
    // Using the MPC notation from ME231A
    Sx.block(0, 0, nx, nx).setIdentity();
    Su.block(0, 0, nx, N * nu).setZero();

    Xlb.block(0, 0, nx, 1) = state.lb;
    Xub.block(0, 0, nx, 1) = state.ub;
    for (int i = 0; i < N; ++i) {
        Sx.block(nx * (i + 1), 0, nx, nx) =
            Sx.block(nx * i, 0, nx, nx) * problem_.A;
        if (i == 0)
            Su.block(nx * (i + 1), 0, nx, nu) = problem_.B;
        else {
            Su.block(nx * (i + 1), 0, nx, nu) = problem_.A * Su.block(nx * i, 0, nx, nu);
            Su.block(nx * (i + 1), nu, nx, nu * (N - 1)) = Su.block(nx * i, 0, nx, nu * (N - 1));
        }
        Qbar.block(nx * i, nx * i, nx, nx) = gains.Q;
        Rbar.block(nu * i, nu * i, nu, nu) = gains.R;

        Ulb.block(nu * i, 0, nu, 1) = input.lb;
        Uub.block(nu * i, 0, nu, 1) = input.ub;

        Xlb.block(nx * (i + 1), 0, nx, 1) = state.lb;
        Xub.block(nx * (i + 1), 0, nx, 1) = state.ub;
    }
    Qbar.block(nx * N, nx * N, nx, nx) = gains.P;
    // Cost function
    H = Su.transpose() * Qbar * Su + Rbar;
    F = Sx.transpose() * Qbar * Su;

    debug_print();
}

Eigen::Matrix<double, Eigen::Dynamic, 1> LinearMPCBase::update(const bool verbose, const Eigen::Matrix<double, Eigen::Dynamic, 1> x0) {
    qp->data_.H = H;
    qp->data_.lb = Ulb;
    qp->data_.ub = Uub;
    qp->data_.A = Su;
    qp->data_.lbA = Xlb - Sx * x0;
    qp->data_.ubA = Xub - Sx * x0;
    qp->data_.g = 2 * F.transpose() * x0;

    if (verbose) {
        qp->print();
        debug_print();
        std::cout << "Sx:\n" << Sx << std::endl;
        std::cout << "x0: " << x0.transpose() << std::endl;
        std::cout << "Sx*x0: " << (Sx * x0).transpose() << std::endl;
        std::cout << "Xlb: " << Xlb.transpose() << "\nXub: " << Xub.transpose() << std::endl;
        std::cout << "Ulb: " << Ulb.transpose() << "\nUub: " << Uub.transpose() << std::endl;
    }

    qp->solve();
    return qp->getOptimizer();
}

void LinearMPCBase::debug_print() {
    std::cout << "LTI Dynamics: \nA: \n"
              << problem_.A << "\nB:\n"
              << problem_.B << std::endl;

    std::cout << "nVars: " << nVars << " nCons: " << nCons << std::endl;
    std::cout << "Sx: " << Sx.rows() << " x " << Sx.cols() << std::endl;
    std::cout << "Su: " << Su.rows() << " x " << Su.cols() << std::endl;
    std::cout << "Qbar: " << Qbar.rows() << " x " << Qbar.cols() << std::endl;
    std::cout << "Rbar: " << Rbar.rows() << " x " << Rbar.cols() << std::endl;
    std::cout << "H: " << H.rows() << " x " << H.cols() << std::endl;
    std::cout << "F: " << F.rows() << " x " << F.cols() << std::endl;
}

// template <typename T>
// void LinearMPCBase<T>::update() {

// }

}  // namespace nonlinear_control
