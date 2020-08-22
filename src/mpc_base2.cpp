#include "mpc_base2.h"
namespace nonlinear_control {

LinearMPCBase2::LinearMPCBase2(const bool _isLTI, const int _N, const int _nx, const int _nu) : isLTI(_isLTI), nx(_nx), nu(_nu), N(_N) {

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

    problem_.A.resize(nx, nx);
    problem_.B.resize(nx, nu);
    state.lb.resize(nx, 1);
    state.ub.resize(nu, 1);
    input.lb.resize(nu, 1);
    input.ub.resize(nu, 1);

    gains.Q.resize(nx, nx);
    gains.Q.setIdentity();
    gains.R.resize(nu, nu);
    gains.R.setIdentity();
    gains.P.resize(nx, nx);
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
    Xlb.resize(nx * N, 1);
    Xlb = -INFINITY * Eigen::MatrixXd(nx*N, 1);
    Xub.resize(nx * N, 1);
    Xub.setOnes();
    Xub = INFINITY * Xub;
}

LinearMPCBase2::~LinearMPCBase2() {}

void LinearMPCBase2::set_gains(Eigen::MatrixXd Q, Eigen::MatrixXd P, Eigen::MatrixXd R) {
    assert(Q.rows()==nx && Q.cols()==nx);
    assert(P.rows()==nx && P.cols()==nx);
    assert(R.rows()==nu && R.cols()==nu);

    for (int i = 0; i < N; ++i) {
        H.block(i*nx, i*nx, nx, nx) = Q;
        if (i==N-1) {
            H.block(i*nx, i*nx, nx, nx)= P;
        }
        H.block(N*nx+i*nu, N*nx+i*nu, nu, nu) = R;
    }
    QP->data_.H = H;
}

// void LinearMPCBase::construct() {
//     // generating the data
//     // Using the MPC notation from ME231A
//     Sx.block(0, 0, nx, nx).setIdentity();
//     Su.block(0, 0, nx, N * nu).setZero();

//     Xlb.block(0, 0, nx, 1) = state.lb;
//     Xub.block(0, 0, nx, 1) = state.ub;
//     for (int i = 0; i < N; ++i) {
//         Sx.block(nx * (i + 1), 0, nx, nx) =
//             Sx.block(nx * i, 0, nx, nx) * problem_.A;
//         if (i == 0)
//             Su.block(nx * (i + 1), 0, nx, nu) = problem_.B;
//         else {
//             Su.block(nx * (i + 1), 0, nx, nu) = problem_.A * Su.block(nx * i, 0, nx, nu);
//             Su.block(nx * (i + 1), nu, nx, nu * (N - 1)) = Su.block(nx * i, 0, nx, nu * (N - 1));
//         }
//         Qbar.block(nx * i, nx * i, nx, nx) = gains.Q;
//         Rbar.block(nu * i, nu * i, nu, nu) = gains.R;

//         Ulb.block(nu * i, 0, nu, 1) = input.lb;
//         Uub.block(nu * i, 0, nu, 1) = input.ub;

//         Xlb.block(nx * (i + 1), 0, nx, 1) = state.lb;
//         Xub.block(nx * (i + 1), 0, nx, 1) = state.ub;
//     }
//     Qbar.block(nx * N, nx * N, nx, nx) = gains.P;
//     // Cost function
//     H = Su.transpose() * Qbar * Su + Rbar;
//     F = Sx.transpose() * Qbar * Su;

//     // debug_print();
// }

// Eigen::Matrix<double, Eigen::Dynamic, 1> LinearMPCBase::update(const bool verbose, const Eigen::Matrix<double, Eigen::Dynamic, 1> x0) {
//     QP->data_.H = H;
//     QP->data_.lb = Ulb;
//     QP->data_.ub = Uub;
//     QP->data_.A = Su;
//     QP->data_.lbA = Xlb - Sx * x0;
//     QP->data_.ubA = Xub - Sx * x0;
//     QP->data_.g = 2 * F.transpose() * x0;

//     if (verbose) {
//         QP->print();
//         debug_print();
//         std::cout << "Sx:\n" << Sx << std::endl;
//         std::cout << "x0: " << x0.transpose() << std::endl;
//         std::cout << "Sx*x0: " << (Sx * x0).transpose() << std::endl;
//         std::cout << "Xlb: " << Xlb.transpose() << "\nXub: " << Xub.transpose() << std::endl;
//         std::cout << "Ulb: " << Ulb.transpose() << "\nUub: " << Uub.transpose() << std::endl;
//     }

//     QP->solve();
//     return QP->getOptimizer();
// }

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

}  // namespace nonlinear_control
