#include "mpc_base.h"
namespace nonlinear_control {

template <typename T> 
MPCBase<T>::MPCBase(const bool _isLTI, const int _N, const int _nx, const int _nu) : isLTI(_isLTI), nx(_nx), nu(_nu), N(_N){
    problem_.N = N;
    problem_.nx = nx;
    problem_.nu = nu;

    problem_.A.resize(nx, nx);
    problem_.B.resize(nx, nu);
    state.lb.resize(nx,1);
    state.ub.resize(nu,1);
    input.lb.resize(nu,1);
    input.ub.resize(nu,1);

    Ax.resize(2*nx, nx);
    Af.resize(2*nx, nx);
    bx.resize(2*nx, 1);
    bf.resize(2*nx, 1);

    Inx.resize(nx, nx); Inx.setIdentity();
    Inu.resize(nu, nu); Inu.setIdentity();
    Onx.resize(nx, nx); Onx.setZero();
    Onu.resize(nu, nu); Onu.setZero();

}

template <typename T> 
MPCBase<T>::~MPCBase() {}


template <typename T> 
void MPCBase<T>::setup() {

    /*********************************/
    /* setting-up the qp for the MPC */
    /*********************************/
    nVars = N*nu+nx;
    nCons = N*nx;

    // generating the data
    // Using the MPC notation from ME231A
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Sx, Su, Qbar, Rbar, H, F;
    Sx.resize(nx*(N+1), nx);  Sx.setZero();
    Su.resize(nx*(N+1), N*nu);Su.setZero();


    Sx.block<nx, nx>(0,0) << Inx;
    // for (int i = 0; i < N; ++i) {
    //     Sx.block<nx, nx>(nx*(i+1), 0) =  Sx.block<nx, nx>(nx*i, 0)*problem_.A;
    // }


    qp = new QPOasesEigen(nVars, nCons);
    qp->options.setToMPC( );
    qp->options.printLevel = qpOASES::PL_LOW;
    qp->setup();

    // // Gains
    // qp->data_.H
}

// template <typename T>
// void MPCBase<T>::update() {

// }

template class MPCBase<float>;
template class MPCBase<double>;

} // namespace nonlinear_control
