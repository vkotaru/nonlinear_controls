#include "controls/linear_mpc.h"
namespace nonlinear_controls {
//////////////////////////////////////////
/// Linear Time Invariant Dynamics MPC
//////////////////////////////////////////
template<typename T>
LinearMPC<T>::LinearMPC(const int &_N, const int &_nx, const int &_nu) : N(_N), nx(_nx), nu(_nu) {
  nVars = N * nu;
  nCons = (N + 1) * nx;

  /// QP solver Setup
  solver = new qpOASES::SQProblem(nVars, nCons);
  options.setToMPC();
  options.printLevel = qpOASES::PL_LOW;
  solver->setOptions(options);
  nWSR = 1e7;

  /// Sizing all Eigen Dynamic Matrices
  A.resize(nx, nx);
  B.resize(nx, nu);

  Sx.resize(nx * (N + 1), nx);
  Sx.setZero();
  Su.resize(nx * (N + 1), N * nu);
  Su.setZero();
  Qbar.resize(nx * (N + 1), nx * (N + 1));
  Qbar.setIdentity();
  Rbar.resize(nu * N, nu * N);
  Rbar.setIdentity();

  Q.resize(nx, nx);
  Q.setIdentity();
  R.resize(nu, nu);
  R.setIdentity();
  P.resize(nx, nx);
  P.setIdentity();

  H.resize(nVars, nVars);
  H.setIdentity();
  F.resize(nx, nVars);
  F.setZero();

  Ulb.resize(nu * N, 1);
  Ulb.setOnes();
  Ulb = -1e6 * Ulb;
  Uub.resize(nu * N, 1);
  Uub.setOnes();
  Uub = 1e6 * Uub;
  Xlb.resize(nx * (N + 1), 1);
  Xlb.setOnes();
  Xlb = -1e6 * Xlb;
  Xub.resize(nx * (N + 1), 1);
  Xub.setOnes();
  Xub = 1e6 * Xub;
}

template<typename T>
LinearMPC<T>::~LinearMPC() = default;

template<typename T>
void LinearMPC<T>::set_mpc_gains(MatrixX<T> Q, MatrixX<T> P, MatrixX<T> R) {
  assert(Q.rows() == nx && Q.cols() == nx);
  assert(P.rows() == nx && P.cols() == nx);
  assert(R.rows() == nu && R.cols() == nu);
  this->Q = Q;
  this->P = P;
  this->R = R;
  for (int i = 0; i < N; ++i) {
    this->Qbar.block(nx * i, nx * i, nx, nx) = Q;
    this->Rbar.block(nu * i, nu * i, nu, nu) = R;
  }
  this->Qbar.block(nx * N, nx * N, nx, nx) = P;
}

template<typename T>
void LinearMPC<T>::set_input_bounds(VectorX<T> lb, VectorX<T> ub) {
  assert(lb.rows() == nu && lb.cols() == 1);
  assert(ub.rows() == nu && ub.cols() == 1);
  Ulb = lb.replicate(N, 1);
  Uub = ub.replicate(N, 1);

  /// set only once
  data_.lb = Ulb.template cast<double>();
  data_.ub = Uub.template cast<double>();
}

template<typename T>
void LinearMPC<T>::set_state_bounds(VectorX<T> lb, VectorX<T> ub) {
  assert(lb.rows() == nx && lb.cols() == 1);
  assert(ub.rows() == nx && ub.cols() == 1);
  Xlb = lb.replicate(N + 1, 1);
  Xub = ub.replicate(N + 1, 1);
}

template<typename T>
void LinearMPC<T>::init_dynamics(MatrixX<T> A, MatrixX<T> B) {
  LinearMPC<T>::update_dynamics(A, B);
}

template<typename T>
void LinearMPC<T>::update_dynamics(MatrixX<T> A, MatrixX<T> B) {
  this->A = A;
  this->B = B;
}

template<typename T>
void LinearMPC<T>::construct() {
  ////////////////////////////////////////
  /// Note: remember to call this function
  /// if dynamics, input bounds or
  /// MPC gains are updated
  ////////////////////////////////////////
  Sx.block(0, 0, nx, nx).setIdentity();
  Su.block(0, 0, nx, N * nu).setZero();

  for (int i = 0; i < N; ++i) {
    Sx.block(nx * (i + 1), 0, nx, nx) =
        Sx.block(nx * i, 0, nx, nx) * A;
    if (i == 0)
      Su.block(nx * (i + 1), 0, nx, nu) = B;
    else {
      Su.block(nx * (i + 1), 0, nx, nu) = A * Su.block(nx * i, 0, nx, nu);
      Su.block(nx * (i + 1), nu, nx, nu * (N - 1)) = Su.block(nx * i, 0, nx, nu * (N - 1));
    }
  }
  // Cost function
  H = Su.transpose() * Qbar * Su + Rbar;
  F = Sx.transpose() * Qbar * Su;

  /// set only once
  data_.H = H.template cast<double>();
  data_.lb = Ulb.template cast<double>();
  data_.ub = Uub.template cast<double>();
  data_.A = Su.template cast<double>();
  QP_INITIALIZED = false;
}


template<typename T>
void LinearMPC<T>::print() {
  std::cout << "--------------------------------------------------"
            << std::endl;
  std::cout << "*           QP setup       *" << std::endl;
  std::cout << "--------------------------------------------------"
            << std::endl;
  std::cout << "H: \n"
            << data_.H << "\nf: \n"
            << data_.g.transpose() << std::endl;
  std::cout << "\nlbA: \n"
            << data_.lbA.transpose() << "\nA: \n"
            << data_.A << "\nubA: \n"
            << data_.ubA.transpose() << std::endl;
  std::cout << "\nlb: " << data_.lb.transpose()
            << "\nub: " << data_.ub.transpose() << std::endl;
  std::cout << "-------------------*****---------------------------"
            << std::endl;
}


template<typename T>
void LinearMPC<T>::updateCArrays(const VectorX<T> &x0) {
  data_.lbA = (Xlb - Sx * x0).template cast<double>();
  data_.ubA = (Xub - Sx * x0).template cast<double>();
  data_.g = (2 * F.transpose() * x0).template cast<double>();

  print();

  Hc = data_.H.transpose().data();
  Ac = data_.A.transpose().data();
  gc = data_.g.transpose().data();
  lbc = data_.lb.transpose().data();
  ubc = data_.ub.transpose().data();
  lbAc = data_.lbA.transpose().data();
  ubAc = data_.ubA.transpose().data();
}

template<typename T>
VectorX<T> LinearMPC<T>::init(const VectorX<T> &x0) {
  this->updateCArrays(x0);
  solver->init(Hc, gc, Ac, lbc, ubc, lbAc, ubAc, nWSR, 0);
}

template<typename T>
VectorX<T> LinearMPC<T>::run(const VectorX<T> &x0) {

  this->updateCArrays(x0);

//  double *cpu_time = 0.01;
  qpOASES::returnValue sol = solver->hotstart(Hc, gc, Ac, lbc, ubc, lbAc, ubAc, nWSR, 0);

  qpOASES::real_t xOpt[nVars];
  solver->getPrimalSolution(xOpt);
  if (sol == qpOASES::SUCCESSFUL_RETURN) {
    std::cout << "SUCCESSFUL_RETURN" << std::endl;
  } else if (sol == qpOASES::RET_MAX_NWSR_REACHED) {
    std::cout << "RET_MAX_NWSR_REACHED" << std::endl;
  } else if (sol == qpOASES::RET_INIT_FAILED) {
    std::cout << "RET_INIT_FAILED" << std::endl;
  } else {
    std::cout << "Something else" << std::endl;
  }

  Eigen::Matrix<T, Eigen::Dynamic, 1> xOptVec;
  xOptVec =
      Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>>(xOpt, nVars, 1).template cast<T>();
  return xOptVec;
}
template<typename T>
VectorX<T> LinearMPC<T>::run(const VectorX<T> &x0, const MatrixX<T> &A, const MatrixX<T> &B) {
  return LinearMPC<T>::run(x0);
}

//////////////////////////////////////////
/// Linear Time Varying Dynamics MPC
/////////////////////////////////////////
template<typename T>
LinearMPCt<T>::LinearMPCt(const int &_N, const int &_nx, const int &_nu) : LinearMPC<T>(_N, _nx, _nu) {

}

template<typename T>
LinearMPCt<T>::~LinearMPCt() = default;

template<typename T>
void LinearMPCt<T>::init_dynamics(MatrixX<T> A, MatrixX<T> B) {
  Astorage.push_back(A);
  Bstorage.push_back(B);
}

template<typename T>
void LinearMPCt<T>::update_dynamics(MatrixX<T> A, MatrixX<T> B) {
  Astorage.pop_front(); // removing A_{k} or is it k-1? (whatever!)
  Astorage.push_back(A); // adding A_{k+N}
  Bstorage.pop_front();
  Bstorage.push_back(B);
}

template<typename T>
void LinearMPCt<T>::construct() {
  this->Sx.block(0, 0, this->nx, this->nx).setIdentity();
  this->Su.block(0, 0, this->nx, this->N * this->nu).setZero();

  for (int i = 0; i < this->N; ++i) {
    this->Sx.block(this->nx * (i + 1), 0, this->nx, this->nx) = Astorage.at(i) * this->Sx.block(this->nx * i, 0, this->nx, this->nx);
    this->Su.block(this->nx * (i + 1), this->nu * i, this->nx, this->nu) = Bstorage.at(i);
    if (i != 0) {
      this->Su.block(this->nx * (i + 1), 0, this->nx, this->nu * i) = Astorage.at(i) * this->Su.block(this->nx * i, 0, this->nx, this->nu * i);
    }
  }
  // Cost function
  this->H = this->Su.transpose() * this->Qbar * this->Su + this->Rbar;
  this->F = this->Sx.transpose() * this->Qbar * this->Su;

  /// set only once
  this->data_.H = this->H.template cast<double>();
  this->data_.lb = this->Ulb.template cast<double>();
  this->data_.ub = this->Uub.template cast<double>();
  this->data_.A = this->Su.template cast<double>();
  this->QP_INITIALIZED = false;
}

template<typename T>
VectorX<T> LinearMPCt<T>::run(const VectorX<T> &x0, const MatrixX<T> &A, const MatrixX<T> &B) {
  /// update:
  ///  constraints bounds
  ///  linear cost
  this->data_.lbA = (this->Xlb - this->Sx * x0).template cast<double>();
  this->data_.ubA = (this->Xub - this->Sx * x0).template cast<double>();
  this->data_.g = (2 * this->F.transpose() * x0).template cast<double>();

  /// solve the QP
//  this->solve();

  /// update dynamics with A_{k+N}, B_{k+N}
  update_dynamics(A, B);
  /// reconstruct the matrices
  construct();

  /// return
//  return this->QP->getOptimizer().template cast<T>();
}

///////////////////////////////////
/// templates
///////////////////////////////////
template class LinearMPC<float>;
template class LinearMPC<double>;
template class LinearMPCt<float>;
template class LinearMPCt<double>;

}  // namespace nonlinear_controls

