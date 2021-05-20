#include "controls/SO3_vblmpc.h"

namespace nonlinear_controls {

template <typename T> SO3VblMPC<T>::SO3VblMPC(bool islti, int _N, double _dt) {
  // temporary inertia matrix
  Eigen::Matrix<T, 3, 3> J;
  J << 112533, 0, 0, 0, 36203, 0, 0, 0, 42673;
  J = J * 1e-6;
  SO3VblMPC(islti, _N, _dt, J);
}

template <typename T>
SO3VblMPC<T>::SO3VblMPC(bool islti, int _N, double _dt,
                        Eigen::Matrix<T, 3, 3> _J)
    : IS_LTI(islti), dt(_dt), N(_N) {
  nx = 6;
  nu = 3;
  inertia_ = _J;
  inertia_inv_ = inertia_.inverse();

  if (!IS_LTI) {
    std::cout << "IS LTV is not yet implemented" << std::endl;
  }

  if (verbose) {
    std::cout << "//////////////////////////////" << std::endl;
    std::cout << "SO3VblMPC constructor" << std::endl;
    std::cout << "inertia\n" << inertia_ << std::endl;
    std::cout << "inertia_inv\n" << inertia_inv_ << std::endl;
  }

  if (IS_LTI) {
    mpcSolver = new LinearMPC<T>(N, nx, nu);
  } else {
    mpcSolver = new LinearMPCt<T>(N, nx, nu);
  }

  // initialize dynamics
  Eigen::Matrix<T, 3, 1> Omd;
  Omd.setZero();
  A.resize(nx, nx);
  B.resize(nx, nu);
  A.setIdentity();
  generate_dynamics(Omd, A);
  B << Eigen::Matrix<T, 3, 3>::Zero(), inertia_inv_;
  B = B * dt;
  if (IS_LTI) {
    mpcSolver->init_dynamics(A, B);
  }
  if (verbose) {
    std::cout << "//////////////////////////////" << std::endl;
    std::cout << "dt " << dt << std::endl;
    std::cout << "Initializing dynamics ... " << std::endl;
    std::cout << "A: \n" << A << std::endl;
    std::cout << "B: \n" << B << std::endl;
  }
  /// gains
  Q.resize(nx, nx);
  P.resize(nx, nx);
  R.resize(nu, nu);
  Q << 1000 * Eigen::Matrix<T, 3, 3>::Identity(),
      Eigen::Matrix<T, 3, 3>::Zero(), Eigen::Matrix<T, 3, 3>::Zero(),
      100 * Eigen::Matrix<T, 3, 3>::Identity();
  P << 6.6514, -0.0000, -0.0000, 0.0894, -0.0000, -0.0000, -0.0000, 6.5123,
      -0.0000, -0.0000, 0.0440, -0.0000, -0.0000, -0.0000, 6.5231, -0.0000,
      -0.0000, 0.0475, 0.0894, -0.0000, -0.0000, 0.0293, -0.0000, -0.0000,
      -0.0000, 0.0440, -0.0000, -0.0000, 0.0141, -0.0000, -0.0000, -0.0000,
      0.0475, -0.0000, -0.0000, 0.0152;
  P = 1e4 * P;
  R = 1 * Eigen::Matrix<T, 3, 3>::Identity();
  mpcSolver->set_mpc_gains(Q, P, R);

  /// bounds
  input_lb.resize(nu, 1);
  input_ub.resize(nu, 1);
  state_lb.resize(nx, 1);
  state_ub.resize(nx, 1);
  input_lb << -10, -10, -10;
  input_ub << 10, 10, 10;
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

template <typename T> SO3VblMPC<T>::~SO3VblMPC() = default;

template <typename T>
void SO3VblMPC<T>::generate_dynamics(Eigen::Matrix<T, 3, 1> Omd,
                                     MatrixX<T> &A) {
  A << -utils::hat<T>(Omd), Eigen::Matrix<T, 3, 3>::Identity(),
      Eigen::Matrix<T, 3, 3>::Zero(),
      inertia_inv_ *
          (utils::hat<T>(inertia_ * Omd) - utils::hat<T>(Omd) * inertia_);
  A = Eigen::Matrix<T, 6, 6>::Identity() + A * dt;
}

template <typename T>
void SO3VblMPC<T>::init_dynamics(Eigen::Matrix<T, 3, 1> Omd) {
  generate_dynamics(Omd, A);
  mpcSolver->init_dynamics(A, B);
  if (verbose) {
    std::cout << "////////////////////////////////////\n"
                 "Omega "
              << Omd.transpose() << std::endl;
    std::cout << "A\n" << A << std::endl;
    std::cout << "B\n" << B << std::endl;
  }
}
template <typename T> void SO3VblMPC<T>::init_dynamics(TSO3<T> attd) {
  generate_dynamics(attd.Omega, A);
  mpcSolver->init_dynamics(A, B);
}

template <typename T>
void SO3VblMPC<T>::run(T dt, TSO3<T> x, TSO3<T> xd, Eigen::Matrix<T, 3, 1> &u) {
  // time-varying trajectory with only operating point used for the dynamics
  // for the next N steps
  auto x0 = x - xd;
  if (verbose) {
    std::cout << "x\n" << std::endl;
    x.print();
    std::cout << "xd\n" << std::endl;
    xd.print();
    std::cout << "x0: " << x0.transpose() << std::endl;
  }

  if (IS_LTI) {
    mpcSolver->print();
    uOpt = mpcSolver->run(x0);

    if (verbose) {
      std::cout << "uOpt " << uOpt.transpose() << std::endl;
    }
  } else {
    generate_dynamics(xd.Omega, A);
    uOpt = mpcSolver->run(x - xd, A, B);
  }
  u = uOpt.block(0, 0, 3, 1);
  u += x.Omega.cross(inertia_ * x.Omega);
  u += -inertia_ * (x.Omega.cross(x.R.transpose() * xd.R * xd.Omega) -
                    x.R.transpose() * xd.R * xd.dOmega);
}

template <typename T>
void SO3VblMPC<T>::set_gains(const MatrixX<T> &_Q, const MatrixX<T> &_P,
                             const MatrixX<T> &_R) {
  this->Q = Q;
  this->P = P;
  this->R = R;

  std::cout << "SO3VblMPC: setting gains\n Q:\n" << std::endl;
  std::cout << _Q << std::endl;
  std::cout << "P:\n";
  std::cout << _P << std::endl;
  std::cout << "R:\n" << std::endl;
  std::cout << _R << std::endl;

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
template class SO3VblMPC<float>;
template class SO3VblMPC<double>;
} // namespace nonlinear_controls
