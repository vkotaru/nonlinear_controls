#include <eigen3/Eigen/Dense>
#include <iostream>
#include <vector>

#include "nonlinear_controls.h"
#include "matplotlibcpp.h"
#include "epigraph.hpp"

namespace plt = matplotlibcpp;
namespace nlc = nonlinear_controls;

enum class CTRL {
  QP_OASES,
  EPIGRAPH_OSQP
};

class EpigraphMPC {
protected:

  int N, nx, nu;
  int nVars, nCons;

  Eigen::MatrixXd Q_, P_, R_; // mpc gains
  Eigen::MatrixXd A_, B_;    // dynamics
  Eigen::MatrixXd Sx, Su, Qbar, Rbar, H, F;
  Eigen::MatrixXd Ulb, Uub, Xlb, Xub; // bounds

  // CVX problem setup
  cvx::OptimizationProblem qp;
  cvx::MatrixX x, u;

  Eigen::VectorXd x0;
  Eigen::VectorXd xd; // goal state

  cvx::osqp::OSQPSolver *solver;

public:
  EpigraphMPC(const int &_N, const int &_nx, const int &_nu) : N(_N), nx(_nx), nu(_nu) {
    nVars = N * nu;
    nCons = (N + 1) * nx;

    // Create variables
    x = qp.addVariable("x", nx, N + 1);
    u = qp.addVariable("u", nu, N);
    x0.resize(nx);
    xd.resize(nx);
  }

  ~EpigraphMPC() = default;

  void set_mpc_gains(const Eigen::MatrixXd &Q, const Eigen::MatrixXd &P, const Eigen::MatrixXd &R) {
    assert(Q.rows() == nx && Q.cols() == nx);
    assert(P.rows() == nx && P.cols() == nx);
    assert(R.rows() == nu && R.cols() == nu);
    this->Q_ = Q;
    this->P_ = P;
    this->R_ = R;
  }

  void set_input_bounds(const Eigen::VectorXd &lb, const Eigen::VectorXd &ub) {
    assert(lb.rows() == nu && lb.cols() == 1);
    assert(ub.rows() == nu && ub.cols() == 1);
    Ulb = lb;
    Uub = ub;
  }

  void set_state_bounds(const Eigen::VectorXd &lb, const Eigen::VectorXd &ub) {
    assert(lb.rows() == nx && lb.cols() == 1);
    assert(ub.rows() == nx && ub.cols() == 1);
    Xlb = lb;
    Xub = ub;
  }

  void init_dynamics(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B) {
    update_dynamics(A, B);
  }

  void update_dynamics(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B) {
    this->A_ = A;
    this->B_ = B;
  }

  void construct() {
    // Dynamics
    for (size_t t = 0; t < N; t++) {
      // dynamics
      qp.addConstraint(cvx::equalTo(x.col(t + 1), cvx::par(A_) * x.col(t) + cvx::par(B_) * u.col(t)));

      // cost
      qp.addCostTerm((x.col(t) - cvx::par(xd)).transpose() * cvx::par(Q_) * (x.col(t) - cvx::par(xd)));
      qp.addCostTerm(u.col(t).transpose() * cvx::par(R_) * u.col(t));


      // State and control limits
//      qp.addConstraint(cvx::box(cvx::par(Xlb), x.col(t), cvx::par(Xub)));
//      qp.addConstraint(cvx::box(cvx::par(Ulb), u.col(t), cvx::par(Uub)));

//      qp.addConstraint(cvx::greaterThan(x.col(t), cvx::par(Xlb)));
//      qp.addConstraint(cvx::lessThan(x.col(t), cvx::par(Xub)));
//      qp.addConstraint(cvx::greaterThan(u.col(t), cvx::par(Ulb)));
//      qp.addConstraint(cvx::lessThan(u.col(t), cvx::par(Uub)));
    }
    qp.addCostTerm((x.col(N) - cvx::par(xd)).transpose() * cvx::par(P_) * (x.col(N) - cvx::par(xd)));
//    qp.addConstraint(cvx::box(cvx::par(Xlb), x.col(N), cvx::par(Xub)));


    // Boundary constraints
    qp.addConstraint(cvx::equalTo(x.col(0), cvx::dynpar(x0)));
    // qp.addConstraint(equalTo(x.col(T), 0.)); // to be free

    solver = new cvx::osqp::OSQPSolver(qp);
//    solver->setAlpha(1.0);
  }

  Eigen::MatrixXd run(const Eigen::VectorXd &x0, const Eigen::VectorXd &xd) {

    this->x0 = x0;
    this->xd = xd;
    solver->solve(false);
    Eigen::MatrixXd x_sol = eval(x);
    Eigen::MatrixXd u_sol = eval(u);

    return u_sol;
  }

};

void run_simulation() {

  nlc::Logger::INFO("Testing 3D position MPC with LTI Dynamics... ");

  //
  // simulation parameters
  //
  double h = 1. / 200.;
  const double simulate_for_seconds = 10;
  int MAX_ITER_STEPS = int(simulate_for_seconds / h);

  //
  // dynamics
  //
  nlc::PointMass3D dynamics_{h};

  //
  // controller
  //
  int nx = 6, nu = 3, N = 10;
//
  CTRL ctrl_type = CTRL::QP_OASES;
  nlc::LinearMPC controller_{N, nx, nu}; // 0.000142

//  CTRL ctrl_type = CTRL::EPIGRAPH_OSQP;
//  EpigraphMPC controller_{N, nx, nu}; // 0.0004

  controller_.init_dynamics(dynamics_.A(), dynamics_.B());
  // gains
  nlc::MatrixXd Q, P, R;
  Q.resize(nx, nx);
  P.resize(nx, nx);
  R.resize(nu, nu);
  Q << 100 * Eigen::Matrix3d::Identity(), Eigen::Matrix3d::Zero(),
      Eigen::Matrix3d::Zero(), 10 * Eigen::Matrix3d::Identity();
  P << 8.1314, 0.0000, -0.0000, 0.6327, -0.0000, -0.0000, 0.0000, 8.1314,
      0.0000, 0.0000, 0.6327, 0.0000, -0.0000, 0.0000, 8.1314, -0.0000, 0.0000,
      0.6327, 0.6327, 0.0000, -0.0000, 0.2606, -0.0000, -0.0000, -0.0000,
      0.6327, 0.0000, -0.0000, 0.2606, 0.0000, -0.0000, 0.0000, 0.6327, -0.0000,
      0.0000, 0.2606;
  P = 1e2 * P;
  R = 1 * Eigen::Matrix<double, 3, 3>::Identity();
  controller_.set_mpc_gains(Q, P, R);
  // state & input bounds
  Eigen::Matrix<double, 6, 1> state_lb, state_ub;
  Eigen::Matrix<double, 3, 1> input_lb, input_ub;
  state_lb << -2., -2.5, 0., -2, -2, -2;
  state_ub << 2., 2.5, 4., 2., 2., 2;
  input_lb = -20. * Eigen::Matrix<double, 3, 1>::Ones();
  input_ub = -input_lb;

  controller_.set_input_bounds(input_lb, input_ub);
  controller_.set_state_bounds(state_lb, state_ub);
//  controller_.set_cpu_time_limit(h);

  // construct the mpc
  controller_.construct();

  //
  // initialize simulation
  //
  Eigen::Matrix<double, 6, 1> goal_state{}, state{};
  goal_state << 1.0, 2.0, 1.5, 0.0, 0.0, 0.0;

  std::cout << "\n X0: " << state.transpose() << std::endl;
  std::cout << "\n Xgoal: " << goal_state.transpose() << std::endl;

  nlc::MatrixXd K, zOpt, uOpt;
  uOpt.resize(nu, 1);
  zOpt.resize(N * nu, 1);
  K.resize(nu, nx);
  K << 30.6287 * Eigen::Matrix3d::Identity(),
      12.4527 * Eigen::Matrix3d::Identity();

  double freq = 0;

  dynamics_.clear_log_vars();
  dynamics_.init(state, true, 0);

  sleep(4);

  double avg_time_elapsed = 0;
  for (int j = 0; j < MAX_ITER_STEPS; ++j) {

//    if (ctrl_type == CTRL::QP_OASES)
    controller_.set_state_bounds(state_lb - goal_state, state_ub - goal_state);

    /// LQR (for sanity check)
    auto uPD = -K * (state - goal_state);

    /// running mpc controller
    auto start = nlc::utils::get_current_time();
//    if (ctrl_type == CTRL::QP_OASES) {
    zOpt = controller_.run((state - goal_state));
//    }
//    else
//    {
//      zOpt = controller_.run(state, goal_state);
//    }
    uOpt = zOpt.block(0, 0, nu, 1);
    auto end = nlc::utils::get_current_time();
    auto time_elapsed = (double(end - start) * 1e-6);
    avg_time_elapsed += time_elapsed;

    // integrating dynamics using the input at N = 0
    std::cout << "PD Input:>>> " << uPD.transpose() << std::endl;
    std::cout << "MPC Input:>>> " << uOpt.transpose() << std::endl;

    uOpt += GRAVITY_VECTOR; // TODO multiply with mass
    dynamics_.step(uOpt, true, time_elapsed);
    state = dynamics_.state();
  }
  avg_time_elapsed = avg_time_elapsed / MAX_ITER_STEPS;
  nlc::Logger::WARN("Average time elapsed: " + std::to_string(avg_time_elapsed));

  dynamics_.plot();

}

int main(int argc, char *argv[]) {
  run_simulation();

  //
  // Incomplete solutions using Epigraph
  //
}
