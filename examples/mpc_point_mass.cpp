#include <eigen3/Eigen/Dense>
#include <iostream>
#include <vector>

#include "nonlinear_controls.h"
// #include "matplotlibcpp.h"
#include "epigraph.hpp"

// namespace plt = matplotlibcpp;
namespace nlc = nonlinear_controls;

enum class SOLVER_TYPE {
  EPIGRAPH = 0,
  QPOASES,
  COUNT
};

void run_simulation(const int horizon, SOLVER_TYPE solve_type_) {

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
  int nx = 6, nu = 3, N = horizon;

  nlc::LinearMPCBase *controller_;

  if (solve_type_ == SOLVER_TYPE::QPOASES) {
    nlc::Logger::STATUS("Initializing qpOASES mpc solver");
    controller_ = new nlc::MPCQPOases(N, nx, nu); // 0.000142
  } else if (solve_type_ == SOLVER_TYPE::EPIGRAPH) {
    nlc::Logger::STATUS("Initializing Epigraph mpc solver");
    controller_ = new nlc::MPCEpigraph(N, nx, nu);
  } else {
    nlc::Logger::ERROR("Unknown solver type! exiting!");
    return;
  }

  controller_->init_dynamics(dynamics_.A(), dynamics_.B());

  nlc::VectorXd q, r;
  q.resize(nx);
  r.resize(nu);
  q << 20, 20, 10, 20, 20, 20;
  r << 1, 1, 1;
  // gains
  nlc::MatrixXd Q, P, R;
  Q.resize(nx, nx);
  P.resize(nx, nx);
  R.resize(nu, nu);
  Q = q.asDiagonal();
  P << 4.8221, 0, 0, 0.8945, 0, 0,
      0, 4.8221, 0, 0, 0.8945, 0,
      0, 0, 3.2500, 0, 0, 0.8945,
      0.8945, 0, 0, 1.0861, 0, 0,
      0, 0.8945, 0, 0, 1.0861, 0,
      0, 0, 0.8945, 0, 0, 1.0861;
  P = 1e3 * P;
  R = r.asDiagonal();
  controller_->set_gains(Q, P, R);
  // state & input bounds
  Eigen::Matrix<double, 6, 1> state_lb, state_ub;
  Eigen::Matrix<double, 3, 1> input_lb, input_ub;

  state_lb << -2., -2.5, 0., -2, -2, -2;
  state_ub << 2., 2.5, 4., 2., 2., 2;

  input_lb = -G_SI * Eigen::Matrix<double, 3, 1>::Ones();
  input_ub = -input_lb;

//  state_lb = -5*Eigen::Matrix<double, 6, 1>::Ones();
//  state_ub = 5*Eigen::Matrix<double, 6, 1>::Ones();

//  controller_->set_input_bounds(input_lb, input_ub);
  controller_->set_state_bounds(state_lb, state_ub);
//  controller_->set_cpu_time_limit(h);

  // construct the mpc
  controller_->construct();

  nlc::MatrixXd K, uOpt;
  uOpt.resize(nu, 1);
//  zOpt.resize(N * nu, 1);
  K.resize(nu, nx);
  K << 30.6287 * Eigen::Matrix3d::Identity(),
      12.4527 * Eigen::Matrix3d::Identity();

  Eigen::Matrix<double, 6, 10> x0s;
  x0s.col(0) << -1, -1, -1, 0, 0, 0;
  x0s.col(1) << -0.012834, 0.94555, 0.414966, 0.542715, 0.05349, 0.539828;
  for (int i = 0; i < 7; i++) {
    x0s.col(i + 2) = Eigen::Matrix<double, 6, 1>::Random();
    x0s(2, i + 2) = std::abs(x0s(2, i + 2));
  }
  x0s.col(9) << 0, 0, -1, 0, 0, 0;

  for (int k = 0; k < 10; k++) {
    //
    // initialize simulation
    //
    Eigen::Matrix<double, 6, 1> goal_state{}, state{};

    state = x0s.col(k);
    goal_state << 0., 0.0, 0, 0.0, 0.0, 0.0;

    /// qpOASES fails for this condition
    ///  state << -0.012834  , 0.94555 ,-0.414966, 0., 0., 0;

    std::cout << "x0: " << state.transpose() << std::endl;
    std::cout << "xgoal: " << goal_state.transpose() << std::endl;

    dynamics_.clear_log_vars();
    dynamics_.init(state, true, 0);

    double avg_time_elapsed = 0;
    double max_time_elapsed = -INFY;
    for (int j = 0; j < MAX_ITER_STEPS; ++j) {
      /// LQR (for sanity check)
      auto uPD = -K * (state - goal_state);

      /// running mpc controller
      auto start = nlc::utils::get_current_time();
      auto zOpt = controller_->run(state, goal_state);
      if (zOpt.has_value()) {
        uOpt = zOpt.value().block(0, 0, nu, 1);
      } else {
        break;
      }
      auto end = nlc::utils::get_current_time();
      auto time_elapsed = (double(end - start) * 1e-6);
      avg_time_elapsed += time_elapsed;
      max_time_elapsed = std::max(time_elapsed, max_time_elapsed);

      // integrating dynamics using the input at N = 0
//    std::cout << "PD Input:>>> " << uPD.transpose() << std::endl;
//    std::cout << "MPC Input:>>> " << uOpt.transpose() << std::endl;
//    std::cout << get_predicted_path(state, goal_state, zOpt) << std::endl;

      uOpt += GRAVITY_VECTOR; // TODO multiply with mass
      dynamics_.step(uOpt, true, time_elapsed);
      state = dynamics_.state();
    }
    avg_time_elapsed = avg_time_elapsed / MAX_ITER_STEPS;
    nlc::Logger::WARN(
        "Average time elapsed: " + std::to_string(avg_time_elapsed) + " max time elapsed: "
            + std::to_string(max_time_elapsed));
    dynamics_.plot();
    controller_->reset();
  }

}

int main(int argc, char *argv[]) {

  int N{5};
  SOLVER_TYPE solver_type_{SOLVER_TYPE::EPIGRAPH};

  auto get_option = [&](const int i) {
    if (std::string("N") == argv[i]) {
      N = std::stoi(argv[i + 1]);
    }
    if (std::string("solver") == argv[i]) {
      switch (std::stoi(argv[i + 1])) {
      case 0:solver_type_ = SOLVER_TYPE::EPIGRAPH;
        break;
      case 1:solver_type_ = SOLVER_TYPE::QPOASES;
        break;
      default:solver_type_ = SOLVER_TYPE::EPIGRAPH;
        break;
      }
    }
  };

  if (argc == 3) {
    get_option(1);
  } else if (argc == 5) {
    get_option(1);
    get_option(3);
  }

  std::cout << "N: " << N << " solver type " << int(solver_type_) << std::endl;
  run_simulation(N, solver_type_);
}
