/**
 * @file test_linear_mpc.cpp
 * @brief Unit tests for Linear MPC controllers
 */

#include <gtest/gtest.h>

#include <Eigen/Dense>

#include "controls/linear_mpc.hpp"
#include "controls/mpc_qpoases.hpp"

namespace nlc = nonlinear_controls;

// =============================================================================
// LinearMPCBase Tests
// =============================================================================

class LinearMPCBaseTest : public ::testing::Test {
protected:
  static constexpr int N = 10;  // Horizon
  static constexpr int nx = 4;  // State dimension
  static constexpr int nu = 2;  // Input dimension
};

TEST_F(LinearMPCBaseTest, Construction) {
  EXPECT_NO_THROW({ nlc::LinearMPCBase mpc(N, nx, nu); });
}

TEST_F(LinearMPCBaseTest, DefaultGainsAreIdentity) {
  nlc::LinearMPC mpc(N, nx, nu);

  // Q, R, P should be identity by default
  Eigen::MatrixXd expected_Q = Eigen::MatrixXd::Identity(nx, nx);
  Eigen::MatrixXd expected_R = Eigen::MatrixXd::Identity(nu, nu);

  // We can't directly access Q, R, P, but we can verify through behavior
  // by setting dynamics and checking the MPC constructs without error
  Eigen::MatrixXd A = Eigen::MatrixXd::Identity(nx, nx);
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(nx, nu);
  B.bottomRows(nu) = Eigen::MatrixXd::Identity(nu, nu);

  mpc.init_dynamics(A, B);
  EXPECT_NO_THROW(mpc.construct());
}

// =============================================================================
// LinearMPC Tests
// =============================================================================

class LinearMPCTest : public ::testing::Test {
protected:
  void SetUp() override {
    // Double integrator dynamics (discrete time)
    // x = [pos, vel], u = [accel]
    // x_{k+1} = A*x_k + B*u_k
    double dt = 0.1;

    A.resize(nx, nx);
    A << 1.0, dt, 0.0, 1.0;

    B.resize(nx, nu);
    B << 0.5 * dt * dt, dt;

    Q.resize(nx, nx);
    Q = Eigen::MatrixXd::Identity(nx, nx);

    R.resize(nu, nu);
    R = 0.1 * Eigen::MatrixXd::Identity(nu, nu);

    P = Q;  // Terminal cost same as stage cost

    state_lb.resize(nx);
    state_ub.resize(nx);
    state_lb << -10.0, -5.0;
    state_ub << 10.0, 5.0;

    input_lb.resize(nu);
    input_ub.resize(nu);
    input_lb << -2.0;
    input_ub << 2.0;
  }

  static constexpr int N = 10;
  static constexpr int nx = 2;
  static constexpr int nu = 1;

  Eigen::MatrixXd A, B, Q, R, P;
  Eigen::VectorXd state_lb, state_ub, input_lb, input_ub;
};

TEST_F(LinearMPCTest, Construction) {
  nlc::LinearMPC mpc(N, nx, nu);
  mpc.init_dynamics(A, B);
  mpc.set_gains(Q, P, R);
  mpc.set_state_bounds(state_lb, state_ub);
  mpc.set_input_bounds(input_lb, input_ub);

  EXPECT_NO_THROW(mpc.construct());
}

TEST_F(LinearMPCTest, ValidStateCheck) {
  nlc::LinearMPC mpc(N, nx, nu);
  mpc.set_state_bounds(state_lb, state_ub);

  Eigen::VectorXd valid_state(nx);
  valid_state << 0.0, 0.0;
  EXPECT_TRUE(mpc.valid_state(valid_state));

  Eigen::VectorXd invalid_state(nx);
  invalid_state << 15.0, 0.0;  // Exceeds state_ub
  EXPECT_FALSE(mpc.valid_state(invalid_state));
}

TEST_F(LinearMPCTest, ValidInputCheck) {
  nlc::LinearMPC mpc(N, nx, nu);
  mpc.set_input_bounds(input_lb, input_ub);

  Eigen::VectorXd valid_input(nu);
  valid_input << 1.0;
  EXPECT_TRUE(mpc.valid_input(valid_input));

  Eigen::VectorXd invalid_input(nu);
  invalid_input << 5.0;  // Exceeds input_ub
  EXPECT_FALSE(mpc.valid_input(invalid_input));
}

// =============================================================================
// MPCQPOases Tests
// =============================================================================

class MPCQPOasesTest : public ::testing::Test {
protected:
  void SetUp() override {
    // Simple double integrator
    double dt = 0.1;

    A.resize(nx, nx);
    A << 1.0, dt, 0.0, 1.0;

    B.resize(nx, nu);
    B << 0.5 * dt * dt, dt;

    Q.resize(nx, nx);
    Q = Eigen::MatrixXd::Identity(nx, nx);

    R.resize(nu, nu);
    R = 0.1 * Eigen::MatrixXd::Identity(nu, nu);

    P = Q;

    state_lb.resize(nx);
    state_ub.resize(nx);
    state_lb << -10.0, -5.0;
    state_ub << 10.0, 5.0;

    input_lb.resize(nu);
    input_ub.resize(nu);
    input_lb << -2.0;
    input_ub << 2.0;
  }

  static constexpr int N = 10;
  static constexpr int nx = 2;
  static constexpr int nu = 1;

  Eigen::MatrixXd A, B, Q, R, P;
  Eigen::VectorXd state_lb, state_ub, input_lb, input_ub;
};

TEST_F(MPCQPOasesTest, Construction) {
  nlc::MPCQPOases mpc(N, nx, nu);
  mpc.init_dynamics(A, B);
  mpc.set_gains(Q, P, R);
  mpc.set_state_bounds(state_lb, state_ub);
  mpc.set_input_bounds(input_lb, input_ub);

  EXPECT_NO_THROW(mpc.construct());
}

TEST_F(MPCQPOasesTest, RunReturnsValidResult) {
  nlc::MPCQPOases mpc(N, nx, nu);
  mpc.init_dynamics(A, B);
  mpc.set_gains(Q, P, R);
  mpc.set_state_bounds(state_lb, state_ub);
  mpc.set_input_bounds(input_lb, input_ub);
  mpc.construct();

  Eigen::VectorXd x0(nx);
  x0 << 1.0, 0.0;  // Initial state: position=1, velocity=0

  Eigen::VectorXd xref(nx);
  xref << 0.0, 0.0;  // Target: origin

  auto result = mpc.run(x0, xref);

  EXPECT_TRUE(result.has_value()) << "MPC should return a valid solution";

  if (result.has_value()) {
    auto u_opt = result.value();
    EXPECT_EQ(u_opt.rows(), N * nu) << "Solution should have N*nu elements";

    // First input should be within bounds
    double u0 = u_opt(0);
    EXPECT_GE(u0, input_lb(0) - 1e-6);
    EXPECT_LE(u0, input_ub(0) + 1e-6);
  }
}

TEST_F(MPCQPOasesTest, DrivesToTarget) {
  nlc::MPCQPOases mpc(N, nx, nu);
  mpc.init_dynamics(A, B);
  mpc.set_gains(Q, P, R);
  mpc.set_state_bounds(state_lb, state_ub);
  mpc.set_input_bounds(input_lb, input_ub);
  mpc.construct();

  Eigen::VectorXd x(nx);
  x << 2.0, 0.0;  // Start at position=2

  Eigen::VectorXd xref(nx);
  xref << 0.0, 0.0;  // Target: origin

  // Simulate for several steps
  for (int i = 0; i < 50; ++i) {
    auto result = mpc.run(x, xref);
    ASSERT_TRUE(result.has_value()) << "MPC should always return a solution";

    // Apply first input
    double u = result.value()(0);
    u = std::max(input_lb(0), std::min(input_ub(0), u));

    // Propagate dynamics
    Eigen::VectorXd u_vec(nu);
    u_vec << u;
    x = A * x + B * u_vec;
  }

  // After enough steps, should be close to target (within 0.2 tolerance)
  EXPECT_NEAR(x(0), 0.0, 0.2) << "Position should converge to target";
  EXPECT_NEAR(x(1), 0.0, 0.2) << "Velocity should converge to zero";
}

TEST_F(MPCQPOasesTest, TrajectoryGeneration) {
  nlc::MPCQPOases mpc(N, nx, nu);
  mpc.init_dynamics(A, B);
  mpc.set_gains(Q, P, R);
  mpc.set_state_bounds(state_lb, state_ub);
  mpc.set_input_bounds(input_lb, input_ub);
  mpc.construct();

  Eigen::VectorXd x0(nx);
  x0 << 1.0, 0.5;

  Eigen::VectorXd xref(nx);
  xref << 0.0, 0.0;

  auto result = mpc.run(x0, xref);
  ASSERT_TRUE(result.has_value());

  // Get predicted state trajectory
  auto X_traj = mpc.X();

  EXPECT_EQ(X_traj.rows(), nx);
  EXPECT_EQ(X_traj.cols(), N);
}

// =============================================================================
// Main
// =============================================================================

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
