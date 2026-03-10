/**
 * @file test_qp_solvers.cpp
 * @brief Unit tests for QP solvers (qpOASES, QPSwift)
 */

#include <gtest/gtest.h>
#include <Eigen/Dense>

#include "common/qpoases_eigen.hpp"
#include "quadprog/qpswift_eigen.h"

namespace nlc = nonlinear_controls;

// =============================================================================
// QPOasesEigen Tests
// =============================================================================

class QPOasesEigenTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Simple QP: min 0.5 x^T H x + g^T x
        // s.t. lb <= x <= ub
        //      lbA <= Ax <= ubA
    }
};

TEST_F(QPOasesEigenTest, SimpleUnconstrainedQP) {
    // min 0.5 * x^T * I * x + 0^T * x = 0.5 * ||x||^2
    // Solution: x = 0
    nlc::QPOasesEigen qp(2, 0);
    qp.data_.H = Eigen::Matrix2d::Identity();
    qp.data_.g = Eigen::Vector2d::Zero();
    qp.data_.lb = -100.0 * Eigen::Vector2d::Ones();
    qp.data_.ub = 100.0 * Eigen::Vector2d::Ones();

    qp.setup();
    auto status = qp.solve();

    EXPECT_EQ(status, qpOASES::SUCCESSFUL_RETURN);

    auto x_opt = qp.getOptimizer();
    EXPECT_NEAR(x_opt(0), 0.0, 1e-6);
    EXPECT_NEAR(x_opt(1), 0.0, 1e-6);
}

TEST_F(QPOasesEigenTest, SimpleLinearCostQP) {
    // min 0.5 * x^T * I * x + [-1, -1]^T * x
    // s.t. -100 <= x <= 100
    // Solution: x = [1, 1]
    nlc::QPOasesEigen qp(2, 0);
    qp.data_.H = Eigen::Matrix2d::Identity();
    qp.data_.g << -1.0, -1.0;
    qp.data_.lb = -100.0 * Eigen::Vector2d::Ones();
    qp.data_.ub = 100.0 * Eigen::Vector2d::Ones();

    qp.setup();
    qp.solve();

    auto x_opt = qp.getOptimizer();
    EXPECT_NEAR(x_opt(0), 1.0, 1e-6);
    EXPECT_NEAR(x_opt(1), 1.0, 1e-6);
}

TEST_F(QPOasesEigenTest, QPWithLinearConstraints) {
    // min 0.5 * x^T * H * x + g^T * x
    // H = [1 0; 0 0.5], g = [1.5; 1.0]
    // s.t. 0.5 <= x1 <= 5.0, -2.0 <= x2 <= 2.0
    //      -1.0 <= x1 + x2 <= 2.0
    // Expected solution from qpOASES example
    nlc::QPOasesEigen qp(2, 1);
    qp.data_.H << 1.0, 0.0, 0.0, 0.5;
    qp.data_.g << 1.5, 1.0;
    qp.data_.A << 1.0, 1.0;
    qp.data_.lbA << -1.0;
    qp.data_.ubA << 2.0;
    qp.data_.lb << 0.5, -2.0;
    qp.data_.ub << 5.0, 2.0;

    qp.options.setToFast();
    qp.options.printLevel = qpOASES::PL_NONE;
    qp.setup();
    auto status = qp.solve();

    EXPECT_EQ(status, qpOASES::SUCCESSFUL_RETURN);

    auto x_opt = qp.getOptimizer();
    // Verify solution satisfies constraints
    EXPECT_GE(x_opt(0), 0.5 - 1e-6);
    EXPECT_LE(x_opt(0), 5.0 + 1e-6);
    EXPECT_GE(x_opt(1), -2.0 - 1e-6);
    EXPECT_LE(x_opt(1), 2.0 + 1e-6);

    double Ax = x_opt(0) + x_opt(1);
    EXPECT_GE(Ax, -1.0 - 1e-6);
    EXPECT_LE(Ax, 2.0 + 1e-6);
}

TEST_F(QPOasesEigenTest, ThreeVariableQP) {
    // min 0.5 * x^T * 2I * x
    // s.t. x1 + 3*x2 + 4*x3 <= -3
    //      -100 <= x <= 100
    nlc::QPOasesEigen qp(3, 1);
    qp.data_.H = 2.0 * Eigen::Matrix3d::Identity();
    qp.data_.g = Eigen::Vector3d::Zero();
    qp.data_.A << 1.0, 3.0, 4.0;
    qp.data_.lbA << -100.0;
    qp.data_.ubA << -3.0;
    qp.data_.lb = -100.0 * Eigen::Vector3d::Ones();
    qp.data_.ub = 100.0 * Eigen::Vector3d::Ones();

    qp.setup();
    auto status = qp.solve();

    EXPECT_EQ(status, qpOASES::SUCCESSFUL_RETURN);

    auto x_opt = qp.getOptimizer();
    // Verify constraint is satisfied
    double constraint_val = x_opt(0) + 3.0 * x_opt(1) + 4.0 * x_opt(2);
    EXPECT_LE(constraint_val, -3.0 + 1e-6);
}

// =============================================================================
// QPSwiftEigen Tests
// Note: QPSwift solver is not fully implemented, so these tests are disabled
// =============================================================================

class QPSwiftEigenTest : public ::testing::Test {
protected:
    void SetUp() override {
        // QPSwift uses inequality constraints: Ax <= b
    }
};

TEST_F(QPSwiftEigenTest, DISABLED_SimpleQPWithBoxConstraints) {
    // Note: QPSwift solver implementation is incomplete
    // This test is disabled until the solver is properly implemented
    nlc::QPSwiftEigen qp(2, 6, 0);
    qp.H << 1.0, 0.0, 0.0, 0.5;
    qp.f << 1.5, 1.0;
    qp.c = 0;

    qp.A << 1.0, 0.0,
           -1.0, 0.0,
            0.0, 1.0,
            0.0, -1.0,
            1.0, 1.0,
           -1.0, -1.0;

    qp.b << 5.0, -0.5, 2.0, 2.0, 2.0, 1.0;

    int status = qp.solve();
    EXPECT_EQ(status, 0) << "QPSwift should return 0 for successful solve";

    // Result is stored in xOpt member
    EXPECT_GE(qp.xOpt(0), 0.5 - 1e-4);
    EXPECT_LE(qp.xOpt(0), 5.0 + 1e-4);
}

// =============================================================================
// Main
// =============================================================================

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
