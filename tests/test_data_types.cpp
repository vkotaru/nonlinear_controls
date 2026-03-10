/**
 * @file test_data_types.cpp
 * @brief Unit tests for data types and manifolds (SO3, SE3, TSE3)
 */

#include <gtest/gtest.h>
#include <Eigen/Dense>

#include "data_types/manifolds.hpp"
#include "data_types/SO3.hpp"
#include "data_types/SE3.hpp"

using namespace nonlinear_controls;
using namespace manifolds;

// =============================================================================
// SpecialOrthogonal Tests
// =============================================================================

TEST(SpecialOrthogonalTest, IdentityInverse) {
    SpecialOrthogonal<float, 4> R4 = SpecialOrthogonal<float, 4>::Identity();
    auto R4_inv = R4.inverse();

    double diff_norm = (R4 - R4_inv).norm();
    EXPECT_LT(diff_norm, 1e-6) << "Identity matrix should equal its inverse";
}

TEST(SpecialOrthogonalTest, OrthogonalityProperty) {
    SpecialOrthogonal<double, 3> R = SpecialOrthogonal<double, 3>::Identity();
    Eigen::Matrix3d product = R * R.transpose();
    Eigen::Matrix3d identity = Eigen::Matrix3d::Identity();

    EXPECT_LT((product - identity).norm(), 1e-10)
        << "R * R^T should equal identity for orthogonal matrices";
}

// =============================================================================
// SO3 Tests
// =============================================================================

class SO3Test : public ::testing::Test {
protected:
    void SetUp() override {
        // Sample rotation matrices from original test
        R << 0.8073, -0.5176, 0.2835,
             0.5229,  0.8501, 0.0630,
            -0.2736,  0.0974, 0.9569;

        Rd << 0.4569, -0.0903, 0.8849,
              0.7532,  0.5684, -0.3309,
             -0.4731,  0.8178, 0.3278;
    }

    SO3 R, Rd;
};

TEST_F(SO3Test, IdentityMatrix) {
    SO3 R_identity = SO3::Identity();
    Eigen::Matrix3d expected = Eigen::Matrix3d::Identity();

    EXPECT_LT((R_identity - expected).norm(), 1e-10);
}

TEST_F(SO3Test, InverseEqualsTranspose) {
    SO3 R_identity = SO3::Identity();
    auto R_inv = R_identity.inverse();
    auto R_trans = R_identity.transpose();

    EXPECT_LT((R_inv - R_trans).norm(), 1e-10)
        << "For SO3, inverse should equal transpose";
}

TEST_F(SO3Test, ErrorComputationConsistency) {
    // Compute error using member function
    auto error_member = R.error(Rd);

    // Compute error using static function
    auto error_static = SO3::error(R, Rd);

    EXPECT_LT((error_member - error_static).norm(), 1e-10)
        << "Member and static error functions should give same result";
}

TEST_F(SO3Test, ErrorComputationManual) {
    // SO3::error uses: vee(0.5 * (Rd^T * R - R^T * Rd))
    // Note the order: Rd^T * R, not R^T * Rd
    auto eR_hat = 0.5 * (Rd.transpose() * R - R.transpose() * Rd);
    Eigen::Vector3d eR_manual;
    eR_manual << eR_hat(2, 1), eR_hat(0, 2), eR_hat(1, 0);

    auto eR_computed = R.error(Rd);

    EXPECT_LT((eR_manual - eR_computed).norm(), 1e-10)
        << "Computed error should match manual calculation";
}

// =============================================================================
// TSE3 Tests
// =============================================================================

class TSE3Test : public ::testing::Test {
protected:
    void SetUp() override {
        // Set up test state
        x.position = Eigen::Vector3d::Random();
        x.velocity = Eigen::Vector3d::Random();
        x.Omega = Eigen::Vector3d::Random();
        x.R << 0.8073, -0.5176, 0.2835,
               0.5229,  0.8501, 0.0630,
              -0.2736,  0.0974, 0.9569;

        // Set up desired state
        xd.position = Eigen::Vector3d::Ones();
        xd.velocity = Eigen::Vector3d::Ones();
        xd.Omega = Eigen::Vector3d::Zero();
        xd.R = Eigen::Matrix3d::Identity();
    }

    TSE3 x, xd;
};

TEST_F(TSE3Test, ErrorOperatorConsistency) {
    // Error using operator-
    auto error_operator = x - xd;

    // Error using member function
    auto error_member = x.error(xd);

    EXPECT_LT((error_operator - error_member).norm(), 1e-10)
        << "operator- and error() should give same result";
}

TEST_F(TSE3Test, ErrorComputationManual) {
    // Manual computation of error
    auto eR_hat = 0.5 * (xd.R.transpose() * x.R - x.R.transpose() * xd.R);
    Eigen::Vector3d eR;
    eR << eR_hat(2, 1), eR_hat(0, 2), eR_hat(1, 0);

    Eigen::Matrix<double, 12, 1> manual_error;
    manual_error << (x.position - xd.position),
                    (x.velocity - xd.velocity),
                    eR,
                    x.Omega - x.R.transpose() * xd.R * xd.Omega;

    auto computed_error = x - xd;

    EXPECT_LT((computed_error - manual_error).norm(), 1e-10)
        << "Computed error should match manual calculation";
}

TEST_F(TSE3Test, ZeroErrorForIdenticalStates) {
    TSE3 x1, x2;
    x1.position = Eigen::Vector3d::Ones();
    x1.velocity = Eigen::Vector3d::Zero();
    x1.R = Eigen::Matrix3d::Identity();
    x1.Omega = Eigen::Vector3d::Zero();

    x2 = x1;  // Copy

    auto error = x1 - x2;
    EXPECT_LT(error.norm(), 1e-10) << "Error between identical states should be zero";
}

// =============================================================================
// Main
// =============================================================================

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
