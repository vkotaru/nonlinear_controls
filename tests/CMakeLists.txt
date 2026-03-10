# ------------------------------------------------------------------------------
# Unit Tests for nonlinear_controls
# ------------------------------------------------------------------------------

include(GoogleTest)

# Helper function to create test executables
function(nlc_add_test TEST_NAME)
    add_executable(${TEST_NAME} ${ARGN})
    target_link_libraries(${TEST_NAME}
        PRIVATE
            nonlinear_controls
            GTest::gtest
            GTest::gtest_main
    )
    target_include_directories(${TEST_NAME}
        PRIVATE
            ${CMAKE_SOURCE_DIR}/include
    )
    gtest_discover_tests(${TEST_NAME})
endfunction()

# ------------------------------------------------------------------------------
# Test Executables
# ------------------------------------------------------------------------------

# Data types and manifolds tests
nlc_add_test(test_data_types test_data_types.cpp)

# MPC solver tests
nlc_add_test(test_linear_mpc test_linear_mpc.cpp)

# QP solver tests
nlc_add_test(test_qp_solvers test_qp_solvers.cpp)

# Note: Legacy test files (test_nlc.cpp, test_qp.cpp, test_mpc.cpp, test_osqp_mpc.cpp)
# are kept for reference but not built. They use manual cout-based assertions.
# See the new test files for Google Test equivalents.
