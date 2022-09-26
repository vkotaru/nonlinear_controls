cmake_minimum_required(VERSION 3.3.2)
project(nonlinear_controls)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-O3")
#set(CMAKE_CXX_FLAGS "-Wall -Os")
SET(CMAKE_BUILD_TYPE Release)

find_package(OsqpEigen)
find_package(Eigen3)

include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})
include_directories(
        include
        third_party/cvxgen
        third_party/qpOASES
        third_party/qpOASES/include
)

set(INCLUDE_DIRS
        include
        third_party/cvxgen
        third_party/qpOASES
        third_party/qpOASES/include
        )

set(QUADPROG_SRCS
        src/quadprog/quadprog.cpp
        src/quadprog/qpswift_eigen.cpp
        )

set(CONTROL_SRCS
        src/SO3_control.cpp
        src/SE3_control.cpp
        src/qp_interface.cpp
        src/SO3_clf.cpp
        src/tracking_mpc.cpp
        include/common/constants.h
        include/common/log.hpp
        examples/mpc_point_mass.cpp
        include/controls/mpc_epigraph.hpp
        include/controls/linear_mpc.hpp)

add_subdirectory(third_party)

set(CMAKE_BUILD_TYPE Debug)

add_library(nonlinear_controls ${CONTROL_SRCS} ${QUADPROG_SRCS})
target_link_libraries(nonlinear_controls qpOASES_cpp epigraph  OsqpEigen::OsqpEigen)
target_include_directories(nonlinear_controls PUBLIC . ${INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/include)

# tests
add_executable(test_main tests/test_nlc.cpp)
target_link_libraries(test_main nonlinear_controls)
add_executable(test_qp tests/test_qp.cpp)
target_link_libraries(test_qp nonlinear_controls)
#add_executable(test_mpc tests/test_mpc.cpp)
#target_link_libraries(test_mpc nonlinear_controls)
add_executable(osqpMPCExample tests/test_osqp_mpc.cpp)
target_link_libraries(osqpMPCExample nonlinear_controls)

# examples
add_executable(example_point_mass examples/mpc_point_mass.cpp)
target_link_libraries(example_point_mass nonlinear_controls)


