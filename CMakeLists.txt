cmake_minimum_required(VERSION 2.8)
project(nonlinear_control)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-Wall -Os")
SET( CMAKE_BUILD_TYPE Debug)

find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
add_definitions(-DMATPLOTLIBCPP_PYTHON_HEADER=Python.h)

include_directories(
include
third_party/cvxgen
third_party/qpOASES
third_party/qpOASES/include 
) # doesn't seem to be including qpOASES.hpp

set(INCLUDE_DIRS 
include
third_party/cvxgen
third_party/qpOASES
third_party/qpOASES/include
third_party/matplotlib-cpp # TODO learn why this is working and not include_directories
)

set(CONTROL_HEADERS
include/manifolds.hpp
include/data_types.hpp 
include/base_control.h 
include/geometric_control.h
include/SO3_control.h
include/SE3_control.h
include/qp_interface.h
include/clf_qp.h
include/mpc_base.h
include/mpc_base2.h
include/SO3_vblmpc.h
)


set(CONTROL_SRCS
src/base_control.cpp
src/geometric_control.cpp
src/SO3_control.cpp
src/SE3_control.cpp
src/qp_interface.cpp
src/clf_qp.cpp
src/mpc_base.cpp
src/mpc_base2.cpp
src/SO3_vblmpc.cpp
)

add_subdirectory(third_party)

add_library(nonlinear_controls ${CONTROL_SRCS})
target_link_libraries(nonlinear_controls qpOASES ${PYTHON_LIBRARIES})
target_include_directories(nonlinear_controls PUBLIC ${INCLUDE_DIRS})

add_executable(test_main tests/test_nlc.cpp)
target_link_libraries(test_main nonlinear_controls)
add_executable(test_qp tests/test_qp.cpp)
target_link_libraries(test_qp nonlinear_controls)
add_executable(test_mpc tests/test_mpc.cpp)
target_link_libraries(test_mpc nonlinear_controls)
add_executable(test_mpc2 tests/test_mpc2.cpp)
target_link_libraries(test_mpc2 nonlinear_controls)
