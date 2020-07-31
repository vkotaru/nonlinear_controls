cmake_minimum_required(VERSION 2.8)
project(nonlinear_control)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-Wall -Os")


include_directories(
include
third_party/cvxgen
)

set(CONTROL_HEADERS
        include/manifolds.hpp
        include/data_types.hpp 
        include/base_control.h 
        include/geometric_control.h
	include/SO3_control.h
	include/SE3_control.h)


set(CONTROL_SRCS
        src/base_control.cpp
        src/geometric_control.cpp
	src/SO3_control.cpp
	src/SE3_control.cpp)

add_library(nonlinear_controls ${CONTROL_HEADERS} ${CONTROL_SRCS})
target_link_libraries(nonlinear_controls m)
target_include_directories(nonlinear_controls PUBLIC "${CMAKE_CURRENT_SOURCE_DIR/include}")

add_executable(test_main src/test_nlc.cpp)
target_link_libraries(test_main nonlinear_controls)
