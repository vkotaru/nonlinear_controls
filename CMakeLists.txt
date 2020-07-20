cmake_minimum_required(VERSION 2.8)
project(nonlinear_control)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(include)

set(CONTROL_SRCS
        src/base_control.cpp
        src/geometric_control.cpp
        src/SO3_control.cpp)

add_library(libcontrol ${CONTROL_SRCS})
target_include_directories(libcontrol PUBLIC .)

add_executable(test_main src/main.cpp ${CONTROL_SRCS})
target_link_libraries(test_main libcontrol)