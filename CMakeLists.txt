cmake_minimum_required(VERSION 2.8)
project(nonlinear_control)

include_directories(include)

set(CONTROL_SRCS
        src/base_controller.cpp)

add_library(control src/main.cpp ${CONTROL_SRCS})

