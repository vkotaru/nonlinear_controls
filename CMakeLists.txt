cmake_minimum_required(VERSION 3.20)

project(nonlinear_controls
    VERSION 1.0.0
    LANGUAGES CXX
    DESCRIPTION "C++ library for nonlinear control systems"
)

# ------------------------------------------------------------------------------
# Build Configuration
# ------------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set default build type if not specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "RelWithDebInfo" "MinSizeRel")
endif()

# Export compile commands for IDE support
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ------------------------------------------------------------------------------
# Options
# ------------------------------------------------------------------------------
option(NLC_BUILD_TESTS "Build unit tests" ON)
option(NLC_BUILD_EXAMPLES "Build example applications" ON)
option(NLC_ENABLE_MATPLOTLIB "Enable matplotlib-cpp for plotting (requires Python)" OFF)

# ------------------------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------------------------
find_package(Eigen3 REQUIRED)
find_package(OsqpEigen QUIET)

if(NOT OsqpEigen_FOUND)
    message(STATUS "OsqpEigen not found. OSQP solver features unavailable.")
endif()

# Optional matplotlib support
if(NLC_ENABLE_MATPLOTLIB)
    find_package(Python3 COMPONENTS Development REQUIRED)
endif()

# ------------------------------------------------------------------------------
# Third-party libraries
# ------------------------------------------------------------------------------
add_subdirectory(third_party)

# ------------------------------------------------------------------------------
# Library Sources
# ------------------------------------------------------------------------------
set(NLC_SOURCES
    src/SO3_control.cpp
    src/SE3_control.cpp
    src/qp_interface.cpp
    src/SO3_clf.cpp
    src/tracking_mpc.cpp
    src/quadprog/quadprog.cpp
)

# ------------------------------------------------------------------------------
# Main Library Target
# ------------------------------------------------------------------------------
add_library(nonlinear_controls ${NLC_SOURCES})
add_library(nonlinear_controls::nonlinear_controls ALIAS nonlinear_controls)

target_include_directories(nonlinear_controls
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/cvxgen
)

target_link_libraries(nonlinear_controls
    PUBLIC
        Eigen3::Eigen
    PRIVATE
        qpOASES
        $<$<BOOL:${OsqpEigen_FOUND}>:OsqpEigen::OsqpEigen>
)

# Need qpOASES headers for public API (mpc_qpoases.hpp)
target_include_directories(nonlinear_controls
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third_party/qpOASES/include>
)

# Matplotlib support (optional)
if(NLC_ENABLE_MATPLOTLIB)
    target_include_directories(nonlinear_controls PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/matplotlib-cpp)
    target_link_libraries(nonlinear_controls PRIVATE Python3::Python)
    target_compile_definitions(nonlinear_controls PRIVATE
        MATPLOTLIBCPP_PYTHON_HEADER=Python.h
        NLC_HAS_MATPLOTLIB)
endif()

# Compiler warnings for our code
target_compile_options(nonlinear_controls PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>
    $<$<CXX_COMPILER_ID:MSVC>:/W4>
)

# ------------------------------------------------------------------------------
# Tests
# ------------------------------------------------------------------------------
if(NLC_BUILD_TESTS)
    enable_testing()
    include(FetchContent)

    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.14.0
    )
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)

    add_subdirectory(tests)
endif()

# ------------------------------------------------------------------------------
# Examples
# ------------------------------------------------------------------------------
if(NLC_BUILD_EXAMPLES)
    if(NLC_ENABLE_MATPLOTLIB)
        add_executable(example_point_mass examples/mpc_point_mass.cpp)
        target_link_libraries(example_point_mass PRIVATE nonlinear_controls)
        target_include_directories(example_point_mass PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/third_party/matplotlib-cpp)
    else()
        message(STATUS "Skipping example_point_mass (requires NLC_ENABLE_MATPLOTLIB)")
    endif()
endif()

# ------------------------------------------------------------------------------
# Installation
# ------------------------------------------------------------------------------
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Install library
install(TARGETS nonlinear_controls qpOASES
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Install headers
install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING
        PATTERN "*.h"
        PATTERN "*.hpp"
)

# Install third-party headers needed by public API
install(DIRECTORY third_party/qpOASES/include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Generate version file
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/nonlinear_controls-config-version.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

# Generate config file
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/nonlinear_controls-config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/nonlinear_controls-config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/nonlinear_controls
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/nonlinear_controls-config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/nonlinear_controls-config-version.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/nonlinear_controls
)

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------
message(STATUS "")
message(STATUS "nonlinear_controls ${PROJECT_VERSION}")
message(STATUS "  Build type:     ${CMAKE_BUILD_TYPE}")
message(STATUS "  Tests:          ${NLC_BUILD_TESTS}")
message(STATUS "  Examples:       ${NLC_BUILD_EXAMPLES}")
message(STATUS "  Matplotlib:     ${NLC_ENABLE_MATPLOTLIB}")
message(STATUS "  OsqpEigen:      ${OsqpEigen_FOUND}")
message(STATUS "")
