# nonlinear_controls

[![CI](https://github.com/vkotaru/nonlinear_controls/actions/workflows/ci.yml/badge.svg)](https://github.com/vkotaru/nonlinear_controls/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](LICENSE)

A modern C++ library for nonlinear control systems, including Model Predictive Control (MPC) and geometric controllers for robotics applications.

## Features

- **Linear MPC**: Time-invariant linear MPC with qpOASES solver (bundled)
- **Geometric Controllers**: SO(3) and SE(3) pose controllers
- **Manifold Data Types**: SO3, SE3, S2 with proper error computations

## Requirements

- C++17 compiler (GCC 9+, Clang 10+)
- CMake 3.20+
- Eigen3

### Optional Dependencies
- OsqpEigen (for OSQP solver support)
- Python3 (for matplotlib plotting in examples)

## Installation

### From Source

```bash
git clone --recursive https://github.com/vkotaru/nonlinear_controls.git
cd nonlinear_controls
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
sudo cmake --install build
```

### CMake Options

| Option | Default | Description |
|--------|---------|-------------|
| `NLC_BUILD_TESTS` | ON | Build unit tests |
| `NLC_BUILD_EXAMPLES` | ON | Build example applications |
| `NLC_ENABLE_MATPLOTLIB` | OFF | Enable matplotlib plotting (requires Python) |

### Using in Your Project

After installation, use `find_package`:

```cmake
find_package(nonlinear_controls REQUIRED)
target_link_libraries(your_target PRIVATE nonlinear_controls::nonlinear_controls)
```

Or add as a subdirectory:

```cmake
add_subdirectory(third_party/nonlinear_controls)
target_link_libraries(your_target PRIVATE nonlinear_controls)
```

## Quick Start

### Linear MPC Example

```cpp
#include <nonlinear_controls.h>

namespace nlc = nonlinear_controls;

int main() {
    const int N = 10, nx = 6, nu = 3;
    nlc::MPCQPOases mpc(N, nx, nu);

    // Setup dynamics: x_{k+1} = A*x_k + B*u_k
    Eigen::MatrixXd A = Eigen::MatrixXd::Identity(nx, nx);
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(nx, nu);
    B.bottomRows(nu).setIdentity();

    // Configure MPC
    mpc.init_dynamics(A, B);
    mpc.set_gains(Q, P, R);  // State cost, terminal cost, input cost
    mpc.set_state_bounds(state_lb, state_ub);
    mpc.set_input_bounds(input_lb, input_ub);
    mpc.construct();

    // Run controller
    Eigen::VectorXd current_state = /* ... */;
    Eigen::VectorXd goal_state = /* ... */;

    auto result = mpc.run(current_state, goal_state);
    if (result.has_value()) {
        auto u_optimal = result.value().head(nu);
        // Apply control input
    }
}
```

### Geometric Controller Example

```cpp
#include <nonlinear_controls.h>

using namespace nonlinear_controls;

// Create SO3 rotation matrices
SO3 R_current, R_desired;
R_current = SO3::Identity();
// ... set rotation values

// Compute rotation error
Eigen::Vector3d rotation_error = R_current.error(R_desired);
```

## Running Tests

```bash
cd build
ctest --output-on-failure
```

## Documentation

Generate documentation with Doxygen:

```bash
doxygen Doxyfile
open docs/html/index.html
```

## Project Structure

```
nonlinear_controls/
├── include/
│   ├── controls/       # MPC and controller implementations
│   ├── common/         # Utilities, QP interfaces
│   ├── data_types/     # SO3, SE3, manifolds
│   └── dynamics/       # Dynamics models
├── src/                # Implementation files
├── tests/              # Unit tests (Google Test)
├── examples/           # Example applications
└── third_party/        # Bundled dependencies (qpOASES)
```

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

This project is licensed under the BSD 3-Clause License - see [LICENSE](LICENSE) for details.

## Acknowledgments

- [qpOASES](https://github.com/coin-or/qpOASES) - Quadratic Programming solver
- [OsqpEigen](https://github.com/robotology/osqp-eigen) - OSQP C++ interface
