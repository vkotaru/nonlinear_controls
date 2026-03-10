# Contributing to nonlinear_controls

Thank you for your interest in contributing! This document provides guidelines for contributing to the project.

## Development Setup

1. Clone the repository with submodules:
   ```bash
   git clone --recursive https://github.com/vkotaru/nonlinear_controls.git
   cd nonlinear_controls
   ```

2. Install dependencies:
   ```bash
   # Ubuntu/Debian
   sudo apt-get install libeigen3-dev cmake

   # macOS
   brew install eigen cmake
   ```

3. Build with tests enabled:
   ```bash
   cmake -B build -DCMAKE_BUILD_TYPE=Debug -DNLC_BUILD_TESTS=ON
   cmake --build build -j$(nproc)
   ```

4. Run tests:
   ```bash
   ctest --test-dir build --output-on-failure
   ```

## Code Style

This project uses clang-format for code formatting. The configuration is in `.clang-format`.

### Before committing:

```bash
# Format all source files
find src include -name '*.cpp' -o -name '*.hpp' -o -name '*.h' | xargs clang-format -i

# Check formatting (CI will fail if this produces changes)
find src include -name '*.cpp' -o -name '*.hpp' -o -name '*.h' | xargs clang-format --dry-run --Werror
```

### Style Guidelines

- Follow the Google C++ Style Guide with modifications in `.clang-format`
- Use `snake_case` for functions and variables
- Use `CamelCase` for class names
- Use `UPPER_CASE` for constants
- Add Doxygen documentation for public APIs

## Testing

- Write tests for all new features
- Run the full test suite before submitting a PR
- Maintain test coverage for critical functionality

### Adding Tests

Tests use Google Test framework. Add new test files to `tests/` and update `tests/CMakeLists.txt`:

```cpp
#include <gtest/gtest.h>
#include "your_header.hpp"

TEST(YourTestSuite, TestName) {
    // Test implementation
    EXPECT_EQ(expected, actual);
}
```

## Pull Request Process

1. Create a feature branch from `master`:
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. Make your changes and ensure:
   - All tests pass
   - Code is formatted with clang-format
   - New features have tests
   - Documentation is updated if needed

3. Commit with clear messages:
   ```bash
   git commit -m "Add feature X that does Y"
   ```

4. Push and create a Pull Request:
   ```bash
   git push origin feature/your-feature-name
   ```

5. Wait for CI to pass and request review

## Reporting Issues

When reporting issues, please include:

- A clear description of the problem
- Steps to reproduce
- Expected vs actual behavior
- Your environment (OS, compiler, CMake version)
- Relevant code snippets or error messages

## Questions?

Feel free to open an issue for questions or discussions about the project.
