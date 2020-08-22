#ifndef NLC_UTILS_H
#define NLC_UTILS_H
#include "eigen3/Eigen/Dense"
namespace nonlinear_control {
namespace utils {

template <typename T>
inline Eigen::Matrix<T, 3, 3> hat(const Eigen::Matrix<T, 3, 1> vector_) {
    return (Eigen::Matrix<T, 3, 3>() << 0.0, -vector_(2), vector_(1),
            vector_(2), 0.0, -vector_(0),
            -vector_(1), vector_(0), 0.0)
        .finished();
}

template <typename T>
inline Eigen::Matrix<T, 3, 1> vee(const Eigen::Matrix<T, 3, 3> mat_) {
    return (Eigen::Matrix<T, 3, 1>() << mat_(2, 1), mat_(0, 2), mat_(1, 0)).finished();
}

#define hatf hat<float>
#define hatd hat<double>

#define veef vee<float>
#define veed vee<double>

}  // namespace utils
}  // namespace nonlinear_control

#endif  // NLC_UTILS_H