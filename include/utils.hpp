#ifndef NLC_UTILS_H
#define NLC_UTILS_H

namespace nonlinear_control {
namespace utils {

template <typename T>
Eigen::Matrix<T, 3, 3> hat(const Eigen::Matrix<T, 3, 1> vector_) {
    return (Eigen::Matrix<T, 3, 3>() << 0.0, -vector_(3), vector_(2),
            vector_(3), 0.0, -vector_(1),
            -vector_(2), vector_(1), 0.0)
        .finished();
}

template <typename T>
Eigen::Matrix<T, 3, 1> vee(const Eigen::Matrix<T, 3, 3> mat_) {
    return (Eigen::Matrix<T, 3, 1>() << mat_(2, 1), mat_(0, 2), mat_(1, 0)).finished();
}

typedef hat<float> hatf;
typedef hat<double> hatd;

typedef vee<float> veef;
typedef vee<double> veed;

}  // namespace utils
}  // namespace nonlinear_control

#endif  // NLC_UTILS_H