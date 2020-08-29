#ifndef NLC_UTILS_H
#define NLC_UTILS_H
#include "eigen3/Eigen/Dense"
#include <sys/time.h>
namespace nonlinear_control {
namespace utils {

static struct timeval tv;
inline unsigned long get_current_time() {
    gettimeofday(&tv, NULL);
    return 1000000 * tv.tv_sec + tv.tv_usec;
}

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

//template <typename T>
//inline Eigen::Matrix<T, 3, 1> vee(const Eigen::Matrix<T, 3, 3> mat_) {
//    return (Eigen::Matrix<T, 3, 1>() << mat_(2, 1), mat_(0, 2), mat_(1, 0)).finished();
//}


///////////////////////////////////////
// TODO Clean and verify the following
///////////////////////////////////////
inline Eigen::Matrix<double, 3, 3> rotmX(const double roll) {
    Eigen::Matrix<double, 3, 3> R;
    R.row(0) << 1, 0, 0;
    R.row(1) << 0, cos(roll), -sin(roll);
    R.row(2) << 0, sin(roll), cos(roll);
    return R;
}
inline Eigen::Matrix<double, 3, 3> rotmY(const double pitch) {
    Eigen::Matrix<double, 3, 3> R;
    R.row(0) << cos(pitch), 0, sin(pitch);
    R.row(1) << 0, 1, 0;
    R.row(2) << -sin(pitch), 0, cos(pitch);
    return R;
}
inline Eigen::Matrix<double, 3, 3> rotmZ(const double yaw) {
    Eigen::Matrix<double, 3, 3> R;
    R.row(0) << cos(yaw), -sin(yaw), 0;
    R.row(1) << sin(yaw), cos(yaw), 0;
    R.row(2) << 0, 0, 1;
    return R;
}
inline Eigen::Matrix<double, 3, 3> rotmZYX(const Eigen::Vector3d eulerZYX) {
    Eigen::Matrix<double, 3, 3> R;
    R = rotmZ(eulerZYX[0]) * rotmY(eulerZYX[1]) * rotmX(eulerZYX[2]);
    return R;
}
inline Eigen::Matrix<double, 3, 3> rotmXYZ(const Eigen::Vector3d eulerXYZ) {
    Eigen::Matrix<double, 3, 3> R;
    R = rotmX(eulerXYZ[0]) * rotmY(eulerXYZ[1]) * rotmZ(eulerXYZ[2]);
    return R;
}
inline Eigen::Vector3d rotmToZYXEulerAngles(const Eigen::Matrix<double, 3, 3>& R) {
    double thetaY = asin(-R(2, 0));
    double thetaZ = atan2(R(1, 0), R(0, 0));
    double thetaX = atan2(R(2, 1), R(2, 2));
    Eigen::Vector3d v;
    v << thetaZ, thetaY, thetaX;
    return v;
}
inline Eigen::Matrix<double, 4, 4> getTransformMat(const Eigen::Matrix<double, 3, 3>& R, const Eigen::Vector3d& pos) {
    Eigen::Matrix<double, 4, 4> T = Eigen::MatrixXd::Zero(4, 4);
    T.block<3, 3>(0, 0) << R;
    T.col(3).head(3) << pos;
    T(4, 4) = 1;
    return T;
}
inline Eigen::Vector3d eulerAngleXYZToZYX(const Eigen::Vector3d& rpy) {
    Eigen::Vector3d ypr;
    //compute rotation matrx
    Eigen::Matrix<double, 3, 3> R = rotmXYZ(rpy);
    ypr = rotmToZYXEulerAngles(R);
}
///////////////////////////////////////////

//////////////////////////////////////////
/// templates
//////////////////////////////////////////
#define hatf hat<float>
#define hatd hat<double>
#define veef vee<float>
#define veed vee<double>
//////////////////////////////////////////
}  // namespace utils
}  // namespace nonlinear_control

#endif  // NLC_UTILS_H
