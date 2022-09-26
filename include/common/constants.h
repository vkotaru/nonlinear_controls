#ifndef NONLINEAR_CONTROLS_COMMON_CONSTANTS_H
#define NONLINEAR_CONTROLS_COMMON_CONSTANTS_H

#include <eigen3/Eigen/Dense>

namespace nonlinear_controls {
#define G_SI 9.81
#define GRAVITY_VECTOR Eigen::Vector3d(0, 0., 9.81)
#define ERROR_TOL 1.e-3
#define SECOND_TO_MILLISECOND 1000

#define E1 Eigen::Vector3d(1., 0., 0)
#define E2 Eigen::Vector3d(0., 1., 0)
#define E3 Eigen::Vector3d(0., 0., 1.)
}

#endif //NONLINEAR_CONTROLS_COMMON_CONSTANTS_H
