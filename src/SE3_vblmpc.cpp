#include "controls/SE3_vblmpc.h"
namespace nonlinear_controls {

SE3VblMPC::SE3VblMPC(bool islti, const int N, const double dt, const double m,
                     const Eigen::Matrix3d &J)
    : GeometricController(), IS_LTI(islti), N(N), mass_(m), inertia_(J),
      pos_mpc_(N, dt), att_mpc_(islti, N, dt, J) {

}

SE3VblMPC::~SE3VblMPC() = default;

void SE3VblMPC::init() {
  att_mpc_.init_dynamics(Eigen::Vector3d::Zero());
}

void SE3VblMPC::run(double dt, TSE3 x, TSE3 xd, Wrench &u) {
  // translational force
  Eigen::Matrix<double, 12, 1> errors = x - xd;
  u.force = pos_mpc_.run(errors.head(6));
  u.force += mass_ * (Eigen::Vector3d() << 0.0, 0.0, 9.81).finished();

  // rotational torque
  Eigen::Vector3d torque;
  att_mpc_.run(dt, x.extractTSO3(), xd.extractTSO3(), torque);
  u.torque = torque;
}

} // namespace nonlinear_controls
