#include "controls/SE3_vblmpc.h"
namespace nonlinear_controls {
template <typename T>
SE3VblMPC<T>::SE3VblMPC(bool islti, const int &N, const T &dt, const T &m,
                        const Eigen::Matrix<T, 3, 3> &J)
    : GeometricController<T>(), IS_LTI(islti), N(N), mass_(m), inertia_(J),
      pos_mpc_(N, dt), att_mpc_(islti, N, dt, J) {



}

template <typename T> SE3VblMPC<T>::~SE3VblMPC() = default;

template <typename T> void SE3VblMPC<T>::init() {
  att_mpc_.init_dynamics(Eigen::Matrix<T, 3, 1>::Zero());
}

template <typename T>
void SE3VblMPC<T>::run(T dt, TSE3<T> x, TSE3<T> xd, Wrench<T> &u) {
  // translational force
  Eigen::Matrix<T, 12, 1> errors = x - xd;
  u.force = pos_mpc_.run(errors.head(6));
  u.force += mass_ * (Eigen::Matrix<T, 3, 1>() << 0.0, 0.0, 9.81).finished();

  // rotational torque
  Eigen::Matrix<T, 3, 1> torque;
  att_mpc_.run(dt, x.extractTSO3(), xd.extractTSO3(), torque);
  u.torque = torque;
}

template class SE3VblMPC<float>;
template class SE3VblMPC<double>;

} // namespace nonlinear_controls
