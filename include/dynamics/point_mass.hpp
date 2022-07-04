#include <eigen3/Eigen/Dense>
#include <iostream>
#include <vector>

#include "nonlinear_controls.h"
#include "matplotlibcpp.h"

namespace nonlinear_controls {
namespace plt = matplotlibcpp;

class PointMass3D {
public:
  struct Vars {
    std::vector<double> t;
    std::vector<double> x, y, z;
    std::vector<double> ux, uy, uz;
    std::vector<double> vx, vy, vz;
    std::vector<double> time_elapsed;

    void clear() {
      t.clear();
      x.clear();
      y.clear();
      z.clear();
      ux.clear();
      uy.clear();
      uz.clear();
      vx.clear();
      vy.clear();
      vz.clear();
      time_elapsed.clear();
    }
  };

protected:
  // log variables for plots
  Vars log_vars{};

  // variables
  Eigen::Matrix<double, 6, 1> state_{};
  double t{0.};

  // dynamics
  Eigen::Matrix<double, 6, 6> A_;
  Eigen::Matrix<double, 6, 3> B_;

  // parameters
  double h{1. / 200.};
  double mass{1.};

public:
  PointMass3D() : PointMass3D(0.005) {}
  explicit PointMass3D(const double dt) : h(dt) {
    A_.setIdentity();
    A_.topRightCorner(3, 3) += dt * Eigen::Matrix3d::Identity();
    B_ << 0.5 * dt * dt * Eigen::Matrix3d::Identity(),
        dt * Eigen::Matrix3d::Identity();
  }

  ~PointMass3D() = default;

  void init(const Eigen::Matrix<double, 6, 1> &_state, bool do_log = false, const double time_elapsed = 0) {
    state_ = _state;
    t = 0;
    if (do_log)
      log(0, state_, Eigen::Vector3d::Zero(), time_elapsed); // Note, dummy input for convenience
  }

  void step(const Eigen::Vector3d &u, bool do_log = false, const double time_elapsed = 0) {
    Eigen::Vector3d net_accel = (u - GRAVITY_VECTOR) / mass;
    state_ = A_ * state_ + B_ * net_accel;
    t += h;
    if (do_log)
      log(t, state_, u, time_elapsed);
  }

  void log(const double _t,
           const Eigen::Matrix<double, 6, 1> &st,
           const Eigen::Vector3d &u,
           const double time_elapsed = 0) {
    log_vars.t.push_back(_t);
    log_vars.x.push_back(st(0));
    log_vars.y.push_back(st(1));
    log_vars.z.push_back(st(2));

    log_vars.vx.push_back(st(3));
    log_vars.vy.push_back(st(4));
    log_vars.vz.push_back(st(5));

    log_vars.ux.push_back(u(0));
    log_vars.uy.push_back(u(1));
    log_vars.uz.push_back(u(2));
    log_vars.time_elapsed.push_back(time_elapsed);
  }

  void plot() const {
    /// plots
    plt::suptitle("Position MPC LTI");
    plt::subplot(2, 2, 1);
    plt::plot(log_vars.t, log_vars.x, "r-");
    plt::plot(log_vars.t, log_vars.y, "g--");
    plt::plot(log_vars.t, log_vars.z, "b-.");
    plt::grid(true);

    plt::subplot(2, 2, 2);
    plt::plot(log_vars.t, log_vars.ux, "r-");
    plt::plot(log_vars.t, log_vars.uy, "g--");
    plt::plot(log_vars.t, log_vars.uz, "b-.");
    plt::grid(true);

    plt::subplot(2, 2, 3);
    plt::plot(log_vars.t, log_vars.vx, "r-");
    plt::plot(log_vars.t, log_vars.vy, "g--");
    plt::plot(log_vars.t, log_vars.vz, "b-.");
    plt::grid(true);

    plt::subplot(2, 2, 4);
    plt::plot(log_vars.t, log_vars.time_elapsed, "r-");
    plt::grid(true);
    plt::show();
  }
  void clear_log_vars() {
    log_vars.clear();
  }

  inline const Eigen::Matrix<double, 6, 6> A() const { return A_; }
  inline const Eigen::Matrix<double, 6, 3> B() const { return B_; }
  inline const Eigen::Matrix<double, 6, 1> state() const { return state_; }

};

}
