#include "eigen3/unsupported/Eigen/MatrixFunctions"
#include <eigen3/Eigen/Dense>

#include <cmath>
#include <fstream>
#include <iostream>
#include <qpOASES.hpp>

#include "common/utils.hpp"
#include "controls/SE3_vblmpc.h"
#include "controls/SO3_clf.h"
#include "controls/double_int_mpc.hpp"
#include "controls/linear_mpc.h"
#include "data_types/data_types.hpp"

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
namespace nlc = nonlinear_controls;

#define nlc_real double

/////////////////////////////////////////////
/// global variables
////////////////////////////////////////////
int N = 5;
double T = 20;
double dt = double(1 / 200.0); // 200 Hz
int nx = 6, nu = 3;
unsigned long startime, previoustime, currenttime;

// storage variable for the plots
std::vector<double> t, x, y, z, ux, uy, uz, eRx, eRy, eRz, Mx, My, Mz;

void clear_plot_vars() {
  t.clear();
  x.clear();
  y.clear();
  z.clear();
  ux.clear();
  uy.clear();
  uz.clear();
  eRx.clear();
  eRy.clear();
  eRz.clear();
  Mx.clear();
  My.clear();
  Mz.clear();
}

template<typename T> struct FlatOutputs {
  Eigen::Matrix<T, 3, 1> x, dx, d2x, u;
};

/////////////////////////////////////////////
/// trajectories
////////////////////////////////////////////
void sinusoidalTraj(double &yd, double &dyd, double &d2yd, const float t,
                    const double freq, const double amp,
                    const double ampPercent, double ampOffset,
                    double freqOffset) {
  yd =
      amp * (1 - ampPercent + ampPercent * pow(cos(freq * t + freqOffset), 2)) +
          ampOffset;
  //    dyd = -amp * ampPercent * freq * sin(2 * freq * t + freqOffset);
  dyd = -2 * amp * ampPercent * freq * cos(freqOffset + freq * t) *
      sin(freqOffset + freq * t);
  //    d2yd =  -2 * amp * ampPercent * pow(freq, 2) * cos(2 * freq * t +
  //    freqOffset);
  d2yd =
      2 * amp * ampPercent * pow(freq, 2) * pow(sin(freqOffset + freq * t), 2) -
          2 * amp * ampPercent * pow(freq, 2) * pow(cos(freqOffset + freq * t), 2);
}

template<typename T> void circularTraj(T t, FlatOutputs<T> &flats) {
  T x0 = 0., y0 = 0., z0 = 0, r = 1, f = 0.5;
  T w = 2 * f * M_PI;

  flats.x << x0 + r * sin(t * w), y0 + r * cos(t * w), z0;
  flats.dx << r * w * cos(t * w), -r * w * sin(t * w), 0;
  flats.d2x << -r * w * w * sin(t * w), -r * w * w * cos(t * w), 0;
  flats.u = flats.d2x;
}

/////////////////////////////////////////////
/// individual tests
////////////////////////////////////////////
void test_position_mpc() {
  /////////////////////////////////////////////
  /// position mpc in 3D with LTI dynamics
  ////////////////////////////////////////////
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Testing 3D position MPC with LTI Dynamics... " << std::endl;
  nx = 6;
  nu = 3;
  N = 1;
  nlc::LinearMPC pos_mpc_(N, nx, nu);

  // dynamics
  nlc::MatrixXd A, B;
  A.resize(nx, nx);
  B.resize(nx, nu);
  A.setIdentity();
  A.topRightCorner(3, 3) += dt * Eigen::Matrix<double, 3, 3>::Identity();
  B << 0.5 * dt * dt * Eigen::Matrix<double, 3, 3>::Identity(),
      dt * Eigen::Matrix<double, 3, 3>::Identity();
  std::cout << "Initializing dynamics ... " << std::endl;
  std::cout << "A: \n" << A << std::endl;
  std::cout << "B: \n" << B << std::endl;
  pos_mpc_.init_dynamics(A, B);

  // gains
  nlc::MatrixXd Q, P, R;
  Q.resize(nx, nx);
  P.resize(nx, nx);
  R.resize(nu, nu);
  Q << 1000 * Eigen::Matrix3d::Identity(), Eigen::Matrix3d::Zero(),
      Eigen::Matrix3d::Zero(), 100 * Eigen::Matrix3d::Identity();
  P << 8.1314, 0.0000, -0.0000, 0.6327, -0.0000, -0.0000, 0.0000, 8.1314,
      0.0000, 0.0000, 0.6327, 0.0000, -0.0000, 0.0000, 8.1314, -0.0000, 0.0000,
      0.6327, 0.6327, 0.0000, -0.0000, 0.2606, -0.0000, -0.0000, -0.0000,
      0.6327, 0.0000, -0.0000, 0.2606, 0.0000, -0.0000, 0.0000, 0.6327, -0.0000,
      0.0000, 0.2606;
  P = 1e4 * P;
  R = 1 * Eigen::Matrix<double, 3, 3>::Identity();
  pos_mpc_.set_mpc_gains(Q, P, R);

  // state & input bounds
  Eigen::Matrix<double, 6, 1> state_lb, state_ub;
  Eigen::Matrix<double, 3, 1> input_lb, input_ub;
  state_lb = -100 * Eigen::Matrix<double, 6, 1>::Ones();
  state_ub = -state_lb;
  input_lb = -100 * Eigen::Matrix<double, 3, 1>::Ones();
  input_ub = -input_lb;
  pos_mpc_.set_input_bounds(input_lb, input_ub);
  pos_mpc_.set_state_bounds(state_lb, state_ub);
  pos_mpc_.set_cpu_time_limit(dt);

  // construct the mpc
  pos_mpc_.construct();

  // initialize simulation
  Eigen::Matrix<double, 6, 1> goal_state, state;
  goal_state << 1.0, 2.0, 3.0, 0.0, 0.0, 0.0;
  state.setZero();
  std::cout << "\n X0: " << state.transpose() << std::endl;
  std::cout << "\n Xgoal: " << goal_state.transpose() << std::endl;
  nlc::MatrixXd K, zOpt, uOpt;
  uOpt.resize(nu, 1);
  zOpt.resize(N * nu, 1);
  K.resize(nu, nx);
  K << 30.6287 * Eigen::Matrix3d::Identity(),
      12.4527 * Eigen::Matrix3d::Identity();

  bool mpc_initialized{false};

  clear_plot_vars();
  for (int j = 0; j < int(2 / dt); ++j) { // only for 2 seconds

    /// LQR (for sanity check)
    auto uPD = -K * (state - goal_state);

    /// running mpc controller
    zOpt = pos_mpc_.run((state - goal_state));
    uOpt = zOpt.block(0, 0, nu, 1);

    // integrating dynamics using the input at N = 0
    std::cout << "PD Input:>>> " << uPD.transpose() << std::endl;
    std::cout << "MPC Input:>>> " << uOpt.transpose() << std::endl;
    state = A * state + B * uOpt;

    t.push_back(dt * j);
    x.push_back(state(0));
    y.push_back(state(1));
    z.push_back(state(2));
    ux.push_back(uOpt(0));
    uy.push_back(uOpt(1));
    uz.push_back(uOpt(2));
  }

  /// plots
  plt::suptitle("Position MPC LTI");
  plt::subplot(1, 2, 1);
  plt::plot(t, x, "r-");
  plt::plot(t, y, "g--");
  plt::plot(t, z, "b-.");
  plt::grid(true);

  plt::subplot(1, 2, 2);
  plt::plot(t, ux, "r-");
  plt::plot(t, uy, "g--");
  plt::plot(t, uz, "b-.");
  plt::grid(true);
  plt::show();
  std::cout << "------------------------------------------" << std::endl;
}

void test_pos_traj_mpc() {
  /////////////////////////////////////////////
  /// position mpc in 3D with LTV dynamics
  ////////////////////////////////////////////
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Testing 3D position MPC with LTV Dynamics... " << std::endl;
  nx = 6;
  nu = 3;
  N = 5;
  T = 5;
  nlc::LinearMPCt pos_mpc_(N, nx, nu);

  // dynamics
  nlc::MatrixXd A, B;
  A.resize(nx, nx);
  B.resize(nx, nu);
  A.setIdentity();
  A.topRightCorner(3, 3) += dt * Eigen::Matrix<double, 3, 3>::Identity();
  B << 0.5 * dt * dt * Eigen::Matrix<double, 3, 3>::Identity(),
      dt * Eigen::Matrix<double, 3, 3>::Identity();
  std::cout << "Initializing dynamics ... " << std::endl;
  std::cout << "A: \n" << A << std::endl;
  std::cout << "B: \n" << B << std::endl;
  for (int i = 0; i < N; ++i) {
    pos_mpc_.init_dynamics(A, B);
  }

  // gains
  nlc::MatrixXd Q, P, R;
  Q.resize(nx, nx);
  P.resize(nx, nx);
  R.resize(nu, nu);
  Q << 1000 * Eigen::Matrix3d::Identity(), Eigen::Matrix3d::Zero(),
      Eigen::Matrix3d::Zero(), 100 * Eigen::Matrix3d::Identity();
  P << 8.1314, 0.0000, -0.0000, 0.6327, -0.0000, -0.0000, 0.0000, 8.1314,
      0.0000, 0.0000, 0.6327, 0.0000, -0.0000, 0.0000, 8.1314, -0.0000, 0.0000,
      0.6327, 0.6327, 0.0000, -0.0000, 0.2606, -0.0000, -0.0000, -0.0000,
      0.6327, 0.0000, -0.0000, 0.2606, 0.0000, -0.0000, 0.0000, 0.6327, -0.0000,
      0.0000, 0.2606;
  P = 1e4 * P;
  R = 1 * Eigen::Matrix<double, 3, 3>::Identity();
  pos_mpc_.set_mpc_gains(Q, P, R);

  // state & input bounds
  Eigen::Matrix<double, 6, 1> state_lb, state_ub;
  Eigen::Matrix<double, 3, 1> input_lb, input_ub;
  state_lb = -100 * Eigen::Matrix<double, 6, 1>::Ones();
  state_ub = -state_lb;
  input_lb = -100 * Eigen::Matrix<double, 3, 1>::Ones();
  input_ub = -input_lb;
  pos_mpc_.set_input_bounds(input_lb, input_ub);
  pos_mpc_.set_state_bounds(state_lb, state_ub);

  // construct the mpc
  pos_mpc_.construct();

  // initialize simulation
  Eigen::Matrix<double, 6, 1> goal_state, state, state_ref;
  goal_state << 1.0, 2.0, 3.0, 0.0, 0.0, 0.0;
  state.setZero();
  state_ref.setZero();
  std::cout << "\n X0: " << state.transpose() << std::endl;
  std::cout << "\n Xgoal: " << goal_state.transpose() << std::endl;
  nlc::MatrixXd K, zOpt, uOpt;
  uOpt.resize(nu, 1);
  zOpt.resize(N * nu, 1);
  K.resize(nu, nx);
  K << 30.6287 * Eigen::Matrix3d::Identity(),
      12.4527 * Eigen::Matrix3d::Identity();
  FlatOutputs<double> flats;

  clear_plot_vars();
  for (int j = 0; j < int(T / dt); ++j) { // only for 2 seconds
    circularTraj<double>(j * dt, flats);
    state_ref << flats.x, flats.dx;

    /// LQR (for sanity check)
    //        uOpt = -K * (state - goal_state);
    /// feed-forward
    //        uOpt = flats.u;

    /// running mpc controller
    zOpt = pos_mpc_.run((state - state_ref), A, B);
    uOpt = zOpt.block(0, 0, nu, 1);

    // integrating dynamics using the input at N = 0
    //        std::cout << "Input:>>> " << uOpt.transpose() << std::endl;
    state = A * state + B * (uOpt + flats.u);

    t.push_back(dt * j);
    x.push_back(state(0));
    y.push_back(state(1));
    z.push_back(state(2));
    ux.push_back(uOpt(0));
    uy.push_back(uOpt(1));
    uz.push_back(uOpt(2));
  }

  /// plots
  plt::suptitle("Position MPC LTV");
  plt::subplot(1, 2, 1);
  plt::plot(t, x, "r-");
  plt::plot(t, y, "g--");
  plt::plot(t, z, "b-.");
  plt::grid(true);

  plt::subplot(1, 2, 2);
  plt::plot(t, ux, "r-");
  plt::plot(t, uy, "g--");
  plt::plot(t, uz, "b-.");
  plt::grid(true);

  std::map<std::string, std::string> keywords;
  keywords.insert(
      std::pair<std::string, std::string>("label", "parametric curve"));
  plt::plot3(x, y, z, keywords);
  plt::grid(true);
  plt::show();
  std::cout << "------------------------------------------" << std::endl;
}

void test_attitude_mpc() {
  /////////////////////////////////////////////
  /// VBL attitude mpc with LTI dynamics
  ////////////////////////////////////////////
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Testing VBL attitude mpc with LTI dynamics... " << std::endl;
  nx = 6;
  nu = 3;
  N = 5;
  nlc::LinearMPC att_mpc_(N, nx, nu);

  // dynamics
  Eigen::Matrix3d inertia, inertia_inv_;
  inertia << 112533, 0, 0, 0, 36203, 0, 0, 0, 42673;
  inertia = inertia * 1e-6;
  inertia_inv_ = inertia.inverse();
  auto generate_dynamics = [&](Eigen::Vector3d Omd) {
    Eigen::Matrix<double, 6, 6> A;
    A << -nlc::utils::hat(Omd), Eigen::Matrix<double, 3, 3>::Identity(),
        Eigen::Matrix<double, 3, 3>::Zero(),
        inertia_inv_ *
            (nlc::utils::hat(inertia * Omd) - nlc::utils::hat(Omd) * inertia);
    A = Eigen::Matrix<double, 6, 6>::Identity() + A * dt;
    return A;
  };
  Eigen::Vector3d Omd = Eigen::Vector3d::Zero();
  nlc::MatrixXd A, B;
  A.resize(nx, nx);
  B.resize(nx, nu);
  A = generate_dynamics(Omd);
  B << Eigen::Matrix<double, 3, 3>::Zero(), inertia_inv_;
  B = B * dt;
  std::cout << "Initializing dynamics ... " << std::endl;
  std::cout << "A: \n" << A << std::endl;
  std::cout << "B: \n" << B << std::endl;
  att_mpc_.init_dynamics(A, B);

  // gains
  nlc::MatrixXd Q, P, R;
  Q.resize(nx, nx);
  P.resize(nx, nx);
  R.resize(nu, nu);
  Q << 1000 * Eigen::Matrix3d::Identity(), Eigen::Matrix3d::Zero(),
      Eigen::Matrix3d::Zero(), 100 * Eigen::Matrix3d::Identity();
  P << 6.6514, -0.0000, -0.0000, 0.0894, -0.0000, -0.0000, -0.0000, 6.5123,
      -0.0000, -0.0000, 0.0440, -0.0000, -0.0000, -0.0000, 6.5231, -0.0000,
      -0.0000, 0.0475, 0.0894, -0.0000, -0.0000, 0.0293, -0.0000, -0.0000,
      -0.0000, 0.0440, -0.0000, -0.0000, 0.0141, -0.0000, -0.0000, -0.0000,
      0.0475, -0.0000, -0.0000, 0.0152;
  P = P * 1e4;
  R = 1 * Eigen::Matrix<double, 3, 3>::Identity();
  att_mpc_.set_mpc_gains(Q, P, R);

  // state & input bounds
  Eigen::Matrix<double, 6, 1> state_lb, state_ub;
  Eigen::Matrix<double, 3, 1> input_lb, input_ub;
  state_lb = -100 * Eigen::Matrix<double, 6, 1>::Ones();
  state_ub = -state_lb;
  input_lb = -100 * Eigen::Matrix<double, 3, 1>::Ones();
  input_ub = -input_lb;
  att_mpc_.set_input_bounds(input_lb, input_ub);
  att_mpc_.set_state_bounds(state_lb, state_ub);

  // construct the mpc
  att_mpc_.construct();

  // initialize simulation
  nlc::TSE3 att, attd;
  att.R.setZero();
  att.R(0, 2) = 1;
  att.R(1, 1) = 1;
  att.R(2, 0) = -1;
  att.print();
  attd.print();
  std::cout << att - attd << std::endl;
  nlc::MatrixXd K, zOpt, uOpt, err;
  uOpt.resize(nu, 1);
  err.resize(nx, 1);
  zOpt.resize(N * nu, 1);
  K.resize(nu, nx);
  K.setZero();
  K << 25.1721, 0., 0., 8.3715, 0., 0., 0., 16.4628, 0., 0., 5.3605, 0., 0., 0.,
      17.9819, 0., 0., 5.8649;

  clear_plot_vars();
  T = 2;
  for (int j = 0; j < int(T / dt); ++j) {
    err = att.extractTSO3() - attd.extractTSO3();
    /// LQR (for sanity check)
    //        uOpt = -K * (err);

    /// running mpc controller
    zOpt = att_mpc_.run(err);
    uOpt = zOpt.block(0, 0, nu, 1);

    /// computing the geometric-controller
    uOpt += att.Omega.cross(inertia * att.Omega);
    uOpt +=
        -inertia * (att.Omega.cross(att.R.transpose() * attd.R * attd.Omega) -
            att.R.transpose() * attd.R * attd.dOmega);

    /// integrating the full nonlinear-dynamics using the input at N = 0
    Eigen::Matrix<double, 3, 3> hat_Om = nlc::utils::hat(att.Omega);
    att.R = att.R * (hat_Om * dt).exp();
    Eigen::Vector3d dOm;
    dOm = inertia_inv_ * (uOpt - att.Omega.cross(inertia * att.Omega));
    att.Omega += dOm * dt;

    t.push_back(dt * j);
    x.push_back(err(0));
    y.push_back(err(1));
    z.push_back(err(2));
    ux.push_back(uOpt(0));
    uy.push_back(uOpt(1));
    uz.push_back(uOpt(2));
  }
  /// plots
  plt::suptitle("Attitude MPC LTI");
  plt::subplot(1, 2, 1);
  plt::plot(t, x, "r-");
  plt::plot(t, y, "g--");
  plt::plot(t, z, "b-.");
  plt::grid(true);

  plt::subplot(1, 2, 2);
  plt::plot(t, ux, "r-");
  plt::plot(t, uy, "g--");
  plt::plot(t, uz, "b-.");
  plt::grid(true);
  plt::show();
  std::cout << "------------------------------------------" << std::endl;
}

void test_att_traj_mpc() {
  /////////////////////////////////////////////
  /// VBL attitude mpc with LTV dynamics
  ////////////////////////////////////////////
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Testing VBL attitude mpc with LTV dynamics... " << std::endl;
  nx = 6;
  nu = 3;
  N = 10;
  nlc::LinearMPCt att_mpc_(N, nx, nu);

  // dynamics
  Eigen::Matrix3d inertia, inertia_inv_;
  inertia << 112533, 0, 0, 0, 36203, 0, 0, 0, 42673;
  inertia = inertia * 1e-6;
  inertia_inv_ = inertia.inverse();
  auto generate_dynamics = [&](Eigen::Vector3d Omd) {
    Eigen::Matrix<double, 6, 6> A;
    A << -nlc::utils::hat(Omd), Eigen::Matrix<double, 3, 3>::Identity(),
        Eigen::Matrix<double, 3, 3>::Zero(),
        inertia_inv_ *
            (nlc::utils::hat(inertia * Omd) - nlc::utils::hat(Omd) * inertia);
    A = Eigen::Matrix<double, 6, 6>::Identity() + A * dt;
    return A;
  };
  Eigen::Vector3d Omd = Eigen::Vector3d::Zero();
  nlc::MatrixXd A, B;
  A.resize(nx, nx);
  B.resize(nx, nu);
  A = generate_dynamics(Omd);
  B << Eigen::Matrix<double, 3, 3>::Zero(), inertia_inv_;
  B = B * dt;

  // gains
  nlc::MatrixXd Q, P, R;
  Q.resize(nx, nx);
  P.resize(nx, nx);
  R.resize(nu, nu);
  Q << 1000 * Eigen::Matrix3d::Identity(), Eigen::Matrix3d::Zero(),
      Eigen::Matrix3d::Zero(), 100 * Eigen::Matrix3d::Identity();
  P << 6.6514, -0.0000, -0.0000, 0.0894, -0.0000, -0.0000, -0.0000, 6.5123,
      -0.0000, -0.0000, 0.0440, -0.0000, -0.0000, -0.0000, 6.5231, -0.0000,
      -0.0000, 0.0475, 0.0894, -0.0000, -0.0000, 0.0293, -0.0000, -0.0000,
      -0.0000, 0.0440, -0.0000, -0.0000, 0.0141, -0.0000, -0.0000, -0.0000,
      0.0475, -0.0000, -0.0000, 0.0152;
  P = P * 1e4;
  R = 1 * Eigen::Matrix<double, 3, 3>::Identity();
  att_mpc_.set_mpc_gains(Q, P, R);

  // state & input bounds
  Eigen::Matrix<double, 6, 1> state_lb, state_ub;
  Eigen::Matrix<double, 3, 1> input_lb, input_ub;
  state_lb = -500 * Eigen::Matrix<double, 6, 1>::Ones();
  state_ub = -state_lb;
  input_lb = -5 * Eigen::Matrix<double, 3, 1>::Ones();
  input_ub = -input_lb;
  att_mpc_.set_input_bounds(input_lb, input_ub);
  att_mpc_.set_state_bounds(state_lb, state_ub);

  // initialize simulation
  nlc::TSE3 att, attd;
  att.R.setZero();
  att.R(0, 2) = 1;
  att.R(1, 1) = 1;
  att.R(2, 0) = -1;
  att.print();
  attd.print();
  std::cout << att - attd << std::endl;
  nlc::MatrixXd K, zOpt, uOpt, uRef, err;
  uOpt.resize(nu, 1);
  uRef.resize(nu, 1);
  err.resize(nx, 1);
  zOpt.resize(N * nu, 1);
  K.resize(nu, nx);
  K.setZero();
  K << 25.1721, 0., 0., 8.3715, 0., 0., 0., 16.4628, 0., 0., 5.3605, 0., 0., 0.,
      17.9819, 0., 0., 5.8649;

  double ydp;
  double dydp;
  double d2ydp;
  double freqp = 0.5;
  for (int j = 0; j < N; ++j) {
    sinusoidalTraj(ydp, dydp, d2ydp, dt * j, freqp, 0.25, 1, -0.125, 3.141 / 4);
    attd.R = nlc::utils::rotmZ(0) * nlc::utils::rotmY(ydp);
    attd.Omega(1) = dydp;
    attd.dOmega(1) = d2ydp;

    if (j == 0) {
      att = attd;
    }

    A = generate_dynamics(attd.Omega);
    att_mpc_.init_dynamics(A, B);
  }

  // construct the mpc
  att_mpc_.construct();

  clear_plot_vars();
  T = 20;
  for (int j = 0; j < int(T / dt); ++j) {
    sinusoidalTraj(ydp, dydp, d2ydp, dt * j, freqp, 0.25, 1, -0.125, 3.141 / 4);
    attd.R = nlc::utils::rotmZ(0) * nlc::utils::rotmY(ydp);
    attd.Omega(1) = dydp;
    attd.dOmega(1) = d2ydp;
    uRef = inertia * attd.dOmega + attd.Omega.cross(inertia * attd.Omega);
    err = att.extractTSO3() - attd.extractTSO3();
    uOpt = uRef;

    /// LQR (for sanity check)
    //        uOpt = -K * (err);

    /// running mpc controller
    //        sinusoidalTraj(ydp, dydp, d2ydp, dt * (j + N), freqp, 0.25, 1,
    //        -0.125, 3.141 / 4); attd.R = nlc::utils::rotmZ(0) *
    //        nlc::utils::rotmY(ydp); attd.Omega(1) = dydp; attd.dOmega(1) =
    //        d2ydp; zOpt = att_mpc_.run(err, generate_dynamics(attd.Omega), B);
    //        uOpt = zOpt.block(0, 0, nu, 1);

    /// computing the geometric-controller
    //        uOpt += att.Omega.cross(inertia * att.Omega);
    //        uOpt += -inertia * (att.Omega.cross(att.R.transpose() * attd.R *
    //        attd.Omega) - att.R.transpose() * attd.R * attd.dOmega); uOpt +=
    //        uRef;

    /// integrating the full nonlinear-dynamics using the input at N = 0
    Eigen::Matrix<double, 3, 3> hat_Om = nlc::utils::hat(att.Omega);
    att.R = att.R * (hat_Om * dt).exp();
    Eigen::Vector3d dOm;
    dOm = inertia_inv_ * (uOpt - att.Omega.cross(inertia * att.Omega));
    att.Omega += dOm * dt;

    t.push_back(dt * j);
    x.push_back(err(0));
    y.push_back(err(1));
    z.push_back(err(2));
    ux.push_back(uOpt(0));
    uy.push_back(uOpt(1));
    uz.push_back(uOpt(2));
  }
  /// plots
  plt::suptitle("Attitude MPC LTV");
  plt::subplot(1, 2, 1);
  plt::plot(t, x, "r-");
  plt::plot(t, y, "g--");
  plt::plot(t, z, "b-.");
  plt::grid(true);

  plt::subplot(1, 2, 2);
  plt::plot(t, ux, "r-");
  plt::plot(t, uy, "g--");
  plt::plot(t, uz, "b-.");
  plt::grid(true);
  plt::show();
  std::cout << "------------------------------------------" << std::endl;
}

void test_double_int_mpc() {
  /////////////////////////////////////////////
  /// Position MPC
  ////////////////////////////////////////////
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Testing Position MPC class... " << std::endl;
  N = 5;
  nlc::DoubleIntMPC<3> pos_mpc_(N, dt);
  // initialize simulation
  Eigen::Matrix<double, 6, 1> goal_state, state;
  goal_state << 1.0, 2.0, 3.0, 0.0, 0.0, 0.0;
  state.setZero();
  std::cout << "\n X0: " << state.transpose() << std::endl;
  std::cout << "\n Xgoal: " << goal_state.transpose() << std::endl;
  nlc::MatrixXd K, zOpt, uOpt;
  uOpt.resize(nu, 1);
  zOpt.resize(N * nu, 1);
  K.resize(nu, nx);
  K << 30.6287 * Eigen::Matrix3d::Identity(),
      12.4527 * Eigen::Matrix3d::Identity();

  clear_plot_vars();
  for (int j = 0; j < int(2 / dt); ++j) { // only for 2 seconds

    /// LQR (for sanity check)
    //        uOpt = -K * (state - goal_state);

    /// running mpc controller
    uOpt = pos_mpc_.run((state - goal_state));

    // integrating dynamics
    state = pos_mpc_.updateState(state, uOpt);

    t.push_back(dt * j);
    x.push_back(state(0));
    y.push_back(state(1));
    z.push_back(state(2));
    ux.push_back(uOpt(0));
    uy.push_back(uOpt(1));
    uz.push_back(uOpt(2));
  }

  /// plots
  plt::suptitle("Position MPC LTI");
  plt::subplot(1, 2, 1);
  plt::plot(t, x, "r-");
  plt::plot(t, y, "g--");
  plt::plot(t, z, "b-.");
  plt::grid(true);

  plt::subplot(1, 2, 2);
  plt::plot(t, ux, "r-");
  plt::plot(t, uy, "g--");
  plt::plot(t, uz, "b-.");
  plt::grid(true);
  plt::show();
  std::cout << "------------------------------------------" << std::endl;
}

void test_se3_vbl_mpc() {
  /////////////////////////////////////////////
  /// SE3 VBL MPC
  ////////////////////////////////////////////
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Testing SE3 MPC class... " << std::endl;
  N = 5;
  Eigen::Matrix3d inertia, inertia_inv_;
  inertia << 112533, 0, 0, 0, 36203, 0, 0, 0, 42673;
  inertia = inertia * 1e-6;
  inertia_inv_ = inertia.inverse();
  double mass = 3;
  nlc::SE3VblMPC mpc_(false, N, dt, mass, inertia);

  // initialize simulation
  nlc::VectorXd err;
  err.resize(12, 1);
  nlc::TSE3 xd, state;
  nlc::Wrench input;
  state.position << 1.0, 2.0, 3.0;
  state.R.setZero();
  state.R(0, 2) = 1;
  state.R(1, 1) = 1;
  state.R(2, 0) = -1;

  std::cout << "X0: " << std::endl;
  state.print();
  std::cout << "Xgoal: " << std::endl;
  xd.print();

  double ydp;
  double dydp;
  double d2ydp;
  double freqp = 0.5;
  for (int j = 0; j < N; ++j) {
    sinusoidalTraj(ydp, dydp, d2ydp, dt * j, freqp, 0.25, 1, -0.125, 3.141 / 4);
    xd.R = nlc::utils::rotmZ(0) * nlc::utils::rotmY(ydp);
    xd.Omega(1) = dydp;
    xd.dOmega(1) = d2ydp;
    mpc_.att_mpc_.init_dynamics(xd.Omega);
  }
  mpc_.att_mpc_.construct();

  clear_plot_vars();
  T = 10;
  startime = nlc::utils::get_current_time();
  for (int j = 0; j < int(T / dt); ++j) { // only for 2 seconds
    sinusoidalTraj(ydp, dydp, d2ydp, dt * j, freqp, 0.25, 1, -0.125, 3.141 / 4);
    xd.R = nlc::utils::rotmZ(0) * nlc::utils::rotmY(ydp);
    xd.Omega(1) = dydp;
    xd.dOmega(1) = d2ydp;
    err = state - xd;

    /// running mpc controller
    mpc_.run(dt, state, xd, input);

    // integrating dynamics (TODO: add this function to SE3VblMPC)
    //
    state.position +=
        dt * state.velocity +
            0.5 * dt * dt *
                (input.force -
                    mass * (Eigen::Vector3d() << 0.0, 0.0, 9.81).finished());
    state.velocity +=
        dt *
            (input.force - mass * (Eigen::Vector3d() << 0.0, 0.0, 9.81).finished());
    state.R = state.R * (nlc::utils::hat(state.Omega * dt)).exp();
    state.Omega += (inertia_inv_ *
        (input.torque - state.Omega.cross(inertia * state.Omega))) *
        dt;

    t.push_back(dt * j);
    x.push_back(state.position(0));
    y.push_back(state.position(1));
    z.push_back(state.position(2));
    ux.push_back(input.force(0));
    uy.push_back(input.force(1));
    uz.push_back(input.force(2));
    eRx.push_back(err(6));
    eRy.push_back(err(7));
    eRz.push_back(err(8));
    Mx.push_back(input.torque(0));
    My.push_back(input.torque(1));
    Mz.push_back(input.torque(2));
  }
  currenttime = nlc::utils::get_current_time();
  printf("Time taken to simulate %.2fs (i.e., %d iter) is: %f\n", T,
         int(T / dt), double(currenttime - startime) / 1000000.0);

  /// plots
  plt::suptitle("SE3 MPC LTV");
  plt::subplot(2, 2, 1);
  plt::plot(t, x, "r-");
  plt::plot(t, y, "g--");
  plt::plot(t, z, "b-.");
  plt::grid(true);

  plt::subplot(2, 2, 2);
  plt::plot(t, ux, "r-");
  plt::plot(t, uy, "g--");
  plt::plot(t, uz, "b-.");
  plt::grid(true);

  plt::subplot(2, 2, 3);
  plt::plot(t, eRx, "r-");
  plt::plot(t, eRy, "g--");
  plt::plot(t, eRz, "b-.");
  plt::grid(true);

  plt::subplot(2, 2, 4);
  plt::plot(t, Mx, "r-");
  plt::plot(t, My, "g--");
  plt::plot(t, Mz, "b-.");
  plt::grid(true);
  plt::show();
  std::cout << "------------------------------------------" << std::endl;
}

void test_att_traj_clf() {
  /////////////////////////////////////////////
  /// VBL attitude mpc with LTV dynamics
  ////////////////////////////////////////////
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Testing VBL attitude CLF!!!" << std::endl;
  nx = 6;
  nu = 3;
  N = 10;

  // dynamics
  Eigen::Matrix3d inertia, inertia_inv_;
  inertia << 0.0049, 0.0000055, 0.0000054, 0.0000055, 0.0053, 0.000021,
      0.0000054, 0.000021, 0.0098;
  inertia_inv_ = inertia.inverse();
  nlc::SO3Clf att_clf_(inertia);
  att_clf_.init();

  // initialize simulation
  nlc::TSO3 att, attd;
  att.R.setZero();
  att.R(0, 2) = 1;
  att.R(1, 1) = 1;
  att.R(2, 0) = -1;
  att.print();
  attd.print();
  std::cout << att - attd << std::endl;

  double ydp;
  double dydp;
  double d2ydp;
  double freqp = 0.5;

  Eigen::Vector3d uRef, uOpt;
  Eigen::Matrix<double, 6, 1> err;
  clear_plot_vars();
  T = 20;
  for (int j = 0; j < int(T / dt); ++j) {
    sinusoidalTraj(ydp, dydp, d2ydp, dt * j, freqp, 0.25, 1, -0.125, 3.141 / 4);
    attd.R = nlc::utils::rotmZ(0) * nlc::utils::rotmY(ydp);
    attd.Omega(1) = dydp;
    attd.dOmega(1) = d2ydp;
    uRef = inertia * attd.dOmega + attd.Omega.cross(inertia * attd.Omega);
    err = att - attd;
    uOpt = uRef;

    /// CLF controller
    att_clf_.run(dt, att, attd, uOpt);

    /// integrating the full nonlinear-dynamics using the input at N = 0
    Eigen::Matrix<double, 3, 3> hat_Om = nlc::utils::hat(att.Omega);
    att.R = att.R * (hat_Om * dt).exp();
    Eigen::Vector3d dOm;
    dOm = inertia_inv_ * (uOpt - att.Omega.cross(inertia * att.Omega));
    att.Omega += dOm * dt;

    t.push_back(dt * j);
    x.push_back(err(0));
    y.push_back(err(1));
    z.push_back(err(2));
    ux.push_back(uOpt(0));
    uy.push_back(uOpt(1));
    uz.push_back(uOpt(2));
  }
  /// plots
  plt::suptitle("Attitude CLF");
  plt::subplot(1, 2, 1);
  plt::plot(t, x, "r-");
  plt::plot(t, y, "g--");
  plt::plot(t, z, "b-.");
  plt::grid(true);

  plt::subplot(1, 2, 2);
  plt::plot(t, ux, "r-");
  plt::plot(t, uy, "g--");
  plt::plot(t, uz, "b-.");
  plt::grid(true);
  plt::show();
  std::cout << "------------------------------------------" << std::endl;
}

/////////////////////////////////////////////
/// main
////////////////////////////////////////////
int main() {

  /// test position (time-invariant) mpc
  test_position_mpc();

  /// test position (trajectory) mpc
  //    test_pos_traj_mpc();

  /// test attitude (time-invariant) mpc
  //  test_attitude_mpc();

  /// test attitude (trajectory) mpc
  //        test_att_traj_mpc();

  /// test position mpc class
  //    test_double_int_mpc();

  /// test SE3 vbl mpc class
  //    test_se3_vbl_mpc();

  /// test attitude (trajectory) clf
  //    test_att_traj_clf();

  return 0;
}
