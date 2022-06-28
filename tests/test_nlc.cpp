#include <eigen3/Eigen/Dense>
#include <iostream>

#include "controls/SO3_control.h"
#include "common/cvxgen_interface.h"
using namespace nonlinear_controls;
using namespace manifolds;

#define GREEN "\033[01;32m"
#define RED "\033[01;31m"
#define NC "\033[0m"

int main() {
  srand((unsigned int) time(0));
  SpecialOrthogonal<float, 4> R4 = SpecialOrthogonal<float, 4>::Identity();
  std::cout << R4 << "\n"
            << R4.inverse() << std::endl;
  if (abs((R4 - R4.inverse()).norm()) < 1e-10) {
    std::cout << GREEN << "SpecialOrthogonal base case tested" << NC << std::endl;
  }

  Eigen::Matrix3d R0d;
  R0d << 3, 5, 7, 2, 5, 7, 6, 5, 4;
  Eigen::Matrix3d b = R0d.inverse();
  std::cout << R0d << "\n"
            << R0d.transpose() << "\n"
            << b << std::endl;

  Eigen::VectorXf v2;

  SO3 Rtest = SO3::Identity();
  std::cout << Rtest << "\n"
            << Rtest.inverse() << std::endl;

  std::cout << "***************************************" << std::endl;
  SO3 R, Rd;
  R << 0.8073, -0.5176, 0.2835, 0.5229, 0.8501, 0.0630, -0.2736, 0.0974,
      0.9569;
  Rd << 0.4569, -0.0903, 0.8849, 0.7532, 0.5684, -0.3309, -0.4731, 0.8178,
      0.3278;

  auto eR = 0.5 * (R.transpose() * Rd - Rd.transpose() * R);
  std::cout << "hat(eR) = " << eR << std::endl;
  std::cout << "eR = " << eR(2, 1) << " " << eR(0, 2) << " " << eR(1, 0) << std::endl;
  std::cout << "R.error(Rd) = [ " << R.error(Rd).transpose() << " ]^T" << std::endl;
  std::cout << "SO3d::error(R, Rd) = [ " << SO3::error(R, Rd).transpose() << " ]^T" << std::endl;

  std::cout << "***************************************" << std::endl;
  TSE3 x, xd;

  const Eigen::Vector3d p = Eigen::Vector3d::Random();
  const Eigen::Vector3d v = Eigen::Vector3d::Random();
  const Eigen::Vector3d Om = Eigen::Vector3d::Random();
  const auto pd = Eigen::Vector3d::Ones();
  const auto vd = Eigen::Vector3d::Ones();
  const auto Omd = Eigen::Vector3d::Zero();
  // Eigen::Matrix<float, 3, 1> p, v, Om;
  // p << 0.8, -0.4, 7.4;
  // v << 0.1, 21, -0.9;
  // Om << 1, 2, 3;

  std::cout << "p: " << p.transpose() << "    v:" << v.transpose() << std::endl;

  x.print();
  x.position << p;
  x.print();
  x.velocity = v;
  x.R = R;
  x.Omega = Om;

  xd.position = pd;
  xd.velocity = vd;
  xd.Omega = Omd;
  xd.R = Eigen::Matrix<double, 3, 3>::Identity();

  auto eR2_hat = 0.5 * (xd.R.transpose() * x.R - x.R.transpose() * xd.R);
  Eigen::Matrix<double, 3, 1> eR2;
  eR2 << eR2_hat(2, 1), eR2_hat(0, 2), eR2_hat(1, 0);

  Eigen::Matrix<double, 12, 1> outside_error;
  outside_error << (p - pd), (v - vd), eR2, x.Omega - x.R.transpose() * xd.R * xd.Omega;

  auto error = x - xd;
  auto error2 = x.error(xd);

  std::cout << "error manual:  [" << outside_error.transpose() << " ]^T" << std::endl;
  std::cout << "error2 :  [" << error2.transpose() << " ]^T" << std::endl;
  std::cout << "error using -: [" << error.transpose() << " ]^T" << std::endl;
  std::cout << "norm of the difference = " << (error - outside_error).norm() << std::endl;
  std::cout << "\n" << "testing done" << std::endl;
  // std::cout << R.config_error(Rd) << std::endl;
  // std::cout << SO3d::config_error(R, Rd) << std::endl;

  // std::cout << "***********************" << std::endl;
  // SO3Controller<float> geo_;

  // std::cout << geo_.state().R << "\n"
  //           << geo_.state().Omega << std::endl;
}
