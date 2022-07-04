
#include "OsqpEigen/OsqpEigen.h"
#include <Eigen/Dense>
#include <iostream>

#include "dynamics/point_mass.hpp"
#include "matplotlibcpp.h"
#include "common/utils.hpp"
#include "common/log.hpp"

namespace plt = matplotlibcpp;
namespace nlc = nonlinear_controls;

void setDynamicsMatrices(Eigen::Matrix<double, 6, 6> &a, Eigen::Matrix<double, 6, 3> &b) {
//  a << 1., 0., 0., 0., 0., 0., 0.1, 0., 0., 0., 0., 0.,
//      0., 1., 0., 0., 0., 0., 0., 0.1, 0., 0., 0., 0.,
//      0., 0., 1., 0., 0., 0., 0., 0., 0.1, 0., 0., 0.,
//      0.0488, 0., 0., 1., 0., 0., 0.0016, 0., 0., 0.0992, 0., 0.,
//      0., -0.0488, 0., 0., 1., 0., 0., -0.0016, 0., 0., 0.0992, 0.,
//      0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0.0992,
//      0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.,
//      0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0.,
//      0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0.,
//      0.9734, 0., 0., 0., 0., 0., 0.0488, 0., 0., 0.9846, 0., 0.,
//      0., -0.9734, 0., 0., 0., 0., 0., -0.0488, 0., 0., 0.9846, 0.,
//      0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.9846;
//
//  b << 0., -0.0726, 0., 0.0726,
//      -0.0726, 0., 0.0726, 0.,
//      -0.0152, 0.0152, -0.0152, 0.0152,
//      -0., -0.0006, -0., 0.0006,
//      0.0006, 0., -0.0006, 0.0000,
//      0.0106, 0.0106, 0.0106, 0.0106,
//      0, -1.456, 0., 1.456,
//      -1.456, 0., 1.456, 0.,
//      -0.3049, 0.3049, -0.3049, 0.3049,
//      -0., -0.0236, 0., 0.0236,
//      0.0236, 0., -0.0236, 0.,
//      0.2107, 0.2107, 0.2107, 0.2107;
}

void setInequalityConstraints(Eigen::Matrix<double, 6, 1> &xMax, Eigen::Matrix<double, 6, 1> &xMin,
                              Eigen::Matrix<double, 3, 1> &uMax, Eigen::Matrix<double, 3, 1> &uMin) {
  double u0 = 10.5916;

  // input inequality constraints
  uMin << 9.6 - u0,
      9.6 - u0,
      9.6 - u0,
      9.6 - u0;

  uMax << 13 - u0,
      13 - u0,
      13 - u0,
      13 - u0;

  // state inequality constraints
  xMin << -M_PI / 6, -M_PI / 6, -OsqpEigen::INFTY, -OsqpEigen::INFTY, -OsqpEigen::INFTY, -1.,
      -OsqpEigen::INFTY, -OsqpEigen::INFTY, -OsqpEigen::INFTY, -OsqpEigen::INFTY,
      -OsqpEigen::INFTY, -OsqpEigen::INFTY;

  xMax << M_PI / 6, M_PI / 6, OsqpEigen::INFTY, OsqpEigen::INFTY, OsqpEigen::INFTY,
      OsqpEigen::INFTY, OsqpEigen::INFTY, OsqpEigen::INFTY, OsqpEigen::INFTY,
      OsqpEigen::INFTY, OsqpEigen::INFTY, OsqpEigen::INFTY;
}

void setWeightMatrices(Eigen::DiagonalMatrix<double, 6> &Q, Eigen::DiagonalMatrix<double, 3> &R) {
  Q.diagonal() << 0, 0, 10., 10., 10., 10., 0, 0, 0, 5., 5., 5.;
  R.diagonal() << 0.1, 0.1, 0.1, 0.1;
}

void castMPCToQPHessian(const Eigen::DiagonalMatrix<double, 6> &Q,
                        const Eigen::DiagonalMatrix<double, 3> &R,
                        int mpcWindow,
                        Eigen::SparseMatrix<double> &hessianMatrix) {

  hessianMatrix.resize(6 * (mpcWindow + 1) + 3 * mpcWindow, 6 * (mpcWindow + 1) + 3 * mpcWindow);

  //populate hessian matrix
  for (int i = 0; i < 6 * (mpcWindow + 1) + 3 * mpcWindow; i++) {
    if (i < 6 * (mpcWindow + 1)) {
      int posQ = i % 6;
      float value = Q.diagonal()[posQ];
      if (value != 0)
        hessianMatrix.insert(i, i) = value;
    } else {
      int posR = i % 3;
      float value = R.diagonal()[posR];
      if (value != 0)
        hessianMatrix.insert(i, i) = value;
    }
  }
}

void castMPCToQPGradient(const Eigen::DiagonalMatrix<double, 6> &Q,
                         const Eigen::Matrix<double, 6, 1> &xRef,
                         int mpcWindow,
                         Eigen::VectorXd &gradient) {

  Eigen::Matrix<double, 6, 1> Qx_ref;
  Qx_ref = Q * (-xRef);

  // populate the gradient vector
  gradient = Eigen::VectorXd::Zero(6 * (mpcWindow + 1) + 3 * mpcWindow, 1);
  for (int i = 0; i < 6 * (mpcWindow + 1); i++) {
    int posQ = i % 6;
    float value = Qx_ref(posQ, 0);
    gradient(i, 0) = value;
  }
}

void castMPCToQPConstraintMatrix(const Eigen::Matrix<double, 6, 6> &dynamicMatrix,
                                 const Eigen::Matrix<double, 6, 3> &controlMatrix,
                                 int mpcWindow,
                                 Eigen::SparseMatrix<double> &constraintMatrix) {
  constraintMatrix.resize(6 * (mpcWindow + 1) + 6 * (mpcWindow + 1) + 3 * mpcWindow,
                          6 * (mpcWindow + 1) + 3 * mpcWindow);

  // populate linear constraint matrix
  for (int i = 0; i < 6 * (mpcWindow + 1); i++) {
    constraintMatrix.insert(i, i) = -1;
  }

  for (int i = 0; i < mpcWindow; i++)
    for (int j = 0; j < 6; j++)
      for (int k = 0; k < 6; k++) {
        float value = dynamicMatrix(j, k);
        if (value != 0) {
          constraintMatrix.insert(6 * (i + 1) + j, 6 * i + k) = value;
        }
      }

  for (int i = 0; i < mpcWindow; i++)
    for (int j = 0; j < 6; j++)
      for (int k = 0; k < 3; k++) {
        float value = controlMatrix(j, k);
        if (value != 0) {
          constraintMatrix.insert(6 * (i + 1) + j, 3 * i + k + 6 * (mpcWindow + 1)) = value;
        }
      }

  for (int i = 0; i < 6 * (mpcWindow + 1) + 3 * mpcWindow; i++) {
    constraintMatrix.insert(i + (mpcWindow + 1) * 6, i) = 1;
  }
}

void castMPCToQPConstraintVectors(const Eigen::Matrix<double, 6, 1> &xMax, const Eigen::Matrix<double, 6, 1> &xMin,
                                  const Eigen::Matrix<double, 3, 1> &uMax, const Eigen::Matrix<double, 3, 1> &uMin,
                                  const Eigen::Matrix<double, 6, 1> &x0,
                                  int mpcWindow, Eigen::VectorXd &lowerBound, Eigen::VectorXd &upperBound) {
  // evaluate the lower and the upper inequality vectors
  Eigen::VectorXd lowerInequality = Eigen::MatrixXd::Zero(6 * (mpcWindow + 1) + 3 * mpcWindow, 1);
  Eigen::VectorXd upperInequality = Eigen::MatrixXd::Zero(6 * (mpcWindow + 1) + 3 * mpcWindow, 1);
  for (int i = 0; i < mpcWindow + 1; i++) {
    lowerInequality.block(6 * i, 0, 6, 1) = xMin;
    upperInequality.block(6 * i, 0, 6, 1) = xMax;
  }
  for (int i = 0; i < mpcWindow; i++) {
    lowerInequality.block(3 * i + 6 * (mpcWindow + 1), 0, 3, 1) = uMin;
    upperInequality.block(3 * i + 6 * (mpcWindow + 1), 0, 3, 1) = uMax;
  }

  // evaluate the lower and the upper equality vectors
  Eigen::VectorXd lowerEquality = Eigen::MatrixXd::Zero(6 * (mpcWindow + 1), 1);
  Eigen::VectorXd upperEquality;
  lowerEquality.block(0, 0, 6, 1) = -x0;
  upperEquality = lowerEquality;
  lowerEquality = lowerEquality;

  // merge inequality and equality vectors
  lowerBound = Eigen::MatrixXd::Zero(2 * 6 * (mpcWindow + 1) + 3 * mpcWindow, 1);
  lowerBound << lowerEquality,
      lowerInequality;

  upperBound = Eigen::MatrixXd::Zero(2 * 6 * (mpcWindow + 1) + 3 * mpcWindow, 1);
  upperBound << upperEquality,
      upperInequality;
}

void updateConstraintVectors(const Eigen::Matrix<double, 6, 1> &x0,
                             Eigen::VectorXd &lowerBound, Eigen::VectorXd &upperBound) {
  lowerBound.block(0, 0, 6, 1) = -x0;
  upperBound.block(0, 0, 6, 1) = -x0;
}

double getErrorNorm(const Eigen::Matrix<double, 6, 1> &x,
                    const Eigen::Matrix<double, 6, 1> &xRef) {
  // evaluate the error
  Eigen::Matrix<double, 6, 1> error = x - xRef;

  // return the norm
  return error.norm();
}

int main() {

  // simulation parameters
  double h = 1. / 200.;
  const double simulate_for_seconds = 10;
  int MAX_ITER_STEPS = int(simulate_for_seconds / h);


  // point-mass dynamics
  nlc::PointMass3D dynamics_{h};


  // set the preview window
  int mpcWindow = 20;
  int nx = 6, nu = 3, N = mpcWindow;

  // allocate the dynamics matrices
  Eigen::Matrix<double, 6, 6> A;
  Eigen::Matrix<double, 6, 3> B;

  // allocate the constraints vector
  Eigen::Matrix<double, 6, 1> xMax;
  Eigen::Matrix<double, 6, 1> xMin;
  Eigen::Matrix<double, 3, 1> uMax;
  Eigen::Matrix<double, 3, 1> uMin;

  // allocate the weight matrices
  Eigen::DiagonalMatrix<double, 6> Q;
  Eigen::DiagonalMatrix<double, 3> R;

  // allocate the initial and the reference state space
  Eigen::Matrix<double, 6, 1> x0;
  Eigen::Matrix<double, 6, 1> xRef;

  // allocate QP problem matrices and vectores
  Eigen::SparseMatrix<double> hessian;
  Eigen::VectorXd gradient;
  Eigen::SparseMatrix<double> linearMatrix;
  Eigen::VectorXd lowerBound;
  Eigen::VectorXd upperBound;

  // set the initial and the desired states
//  x0 << 0, 0, 0, 0, 0, 0;
  x0 <<  -0.012834, 0.94555, -0.414966,0,0,0;
  xRef << 0, 0, 1, 0, 0, 0;

  // set MPC problem quantities
  A = dynamics_.A();
  B = dynamics_.B();
  Q.diagonal() << 20, 20, 10, 20, 20, 20;
  R.diagonal() << 1., 1., 1.;

  xMin << -2., -2.5, 0., -2, -2, -2;
  xMax << 2., 2.5, 4., 2., 2., 2;
  uMin << -20, -20, -20;
  uMax << 20, 20, 20;

  // cast the MPC problem as QP problem
  castMPCToQPHessian(Q, R, mpcWindow, hessian);
  castMPCToQPGradient(Q, xRef, mpcWindow, gradient);
  castMPCToQPConstraintMatrix(A, B, mpcWindow, linearMatrix);
  castMPCToQPConstraintVectors(xMax, xMin, uMax, uMin, x0, mpcWindow, lowerBound, upperBound);

  // instantiate the solver
  OsqpEigen::Solver solver;

  // settings
  //solver.settings()->setVerbosity(false);
  solver.settings()->setWarmStart(true);

  // set the initial data of the QP solver
  solver.data()->setNumberOfVariables(6 * (mpcWindow + 1) + 3 * mpcWindow);
  solver.data()->setNumberOfConstraints(2 * 6 * (mpcWindow + 1) + 3 * mpcWindow);
  if (!solver.data()->setHessianMatrix(hessian))
    return 1;
  if (!solver.data()->setGradient(gradient))
    return 1;
  if (!solver.data()->setLinearConstraintsMatrix(linearMatrix))
    return 1;
  if (!solver.data()->setLowerBound(lowerBound))
    return 1;
  if (!solver.data()->setUpperBound(upperBound))
    return 1;

  // instantiate the solver
  if (!solver.initSolver())
    return 1;

  // controller input and QPSolution vector
  Eigen::Vector3d ctr;
  Eigen::VectorXd QPSolution;

  // number of iteration steps
  int numberOfSteps = 50;

  for (int i = 0; i < numberOfSteps; i++) {

    auto start = nlc::utils::get_current_time();
    // solve the QP problem
    if (solver.solveProblem() != OsqpEigen::ErrorExitFlag::NoError)
      return 1;
    auto end =  nlc::utils::get_current_time();
    nlc::Logger::INFO("time elapsed "+std::to_string(double(end-start)*1.e-6));

    // get the controller input
    QPSolution = solver.getSolution();
    ctr = QPSolution.block(6 * (mpcWindow + 1), 0, 3, 1);

    // save data into file
    auto x0Data = x0.data();

    // propagate the model
    x0 = A * x0 + B * ctr;

    // update the constraint bound
    updateConstraintVectors(x0, lowerBound, upperBound);
    if (!solver.updateBounds(lowerBound, upperBound))
      return 1;
  }
  return 0;
}
