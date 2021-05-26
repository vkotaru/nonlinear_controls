//
// Created by kotaru on 5/26/21.
//

#ifndef NONLINEAR_CONTROLS_QPSWIFT_EIGEN_H
#define NONLINEAR_CONTROLS_QPSWIFT_EIGEN_H
#include "Prime.h"
#include "quadprog/quadprog.h"
#include <iostream>

namespace nonlinear_controls {

class QPSwiftEigen : public QuadProg{
protected:
public:
  QPSwiftEigen(const int &n, const int &m, const int &p);
  ~QPSwiftEigen() = default;

  int solve() override;
};

} // namespace nonlinear_controls

#endif // NONLINEAR_CONTROLS_QPSWIFT_EIGEN_H
