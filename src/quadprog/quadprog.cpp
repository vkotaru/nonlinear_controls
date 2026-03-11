//
// Created by kotaru on 5/26/21.
//
#include "quadprog/quadprog.h"

#include <stdexcept>

namespace nonlinear_controls {

QuadProg::QuadProg(const int& n, const int& m, const int& p) : nv(n), ni(m), ne(p) {
  H.resize(nv, nv);
  f.resize(nv);
  A.resize(ni, nv);
  b.resize(ni);
  Aeq.resize(ne, nv);
  beq.resize(ne);
  lb.resize(nv);
  ub.resize(nv);
  x0.resize(nv);
  x_opt.resize(nv);
}

QuadProg::~QuadProg() = default;

int QuadProg::solve() {
  throw std::logic_error("QuadProg::solve is not implemented");
}

}  // namespace nonlinear_controls
