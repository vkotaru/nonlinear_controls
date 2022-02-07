#ifndef __NLC_CLF_QP_H__
#define __NLC_CLF_QP_H__

#include "epigraph.hpp"
#include <chrono>

namespace nonlinear_controls {

class ClfQP {
protected:
  size_t nVars;
  size_t nCons;
  bool use_penality_term = false;

  Eigen::VectorXd xref;
  Eigen::VectorXd xprev;

  Eigen::MatrixXd Qref;
  Eigen::MatrixXd Qprev;
  double p{10000};

  Eigen::VectorXd xMin, xMax;

  Eigen::MatrixXd Aineq;
  Eigen::VectorXd bineq;

  cvx::OptimizationProblem qp;
  cvx::osqp::OSQPSolver *solver;

  cvx::VectorX xVar;
  cvx::Scalar dVar;

public:
  ClfQP(const size_t NVARS, const size_t NCONS,
        const bool enable_penality = false)
      : nVars(NVARS), nCons(NCONS), use_penality_term(enable_penality) {

    if (nCons < 1) {
      std::cerr << "Alteast one constraint should be provided\n";
    }

    xref.resize(nVars);
    xprev.resize(nVars);

    Qref.resize(nVars, nVars);
    Qprev.resize(nVars, nVars);

    xMin.resize(nVars);
    xMax.resize(nVars);

    // setup the solver
    xVar = qp.addVariable("x", nVars);
    dVar = qp.addVariable("d");
    qp.addCostTerm((xVar - cvx::dynpar(xref)).transpose() * cvx::par(Qref) *
                   (xVar - cvx::dynpar(xref)));
    qp.addCostTerm((xVar - cvx::dynpar(xprev)).transpose() * cvx::par(Qprev) *
                   (xVar - cvx::dynpar(xprev)));
    if (use_penality_term) {
      qp.addCostTerm(cvx::par(p) * dVar * dVar);
    }

    for (size_t i = 0; i < nVars; i++) {
      qp.addConstraint(cvx::box(xMin[i], xVar[i], xMax[i]));
    }
    Aineq.resize(nVars, nCons);
    bineq.resize(nCons);
    qp.addConstraint(
        cvx::greaterThan(cvx::dynpar(Aineq) * xVar, cvx::dynpar(bineq)));
  }

  ~ClfQP() = default;

  std::optional<Eigen::VectorXd> solve() { return std::nullopt; }
};

} // namespace nonlinear_controls
#endif // __NLC_CLF_QP_H__