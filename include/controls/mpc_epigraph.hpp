#ifndef _NONLINEAR_CONTROLS_MPC_EPIGRAPH_HPP_
#define _NONLINEAR_CONTROLS_MPC_EPIGRAPH_HPP_

#include "controls/linear_mpc.hpp"
#include "epigraph.hpp"

namespace nonlinear_controls {

class MPCEpigraph : public LinearMPCBase {
protected:
  // CVX problem setup
  cvx::OptimizationProblem qp;
  cvx::MatrixX x, u;

  cvx::osqp::OSQPSolver *solver{};

public:
  MPCEpigraph(const long &N, const long &nx, const long &nu) : LinearMPCBase(N, nx, nu) {
    // Create variables
    x = qp.addVariable("x", nx, N + 1);
    u = qp.addVariable("u", nu, N);
    x0.resize(nx);
    xref.resize(nx);
  }

  ~MPCEpigraph() = default;

  void construct() override {
    // Dynamics
    for (long t = 0; t < N; t++) {
      // dynamics
      qp.addConstraint(cvx::equalTo(x.col(t + 1), cvx::par(A) * x.col(t) + cvx::par(B) * u.col(t)));

      // cost
      qp.addCostTerm((x.col(t) - cvx::par(xref)).transpose() * cvx::par(Q) * (x.col(t) - cvx::par(xref)));
      qp.addCostTerm(u.col(t).transpose() * cvx::par(R) * u.col(t));


      // State and control limits
//      qp.addConstraint(cvx::box(cvx::par(state_bnds_.lb), x.col(t), cvx::par(state_bnds_.ub)));
//      qp.addConstraint(cvx::box(cvx::par(input_bnds_.lb), u.col(t), cvx::par(input_bnds_.ub)));
    }
    qp.addCostTerm((x.col(N) - cvx::par(xref)).transpose() * cvx::par(P) * (x.col(N) - cvx::par(xref)));
//    qp.addConstraint(cvx::box(cvx::par(Xlb), x.col(N), cvx::par(Xub)));


    // Boundary constraints
    qp.addConstraint(cvx::equalTo(x.col(0), cvx::dynpar(x0)));
    // qp.addConstraint(equalTo(x.col(T), 0.)); // to be free

    solver = new cvx::osqp::OSQPSolver(qp);
//    solver->setAlpha(1.0);
  }

  MatrixXd run(const VectorXd &x0_, const VectorXd &xd_) override {

    this->x0 = x0_;
    this->xref = xd_;
    solver->solve(verbose);
    Eigen::MatrixXd x_sol = eval(x);
    Eigen::MatrixXd u_sol = eval(u);

    return u_sol;
  }
  void reset() override {

  }

};

} // namespace nonlinear_controls
#endif //_NONLINEAR_CONTROLS_MPC_EPIGRAPH_HPP_
