#ifndef NONLINEAR_CONTROLS_CONTROLS_MPC_QPOASES_HPP
#define NONLINEAR_CONTROLS_CONTROLS_MPC_QPOASES_HPP

#include <deque>
#include <memory>
#include <sstream>
#include <vector>

#include "common/log.hpp"
#include "common/qpoases_eigen.hpp"
#include "controls/linear_mpc.hpp"

namespace nonlinear_controls {

class MPCQPOases : public LinearMPC {
protected:
  std::unique_ptr<qpOASES::SQProblem> solver;
  qpOASES::int_t nWSR = 10;
  qpOASES::Options options;
  qpOASES::real_t cpu_time_limit{};
  QPOasesData data_;
  qpOASES::real_t *Hc{}, *Ac{}, *gc{}, *lbc{}, *ubc{}, *lbAc{}, *ubAc{};
  bool solver_initialized{false};

  void create_solver() {
    solver = std::make_unique<qpOASES::SQProblem>(n_vars, n_cons);
    options.setToMPC();
    options.printLevel = qpOASES::PL_LOW;
    solver->setOptions(options);
    nWSR = 1e7;
  }

  std::optional<MatrixXd> solve(const VectorXd& _x0) {
    this->updateCArrays(_x0);
    /**
     * NOTE: the qpOASES solver works only if the
     * "n" is locally defined and passed to the solver
     * otherwise, the solver becomes infeasible
     *
     * Please refer qpOASES GitHub/forums for discussion on this
     * https://github.com/coin-or/qpOASES/issues/83#issuecomment-762053495
     *
     * nWSR is currently used as an input and an output variable.
     * This means that the maximum number of allowed working set changes will decrease monotonically
     * and will eventually trigger the error. Just make sure to set nWSR before each call to
     *hotstart.
     **/
    int n = nWSR;
    qpOASES::returnValue sol;
    if (!solver_initialized) {
      sol = solver->init(Hc, gc, Ac, lbc, ubc, lbAc, ubAc, n, 0);
      solver_initialized = true;
    } else {
      sol = solver->hotstart(Hc, gc, Ac, lbc, ubc, lbAc, ubAc, n, 0);
    }

    qpOASES::real_t x_opt[n_vars];
    solver->getPrimalSolution(x_opt);

    if (verbose) {
      if (sol == qpOASES::SUCCESSFUL_RETURN) {
        Logger::SUCCESS("QP solved: SUCCESSFUL_RETURN");
      } else if (sol == qpOASES::RET_MAX_NWSR_REACHED) {
        Logger::WARN("QP solver: RET_MAX_NWSR_REACHED");
      } else if (sol == qpOASES::RET_INIT_FAILED) {
        Logger::ERROR("QP solver: RET_INIT_FAILED");
      } else {
        Logger::ERROR("QP solver: Unknown return value");
      }
    }

    Eigen::Matrix<double, Eigen::Dynamic, 1> x_opt_vec;
    x_opt_vec = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>>(x_opt, n_vars, 1);

    Uarray = x_opt_vec;  // saved for generating the trajectory
    return std::optional<MatrixXd>{x_opt_vec};
  }

public:
  MPCQPOases(const long& N, const long& nx, const long& nu) : LinearMPC(N, nx, nu) {
    /// QP solver Setup
    create_solver();
  }

  void set_input_bounds(const VectorXd& lb, const VectorXd& ub) override {
    LinearMPC::set_input_bounds(lb, ub);
    /// set only once
    data_.lb = Ulb;
    data_.ub = Uub;
  }
  void construct() override {
    LinearMPC::construct();

    /// set only once
    data_.H = H;
    data_.lb = Ulb;
    data_.ub = Uub;
    data_.A = Su;
    qp_initialized = false;
  }
  void print() override {
    std::ostringstream oss;
    oss << "--------------------------------------------------\n"
        << "*           QP setup       *\n"
        << "--------------------------------------------------\n"
        << "H: \n"
        << data_.H << "\nf: \n"
        << data_.g.transpose() << "\n\nlbA: \n"
        << data_.lbA.transpose() << "\nA: \n"
        << data_.A << "\nubA: \n"
        << data_.ubA.transpose() << "\n\nlb: " << data_.lb.transpose()
        << "\nub: " << data_.ub.transpose()
        << "\n-------------------*****---------------------------";
    Logger::INFO(oss.str());
  }
  void updateCArrays(const VectorXd& _x0) override {
    data_.lbA = (Xlb - Sx * _x0);
    data_.ubA = (Xub - Sx * _x0);
    data_.g = (2 * F.transpose() * _x0);
    //  print();
    Hc = data_.H.transpose().data();
    Ac = data_.A.transpose().data();
    gc = data_.g.transpose().data();
    lbc = data_.lb.transpose().data();
    ubc = data_.ub.transpose().data();
    lbAc = data_.lbA.transpose().data();
    ubAc = data_.ubA.transpose().data();
  }
  std::optional<MatrixXd> run(const VectorXd& _x0, const VectorXd& xd_) override {
    if (!valid_state(_x0)) {
      Logger::ERROR("Infeasible initial condition!");
      return std::nullopt;
    }
    this->x0 = _x0;
    this->xref = xd_;

    Xlb = (state_bnds_.lb - xd_).replicate(N + 1, 1);
    Xub = (state_bnds_.ub - xd_).replicate(N + 1, 1);
    return solve(_x0 - xd_);
  }
  void reset() override { solver_initialized = false; }
  Eigen::MatrixXd X() override {
    // use this function only after the solve function is implemented
    Eigen::MatrixXd x_traj;
    x_traj.resize(nx, N);
    x_traj.setZero();
    Eigen::VectorXd x = x0 - xref;
    for (int i = 0; i < N; ++i) {
      x = A * x + B * Uarray.block(nu * i, 0, nu, 1);
      x_traj.col(i) = x + xref;
    }
    return x_traj;
  }
};

}  // namespace nonlinear_controls
#endif  // NONLINEAR_CONTROLS_CONTROLS_MPC_QPOASES_HPP
