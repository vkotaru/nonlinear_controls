#ifndef NONLINEAR_CONTROLS_MPC_BASE_H
#define NONLINEAR_CONTROLS_MPC_BASE_H
#include "deque"
#include "qpoases_eigen.hpp"
#include "vector"
namespace nonlinear_control {

struct MPCGains {
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Q, R, P;
};

struct MPCBounds {
    Eigen::Matrix<double, Eigen::Dynamic, 1> lb, ub;
};

struct MPCProblem {
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> A, B;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> lb, ub;
    int N, nx, nu;
};

class MPCBase {
   private:
    QPOasesEigen *qp;

    bool isLTI;
    const int nx, nu, N;
    int nVars, nCons;

   public:
    MPCBase(const bool _isLTI, const int _N, const int _nx, const int _nu);
    ~MPCBase();

    MPCProblem problem_;
    MPCBounds state, input;
    MPCGains gains;

    // temporary (move to private after debug is complete)
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Au, Ax, Af, bx, bu, bf;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Sx, Su, Qbar, Rbar, H, F;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Ulb, Uub, Xlb, Xub;

    void setup();
    void construct();
    void debug_print();
    void print() {
      qp->print();
    }
    Eigen::Matrix<double, Eigen::Dynamic, 1> update(const bool verbose, const Eigen::Matrix<double, Eigen::Dynamic, 1> &x0);
};

}  // namespace nonlinear_control

#endif  // NONLINEAR_CONTROLS_MPC_BASE_H