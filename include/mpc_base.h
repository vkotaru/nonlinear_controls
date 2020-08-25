#ifndef NONLINEAR_CONTROLS_MPC_BASE_H
#define NONLINEAR_CONTROLS_MPC_BASE_H
#include "deque"
#include "qpoases_eigen.hpp"
#include "vector"
namespace nonlinear_control {
//template<typename T>

template <typename T>
struct MPCGains {
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Q, R, P;
};

template <typename T>
struct MPCBounds {
    Eigen::Matrix<T, Eigen::Dynamic, 1> lb, ub;
};

template <typename T>
struct MPCProblem {
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A, B;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> lb, ub;
    int N, nx, nu;
};

class LinearMPCBase {
  private:
    QPOasesEigen* qp;

    bool isLTI;
    const int nx, nu, N;
    int nVars, nCons;

  public:
    LinearMPCBase(const bool _isLTI, const int _N, const int _nx, const int _nu);
    ~LinearMPCBase();

    MPCProblem<double> problem_;
    MPCBounds<double> state, input;
    MPCGains<double> gains;

    // temporary (move to private after debug is complete)
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Au, Ax, Af, bx, bu, bf;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Sx, Su, Qbar, Rbar, H, F;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Ulb, Uub, Xlb, Xub;

    void construct();
    void debug_print();
    Eigen::Matrix<double, Eigen::Dynamic, 1> update(const bool verbose, const Eigen::Matrix<double, Eigen::Dynamic, 1> x0);
};

}  // namespace nonlinear_control

#endif  // NONLINEAR_CONTROLS_MPC_BASE_H
