#ifndef NONLINEAR_CONTROLS_MPC_BASE2_H
#define NONLINEAR_CONTROLS_MPC_BASE2_H

#include "deque"
#include "vector"
#include <assert.h>
#include "mpc_base.h"
#include "qpoases_eigen.hpp"


namespace nonlinear_control {

class LinearMPCBase2 {
   // mpc without substitution
   private:
    QPOasesEigen *QP;

    bool isLTI;
    const int nx, nu, N;
    int nVars, nCons;

    MPCProblem problem_;
    MPCBounds state, input;
    MPCGains gains;

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> G0, E0, H, F;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Ulb, Uub, Xlb, Xub, lbA, ubA;
   public:
    LinearMPCBase2(const bool _isLTI, const int _N, const int _nx, const int _nu);
    ~LinearMPCBase2();

    void set_max_iter(const int MAX_ITER) {
        QP->data_.nWSR = MAX_ITER;
    }

    void set_gains(Eigen::MatrixXd Q, Eigen::MatrixXd P, Eigen::MatrixXd R);
    void setup();
    void construct();
    void debug_print();
    Eigen::Matrix<double, Eigen::Dynamic, 1> update(const bool verbose, const Eigen::Matrix<double, Eigen::Dynamic, 1> &x0);
};

}  // namespace nonlinear_control

#endif  // NONLINEAR_CONTROLS_MPC_BASE2_H