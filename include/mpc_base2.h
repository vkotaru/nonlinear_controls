#ifndef NONLINEAR_CONTROLS_MPC_BASE2_H
#define NONLINEAR_CONTROLS_MPC_BASE2_H

#include "deque"
#include "vector"
#include <assert.h>
#include "mpc_base.h"
#include "qpoases_eigen.hpp"


namespace nonlinear_control {
template<typename T>
using MatrixX =  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template<typename T>
using VectorX =  Eigen::Matrix<T, Eigen::Dynamic, 1>;

template <typename T>
class LinearMPCBase2 {
    // mpc without substitution
  private:
    QPOasesEigen* QP;

    bool isLTI;
    const int nx, nu, N;
    int nVars, nCons;

    MPCProblem<T> problem_;
    MPCBounds<T> state, input;
    MPCGains<T> gains;

    MatrixX<T> G0, E0, H, g;
    MatrixX<T> Ulb, Uub, Xlb, Xub, lbA, ubA, tolCons;
    std::deque <MatrixX<T>> Astorage, Bstorage;
  public:
    LinearMPCBase2(const bool _isLTI, const int _N, const int _nx, const int _nu);
    ~LinearMPCBase2();
    void setup();
    void construct();
    void debug_print();
    VectorX<T> update(const VectorX<T> x0);
    VectorX<T> update(const VectorX<T> x0, const MatrixX<T> A, const MatrixX<T> B);

    // Dynamics
    void init_dynamics(MatrixX<T> A, MatrixX<T> B);
    void init_dynamics(const int N, MatrixX<T> A, MatrixX<T> B);
    void update_dynamics(MatrixX<T> A, MatrixX<T> B);

    // setters
    void set_max_iter(const int MAX_ITER) {
        QP->data_.nWSR = MAX_ITER;
    }
    void set_tol(const T tol) {
        tolCons = tol * VectorX<T>::Ones(nCons, 1);
    }
    void set_mpc_cost(MatrixX<T> _H);
    void set_mpc_cost(MatrixX<T> _H, VectorX<T> _g);
    void set_mpc_gains(MatrixX<T> Q, MatrixX<T> P, MatrixX<T> R);
    void set_mpc_gains(MPCGains<T> gains);
    void set_input_bounds(VectorX<T> lb, VectorX<T> ub);
    void set_input_bounds(MPCBounds<T> _input);
    void set_state_bounds(VectorX<T> lb, VectorX<T> ub);
    void set_state_bounds(MPCBounds<T> _state);
    void set_bounds(MPCBounds<T> _state, MPCBounds<T> _input);
    void set_mpc_bounds(); // << main bounds setter

    // getters
    MatrixX<T> getCurrentA() {
        return Astorage.at(0);
    }
    MatrixX<T> getCurrentB() {
        return Bstorage.at(0);
    }
};

}  // namespace nonlinear_control

#endif  // NONLINEAR_CONTROLS_MPC_BASE2_H
