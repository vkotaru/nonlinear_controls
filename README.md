# nonlinear_controls
Nonlinear controls, linear controls, lyapunov based controls, barrier functions, MPC, geometric etc. 


### MPC Usage
Makes use of `Eigen` and `qpOASES` (currently, the applications are very specific)

#### LTI Dynamics: `LinearMPC`
```
namespace nlc = nonlinear_controls;

// solvers
nlc::MPCQPOases
nlc::MPCEpigraph

nlc::MPCQPOases mpc_(N, nx, nu);
mpc_.init_dynamics(A, B);
mpc_.set_gains(Q, P, R);
mpc_.set_input_bounds(input_lb, input_ub);
mpc_.set_state_bounds(state_lb, state_ub);

mpc_.construct(); // constructs the necessary cost and constraint matrices 
for(;;) {
  zOpt = mpc_.run(err_state);
  uOpt = zOpt.block(0, 0, nu, 1);
}
```



