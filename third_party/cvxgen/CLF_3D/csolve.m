% csolve  Solves a custom quadratic program very rapidly.
%
% [vars, status] = csolve(params, settings)
%
% solves the convex optimization problem
%
%   minimize(quad_form(x, Q) + c'*x)
%   subject to
%     A*x <= b
%     xlb <= x
%     x <= xub
%
% with variables
%        x   3 x 1
%
% and parameters
%        A   1 x 3
%        Q   3 x 3    PSD
%        b   1 x 1
%        c   3 x 1
%      xlb   3 x 1
%      xub   3 x 1
%
% Note:
%   - Check status.converged, which will be 1 if optimization succeeded.
%   - You don't have to specify settings if you don't want to.
%   - To hide output, use settings.verbose = 0.
%   - To change iterations, use settings.max_iters = 20.
%   - You may wish to compare with cvxsolve to check the solver is correct.
%
% Specify params.A, ..., params.xub, then run
%   [vars, status] = csolve(params, settings)
% Produced by CVXGEN, 2020-07-30 04:41:53 -0400.
% CVXGEN is Copyright (C) 2006-2017 Jacob Mattingley, jem@cvxgen.com.
% The code in this file is Copyright (C) 2006-2017 Jacob Mattingley.
% CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial
% applications without prior written permission from Jacob Mattingley.

% Filename: csolve.m.
% Description: Help file for the Matlab solver interface.
