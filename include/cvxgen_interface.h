#ifndef NLC_CVXGEN_INTERFACE_H
#define NLC_CVXGEN_INTERFACE_H

namespace nonlinear_control {

#ifndef NLC_CVXGEN_CLF_3D
#define NLC_CVXGEN_CLF_3D
namespace clf3D {
#include "CLF_3D/solver.h"
Workspace work;
Vars vars;
Params params;
Settings settings;
#include "CLF_3D/ldl.c"
#include "CLF_3D/matrix_support.c"
#include "CLF_3D/solver.c"
#include "CLF_3D/util.c"
}  // namespace clf3D
#endif  // NLC_CVXGEN_CLF_3D


#ifndef NLC_CVXGEN_MPC_3D
#define NLC_CVXGEN_MPC_3D
namespace mpc3D {
#include "MPC_3D/solver.h"
Workspace work;
Vars vars;
Params params;
Settings settings;
#include "MPC_3D/ldl.c"
#include "MPC_3D/matrix_support.c"
#include "MPC_3D/solver.c"
#include "MPC_3D/util.c"
}  // namespace mpc3D
#endif  // NLC_CVXGEN_MPC_3D


}  // namespace nonlinear_control

#endif  // NLC_CVXGEN_INTERFACE_H