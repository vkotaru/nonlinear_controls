#ifndef NLC_CLF_QP_H
#define NLC_CLF_QP_H

#include "cvxgen_interface.h"
#include "qp_interface.h"
namespace nonlinear_control {

class ClfQP : public QPInterface {
   protected:
    long num_iters;

   public:
    ClfQP();
    ~ClfQP();

    virtual void setup() override;
    virtual void solve() override;
    Eigen::Vector3d getOptimizer();
};

}  // namespace nonlinear_control
#endif  // NLC_CLF_QP_H