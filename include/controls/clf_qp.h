#ifndef NLC_CLF_QP_H
#define NLC_CLF_QP_H

#include "common/cvxgen_interface.h"
#include "common/qp_interface.h"
namespace nonlinear_controls {

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

}  // namespace nonlinear_controls
#endif  // NLC_CLF_QP_H