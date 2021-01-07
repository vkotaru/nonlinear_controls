#include "controls/clf_qp.h"
namespace nonlinear_controls {

namespace clf3D {
#include "CLF_3D/ldl.c"
#include "CLF_3D/matrix_support.c"
#include "CLF_3D/solver.c"
#include "CLF_3D/util.c"
Workspace work;
Vars vars;
Params params;
Settings settings;
}  // namespace clf3D
ClfQP::ClfQP() : QPInterface() {
    this->problem_.H.resize(3, 3);
    this->problem_.f.resize(3, 1);
    this->problem_.A.resize(1, 3);
    this->problem_.b.resize(1, 1);
    this->problem_.xlb.resize(3, 1);
    this->problem_.xub.resize(3, 1);
}

ClfQP::~ClfQP() {
}

void ClfQP::setup() {
    clf3D::set_defaults();
    clf3D::setup_indexing();
    clf3D::settings.verbose = 1;
    clf3D::vars.x[0] = 0.0;
    clf3D::vars.x[1] = 0.0;
    clf3D::vars.x[2] = 0.0;
}

void ClfQP::solve() {
    // transfer data
    Eigen::Map<Eigen::Matrix<double, 3, 3>>(&clf3D::params.Q[0], 3, 3) = this->problem_.H;
    Eigen::Map<Eigen::Matrix<double, 3, 1>>(&clf3D::params.c[0], 3, 1) = this->problem_.f;
    Eigen::Map<Eigen::Matrix<double, 1, 3>>(&clf3D::params.A[0], 1, 3) = this->problem_.A;
    Eigen::Map<Eigen::Matrix<double, 1, 1>>(&clf3D::params.b[0], 1, 1) = this->problem_.b;
    Eigen::Map<Eigen::Matrix<double, 3, 1>>(&clf3D::params.xlb[0], 3, 1) = this->problem_.xlb;
    Eigen::Map<Eigen::Matrix<double, 3, 1>>(&clf3D::params.xub[0], 3, 1) = this->problem_.xub;

    if (clf3D::settings.verbose > 0) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                printf("Q[%d, %d] = %f ", i, j, clf3D::params.Q[i + j * 3]);
            }
            printf("c[%d]=%f ", i, clf3D::params.c[i]);
            printf("A[%d]=%f ", i, clf3D::params.A[i]);
            printf("xlb[%d]=%f ", i, clf3D::params.xlb[i]);
            printf("xub[%d]=%f\n", i, clf3D::params.xub[i]);
        }
        printf("b=%f\n", clf3D::params.b[0]);
    }

    // solve
    num_iters = clf3D::solve();
}
 
Eigen::Vector3d ClfQP::getOptimizer() {
    return Eigen::Vector3d(clf3D::vars.x[0], clf3D::vars.x[1], clf3D::vars.x[2]);
}

}  // namespace nonlinear_controls
