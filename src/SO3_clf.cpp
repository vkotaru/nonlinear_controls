#include "SO3_clf.h"

namespace nonlinear_control {

template <typename T>
SO3Clf<T>::SO3Clf(const Eigen::Matrix<T, 3, 3>& J) : GeometricController<T>(), inertia_(J) {
    Eigen::EigenSolver<Eigen::Matrix<T, 3, 3> > s(this->inertia_);
    inertia_scaled_ = this->inertia_ * (1 / s.eigenvalues().real().minCoeff() ); // 0.00489992 is min eigenvalue of interia matrix

    /// QP setup
    //    solver = new qpOASES::SQProblem(4, 1);
    options.setToMPC();
    options.printLevel = qpOASES::PL_LOW;
    solver.setOptions(options);

    for (int i = 0; i < 16; i++) {
        H[i] = 0.0;
        H_new[i] = 0.0;
    }
    for (int j = 0; j < 4; j++) {
        H[4 * j + j] = 1.0;
        A[j] = 0.0;
        g[j] = 0.0;
        lb[j] = -1000;
        ub[j] = 1000;

        H_new[4 * j + j] = 1.0;
        A_new[j] = 0.0;
        g_new[j] = 0.0;
        lb_new[j] = -1000;
        ub_new[j] = 1000;
    }
    H[15] = 4e2;
    A[3] = -1;
    lbA[0] = -INFINITY;
    ubA[0] =  INFINITY;
    nWSR = 10;
    H_new[15] = 4e2;
    A_new[3] = -1;
    lbA_new[0] = -INFINITY;
    ubA_new[0] =  INFINITY;
    nWSR_new  = 10;
}

template <typename T>
SO3Clf<T>::~SO3Clf() {}

template <typename T>
void SO3Clf<T>::print_qp_setup() {
    std::cout << "--------------------------------------------------" << std::endl;
    std::cout << "*           QP setup       *" << std::endl;
    std::cout << "--------------------------------------------------" << std::endl;
    printf("Cost function: \n\n");
    for (int i = 0; i < 4; i++) {
        printf("[");
        for (int j = 0; j < 4; j++) {
            printf("%f\t", H[4 * i + j]);
        }
        printf("]\t[%f]\n", g[i]);
    }
    printf("Bounds: \n\n");
    for (int i = 0; i < 4; i++) {
        printf("[%f]\t<x%d<\t[%f]\n", lb[i], (i + 1), ub[i]);
    }
    printf("Constraints: \n\n");
    printf("[%f]<[", lbA[0]);
    for (int i = 0; i < 4; i++) printf("\t%f", A[i]);
    printf("]x<[%f]\n", ubA[0]);
    printf("\n");
    std::cout << "-------------------*****---------------------------" << std::endl;

}

template <typename T>
void SO3Clf<T>::print_qp2_setup() {
    std::cout << "--------------------------------------------------" << std::endl;
    std::cout << "*           QP_new setup       *" << std::endl;
    std::cout << "--------------------------------------------------" << std::endl;
    printf("Cost function: \n\n");
    for (int i = 0; i < 4; i++) {
        printf("[");
        for (int j = 0; j < 4; j++) {
            printf("%f\t", H_new[4 * i + j]);
        }
        printf("]\t[%f]\n", g_new[i]);
    }
    printf("Bounds: \n\n");
    for (int i = 0; i < 4; i++) {
        printf("[%f]\t<x%d<\t[%f]\n", lb_new[i], (i + 1), ub_new[i]);
    }
    printf("Constraints: \n\n");
    printf("[%f]<[", lbA_new[0]);
    for (int i = 0; i < 4; i++) printf("\t%f", A_new[i]);
    printf("]x<[%f]\n", ubA_new[0]);
    printf("\n");
    std::cout << "-------------------*****---------------------------" << std::endl;

}

template <typename T>
void SO3Clf<T>::init() {
    std::cout << "/////////////////////////////////////////" << std::endl;
    std::cout << "///////////////// init //////////////////" << std::endl;
    std::cout << "inertia \n" << this->inertia_ << std::endl;
    solver.init(H, g, A, lb, ub, lbA, ubA, nWSR);
    print_qp_setup();
    std::cout << "/////////////////////////////////////////" << std::endl;
}

template <typename T>
void SO3Clf<T>::run(T dt, TSO3<T> x, TSO3<T> xd, Eigen::Matrix<T, 3, 1>& u) {

    /// computing errors
    Eigen::Matrix<T, 6, 1> error = x - xd;
    eR = error.block(0, 0, 3, 1);
    eOmega = error.block(3, 0, 3, 1);

    dR = x.R * utils::hat(x.Omega);
    dRc = xd.R * utils::hat(xd.Omega);
    Eigen::Matrix<T, 3, 3> m1, m2;
    m1  = (dRc.transpose() * x.R - x.R.transpose() * dRc);
    m2 = (xd.R.transpose() * dR - dR.transpose() * xd.R);
    deR = 0.5 * (utils::vee(m1) + utils::vee(m2));

    /// setting up the CLF-QP
    V2 = (eOmega.transpose() * inertia_scaled_ * eOmega * 0.5
          + epsilon2 * (eR.transpose() * eOmega) + c2 * (eR.transpose() * eR) * 0.5).value();
    LgV2 = (eOmega.transpose() * inertia_scaled_ + epsilon2 * eR.transpose());
    LfV2 = ((epsilon2 * eOmega.transpose() + c2 * eR.transpose()) * deR
            -  LgV2 * ( dR.transpose() * xd.R * xd.Omega + x.R.transpose() * xd.R * xd.dOmega)).value();

    for (int i = 0; i < 3; i ++) {
        A_new[i] = LgV2(0, i);
    }
    ubA_new[0] = -LfV2 - eta2 * V2;

    print_qp2_setup();
    /// solving the QP
    qpOASES::int_t nwsr = 1000;
    qpOASES::real_t cpu_time = 1 / 500; // << modify this >>
    qpOASES::returnValue sol_info = solver.hotstart(H_new, g_new, A_new, lb_new, ub_new, lbA_new, ubA_new, nwsr);
    if (sol_info == qpOASES::SUCCESSFUL_RETURN ) {
        this->pause = false;
        std::cout << "Optimal solution found" << std::endl;
    } else {
        this->pause = true;
        std::cout << "Optimal solution NOT found" << std::endl;
    }
    solver.getPrimalSolution(xOpt);
    Eigen::Matrix<T, 3, 1> dOmega;
    dOmega << xOpt[0], xOpt[1], xOpt[2];
    printf( "\nxOpt = [ %e, %e, %e ];  objVal = %e\n\n", xOpt[0], xOpt[1], xOpt[2], solver.getObjVal() );

    //    while(1);
    /// computing the input
    u = this->inertia_ * dOmega + utils::hat(x.Omega) * this->inertia_ * x.Omega;

}


template class SO3Clf<float>;
template class SO3Clf<double>;

} // namespace nonlinear_control
