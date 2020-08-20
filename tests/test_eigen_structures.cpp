#include <iostream>
#include <eigen3/Eigen/Dense>
#include <deque>
#include <sys/time.h>

int main() {


    Eigen::Matrix<double, 15, 15> bigM = Eigen::Matrix<double, 15, 15>::Random();
    Eigen::Matrix<double, 15, 15> *bigMPtr;
    bigMPtr = &bigM;

    std::cout << "bigM: \n" << bigM << std::endl;
    auto A = bigM.block(0,0,3,3);
    A.setIdentity();
    std::cout << "bigM: \n" << bigM << std::endl;

    for (int it = 0; it < 225; ++it) {
        printf("bigM[%d] = %f = %f\n", it, bigM(it), *(bigMPtr+it));
    }

    // std::deque<Eigen::Matrix<double, 3, 3>> matrixStorage;
    // matrixStorage.push_back(bigM.block(0, 0, 3, 3));
    // matrixStorage.at(0).setZero();

    // std::cout << "bigM: \n" << bigM << std::endl;


    // struct timeval stop, start;
    // printf("Beginning inverse computation...");
    // gettimeofday(&start, NULL);
    // bigM.inverse();
    // gettimeofday(&stop, NULL);
    // printf("took %lu us\n", (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec);
    return 0;
}


