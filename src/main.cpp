#include <iostream>
#include "ukf.h"

int main() {
    UKF ukf;
    VectorXd z = VectorXd(5);
    MatrixXd S = MatrixXd(3,3);
    ukf.PredictRadarMeasurement(&z, &S);
    std::cout << "z  = " << std::endl << z << std::endl;
    std::cout << "S  = " << std::endl << S << std::endl;
    return 0;
}