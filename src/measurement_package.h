//
// Created by Oleg Leyzerov on 29/03/2017.
//

#ifndef KALMAN2D_MEASUREMENT_PACKAGE_H
#define KALMAN2D_MEASUREMENT_PACKAGE_H

#include "Eigen/Dense"

class MeasurementPackage {
public:
    long timestamp_;

    enum SensorType {
        LASER, RADAR
    } sensor_type_;

    Eigen::VectorXd raw_measurements_;
};

#endif //KALMAN2D_MEASUREMENT_PACKAGE_H
