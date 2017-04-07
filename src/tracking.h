//
// Created by Oleg Leyzerov on 29/03/2017.
//

#ifndef KALMAN2D_TRACKING_H
#define KALMAN2D_TRACKING_H

#include "measurement_package.h"
#include "ukf.h"
#include "Eigen/Dense"

using Eigen::VectorXd;

class Tracking {
public:
    Tracking();
    virtual ~Tracking();
    void ProcessMeasurement(const MeasurementPackage &measurement_pack, const VectorXd &ground_truth);
    VectorXd getRMSE();
    //KalmanFilter kf_;
    UKF ukf_;
    bool is_initialized_;
private:

    long previous_timestamp_;

    float noise_ax;
    float noise_ay;
    VectorXd RMSE;


};




#endif //KALMAN2D_TRACKING_H
