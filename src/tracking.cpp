//
// Created by Oleg Leyzerov on 29/03/2017.
//

#include <iostream>
#include "tracking.h"
#include "Eigen/Dense"


using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;

Tracking::Tracking(){
    use_radar_ = true; // true if we should use measurement from RADAR
    use_laser_ = true; // true if we should use measurement from LIDAR
    is_initialized_ = false;
    previous_timestamp_ = 0;
    ukf_.epsilon_ = 0.001;
    ukf_.n_x_ = 5;
    ukf_.n_aug_ = 7;
    ukf_.lambda_ = 3 - ukf_.n_aug_;
    ukf_.n_sigma_ = 2 * ukf_.n_aug_ + 1;
    double lna = ukf_.lambda_ + ukf_.n_aug_;
    double weight = 0.5 / lna;
    ukf_.weights_ = VectorXd(2 * ukf_.n_aug_ + 1);
    ukf_.weights_(0) = ukf_.lambda_ / lna;
    for (int i = 2 * ukf_.n_aug_; i>0; i--) {
        ukf_.weights_(i) = weight;
    }
    ukf_.std_a_ = 3; // 3 m/s^2
    ukf_.std_yawdd_ = 0.4; // 0.9 rad
    ukf_.std_lpx_ = 0.15;
    ukf_.std_lpy_ = 0.15;
    ukf_.std_rro_ = 0.3;
    ukf_.std_rphi_ = 0.03;
    ukf_.std_rrod_ = 0.3;
    ukf_.z_lidar_ = 2;
    ukf_.R_lidar_ = MatrixXd::Zero(2,2);
    ukf_.R_lidar_(0, 0) = ukf_.std_lpx_ * ukf_.std_lpx_;
    ukf_.R_lidar_(1, 1) = ukf_.std_lpy_ * ukf_.std_lpy_;
    ukf_.z_radar_ = 3;
    ukf_.R_radar_ = MatrixXd::Zero(3,3);
    ukf_.R_radar_(0, 0) = ukf_.std_rro_ * ukf_.std_rro_;
    ukf_.R_radar_(1, 1) = ukf_.std_rphi_ * ukf_.std_rphi_;
    ukf_.R_radar_(2, 2) = ukf_.std_rrod_ * ukf_.std_rrod_;

    ukf_.x_ = VectorXd(ukf_.n_x_); // 5D state vector
    ukf_.x_aug_ = VectorXd::Zero(ukf_.n_aug_); // 7D state vector
    ukf_.P_ = MatrixXd(ukf_.n_x_, ukf_.n_x_); // state covariance matrix
    ukf_.P_ <<
            1, 0, 0, 0, 0,
            0, 1, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;

    ukf_.P_aug_ = MatrixXd::Zero(ukf_.n_aug_, ukf_.n_aug_);
    ukf_.P_aug_.topLeftCorner(ukf_.n_x_, ukf_.n_x_) = ukf_.P_;
    ukf_.P_aug_(ukf_.n_x_, ukf_.n_x_) = ukf_.std_a_ * ukf_.std_a_;
    ukf_.P_aug_(ukf_.n_x_ + 1, ukf_.n_x_ + 1) = ukf_.std_yawdd_ * ukf_.std_yawdd_;
    RMSE = VectorXd(4); // root mean squared error
    RMSE << 0, 0, 0, 0;

}
Tracking::~Tracking(){}

void Tracking::ProcessMeasurement(const MeasurementPackage &measurement_pack, const VectorXd &ground_truth) {
    if (use_laser_ == false && measurement_pack.sensor_type_ == MeasurementPackage::LASER) return;
    if (use_radar_ == false && measurement_pack.sensor_type_ == MeasurementPackage::RADAR) return;
    if (!is_initialized_) {

        if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            ukf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0, 0;
        } else if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            // converting rho, phi, rho dot into position and velocity
            double px = measurement_pack.raw_measurements_[0] * cos(measurement_pack.raw_measurements_[1]);
            double py = measurement_pack.raw_measurements_[0] * sin(measurement_pack.raw_measurements_[1]);
            double vx = measurement_pack.raw_measurements_[2] * cos(measurement_pack.raw_measurements_[1]);
            double vy = measurement_pack.raw_measurements_[2] * sin(measurement_pack.raw_measurements_[1]);

            double v = sqrt(vx*vx + vy*vy);
            double yaw = 0;//measurement_pack.raw_measurements_[1];

            if (px < ukf_.epsilon_ && py < ukf_.epsilon_) {
                px = ukf_.epsilon_;
                py = ukf_.epsilon_;
            }
            ukf_.x_ << px, py, v, yaw, 0;
        }

        previous_timestamp_ = measurement_pack.timestamp_;
        is_initialized_ = true;
        return;
    }
    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;

    previous_timestamp_ = measurement_pack.timestamp_;

    while (dt>0.1) {ukf_.Predict(0.01); dt -= 0.01;}
    ukf_.Predict(dt);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        ukf_.UpdateRadar(measurement_pack.raw_measurements_);
    } else {
        ukf_.UpdateLidar(measurement_pack.raw_measurements_);
    }

    VectorXd residual;
    VectorXd ukf_x_cartesian = VectorXd(4);
    double x_estimate_ = ukf_.x_(0);
    double y_estimate_ = ukf_.x_(1);
    double vx_estimate_ = ukf_.x_(2) * cos(ukf_.x_(3));
    double vy_estimate_ = ukf_.x_(2) * sin(ukf_.x_(3));
    ukf_x_cartesian << x_estimate_, y_estimate_, vx_estimate_, vy_estimate_;

    residual = ukf_x_cartesian - ground_truth;
    residual = residual.array() * residual.array();
    RMSE += residual;

    //std::cout << "x_= " << ukf_.x_ << std::endl;
    //std::cout << "P_= " << ukf_.P_ << std::endl;

}

VectorXd Tracking::getRMSE() {
    return RMSE;
}


