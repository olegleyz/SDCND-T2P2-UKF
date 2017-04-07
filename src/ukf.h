//
// Created by Oleg Leyzerov on 04/04/2017.
//
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;


#ifndef SRC_UKF_H
#define SRC_UKF_H

class UKF {
public:
    UKF(); // constructor
    virtual ~UKF(); // destructor

    MatrixXd AugmentedSigmaPoints();
    void SigmaPointPrediction(double delta_t);
    void PredictMeanAndCovariance();
    void Predict(double delta_t);
    void PredictMeasurement(const VectorXd &z, int sensor_type);
    void UpdateState(MatrixXd &Zsig, VectorXd &z_pred, MatrixXd &S, const VectorXd &z, int sensor_type);
    void UpdateRadar(const VectorXd &z);
    void UpdateLidar(const VectorXd &z);
    double NormalizeAngle(double angle);

    double lambda_;
    float epsilon_ = 0.001;
    int n_x_; // state dimension
    int n_aug_; // augmented dimension
    int n_sigma_; // amount of augmented sigma points (2*n_aug+1)
    VectorXd weights_; // weights vector
    VectorXd x_; // state vector
    VectorXd x_aug_; // augmented state vector
    MatrixXd P_; // process covariance matrix
    MatrixXd P_aug_; // process covariance augmented matrix
    double std_a_; // process noise std. dev. longitudinal acceleration
    double std_yawdd_; // process noise std. dev. yaw acceleration
    double std_lpx_; // measurement noise LIDAR position px
    double std_lpy_; // measurement noise LIDAR position py
    double std_rro_; // measurement noise RADAR ro
    double std_rphi_; // measurement noise RADAR phi
    double std_rrod_; // measurement noise RADAR ro dot
    int z_radar_; // measurement dimension for radar (r, phi, r_dot)
    MatrixXd R_radar_; // Radar measurement covariance matrix
    int z_lidar_; // measurement dimension for lidar (px, py)
    MatrixXd R_lidar_; // LIDAR measurement covariance matrix

    MatrixXd Xsig_pred_; //
    void Prediction(double delta_t);

private:
};

#endif //SRC_UKF_H
