//
// Created by Oleg Leyzerov on 04/04/2017.
//

#include <iostream>
#include "ukf.h"
#include "Eigen/Dense"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using namespace std;

UKF::UKF() {

};
UKF::~UKF() {};

MatrixXd UKF::AugmentedSigmaPoints() {

    MatrixXd Xsig_aug = MatrixXd(n_aug_,n_sigma_); // matrix of Sigma Points

    P_aug_.fill(0.0);
    P_aug_.topLeftCorner(n_x_,n_x_) = P_;
    P_aug_(n_x_, n_x_) = std_a_ * std_a_;
    P_aug_(n_x_+1, n_x_+1) = std_yawdd_ * std_yawdd_;

    MatrixXd A = P_aug_.llt().matrixL(); // sqrt(P_aug)
    MatrixXd lA = sqrt(lambda_ + n_aug_) * A; // sqrt(lambda + n_aug) * sqrt(P_aug)
    x_aug_.head(n_x_) = x_; // generating x_aug (px,py,v,yaw,yawd,mean_noise_a,mean_noise_yawd)

    Xsig_aug.col(0) = x_aug_;
    for (int i = n_aug_; i >0 ; i--) {
        Xsig_aug.col(i) = x_aug_ + lA.col(i - 1);
        Xsig_aug.col(i + n_aug_) = x_aug_ - lA.col(i - 1);
    }
    return Xsig_aug;

}
void UKF::SigmaPointPrediction(double delta_t) {
    MatrixXd Xsig_aug = UKF::AugmentedSigmaPoints(); // Generating augmented Sigma Points, dim (7, 15)
    MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1); // Predicted Sigma Points
    for (int i=0; i< n_sigma_; i++) {
        VectorXd x = Xsig_aug.col(i);
        double pxdot, pydot, yawdot, cos_yaw, sin_yaw;
        cos_yaw = cos(x[3]); // cos(yaw)
        sin_yaw = sin(x[3]); // sin(yaw)
        if (fabs(x[4]) < epsilon_) { // check if yawd ~0
            double vdt = x[2] * delta_t; // v * delta_t
            pxdot = vdt * cos_yaw;
            pydot = vdt * sin_yaw;
            yawdot = 0;
        } else {
            double vyawd = x[2] / x[4];
            double angle = x[3] + x[4] * delta_t;
            pxdot = vyawd * (sin(angle) - sin_yaw);
            pydot = vyawd * (-cos(angle) + cos_yaw);
            yawdot = x[4] * delta_t;
        }
        VectorXd xdot = VectorXd(n_x_); // predicted state vector without noise yet
        xdot << pxdot, pydot, 0, yawdot, 0;

        VectorXd noise = VectorXd(n_x_); // process noise vector
        double nu_a = x(5);
        double nu_yawdd = x(6);
        double dtna = delta_t * nu_a;
        double dtny = delta_t * nu_yawdd;
        double half_dt = 0.5 * delta_t;
        double pxn = half_dt * dtna * cos_yaw; // 1/2 * dt^2 * cos(yaw) * nu_a
        double pyn = half_dt * dtna * sin_yaw; // 1/2 * dt^2 * sin(yaw) * nu_a
        double pyawn = half_dt * dtny;

        noise << pxn, pyn, dtna, pyawn, dtny;

        Xsig_pred.col(i) = x.head(n_x_) + xdot + noise; // Matrix of predicted Sigma Points
    }
    Xsig_pred_ = Xsig_pred;
}
void UKF::PredictMeanAndCovariance() {

    VectorXd x = VectorXd::Zero(n_x_); // predicted mean vector
    for (int i = 0; i < n_sigma_; i++) {
        x += weights_(i) * Xsig_pred_.col(i);
    }
    MatrixXd P = MatrixXd::Zero(n_x_, n_x_); // predicted process covariance matrix
    for (int i = 0; i < n_sigma_; i++) {
        VectorXd dif = Xsig_pred_.col(i) - x;
        dif(3) = UKF::NormalizeAngle(dif(3));

        P = P + weights_(i) * dif * dif.transpose();
    }
    x_ = x; // update state vector
    P_ = P; // update covariance matrix
}
void UKF::Predict(double delta_t) {
    UKF::SigmaPointPrediction(delta_t); // Generating and Predicting Sigma Points, dim (5, 15)
    UKF::PredictMeanAndCovariance(); // Predict mean and covariance of predicted Sigma Points
}

void UKF::PredictMeasurement(const VectorXd &z, int sensor_type) {


    if (sensor_type == 0) {
        double radr, radphi, radrd;
        MatrixXd Zsig;
        VectorXd z_pred;
        MatrixXd S;

        Zsig = MatrixXd(z_radar_, n_sigma_); // matrix for sigma points in measurement space
        z_pred = VectorXd::Zero(z_radar_); // mean predicted measurement
        S = MatrixXd::Zero(z_radar_, z_radar_); // measurement covariance matrix


        for (int i = 0; i < n_sigma_; i++) {
            VectorXd x = Xsig_pred_.col(i);
            // check for nonsense data

            if (sensor_type == 0) {
                if (fabs(x[0]) < epsilon_ && fabs(x[1]) < epsilon_) { //if px and py ~ 0
                    x[0] = epsilon_;
                    x[1] = epsilon_;
                } else if (fabs(x[0]) < epsilon_) x[0] = epsilon_; // if px ~ 0
                radr = sqrt(x[0] * x[0] + x[1] * x[1]);
                radphi = atan2(x[1], x[0]);
                radrd = (x[0] * cos(x[3]) * x[2] + x[1] * sin(x[3]) * x[2]) / radr;
                Zsig.col(i) << radr, radphi, radrd;
            }
        }

        for (int i = 0; i < n_sigma_; i++) {
            z_pred += weights_(i) * Zsig.col(i);
        }

        for (int i = 0; i < n_sigma_; i++) {
            VectorXd dif = Zsig.col(i) - z_pred;
            dif(1) = UKF::NormalizeAngle(dif(1));
            S += weights_(i) * dif * dif.transpose();
        }
        S += R_radar_;

        UKF::UpdateState(Zsig, z_pred, S, z, sensor_type);

    } else if (sensor_type==1) {
        MatrixXd H_ = MatrixXd(2, 5);
        H_ << 1, 0, 0, 0, 0,
                0, 1, 0, 0, 0;

        VectorXd z_pred = H_ * x_;
        VectorXd y = z - z_pred;
        MatrixXd Ht = H_.transpose();
        MatrixXd S = H_ * P_ * Ht + R_lidar_;
        MatrixXd K = P_ * Ht * S.inverse();
        x_ = x_ + K * y;
        long x_size = x_.size();
        MatrixXd I_ = MatrixXd::Identity(x_size, x_size);
        P_ = (I_ - K * H_) * P_;
    }

}
double UKF::NormalizeAngle(double angle) {
    if (angle > M_PI) {
        double temp = fmod((angle - M_PI), (2 * M_PI)); // -= 2. * M_PI;
        angle = temp - M_PI;
    } // phi normalization
    if (angle < -M_PI) {
        double temp = fmod((angle + M_PI) ,(2 * M_PI));
        angle = temp + M_PI;
    }
    return angle;
}
void UKF::UpdateState(MatrixXd &Zsig, VectorXd &z_pred, MatrixXd &S, const VectorXd &z, int sensor_type) {
    // Cross correlation between sigma points in state space and measurement space
    MatrixXd Tc;
    Tc = MatrixXd::Zero(n_x_, sensor_type == 0 ? z_radar_ : z_lidar_);

    for (int i = 0; i < n_sigma_; i++) {
        VectorXd dif_x = Xsig_pred_.col(i) - x_;
        dif_x(3) = UKF::NormalizeAngle(dif_x(3));
        VectorXd dif_z = Zsig.col(i) - z_pred;
        dif_z(1) = UKF::NormalizeAngle(dif_z(1));
        Tc += weights_(i) * dif_x * dif_z.transpose();
    }

    MatrixXd K = Tc * S.inverse(); // Kalman Gain
    VectorXd dif_z = z - z_pred; // z - incoming radar measurement
    if (sensor_type == 0) {
        dif_z(1) = UKF::NormalizeAngle(dif_z(1));
    }
    x_ += K * dif_z; // state update
    P_ -= K * S * K.transpose(); // process covariance update
 }

void UKF::UpdateRadar(const VectorXd &z) {
    UKF::PredictMeasurement(z, 0); // predict Radar meas. Sigma, state and covar. matrix
}

void UKF::UpdateLidar(const VectorXd &z) {
    UKF::PredictMeasurement(z, 1); // predict Radar meas. Sigma, state and covar. matrix
}