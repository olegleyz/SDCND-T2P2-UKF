//
// Created by Oleg Leyzerov on 04/04/2017.
//

#include <iostream>
#include "ukf.h"
#include "Eigen/Dense"

using Eigen::VectorXd;
using Eigen::MatrixXd;

UKF::UKF() {};
UKF::~UKF() {};

void UKF::GenerateSigmaPoints(MatrixXd* Xsig_out) {
    int n_x = 5; // state vector dimension
    double lambda = 3 - n_x; // lambda for sigma points generation

    VectorXd x = VectorXd(n_x); // state vector
    x << 5.7441,
        1.3800,
        2.2049,
        0.5015,
        0.3528;

    MatrixXd P = MatrixXd(n_x, n_x); // process covariance matrix

    P << 0.0043, -0.0013, 0.0030, -0.0022, -0.0020,
        -0.0013, 0.0077, 0.0011, 0.0071, 0.0060,
        0.0030, 0.0011, 0.0054, 0.0007, 0.0008,
        -0.0022, 0.0071, 0.0007, 0.0098, 0.0100,
        -0.0020, 0.0060, 0.0008, 0.0100, 0.0123;

    MatrixXd Xsig = MatrixXd (n_x, 2 * n_x + 1); // matrix of Sigma Points

    MatrixXd A = P.llt().matrixL();

    MatrixXd lA = pow ((lambda + n_x), 0.5) * A;

    Xsig.col(0) = x;
    for (int i = 1; i <= n_x; i++) {
        Xsig.col(i) = x + lA.col(i-1);
        Xsig.col(i + n_x) = x - lA.col(i-1);
    }
    *Xsig_out = Xsig;

}

void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out){
    int n_x = 5; // state dimension
    int n_aug = 7; // px,py,v, yaw, yaw dot, nu_a, nu_sai

    double std_a = 0.2; // process noise standard deviation longitudinal acceleration in m/s^2
    double std_yawdd = 0.2; // process noise standard deviation yaw acceleration in rad/s^2

    double lambda = 3 - n_aug;

    VectorXd x = VectorXd(n_x); // state vector
    x << 5.7441,
            1.3800,
            2.2049,
            0.5015,
            0.3528;

    VectorXd x_aug = VectorXd::Zero(n_aug); // augmented mean vector

    x_aug.head(5) = x;

    MatrixXd P = MatrixXd(n_x, n_x); // process covariance matrix

    P << 0.0043, -0.0013, 0.0030, -0.0022, -0.0020,
            -0.0013, 0.0077, 0.0011, 0.0071, 0.0060,
            0.0030, 0.0011, 0.0054, 0.0007, 0.0008,
            -0.0022, 0.0071, 0.0007, 0.0098, 0.0100,
            -0.0020, 0.0060, 0.0008, 0.0100, 0.0123;

    MatrixXd P_aug = MatrixXd::Zero(n_aug, n_aug); // augmented process covariance matrix
    P_aug.topLeftCorner(n_x, n_x) = P;
    P_aug(n_x, n_x) = std_a * std_a;
    P_aug(n_x + 1, n_x + 1) = std_yawdd * std_yawdd;

    MatrixXd Xsig_aug = MatrixXd (n_aug, 2 * n_aug + 1); // matrix of Sigma Points

    MatrixXd A = P_aug.llt().matrixL();

    MatrixXd lA = pow ((lambda + n_aug), 0.5) * A;

    Xsig_aug.col(0) = x_aug;
    for (int i = 1; i <= n_aug; i++) {
        Xsig_aug.col(i) = x_aug + lA.col(i-1);
        Xsig_aug.col(i + n_aug) = x_aug - lA.col(i-1);
    }
    *Xsig_out = Xsig_aug;
}
void UKF::SigmaPointPrediction(MatrixXd* Xsig_out) {
    int n_x = 5; // state dimenstion
    int n_aug = 7; //augmented dimension

    //create example sigma point matrix
    MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);
    Xsig_aug <<
             5.7441,  5.85768,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.63052,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,
            1.38,  1.34566,  1.52806,     1.38,     1.38,     1.38,     1.38,     1.38,   1.41434,  1.23194,     1.38,     1.38,     1.38,     1.38,     1.38,
            2.2049,  2.28414,  2.24557,  2.29582,   2.2049,   2.2049,   2.2049,   2.2049,   2.12566,  2.16423,  2.11398,   2.2049,   2.2049,   2.2049,   2.2049,
            0.5015,  0.44339, 0.631886, 0.516923, 0.595227,   0.5015,   0.5015,   0.5015,   0.55961, 0.371114, 0.486077, 0.407773,   0.5015,   0.5015,   0.5015,
            0.3528, 0.299973, 0.462123, 0.376339,  0.48417, 0.418721,   0.3528,   0.3528,  0.405627, 0.243477, 0.329261,  0.22143, 0.286879,   0.3528,   0.3528,
            0,        0,        0,        0,        0,        0,  0.34641,        0,         0,        0,        0,        0,        0, -0.34641,        0,
            0,        0,        0,        0,        0,        0,        0,  0.34641,         0,        0,        0,        0,        0,        0, -0.34641;

    MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
    double delta_t = 0.1; // time diff in sec

    for (int i = 0; i < (2 * n_aug + 1); i++) {
        VectorXd x = Xsig_aug.col(i).head(n_x);
        double pxdot, pydot, yawdot, cos_yaw, sin_yaw;
        cos_yaw = cos(x[3]);
        sin_yaw = sin(x[3]);
        if (x[4] < epsilon) {
            double xdt = x[2] * delta_t;
            pxdot = xdt * cos_yaw;
            pydot = xdt * sin_yaw;
            yawdot = 0;
        } else {
            double x2x4 = x[2] / x[4];
            double angle = x[3] + x[4] * delta_t;
            pxdot = x2x4 * (sin(angle) - sin_yaw);
            pydot = x2x4 * (-cos(angle) + cos_yaw);
            yawdot = x[4] * delta_t;
        }
        VectorXd xdot = VectorXd(n_x);
        xdot << pxdot, pydot, 0, yawdot, 0;

        VectorXd noise = VectorXd(n_x);
        double nu_a = Xsig_aug(5,i);
        double nu_yawdd = Xsig_aug(6,i);
        double dtna = delta_t * nu_a;
        double dtny = delta_t * nu_yawdd;
        double half_dt = 0.5 * delta_t;
        double pxn = half_dt * dtna * cos_yaw; // 1/2 * dt^2 * cos(yaw) * nu_a
        double pyn = half_dt * dtna * sin_yaw; // 1/2 * dt^2 * sin(yaw) * nu_a
        double pyawn = half_dt * dtny;

        noise << pxn, pyn, dtna, pyawn, dtny;

        Xsig_pred.col(i) = x + xdot + noise;
    }
    *Xsig_out = Xsig_pred;
}
void UKF::PredictMeanAndCovariance(VectorXd* x_out, MatrixXd* P_out) {
    int n_x = 5; // state dimension
    int n_aug = 7; // augmented dimension
    int len = 2 * n_aug + 1;

    double lambda = 3 - n_aug; // spreading parameter
    MatrixXd Xsig_pred = MatrixXd(n_x, n_aug * 2 + 1);
    Xsig_pred <<
            5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
            1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
            2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
            0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
            0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;

    VectorXd weights = VectorXd(len); // weights vector
    weights(0) = lambda / (lambda + n_aug);
    for (int i = len-1; i>0; i--) {
        weights(i) = 1 / (2 * (lambda + n_aug));
    }

    VectorXd x = VectorXd::Zero(n_x); // predicted mean vector
    for (int i = len; i>0; i--) {
        x += weights(i - 1) * Xsig_pred.col(i - 1);
    }

    MatrixXd P = MatrixXd::Zero(n_x, n_x); // predicted process covariance matrix
    for (int i = len; i>0; i--) {
        VectorXd dif = Xsig_pred.col(i - 1) - x;
        if (dif[3] > M_PI) dif[3] -= 2. * M_PI;
        if (dif[3] < -M_PI) dif[3] += 2. * M_PI;
        P += weights(i - 1) * dif * dif.transpose();
    }
    *x_out = x;
    *P_out = P;
}
void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out) {
    int n_x = 5;
    int n_aug = 7;
    int n_z = 3; // measurement dimension for radar (r, phi, r_dot)
    double lambda = 3 - n_aug;

    VectorXd weights = VectorXd(2 * n_aug + 1);
    weights(0) = lambda / (lambda + n_aug);
    double weight = 0.5 / (lambda + n_aug);
    for (int i = 2 * n_aug; i>0; i--) {
        weights(i) = weight;
    }

    // Radar measurement covariance matrix
    double std_radr = 0.3;
    double std_radphi = 0.0175;
    double std_radrd = 0.1;

    MatrixXd R = MatrixXd::Zero(3, 3);
    R(0, 0) = std_radr * std_radr;
    R(1, 1) = std_radphi * std_radphi;
    R(2, 2) = std_radrd * std_radrd;

    MatrixXd Xsig_pred = MatrixXd (n_x, 2 * n_aug + 1);
    Xsig_pred <<
            5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
            1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
            2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
            0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
            0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;

    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1); // matrix for sigma points in measurement space

    double radr, radphi, radrd;
    for (int i = 2 * n_aug + 1; i > 0; i--) {
        VectorXd x = Xsig_pred.col(i - 1);
        if (x[0] < epsilon && x[1] < epsilon) {
            x[0] = epsilon;
            x[1] = epsilon;
        } else if (x[0] < epsilon) x[0] = epsilon;
        radr = sqrt(x[0]*x[0] + x[1]*x[1]);
        radphi = atan2(x[1], x[0]);
        radrd = (x[0] * cos(x[3]) * x[2] + x[1] * sin(x[3]) * x[2]) / radr;
        Zsig.col(i - 1) << radr, radphi, radrd;
    }

    VectorXd z_pred = VectorXd::Zero(n_z); // mean predicted measurement
    for (int i = 2 * n_aug + 1; i>0; i--) {
        z_pred += weights(i - 1) * Zsig.col(i - 1);
    }

    MatrixXd S = MatrixXd::Zero(n_z, n_z); // measurement covariance matrix
    for (int i = 2 * n_aug + 1; i>0; i--) {
        VectorXd dif = Zsig.col(i - 1) - z_pred;
        if (dif[1] > M_PI) dif[1] -= 2. * M_PI;
        if (dif[1] < -M_PI) dif[1] += 2. * M_PI;
        S += weights(i - 1) * dif * dif.transpose();
    }
    S += R;
    *z_out = z_pred;
    *S_out = S;

}
