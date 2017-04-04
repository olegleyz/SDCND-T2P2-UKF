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
    void GenerateSigmaPoints(MatrixXd* Xsig_out);
    void AugmentedSigmaPoints(MatrixXd* Xsig_out);
    void SigmaPointPrediction(MatrixXd* Xsig_out);
    void PredictMeanAndCovariance(VectorXd* x_out, MatrixXd* P_out);
    void PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out);
    int n_x;
    double lambda;
    float epsilon = 0.001;

private:
};

#endif //SRC_UKF_H
