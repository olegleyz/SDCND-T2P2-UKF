#include <iostream>
#include <vector>
#include <fstream>
//#include "ukf.h"
#include "measurement_package.h"
#include "tracking.h"

using namespace std;

int main(int argc, char* argv[]) {

    if (argc != 3) {
        cerr << "Input and/or output files are missing" << endl;
        return 0;
    }
    vector<MeasurementPackage> measurement_pack_list;
    vector<VectorXd> ground_truth_list;

    string in_file_name_ = argv[1];

    ifstream in_file(in_file_name_.c_str(), ifstream::in);
    if (!in_file.is_open()) {
        cerr << "Can't open input file: " << in_file_name_ << endl;
    }

    string out_file_name_ = argv[2];
    ofstream  out_file(out_file_name_.c_str(), ofstream::out);
    if (!out_file.is_open()) {
        cerr << "Can't open output file: " << out_file_name_ << endl;
    }
    string line;
    int i = 0;
    while (getline(in_file, line) && (i < 200000)) {
        MeasurementPackage meas_package;
        VectorXd ground_truth(4);
        istringstream iss(line);
        string sensor_type;
        iss >> sensor_type;
        long timestamp;
        if (sensor_type.compare("L") == 0) {
            meas_package.sensor_type_ = MeasurementPackage::LASER;
            meas_package.raw_measurements_ = VectorXd(2);
            float x;
            float y;
            iss >> x;
            iss >> y;
            meas_package.raw_measurements_ << x, y;
            iss >> timestamp;
            meas_package.timestamp_ = timestamp;
            measurement_pack_list.push_back(meas_package);
        } else if (sensor_type.compare("R") == 0) {
            meas_package.sensor_type_ = MeasurementPackage::RADAR;
            meas_package.raw_measurements_ = VectorXd(3);
            float rho, phi, rhodot;
            iss >> rho;
            iss >> phi;
            iss >> rhodot;
            meas_package.raw_measurements_ << rho, phi, rhodot;
            iss >> timestamp;
            meas_package.timestamp_ = timestamp;
            measurement_pack_list.push_back(meas_package);
        }
        float gtpx, gtpy, gtvx, gtvy;
        iss >> gtpx;
        iss >> gtpy;
        iss >> gtvx;
        iss >> gtvy;
        ground_truth << gtpx, gtpy, gtvx, gtvy;
        ground_truth_list.push_back(ground_truth);
        i++;
    }

    Tracking tracking;
    size_t N = measurement_pack_list.size();
    for (size_t k = 0; k < N; k++) {
        tracking.ProcessMeasurement(measurement_pack_list[k],ground_truth_list[k]);

        if (measurement_pack_list[k].sensor_type_==MeasurementPackage::RADAR) {
            out_file << measurement_pack_list[k].raw_measurements_[0] << "\t";
            out_file << measurement_pack_list[k].raw_measurements_[1] << "\t";
        } else {
            float ro = measurement_pack_list[k].raw_measurements_[0];
            float phi = measurement_pack_list[k].raw_measurements_[1];
            out_file << ro * cos(phi) << "\t";
            out_file << ro * sin(phi) << "\t";
        }
        out_file << ground_truth_list[k][0] << "\t";
        out_file << ground_truth_list[k][1] << "\t";
        out_file << ground_truth_list[k][2] << "\t";
        out_file << ground_truth_list[k][3] << "\n";
    }

    // calculating RMSE
    VectorXd RMSE = tracking.getRMSE();
    RMSE = (RMSE/N).array().sqrt();
    cout << "RMSE" << endl;
    cout << RMSE[0] << endl;
    cout << RMSE[1] << endl;
    cout << RMSE[2] << endl;
    cout << RMSE[3] << endl;

//    cout << tracking.ukf_.x_ << endl;
//    cout << tracking.ukf_.P_ << endl;

    if (in_file.is_open()) {
        in_file.close();
    }

    if (out_file.is_open()) {
        out_file.close();
    }


    /*
    Tracking tracking;
    tracking.is_initialized_ = true;
    tracking.ukf_.x_ <<
                     5.93637,
            1.49035,
            2.20528,
            0.536853,
            0.353577;
    tracking.ukf_.std_a_ = 0.2; // 3 m/s^2
    tracking.ukf_.std_yawdd_ = 0.2;
    tracking.ukf_.P_aug_(tracking.ukf_.n_x_, tracking.ukf_.n_x_) = tracking.ukf_.std_a_ * tracking.ukf_.std_a_;
    tracking.ukf_.P_aug_(tracking.ukf_.n_x_ + 1, tracking.ukf_.n_x_ + 1) = tracking.ukf_.std_yawdd_ * tracking.ukf_.std_yawdd_;


    MatrixXd Xsig_pred = MatrixXd(5,15);
    Xsig_pred <<
              5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
            1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
            2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
            0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
            0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;
    tracking.ukf_.Xsig_pred_ = Xsig_pred;
    //cout << tracking.ukf_.Xsig_pred_ << endl;

    tracking.ukf_.PredictMeanAndCovariance();



    VectorXd z = VectorXd(3);
    z <<
            5.9214,   //rho in m
            0.2187,   //phi in rad
            2.0062;
    tracking.ukf_.PredictMeasurement(z, 0); // predict Radar meas. Sigma, state and covar. matrix
    std::cout << tracking.ukf_.x_ << std::endl;
    std::cout << tracking.ukf_.P_ << std::endl;
*/
/*
    Tracking tracking;
    tracking.is_initialized_ = true;
    tracking.ukf_.x_ <<
             5.93637,
            1.49035,
            2.20528,
            0.536853,
            0.353577;
    tracking.ukf_.SigmaPointPrediction(0);
    cout << tracking.ukf_.Xsig_pred_ << endl;
*/



    return 0;


}