#include <iostream>
#include "measurement_package.h"
#include <vector>
#include <fstream>
#include "tracking.h"

using namespace std;
using Eigen::VectorXd;
using std::vector;


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
    while (getline(in_file, line) && (i >= -1)) {
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

        out_file << tracking.kf_.x_[0] << "\t";
        out_file << tracking.kf_.x_[1] << "\t";
        out_file << tracking.kf_.x_[2] << "\t";
        out_file << tracking.kf_.x_[3] << "\t";
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


    if (in_file.is_open()) {
        in_file.close();
    }

    if (out_file.is_open()) {
        out_file.close();
    }

    return 0;
}