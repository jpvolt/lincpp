#include "lin/mat.hpp"
#include "lin/kalman.hpp"
#include <iostream>
#define WITHOUT_NUMPY
#include "matplotlibcpp.h"
#include <iterator>
#include <random>
#include <vector>


namespace plt = matplotlibcpp;

int main(int argc, char *argv[]){
    double dt = 0.1;
    lin::Mat<double> Xk(1,1);
    Xk = 3;
    lin::Mat<double> A(1,1);
    A = 1;

    lin::Mat<double> B(1,1);
    B = 0;

    lin::Mat<double> Hk(1,1);
    Hk = 1;

    lin::Mat<double> Rk(1,1);
    Rk = 0.1;

    lin::Mat<double> Pk(1,1);
    Pk = 1;

    lin::Mat<double> control(1,1); 
    control = 0;
    lin::Mat<double> Qk(1,1);
    Qk = 0.01;


    // Define random generator with Gaussian distribution
    const double mean = 0.0;
    const double stddev = 0.2;
    std::default_random_engine generator;
    std::normal_distribution<double> dist(mean, stddev);

    kalman::KF kf;

    kf.setState(&Xk);
    kf.setCovariance(&Pk);
    kf.setSensorError(&Rk);
    kf.setStateTransition(&A);
    kf.setSensorTransition(&Hk);
    kf.setControlTransition(&B);
    

    lin::Mat<double> sensor(1,1);
    sensor = 0;

    kf.init();
    std::vector<double> vr;
    std::vector<double> vs;
    std::vector<double> vk;
    std::vector<double> gpsk;

    double gps = 0;
    double v = 4;
    for(int i=0;i<90;i++){

        vr.push_back(v);
        sensor = v + dist(generator);
        vs.push_back(sensor(0,0));

        kf.predict(control, Qk);
        kf.update(sensor);

        vk.push_back(Xk(0,0));

    }

    plt::plot(vr);
    plt::plot(vs);
    plt::plot(vk);
    plt::show();

}

