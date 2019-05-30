#include "lin/mat.hpp"
#include "lin/kalman.hpp"
#include <iostream>
#define WITHOUT_NUMPY
#include "matplotlibcpp.h"
#include <iterator>
#include <random>
#include <vector>
#include <cmath>

#define print std::cout<<"here"<<std::endl;

namespace plt = matplotlibcpp;

lin::Mat<double> F(lin::Mat<double> X, lin::Mat<double> u){
    lin::Mat<double> out(3,1);
    out(0,0) = X(0,0) + X(1,0)*u(0,0);
    out(1,0) = X(1,0);
    out(2,0) = X(2,0);
    return out;
}
lin::Mat<double> JG(lin::Mat<double> X){
    lin::Mat<double> out(1,3);
    double denon = std::sqrt(std::pow(X(0,0),2) + std::pow(X(2,0),2));
    out(0,0) = X(0,0)/denon;
    out(0,1) = 0;
    out(0,2) = X(2,0)/denon;
    return out;
}

lin::Mat<double> G(lin::Mat<double> X){
    lin::Mat<double> out(1,1);
    out = std::sqrt(std::pow(X(0,0),2) + std::pow(X(2,0),2));
    return out;
}

int main(int argc, char *argv[]){

    // Define random generator with Gaussian distribution
    std::default_random_engine generator;
    std::normal_distribution<double> radar(0, 5);
    std::normal_distribution<double> dist(0, 0.2);
    
    lin::Mat<double> X(3,1);
    X = {0, 200, 1100};

    lin::Mat<double> R(1,1);
    R = 25;


    lin::Mat<double> Q(3,3);
    Q = 0;
    Q(2,2) = 0.1;
    Q(0,0) = 1.5625e-07;
    Q(0,1) = 6.2500e-06;
    Q(1,0) = 6.2500e-06;
    Q(1,1) = 2.5000e-04;

    lin::Mat<double> Z(1,1);
    Z = 0;

    lin::Mat<double> control(1,1);
    control = 0.05;

    kalman::EKFin ekf;

    ekf.setState(&X);
    ekf.setSensorNoiseCovariance(R);
    ekf.setProcessNoiseCovariance(Q);
    ekf.setSensor(&Z);
    ekf.setControl(&control);
    ekf.setF(F);
    ekf.setG(G);
    ekf.setGJacobian(JG);
    ekf.init();

    std::vector<double> positionr;
    std::vector<double> velocityr;
    std::vector<double> altituder;
    std::vector<double> positionk;
    std::vector<double> velocityk;
    std::vector<double> altitudek;

    double pos = 0;
    double vel = 100;
    double alt = 1000;

    for(int i=0;i<200;i++){
        
        vel+= 0.01*dist(generator);
        alt+= 0.01*dist(generator);
        pos+= vel*0.1;
        double err = pos + 0.05*dist(generator);
        Z(0,0) = std::sqrt(std::pow(pos,2) + std::pow(alt,2)) + err;
        positionr.push_back(pos);
        velocityr.push_back(vel);
        altituder.push_back(alt);

        ekf.predict();
        ekf.update();
        positionk.push_back(X(0,0));
        velocityk.push_back(X(1,0));
        altitudek.push_back(X(2,0));

        
    }

    plt::named_plot("pos real", positionr);
    plt::named_plot("pos kalman", positionk);
    plt::named_plot("vel real", velocityr);
    plt::named_plot("vel kalman", velocityk);
    plt::named_plot("altitude real", altituder);
    plt::named_plot("altitude kalman", altitudek);
    plt::legend();
    plt::show();

}

