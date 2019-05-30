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
    out(0,0) = X(0,0) + std::cos(X(2,0))*u(0,0);
    out(1,0) = X(1,0) + std::sin(X(2,0))*u(0,0);
    out(2,0) = X(2,0) + 0.1*std::tan(u(1,0))*u(0,0);
    return out;
}

lin::Mat<double> G(lin::Mat<double> X){
    lin::Mat<double> out(4,1);
    out(0,0) = X(0,0);
    out(1,0) = X(1,0);
    out(2,0) = X(0,0);
    out(3,0) = X(1,0);
    return out;
}

int main(int argc, char *argv[]){

    // Define random generator with Gaussian distribution
    std::default_random_engine generator;
    std::normal_distribution<double> gps1(0, 1);
    std::normal_distribution<double> gps2(0, 0.64);
    std::normal_distribution<double> dist(0, 0.2);
    
    lin::Mat<double> X(3,1);
    lin::Mat<double> Xr(3,1);
    Xr =0;
    X = 0;

    lin::Mat<double> R(4,4);
    R = 0;
    R(0,0) = 1;
    R(1,1) = 1;
    R(2,2) = 0;
    R(3,3) = 0;

    lin::Mat<double> P(3,3);

    lin::Mat<double> Q(3,3);
    Q = 0.04;

    lin::Mat<double> Z(4,1);
    Z = 0;

    lin::Mat<double> control(2,1);
    control = 0;
    control(0,0) = 0.1; // dt

    kalman::EKFin ekf;

    ekf.setState(&X);
    ekf.setSensorNoiseCovariance(R);
    ekf.setProcessNoiseCovariance(Q);
    ekf.setSensor(&Z);
    ekf.setControl(&control);
    ekf.setF(F);
    ekf.setG(G);
    ekf.init();

    std::vector<double> gps1x;
    std::vector<double> gps1y;
    std::vector<double> gps2x;
    std::vector<double> gps2y;
    std::vector<double> realx;
    std::vector<double> realy;
    std::vector<double> kx;
    std::vector<double> ky;


    for(int i=0;i<1000;i++){
        
        double s = dist(generator);
        control(1,0) = s;

        Z(0,0) = Xr(0,0) + gps1(generator);
        Z(1,0) = Xr(1,0) + gps1(generator);
        Z(2,0) = Xr(0,0) + gps2(generator);
        Z(3,0) = Xr(1,0) + gps2(generator);
        gps1x.push_back(Z(0,0));
        gps1y.push_back(Z(1,0));
        gps2x.push_back(Z(2,0));
        gps2y.push_back(Z(3,0));

        ekf.predict();
        ekf.update();
        kx.push_back(X(0,0));
        ky.push_back(X(1,0));

        Xr = F(Xr, control);
        realx.push_back(Xr(0,0));
        realy.push_back(Xr(1,0));
        
    }

    plt::named_plot("X gps1", gps1x);
    plt::named_plot("X gps2", gps2x);
    plt::named_plot("X real", realx);
    plt::named_plot("X kalman", kx);
    plt::legend();
    plt::show();

}

