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
double dt = 0.01;

lin::Mat<double> F(lin::Mat<double> X, lin::Mat<double> u){
    return X;
}

lin::Mat<double> G(lin::Mat<double> X){
    lin::Mat<double> out(3,1);
    out(0,0) = X(0,0);
    out(1,0) = X(0,0);
    out(2,0) = X(0,0);
    return out;
}

int main(int argc, char *argv[]){

    // Define random generator with Gaussian distribution
    std::default_random_engine generator;
    std::normal_distribution<double> dist1(0, 1);
    std::normal_distribution<double> dist2(0, 0.64);
    std::normal_distribution<double> dist3(0, 0.32);
    std::normal_distribution<double> dist(0, 0.5);
    
    lin::Mat<double> X(1,1);
    lin::Mat<double> Xr(1,1);
    Xr =0;
    X = 0;

    lin::Mat<double> R(3,3);
    R = 0;
    R(0,0) = 1;
    R(1,1) = 0.4;
    R(2,2) = 0.10;

    lin::Mat<double> P(1,1);

    lin::Mat<double> Q(1,1);
    Q = 0.25;

    lin::Mat<double> Z(3,1);
    Z = 0;

    lin::Mat<double> control(1,1);
    control = 0;

    kalman::EKFin ekf;

    ekf.setState(&X);
    ekf.setSensorError(&R);
    ekf.setSensor(&Z);
    ekf.setCovariance(&P);
    ekf.setControl(&control);
    ekf.setF(F);
    ekf.setG(G);
    ekf.computeJacobians();

    std::vector<double> temp1;
    std::vector<double> temp2;
    std::vector<double> temp3;
    std::vector<double> tempreal;
    std::vector<double> tempkalman;


    for(int i=0;i<1000;i++){



        Z(0,0) = Xr(0,0) + dist1(generator);
        Z(1,0) = Xr(0,0) + dist2(generator);
        Z(2,0) = Xr(0,0) + dist3(generator);
        tempreal(Xr(0,0));
        temp1.push_back(Z(0,0));
        temp2.push_back(Z(1,0));
        temp3.push_back(Z(2,0));

        ekf.predict(Q);
        ekf.update();
        tempkalman.push_back(X(0,0));

        Xr = Xr(0,0) + dist(generator);
        
    }

    plt::named_plot("temp1", temp1);
    plt::named_plot("temp2", temp2);
    plt::named_plot("temp3", temp3);
    plt::named_plot("real temp", tempreal);
    plt::named_plot("kalman temp", tempkalman);
    plt::legend();
    plt::show();

}

