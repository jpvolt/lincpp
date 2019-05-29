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

lin::Mat<double> F(lin::Mat<double> x, lin::Mat<double> u){
    lin::Mat<double> out(3,1);
    out(0,0) = x(1,0)*std::cos(x(2,0))*dt;
    out(1,0) = x(1,0);
    out(2,0) = x(2,0);
    return out;
}

lin::Mat<double> G(lin::Mat<double> x){
    lin::Mat<double> out(4,1);
    out = 0;
    out(0,0) = x(1,0);
    out(1,0) = x(1,0);
    return out;
}

int main(int argc, char *argv[]){

    // Define random generator with Gaussian distribution
    const double mean = 0.0;
    const double stddev = 0.3;
    std::default_random_engine generator;
    std::normal_distribution<double> dist(mean, stddev);
    
    lin::Mat<double> X(3,1);
    X = 100;

    lin::Mat<double> R(4,4);
    R = 0;
    R(0,0) = 0.64;
    R(1,1) = 0.32;
    R(2,2) = 0.20;
    R(3,3) = 0.10;

    lin::Mat<double> P(3,3);
    P = P.I();
    P = P*0.5;

    lin::Mat<double> Q(3,3);
    Q = Q.I()*0.1;

    lin::Mat<double> Z(4,1);
    Z = 0;

    lin::Mat<double> control(4,1);
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
    print


    std::vector<double> vr;
    std::vector<double> vs;
    std::vector<double> vt;
    std::vector<double> vk;
    std::vector<double> xr;
    std::vector<double> xk;


    double x = 0;
    double v = 0;
    double a = 0;
    double  LO = -2;
    double HI = 4;
    print
    for(int i=0;i<90;i++){


        double j =  LO + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(HI-LO)));
        double s =  LO + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(2)));

        std::cout<<"x:"<<x<<", v:"<<v<<", a:"<<a<<", s:"<<s<<std::endl;
        Z(0,0) = v + 6*dist(generator);
        Z(1,0) = v + 3*dist(generator);
        Z(2,0) = a + 2*dist(generator);
        Z(3,0) = s + 1*dist(generator);
        vs.push_back(Z(0,0));
        vt.push_back(Z(1,0));

        ekf.predict(Q);
        ekf.update();
        xk.push_back(X(0,0));
        vk.push_back(X(1,0));

        std::cout<<"xk:"<<X(0,0)<<", vk:"<<X(1,0)<<std::endl;
        // env update
        x = x + v*std::cos(s);
        v = v + a*dt;
        a = a + j*dt;

    }

    plt::named_plot("Real speed",vr);
    plt::named_plot("Sensor 1",vs);
    plt::named_plot("sensor 2",vt);
    plt::named_plot("kalman speed",vk);
    //plt::named_plot("Real positon",gpsr);
    //plt::named_plot("gps sensor",gpss);
    //plt::named_plot("Kalman position",gpsk);
    plt::legend();
    plt::show();

}

