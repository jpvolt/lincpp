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
    lin::Mat<double> out(2,1);
    out(0,0) = x(0,0) + cos(x(1,0))*u(0,0);
    out(1,0) = x(1,0);
    return out;
}

lin::Mat<double> G(lin::Mat<double> x){
    lin::Mat<double> out(2,1);
    out(0,0) = x(0,0);
    out(1,0) = x(1,0);
    return out;
}

int main(int argc, char *argv[]){

    // Define random generator with Gaussian distribution
    const double mean = 0.0;
    const double stddev = 0.3;
    std::default_random_engine generator;
    std::normal_distribution<double> dist(mean, stddev);
    
    lin::Mat<double> X(2,1);
    X = 0;

    lin::Mat<double> R(2,2);
    R = 0;
    R(0,0) = 0.64;
    R(1,1) = 0.32;

    lin::Mat<double> Q(2,2);
    Q = Q.I()*0.001;

    lin::Mat<double> Z(2,1);
    Z = 0;

    lin::Mat<double> control(1,1);
    control = dt;

    kalman::EKFin ekf;

    ekf.setState(&X);
    ekf.setSensorNoiseCovariance(R);
    ekf.setSensor(&Z);
    ekf.setProcessNoiseCovariance(Q);
    ekf.setControl(&control);
    ekf.setF(F);
    ekf.setG(G);
    ekf.init(50);


    std::vector<double> xr;
    std::vector<double> tr;
    std::vector<double> xk;
    std::vector<double> tk;
    std::vector<double> xs;
    std::vector<double> ts;


    double x = 0;
    double t = 0;
    double  LO = -2;
    double HI = 4;
    lin::Mat<double> s(2,1);
    s = 0;
        for(int i=0;i<100;i++){



        s = {X(0,0), X(1,0) + X(1,0)*0.05*dist(generator)};
        s = F(s, control);

        xr.push_back(s(0,0));
        tr.push_back(s(1,0));

        Z(0,0) = s(0,0) + 3*dist(generator);
        Z(1,0) = s(1,0) + 1*dist(generator);
        xs.push_back(Z(0,0));
        ts.push_back(Z(1,0));

        ekf.predict_update();

        xk.push_back(X(0,0));
        tk.push_back(X(1,0));


    }

    plt::named_plot("Real x",xr);
    plt::named_plot("Sensor x",xs);
    plt::named_plot("kalman x",xk);
    plt::legend();
    plt::show();

}

