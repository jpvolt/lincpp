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
    lin::Mat<double> out(5,1);
    out(0,0) = x(0,0) + x(2,0)*std::cos(x(3,0))*u(0,0);
    out(1,0) = x(1,0) + x(2,0)*std::sin(x(3,0))*u(0,0);
    out(2,0) = x(2,0) + u(1,0)*u(0,0);
    out(3,0) = x(3,0);
    out(4,0) = x(4,0) + x(2,0)*std::tan(u(2,0))*u(0,0);
    return out;
}

lin::Mat<double> G(lin::Mat<double> x){
    lin::Mat<double> out(6,1);
    out(0,0) = x(0,0);
    out(1,0) = x(1,0);
    out(2,0) = x(2,0);
    out(3,0) = x(2,0);
    out(4,0) = x(2,0);
    out(5,0) = x(3,0);
    return out;
}

int main(int argc, char *argv[]){

    // Define random generator with Gaussian distribution
    std::default_random_engine generator;
    std::normal_distribution<double> distGPS(0, 1);
    std::normal_distribution<double> distPIR(0, 0.64);
    std::normal_distribution<double> distYasa(0, 0.32);
    std::normal_distribution<double> distAcc(0, 0.5);
    
    lin::Mat<double> X(5,1);
    lin::Mat<double> Xr(5,1);
    Xr =0;
    X = 0;

    lin::Mat<double> R(6,6);
    R = 0;
    R(0,0) = 1; //1.84;
    R(1,1) = 1; //1.84;
    R(3,3) = 0.64*0.64;
    R(5,5) = 0.5*0.5;
    /*
    R(2,2) = 1; //1.00;
    R(4,4) = 0.32*0.32;
    */


    lin::Mat<double> P(5,5);

    lin::Mat<double> Q(5,5);
    Q = Q.I()*0.2;

    lin::Mat<double> Z(6,1);
    Z = 0;

    lin::Mat<double> control(3,1);
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
    std::vector<double> yr;
    std::vector<double> yk;
    std::vector<double> xk;


    double  LO = -2;
    double HI = 4;
    print
    for(int i=0;i<150;i++){


        double a =  LO + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(HI-LO)));
        double s =  LO + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(2)));

        std::cout<<"x:"<<Xr(0,0)<<", v:"<<Xr(2,0)<<", a:"<<a<<", s:"<<s<<std::endl;

        control(0,0) = 0.1;
        control(1,0) = a;
        control(2,0) = s;
        Z(0,0) = Xr(0,0) + distGPS(generator);
        Z(1,0) = Xr(1,0) + distGPS(generator);
        Z(2,0) = Xr(2,0) + distGPS(generator);
        Z(3,0) = Xr(2,0) + distPIR(generator);
        Z(4,0) = Xr(2,0) + distYasa(generator);
        Z(5,0) = a + distAcc(generator);
        vs.push_back(Z(2,0));
        vt.push_back(Z(4,0));

        ekf.predict(Q);
        ekf.update();
        xk.push_back(X(0,0));
        xr.push_back(Xr(0,0));
        yr.push_back(Xr(1,0));
        yk.push_back(X(1,0));
        vk.push_back(X(2,0));
        vr.push_back(Xr(2,0));

        std::cout<<"xk:"<<X(0,0)<<", vk:"<<X(2,0)<<", theta:"<<X(4,0)<<std::endl;
        // env update
        control(1,0) = a + distAcc(generator);
        control(2,0) = s + distAcc(generator);
        Xr = F(Xr, control);
    }

    plt::named_plot("Real speed",vr);
    //plt::named_plot("Sensor 1",vs);
    //plt::named_plot("sensor 2",vt);
    plt::named_plot("kalman speed",vk);
    //plt::named_plot("Real positon",gpsr);
    //plt::named_plot("gps sensor",gpss);
    //plt::named_plot("Kalman position",gpsk);
    plt::named_plot("Kalman position", xk);
    plt::named_plot("Real position", xr);
    plt::legend();
    plt::show();

}

