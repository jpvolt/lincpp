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
    lin::Mat<double> Xk(2,1);
    Xk = 100;
    lin::Mat<double> A(2,2);
    double a[] = {1, dt, 0, 1};
    A = a;
    std::cout<<"A:"<<A;
    lin::Mat<double> B(2,1);
    B(0,0) = dt*dt/2;
    B(1,0) = dt;
    lin::Mat<double> Hk(3,2);
    double h[] = {1 , dt, 0, 1 , 0, 1};
    Hk = h;
    std::cout<<"Hk:"<<Hk;
    lin::Mat<double> Rk(3,3);
    Rk = Rk.I();
    Rk(0,0) = 2;
    Rk(1,1) = 1.64;
    Rk(2,2) = 1.3;
    std::cout<<"Rk:"<<Rk;
    lin::Mat<double> Pk(2,2);
    Pk = Pk.I();
    std::cout<<"Pk:"<<Pk;
    lin::Mat<double> control(1,1);

    lin::Mat<double> Qk(2,2);
    Qk = Qk.I()*0.1;
    std::cout<<"Qk:"<<Qk;
    lin::Mat<double> Zk(3,1);

    // Define random generator with Gaussian distribution
    const double mean = 0.0;
    const double stddev = 0.3;
    std::default_random_engine generator;
    std::normal_distribution<double> dist(mean, stddev);

    kalman::KF kf;

    kf.setState(&Xk);
    kf.setCovariance(&Pk);
    kf.setSensorError(&Rk);
    kf.setStateTransition(&A);
    kf.setSensorTransition(&Hk);
    kf.setControlTransition(&B);
    kf.init();

    std::vector<double> vr;
    std::vector<double> vs;
    std::vector<double> vt;
    std::vector<double> vk;
    std::vector<double> gpsr;
    std::vector<double> gpss;
    std::vector<double> gpsk;

    double gps = 0;
    double v = 0;
    for(int i=0;i<900;i++){

        double  LO = -2;
        double HI = 4;
        double a =  LO + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(HI-LO)));

        //control = a + dist(generator)*0.3;
        control = 0;
        std::cout<<"pos:"<<gps<<", v:"<<v<<", a:"<<a<<std::endl;
        Zk(0,0) = gps + 6*dist(generator);
        Zk(1,0) = v + 3*dist(generator) + 2;
        Zk(2,0) = v + 2*dist(generator) - 1.5;
        gpss.push_back(Zk(0,0));
        vs.push_back(Zk(1,0));
        vt.push_back(Zk(2,0));

        kf.predict(control, Qk);
        kf.update(Zk);

        // env update
        v = v + a*dt;
        gps = gps + v*dt;
        gpsr.push_back(gps);
        vr.push_back(v);

        gpsk.push_back(Xk(0,0));
        vk.push_back(Xk(1,0));

    }

    plt::named_plot("Real speed",vr);
    plt::named_plot("Sensor 1",vs);
    plt::named_plot("sensor 2",vt);
    plt::named_plot("kalman speed",vk);
    plt::named_plot("Real positon",gpsr);
    plt::named_plot("gps sensor",gpss);
    plt::named_plot("Kalman position",gpsk);
    plt::plot(gpsr);
    plt::plot(gpss);
    plt::plot(gpsk);
    plt::legend();
    plt::show();

}

