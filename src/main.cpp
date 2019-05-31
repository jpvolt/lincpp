#include "lin/mat.hpp"
#include "lin/kalman.hpp"
#include <iostream>
#define WITHOUT_NUMPY
#include "matplotlibcpp.h"
#include <iterator>
#include <random>
#include <vector>
#include <cmath>
#include<fstream>
#include "random.hpp"

#define print std::cout<<"here"<<std::endl;

namespace plt = matplotlibcpp;
using Random = effolkronium::random_static;


class RadarSim{
    public:
        double pos, vel, alt, dt;
        RadarSim(double pos, double vel, double alt, double dt):pos(pos), vel(vel), alt(alt), dt(dt){}
        double getRange(double rand){
            vel = vel + 0.1*rand;
            alt = alt + 0.1*rand;
            pos = pos + vel*dt;

            double err = pos * 0.05*rand;
            double slant_dist = std::sqrt(std::pow(pos,2)+std::pow(alt,2));

            return slant_dist + err;
        }
};

    lin::Mat<double> out(3,1);
    lin::Mat<double> A(3,3);
    A = A.I();
    A(0,1) = 0.05;
    out = A*X;
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

    double dt = 0.05;

    RadarSim radar(0,100,1000,dt);
    kalman::EKFin ekf;
    
    lin::Mat<double> X(3,1);
    X = {radar.pos-100, radar.vel+100, radar.alt+1000};

    lin::Mat<double> R(1,1);
    R = 25;

    // copyed from example
    lin::Mat<double> Q(3,3);
    Q = {1.5625e-07, 6.2500e-06, 0.0000e+00,
        6.2500e-06, 2.5000e-04, 0.0000e+00,
        0.0000e+00, 0.0000e+00, 1.0000e-01};


    lin::Mat<double> pkC(X.rows, X.rows);
    ekf.Pk = pkC.I()*50;

    lin::Mat<double> Z(1,1);
    Z = 0;

    lin::Mat<double> control(1,1);
    control = dt;


    ekf.setState(&X);
    ekf.setSensorNoiseCovariance(R);
    ekf.setProcessNoiseCovariance(Q);
    ekf.setSensor(&Z);
    ekf.setControl(&control);
    ekf.setF(F);
    ekf.setG(G);
    ekf.setGJacobian(JG);
    ekf.init(1);

    std::vector<double> positionr;
    std::vector<double> velocityr;
    std::vector<double> altituder;
    std::vector<double> positionk;
    std::vector<double> velocityk;
    std::vector<double> altitudek;
    std::vector<double> time;
    std::vector<double> zs;

    std::ifstream myReadFile;
     myReadFile.open("in.txt");

     std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    for(int i=0;i<400;i++){

        double p,v,a,z;
        
        myReadFile>>p;
        myReadFile>>v;
        myReadFile>>a;
        myReadFile>>z;

     
        Z = radar.getRange(dist(mt));  
        p = radar.pos;
        v = radar.vel;
        a = radar.alt;
   
        positionr.push_back(p);
        velocityr.push_back(v);
        altituder.push_back(a);

        ekf.predict_update();

        positionk.push_back(X(0,0));
        velocityk.push_back(X(1,0));
        altitudek.push_back(X(2,0));


        time.push_back((i*dt));
        
    }

    plt::named_plot("pos real", time, positionr);
    plt::named_plot("pos kalman",time, positionk);
    plt::named_plot("vel real", time, velocityr);
    plt::named_plot("vel kalman",time, velocityk);
    plt::named_plot("alt real", time, altituder);
    plt::named_plot("alt kalman",time, altitudek);
    plt::legend();
    plt::show();

}

