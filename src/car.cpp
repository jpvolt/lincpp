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

lin::Mat<double> F(lin::Mat<double> X, lin::Mat<double> u){
    lin::Mat<double> out(5,1);
    out(0,0) = X(0,0) + X(2,0) * std::cos(X(3,0)) * u(0,0); // x
    out(1,0) = X(1,0) + X(2,0) * std::sin(X(3,0)) * u(0,0); // y
    out(2,0) = X(2,0) + u(1,0)*u(0,0); // v
    out(3,0) = X(3,0) + X(2,0)*std::tan(u(2,0))*u(0,0); // theta
    out(4,0) = u(1,0); // a
    return out;
}
lin::Mat<double> G(lin::Mat<double> X){
    lin::Mat<double> out(8,1);
    // gps 1 and gps 2
    out(0,0) = X(0,0);
    out(1,0) = X(1,0);
    out(2,0) = X(0,0);
    out(3,0) = X(1,0);
    // speed sensors
    out(4,0) = X(2,0);
    out(5,0) = X(2,0);
    out(6,0) = X(2,0);
    // acelerometer
    out(7,0) = X(4,0);
    return out;
}

std::random_device rd;

class CarSim{
    private:
        std::mt19937 mt = std::mt19937(rd());
        std::uniform_real_distribution<double> dist = std::uniform_real_distribution<double>(0.0, 1.0);
        std::uniform_real_distribution<double> distgps1 = std::uniform_real_distribution<double>(0.0, 1.0);
        std::uniform_real_distribution<double> distgps2 = std::uniform_real_distribution<double>(0.0, 0.89);
        std::uniform_real_distribution<double> distpir = std::uniform_real_distribution<double>(0.0, 0.77);
        std::uniform_real_distribution<double> distyasa = std::uniform_real_distribution<double>(0.0, 0.45);
        std::uniform_real_distribution<double> distacc = std::uniform_real_distribution<double>(0.0, 0.54);
    public:
        double dt, delta, a;
        lin::Mat<double> estado;
        lin::Mat<double> u;
        CarSim(lin::Mat<double> X, double dt){
            estado = X;
            u = lin::Mat<double>(3,1);
        }
        lin::Mat<double> getMeasurement(){
            lin::Mat<double> x = F(estado, u);
            // add error to model
            estado(0,0) = x(0,0) + 0.1*dist(mt);
            estado(1,0) = x(1,0) + 0.1*dist(mt);
            estado(2,0) = x(2,0) + 0.1*dist(mt);
            estado(3,0) = x(3,0) + 0.1*dist(mt);
            estado(4,0) = x(4,0) + 0.1*dist(mt);

            lin::Mat<double> m;
            m = G(x);
            // add error to measurment
            m(0,0) = m(0,0) + distgps1(mt);
            m(1,0) = m(1,0) + distgps1(mt);
            m(2,0) = m(2,0) + distgps2(mt);
            m(3,0) = m(3,0) + distgps2(mt);
            m(4,0) = m(4,0) + distgps1(mt);
            m(5,0) = m(5,0) + distpir(mt);
            m(6,0) = m(6,0) + distyasa(mt);
            m(7,0) = m(7,0) + distacc(mt);
            return m;
        }
        lin::Mat<double> getControl(){
            a = 0.1;
            delta = 1;
            u(0,0) = dt;
            u(1,0) = a;
            u(2,0) = delta;
            return u;
        }
};


int main(int argc, char *argv[]){

    double dt = 0.05;

    kalman::EKFin ekf;
    
    lin::Mat<double> X(5,1);
    X = {10, 0, 0, 0, 0};

    CarSim car(X ,0.05);

    lin::Mat<double> R(8,8);
    R = 0;
    R(0,0) = 1;
    R(1,0) = 1;
    R(2,0) = 0.8;
    R(3,0) = 0.8;
    R(4,0) = 1;
    R(5,0) = 0.6;
    R(6,0) = 0.2;
    R(7,0) = 0.3;

    // copyed from example
    lin::Mat<double> Q(5,5);
    Q = Q.I()*0.01;

    lin::Mat<double> Z(8,1);

    lin::Mat<double> control(3,1);

    ekf.setState(&X);
    ekf.setSensorNoiseCovariance(R);
    ekf.setProcessNoiseCovariance(Q);
    ekf.setSensor(&Z);
    ekf.setControl(&control);
    ekf.setF(F);
    ekf.setG(G);
    ekf.init(50);

    for(int i=0;i<400;i++){
        control = car.getControl();
        Z = car.getMeasurement();
        ekf.predict_update();
       }
   /*
    plt::named_plot("pos real", time, positionr);
    plt::named_plot("pos kalman",time, positionk);
    plt::named_plot("vel real", time, velocityr);
    plt::named_plot("vel kalman",time, velocityk);
    plt::named_plot("alt real", time, altituder);
    plt::named_plot("alt kalman",time, altitudek);
    plt::legend();
    plt::show();
*/
}

