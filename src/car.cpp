#include "lin/mat.hpp"
#include "lin/kalman.hpp"
#include <iostream>
#define WITHOUT_NUMPY
#include "matplotlibcpp.h"
#include <iterator>
#include <random>
#include <vector>
#include <cmath>
#include "rapidcsv.h"

#define dtorad 0.0174533

namespace plt = matplotlibcpp;
double dt = 0.1;

lin::Mat<double> F(lin::Mat<double> x, lin::Mat<double> u){
    lin::Mat<double> out(4,1);
    out(0,0) = x(0,0) + x(2,0)*std::cos(x(3,0))*u(0,0); // x = x + v.cos(theta).dt
    out(1,0) = x(1,0) + x(2,0)*std::sin(x(3,0))*u(0,0); // y = y + v.sin(theta).dt
    out(3,0) = x(3,0) + x(2,0)*std::tan(u(1,0))*u(0,0)/u(2,0); // theta = theta + v.tan(delta).dt/l
    out(2,0) = x(2,0) + (u(3,0)/(0.228*250))*u(0,0); // v = v + (T/r*m)*dt
    return out;
}

lin::Mat<double> G(lin::Mat<double> x){
    lin::Mat<double> out(5,1);
    out(0,0) = x(0,0); // gpsx = x
    out(1,0) = x(1,0); // gpsy = y
    out(2,0) = x(2,0); // gpsv = v
    out(3,0) = x(2,0); // vpir = v
    out(4,0) = x(2,0); // vvcu = v
    return out;
}
lin::Mat<double> GPStoCart(double latitude, double longitude){
    lin::Mat<double> out(2,1);
    double R = 6371000;
    out(1,0) =  R * std::cos(latitude*dtorad) * cos(longitude*dtorad);
    out(0,0) =  R * std::cos(latitude*dtorad) * sin(longitude*dtorad);
    return out;
}

int main(int argc, char *argv[]){
    
    lin::Mat<double> X(4,1);
    X = 0;

    lin::Mat<double> R(5,5);
    R = 2;
    R(0,0) = 9.3; // gpsx error
    R(1,1) = 9.3; // gpsy error
    R(2,2) = 8.3; // gpsv error
    R(3,3) = 5; // vpir error
    R(4,4) = 3; // vcu error

    lin::Mat<double> Q(4,4);
    Q = Q.I()*0.005;

    lin::Mat<double> Z(5,1);
    Z = 0;

    lin::Mat<double> control(4,1);
    control(0,0) = dt; // dt
    control(2,0) = 1.5; // l

    kalman::EKFin ekf;

    ekf.setState(&X);
    ekf.setSensorNoiseCovariance(R);
    ekf.setSensor(&Z);
    ekf.setProcessNoiseCovariance(Q);
    ekf.setControl(&control);
    ekf.setF(F);
    ekf.setG(G);
    ekf.init(1);


    // load gps csv
    rapidcsv::Document gpsCSV( "Final_Log_GPS.csv", rapidcsv::LabelParams(0,-1), rapidcsv::SeparatorParams(';'),rapidcsv::ConverterParams(true));
    
    // load gps data into vectors
    std::vector<double> longitude = gpsCSV.GetColumn<double>("Longitude");
    std::vector<double> latitude  = gpsCSV.GetColumn<double>("Latitude");
    std::vector<double> gpsSpeed  = gpsCSV.GetColumn<double>("Speed [km/h]");
    std::vector<double> timeGps  = gpsCSV.GetColumn<double>("Time [s]");
    // compute gps frequency from time data
    double gps_freq = timeGps.size()/(timeGps.back()-timeGps[0]);
    double total_time = timeGps.back()-timeGps[0];

    // load steering csv 
    rapidcsv::Document steerCSV("Final_Log_Steering.csv", rapidcsv::LabelParams(0,-1), rapidcsv::SeparatorParams(';'),rapidcsv::ConverterParams(true));
    // load steering data into vectors
    std::vector<double> steering = steerCSV.GetColumn<double>("Steering Angle [o]");
    std::vector<double> timeSteer  = steerCSV.GetColumn<double>("Time [s]");
    // compute steering sensor frequency from time data
    double steer_freq = timeSteer.size()/(timeSteer.back()-timeSteer[0]);

    // load pir(wheel speed) csv 
    rapidcsv::Document pirCSV("Final_Log_PirF_Filtro.csv", rapidcsv::LabelParams(0,-1), rapidcsv::SeparatorParams(';'),rapidcsv::ConverterParams(true));
    // load pir data into vectors
    std::vector<double> pir = pirCSV.GetColumn<double>("Average Speed (aprox.) [km/h]");
    std::vector<double> timePir  = pirCSV.GetColumn<double>("Time [s]");
    // compute pir sensor frequency from time data
    double pir_freq = timePir.size()/(timePir.back()-timePir[0]);


    // load vcu csv 
    rapidcsv::Document vcuCSV("Final_Log_VCU_Info_1.csv", rapidcsv::LabelParams(0,-1), rapidcsv::SeparatorParams(';'),rapidcsv::ConverterParams(true));
    // load pir data into vectors
    std::vector<double> vcu = vcuCSV.GetColumn<double>("Speed (VCU) [km/h]");
    std::vector<double> vcuTorque = vcuCSV.GetColumn<double>("Reported Torque [Nm]");
    std::vector<double> timeVcu  = vcuCSV.GetColumn<double>("Time [s]");
    // compute pir sensor frequency from time data
    double vcu_freq = timeVcu.size()/(timeVcu.back()-timeVcu[0]);


    double max_freq = vcu_freq;
    int index_gps, index_pir, index_steer, index_vcu;

    int log_ofset = 22000;
    int gpsoffset = -1180;

    index_gps = (gps_freq/max_freq)*(log_ofset+gpsoffset);

    // initialize X 
    lin::Mat<double> cart;
    cart = GPStoCart(latitude[index_gps], longitude[index_gps]);
    X(0,0) = cart(0,0); 
    X(1,0) = cart(1,0); 

    cart = GPStoCart(latitude[index_gps + 5], longitude[index_gps + 5]);
    double dx, dy;
    dx = cart(0,0) - X(0,0); 
    dy = cart(1,0) - X(1,0); 
    std::cout<<"dy/dx:" <<dx/dy<<std::endl;
    X(3,0) = std::atan(dy/dx);
    X(2,0) = pir[0];


    std::vector<double> gpsx, gpsy, gpsv , pirv,vcuv,  x, y, v, theta, delta, dgpsk;


    for(int i = log_ofset;i<36000;i++){
        index_gps = (gps_freq/max_freq)*(i+gpsoffset);
        index_steer = (steer_freq/max_freq)*i;
        index_vcu = (vcu_freq/max_freq)*i;
        index_pir = (pir_freq/max_freq)*i;

        std::cout << "indexvcu "<< index_vcu << " i:" << i << std::endl;
        cart = GPStoCart(latitude[index_gps], longitude[index_gps]);
        Z(0,0) = cart(0,0);
        Z(1,0) = cart(1,0);
        Z(2,0) = gpsSpeed[index_gps];
        Z(3,0) = pir[index_pir];
        Z(4,0) = vcu[index_vcu];
        gpsx.push_back(cart(0,0));
        gpsy.push_back(cart(1,0));
        gpsv.push_back(gpsSpeed[index_gps]);
        pirv.push_back(pir[index_pir]);
        vcuv.push_back(vcu[index_vcu]);

        control(0,0) = 0.18/vcu_freq;
        // codigo precisa dessa constante pra sincronizar os logs ... dt foi determinado manualmente ... NEED FIX
        control(1,0) = steering[index_steer] * dtorad * 0.2727; // multiply by 0.017 to convert from degres to radians
        control(3,0) = vcuTorque[index_vcu];

        delta.push_back(steering[index_steer]);
        ekf.predict_update();

        x.push_back(X(0,0));
        y.push_back(X(1,0));
        v.push_back(X(2,0));
        theta.push_back(X(3,0));

        double d = std::sqrt(std::pow(x.back()-gpsx.back(),2)+std::pow(y.back()-gpsy.back(),2));
        dgpsk.push_back(d);

    }

    plt::subplot(2, 2, 1);
    plt::named_plot("steer",delta);
    plt::named_plot("theta",theta);
    plt::legend();
    //plt::named_plot("pir v", pirv);
    plt::subplot(2, 2, 2);
    plt::named_plot("gps v",gpsv);
    plt::named_plot("pir v", pirv);
    plt::named_plot("vcu v", vcuv);
    plt::named_plot("kalman v",v);
    plt::legend();
    plt::subplot(2, 2, 3);
    plt::named_plot("kalman path", x, y);
    plt::named_plot("gps path", gpsx, gpsy);
    plt::legend();
    plt::subplot(2, 2, 4);
    plt::named_plot("gps and kalman path distance", dgpsk);
    plt::legend();
    plt::show();

}



