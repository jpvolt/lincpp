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


namespace plt = matplotlibcpp;
double dt = 0.01;

lin::Mat<double> F(lin::Mat<double> x, lin::Mat<double> u){
    lin::Mat<double> out(4,1);
    out(0,0) = x(0,0) + x(2,0)*std::cos(x(3,0))*u(0,0); // x = x + v.cos(theta).dt
    out(1,0) = x(1,0) + x(2,0)*std::sin(x(3,0))*u(0,0); // y = y + v.sin(theta).dt
    out(3,0) = x(3,0) + x(2,0)*std::tan(u(1,0))*u(0,0)/u(2,0); // theta = theta + v.tan(delta).dt/l
    out(2,0) = x(2,0); // v = v
    return out;
}

lin::Mat<double> G(lin::Mat<double> x){
    lin::Mat<double> out(4,1);
    out(0,0) = x(0,0); // gpsx = x
    out(1,0) = x(1,0); // gpsy = y
    out(2,0) = x(2,0); // vpir = v
    out(3,0) = x(2,0); // vvcu = v
    return out;
}
lin::Mat<double> GPStoCart(double latitude, double longitude){
    lin::Mat<double> out(2,1);
    double R = 6371000;
    out(0,0) =  R * std::cos(latitude*0.017) * cos(longitude*0.017);
    out(1,0) =  R * std::cos(latitude*0.017) * sin(longitude*0.017);
    return out;
}

int main(int argc, char *argv[]){
    
    lin::Mat<double> X(4,1);
    X = 0;

    lin::Mat<double> R(4,4);
    R = 0;
    R(0,0) = 3.3; // gpsx error
    R(1,1) = 3.3; // gpsy error
    R(2,2) = 2; // vpir error
    R(3,3) = 100; // vcu error

    lin::Mat<double> Q(4,4);
    Q = Q.I();

    lin::Mat<double> Z(4,1);
    Z = 0;

    lin::Mat<double> control(3,1);
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
    ekf.init(100);


    // load gps csv
    rapidcsv::Document gpsCSV( "Final_Log_GPS.csv", rapidcsv::LabelParams(0,-1), rapidcsv::SeparatorParams(';'),rapidcsv::ConverterParams(true));
    
    // load gps data into vectors
    std::vector<double> longitude = gpsCSV.GetColumn<double>("Longitude");
    std::vector<double> latitude  = gpsCSV.GetColumn<double>("Latitude");
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
    rapidcsv::Document pirCSV("Final_Log_PirF.csv", rapidcsv::LabelParams(0,-1), rapidcsv::SeparatorParams(';'),rapidcsv::ConverterParams(true));
    // load pir data into vectors
    std::vector<double> pir = pirCSV.GetColumn<double>("Average Speed (aprox.) [km/h]");
    std::vector<double> timePir  = pirCSV.GetColumn<double>("Time [s]");
    // compute pir sensor frequency from time data
    double pir_freq = timePir.size()/(timePir.back()-timePir[0]);


    // load vcu csv 
    rapidcsv::Document vcuCSV("Final_Log_VCU_Info_1.csv", rapidcsv::LabelParams(0,-1), rapidcsv::SeparatorParams(';'),rapidcsv::ConverterParams(true));
    // load pir data into vectors
    std::vector<double> vcu = vcuCSV.GetColumn<double>("Speed (VCU) [km/h]");
    std::vector<double> timeVcu  = vcuCSV.GetColumn<double>("Time [s]");
    // compute pir sensor frequency from time data
    double vcu_freq = timeVcu.size()/(timeVcu.back()-timeVcu[0]);

    int index_gps, index_pir, index_steer, index_vcu;

    int log_ofset = 31000;

    index_gps = (gps_freq/pir_freq)*log_ofset;

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


    std::vector<double> gpsx, gpsy,  pirv,vcuv,  x, y, v, theta;

    for(int i = log_ofset;i<33000;i++){
        index_gps = (gps_freq/pir_freq)*i;
        index_steer = (steer_freq/pir_freq)*i;
        index_vcu = (vcu_freq/pir_freq)*i;
        index_pir = i;

        cart = GPStoCart(latitude[index_gps], longitude[index_gps]);
        Z(0,0) = cart(0,0);
        gpsx.push_back(Z(0,0));
        Z(1,0) = cart(1,0);
        gpsy.push_back(Z(1,0));
        Z(2,0) = pir[index_pir];
        pirv.push_back(Z(2,0));
        //Z(3,0) = vcu[index_vcu];
        vcuv.push_back(Z(3,0));

        control(1,0) = steering[index_steer] * 0.017 ; // multiply by 0.017 to convert from degres to radians

        ekf.predict_update();

        x.push_back(X(0,0));
        y.push_back(X(1,0));
        v.push_back(X(2,0));
        theta.push_back(X(3,0));

    }


    //plt::named_plot("Gps x", gpsx);
    //plt::named_plot("kalman x",x);
    plt::named_plot("pir v", pirv);
    plt::named_plot("vcu v", vcuv);
    plt::named_plot("kalman v",v);
    plt::legend();
    plt::show();

}



