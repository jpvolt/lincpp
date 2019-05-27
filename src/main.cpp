#include "lin/mat.hpp"
#include "lin/kalman.hpp"
#include <iostream>

int main(int argc, char *argv[]){
    double dt = 0.1;
    lin::Mat<double> Xk(2,1);
    Xk = 0;
    lin::Mat<double> A(2,2);
    double a[] = {1, dt, 0, 1};
    A = a;

    lin::Mat<double> B(2,1);
    double b[] = {dt*dt/2, dt};
    B = b;

    lin::Mat<double> Hk(2,2);
    double h[] = {1, 0, 0, 1};
    Hk = h;

    lin::Mat<double> Rk(2,2);
    double r[] = {5, 0, 0, 2};
    Rk = r;

    lin::Mat<double> Pk(2,2);
    double p[] = {1, 0, 0, 1};
    Pk = p;

    lin::Mat<double> control(1,1); // aceleration
    control = 0;
    lin::Mat<double> Qk(2,2);
    Qk = 0;

    kalman::KF kf;

    kf.setState(&Xk);
    kf.setCovariance(&Pk);
    kf.setSensorError(&Rk);
    kf.setStateTransition(&A);
    kf.setSensorTransition(&Hk);
    kf.setControlTransition(&B);
    

    lin::Mat<double> sensor(2,1);
    sensor = 0;

    kf.init();
    kf.predict(control, Qk);
    kf.update(sensor);

}

