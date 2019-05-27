#include "lin/mat.hpp"
#include "lin/kalman.hpp"
#include <iostream>

int main(int argc, char *argv[]){
    lin::Mat<double> state(3,1);
    lin::Mat<double> transition(3,3);
    lin::Mat<double> control(3,3);
    double a[] = {0,1,2,3,4,5,6,7,8,9};
    control = a;
    std::cout<<control;
    kalman::KF kf;
}

