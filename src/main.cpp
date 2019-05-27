#include "lin/mat.hpp"
#include "lin/kalman.hpp"
#include <iostream>

int main(int argc, char *argv[]){
    lin::Mat<double> state(3,1);
    lin::Mat<double> transition(3,3);
    lin::Mat<double> control(3,3);
    kalman::EKF ekf;
    ekf.setState(&state);
    ekf.setTransition(&transition); 
}

