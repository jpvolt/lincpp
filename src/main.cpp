#include "lin/mat.hpp"
#include "lin/kalman.hpp"
#include <iostream>

int main(int argc, char *argv[]){
    lin::Mat<double> state(3,1);
    lin::Mat<double> transition(3,3);
    lin::Mat<double> control(3,3);
    lin::Mat<double> test(3,3);
    test(0,0) = 3;
    test(0,1) = 0;
    test(0,2) = 2;
    test(1,0) = 2;
    test(1,1) = 0;
    test(1,2) = -2;
    test(2,0) = 0;
    test(2,1) = 1;
    test(2,2) = 1;
    std::cout<<"det:"<<test.determinant()<<std::endl;
    bool s;
    lin::Mat<double> inv(3,3);
    inv = test.inverse(s);
    std::cout<<"inv:"<<inv;
    kalman::EKF ekf;
    ekf.setState(&state);
    ekf.setTransition(&transition); 
}

