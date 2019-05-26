#include "lin/mat.hpp"
#include "lin/kalman.hpp"
#include <iostream>

int main(int argc, char *argv[]){
    lin::Mat<double> state(3,1);
    lin::Mat<double> transition(3,3);
    lin::Mat<double> control(3,3);
    lin::Mat<int> test(4,4);
    test(0,0) = 1;
    test(0,1) = 3;
    test(0,2) = 5;
    test(0,3) = 9;
    test(1,0) = 1;
    test(1,1) = 3;
    test(1,2) = 1;
    test(1,3) = 7;
    test(2,0) = 4;
    test(2,1) = 3;
    test(2,2) = 9;
    test(2,3) = 7;
    test(3,0) = 5;
    test(3,1) = 2;
    test(3,2) = 0;
    test(3,3) = 9;
    std::cout<<"det:"<<test.determinant()<<std::endl;
    kalman::EKF ekf;
    ekf.setState(&state);
    ekf.setTransition(&transition); 
}
