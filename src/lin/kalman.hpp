#pragma once 
#include <string>
#include "mat.hpp"
namespace kalman{

class KF{
    private:
        lin::Mat<std::string>* statesNames;
        lin::Mat<double>* state;
        lin::Mat<double>* stateTransitionMat;
        lin::Mat<double>* controlTransitionMat;
    public:
        KF(){}
        void setStateNames(lin::Mat<std::string> *names){statesNames = names;}
        void setState(lin::Mat<double> *statesValue){state = statesValue;}
        void setStateTransition(lin::Mat<double> *stateTransitionMat){this->stateTransitionMat = stateTransitionMat;}
        void setControlTransition(lin::Mat<double> *controlTransitionMat){this->controlTransitionMat = controlTransitionMat;}
        lin::Mat<std::string> getStateNames(){return *statesNames;}
        lin::Mat<double> getState(){return *state;}
        void estimate(lin::Mat<double> control){
            lin::Mat<double> Xk(state->rows, state->cols);
            Xk = ((*stateTransitionMat)*(*state)) + (*controlTransitionMat)*control;
        }
};
class EKF{
    private:
        lin::Mat<std::string>* statesNames;
        lin::Mat<double>* state;
        lin::Mat<double>* transition;
    public:
        EKF(){}
        void setStateNames(lin::Mat<std::string> *names){statesNames = names;}
        void setState(lin::Mat<double> *statesValue){state = statesValue;}
        void setTransition(lin::Mat<double> *transitionMat){transition = transitionMat;}
        lin::Mat<std::string> getStateNames(){return *statesNames;}
        lin::Mat<double> getState(){return *state;}
        void estimate(lin::Mat<double> control){
            lin::Mat<double> Xk(state->rows, state->cols);
            Xk = ((*transition)*(*state)) + control;
        }
};

}
