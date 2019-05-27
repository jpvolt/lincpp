#pragma once 
#include <string>
#include "mat.hpp"
namespace kalman{

class KF{
    private:
        lin::Mat<std::string>* statesNames;
        lin::Mat<double>* Xk; // state matrix
        lin::Mat<double>* A; // control to control transition matrix state to state
        lin::Mat<double>* B; // control transition matrix control to state 
        lin::Mat<double>* Pk; // process covariance matrix
        lin::Mat<double>* Hk; // sensor transition matrix -  state to sensor 
        lin::Mat<double>* Rk; // sensor covariance matrix
        lin::Mat<double>* K; // kalman gain
    public:
        KF(){}
        void setStateNames(lin::Mat<std::string> *names){statesNames = names;}
        void setState(lin::Mat<double> *statesValue){Xk = statesValue;}
        void setStateTransition(lin::Mat<double> *stateTransitionMat){A = stateTransitionMat;}
        void setControlTransition(lin::Mat<double> *controlTransitionMat){B = controlTransitionMat;}
        void setCovariance(lin::Mat<double> *covarianceMat){Pk = covarianceMat;}
        void setSensorTransition(lin::Mat<double> *sensorTransitionMat){Hk = sensorTransitionMat;}
        void setSensorError(lin::Mat<double> *sensorErrorMat){Rk = sensorErrorMat;}
        lin::Mat<std::string> getStateNames(){return *statesNames;}
        lin::Mat<double> getState(){return *Xk;}
        void init(){
            K = new lin::Mat<double>(2,2);
            bool s;
            *K = *(Pk)*(*Hk).T()*((*Hk)*(*Pk)*(*Hk).T() + (*Rk)).inverse(s); // initialize kalman gain
        }
        void predict(lin::Mat<double> control, lin::Mat<double> Q){
            *Xk = ((*A)*(*Xk)) + (*B)*control; // update prediction
            std::cout<<"Model:"<<*Xk;
            std::cout<<std::endl;
            std::cout<<"Control:"<<(*B)*control;
            std::cout<<std::endl;
            *Pk = (*A)*(*Pk)*(*A).T() + Q; // update covariance
            std::cout<<"process Covariance:"<<(*Pk);
            std::cout<<std::endl;
        }
        void update(lin::Mat<double> sensorData){
            bool s;
            *K = *(Pk)*(*Hk).T()*((*Hk)*(*Pk)*(*Hk).T() + (*Rk)).inverse(s); // update kalman gain
            std::cout<<"inverse:"<<s<<std::endl;
            *Xk = (*Xk) + (*K)*(sensorData - ((*Hk)*(*Xk)));
            std::cout<<"State:"<<*Xk;
            std::cout<<std::endl;
            lin::Mat<double> inter(2,2);
            inter = (*K)*(*Hk);
            *Pk = (inter.I() - inter)*(*Pk);
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
