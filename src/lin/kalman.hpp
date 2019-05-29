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
class EKFin{
    private:
        std::vector<lin::Mat<double>(*)(lin::Mat<double>, lin::Mat<double>)> F;
        lin::Mat<double> JF;
        std::vector<lin::Mat<double>(*)(lin::Mat<double>)> G;
        lin::Mat<double> JG;
        lin::Mat<double> Xk_1; // old state matrix
        lin::Mat<double>* Xk; // state matrix
        lin::Mat<double> Uk_1; //  old control matrix
        lin::Mat<double> Uk; // control matrix
    public:
        EKFin(){}
        void setState(lin::Mat<double> *statesValue){Xk = statesValue;}
        void addToF(lin::Mat<double> (*f)(lin::Mat<double>, lin::Mat<double>)){
            F.push_back(f);
        }
        void addToG(lin::Mat<double>(*g)(lin::Mat<double>)){
            G.push_back(g);
        }
        void computeJacobians(double delta = 0.01){
            int i,j;
            lin::Mat<double> J(Xk->rows, Xk->cols);
            J = 0;
            JF = J;
            JG = J;
            lin::Mat<double> X0, X1, x,u;

            X0 = F[i](Xk_1, Uk_1);
            for(i=0;i<Xk->rows;i++){
                x = Xk_1;
                x(i,0) = x(i,0) + delta;
                X1 = F[i](x, Uk);
                for(j=0;j<Xk->rows;j++){
                    JF(j,i) = (x(i,0) - Xk_1(i,0))/delta;
                }
            }

        }
};

}
