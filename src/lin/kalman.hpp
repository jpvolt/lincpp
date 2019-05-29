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
            *Pk = (*A)*(*Pk)*(*A).T() + Q; // update covariance
        }
        void update(lin::Mat<double> sensorData){
            bool s;
            *K = *(Pk)*(*Hk).T()*((*Hk)*(*Pk)*(*Hk).T() + (*Rk)).inverse(s); // update kalman gain
            *Xk = (*Xk) + (*K)*(sensorData - ((*Hk)*(*Xk)));
            lin::Mat<double> inter;
            inter = (*K)*(*Hk);
            *Pk = (inter.I() - inter)*(*Pk);
        }
};
class EKFin{
    private:
        lin::Mat<double>(*F)(lin::Mat<double>, lin::Mat<double>);
        lin::Mat<double> JF;
        lin::Mat<double>(*G)(lin::Mat<double>);
        lin::Mat<double> JG;
        lin::Mat<double>* Xk; // state matrix
        lin::Mat<double>* Uk; // control matrix
        lin::Mat<double>* Pk; // process covariance matrix
        lin::Mat<double> K; // kalman gain
        lin::Mat<double>* Rk; // sensor covariance matrix
        lin::Mat<double>* Zk; // sensor  matrix
    public:
        EKFin(){}
        void setState(lin::Mat<double> *statesValue){Xk = statesValue;}
        void setCovariance(lin::Mat<double> *covarianceMat){Pk = covarianceMat;}
        void setSensorError(lin::Mat<double> *sensorErrorMat){Rk = sensorErrorMat;}
        void setSensor(lin::Mat<double> *sensorMat){Zk = sensorMat;}
        void setControl(lin::Mat<double> *controlMat){Uk = controlMat;}
        void setF(lin::Mat<double> (*f)(lin::Mat<double>, lin::Mat<double>)){
            F = f;
        }
        void setG(lin::Mat<double>(*g)(lin::Mat<double>)){
            G = g;
        }
        void computeJacobians(double delta = 0.1){
            int i,j;
            lin::Mat<double> J(Xk->rows, Xk->rows);
            lin::Mat<double> J2(Zk->rows, Xk->rows);
            J = 0;
            J2 = 0;
            JF = J;
            JG = J2;
            lin::Mat<double> X0, X1, x, Xf;

            for(i=0;i<Xk->rows;i++){
                x = (*Xk);
                x(i,0) = x(i,0) - delta;
                X0 = F(x, (*Uk));
                x(i,0) = x(i,0) + 2*delta;
                X1 = F(x, (*Uk));
                Xf = X1 - X0;
                for(j=0;j<Xk->rows;j++){
                    JF(j,i) = Xf(i,0)/2*delta;
                }
                x = (*Xk);
                x(i,0) = x(i,0) - delta;
                X0 = G(x);
                x(i,0) = x(i,0) + 2*delta;
                X1 = G(x);
                Xf = X1 - X0;
                for(j=0;j<Zk->rows;j++){
                    JG(j,i) = Xf(i,0)/2*delta;
                }
            }
            std::cout << " JF :" << JF;
            std::cout << " JG :" << JG;
            std::cout << std::endl;
        }
        void predict(lin::Mat<double> Q){
            (*Xk) =  F((*Xk), (*Uk));
            (*Pk) = JF*(*Pk)*JF.T() + Q;
        }
        void update(){
            bool b;
            K = (*Pk)*JG.T()*(JG*(*Pk)*JG.T() + (*Rk)).inverse(b);
            (*Xk) = (*Xk) + K*((*Zk) - G((*Xk)));
            lin::Mat<double> inter;
            inter = K*JG;
            *Pk = (inter.I() - inter)*(*Pk);
            computeJacobians();
        }
};

}
