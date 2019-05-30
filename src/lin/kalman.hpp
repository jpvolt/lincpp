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
        lin::Mat<double>(*F)(lin::Mat<double>, lin::Mat<double>); // f(x,u)
        lin::Mat<double> JF; // f(x,u) jacobian
        lin::Mat<double>(*G)(lin::Mat<double>); // g(x)
        lin::Mat<double> JG; // g(x) jacobian
        lin::Mat<double>* Xk; // state matrix
        lin::Mat<double>* Uk; // control matrix
        lin::Mat<double> Qk; // process noise covariance matrix
        lin::Mat<double> Pk; // process covariance matrix
        lin::Mat<double> Rk; // sensor noise covariance matrix
        lin::Mat<double> K; // kalman gain
        lin::Mat<double>* Zk; // sensor  matrix
        void computeJacobians(double delta = 0.01){
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
                    JF(j,i) = Xf(j,0)/(2*delta);
                }
                x = (*Xk);
                x(i,0) = x(i,0) - delta;
                X0 = G(x);
                x(i,0) = x(i,0) + 2*delta;
                X1 = G(x);
                Xf = X1 - X0;
                for(j=0;j<Zk->rows;j++){
                    JG(j,i) = Xf(j,0)/(2*delta);
                }
            }
            std::cout << " X :" << (*Xk);
            std::cout << " JF :" << JF;
            std::cout << " JG :" << JG;
            std::cout << std::endl;
        }
    public:
        EKFin(){}
        void setState(lin::Mat<double> *statesValue){Xk = statesValue;}
        void setProcessNoiseCovariance(lin::Mat<double> covarianceMat){Qk = covarianceMat;}
        void setSensorNoiseCovariance(lin::Mat<double> covarianceMat){Rk = covarianceMat;}
        void setSensor(lin::Mat<double> *sensorMat){Zk = sensorMat;}
        void setControl(lin::Mat<double> *controlMat){Uk = controlMat;}
        void setF(lin::Mat<double> (*f)(lin::Mat<double>, lin::Mat<double>)){F = f;}
        void setG(lin::Mat<double>(*g)(lin::Mat<double>)){G = g;}
        void init(){
            lin::Mat<double> pkC(Xk->rows, Xk->rows);
            Pk = pkC.I();
        }
        void predict(){
            (*Xk) =  F((*Xk), (*Uk));
            computeJacobians();
            Pk = JF*Pk*JF.T() + Qk;
        }
        void update(){
            lin::Mat<double> Yk, Sk, Ik;
            Yk = (*Zk) - G((*Xk));
            Sk = JG*Pk*JG.T() + Rk;
            bool b;
            K = Pk*JG.T()*Sk.inverse(b);
            (*Xk) = (*Xk) + K*Yk;
            Ik = K*JG;
            Pk = (Ik.I() - Ik)*Pk;
        }
};

}
