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
   protected:
        lin::Mat<double>(*F)(lin::Mat<double>, lin::Mat<double>); // f(x,u)
        lin::Mat<double> JF; // f(x,u) jacobian
        lin::Mat<double>(*G)(lin::Mat<double>); // g(x)
        lin::Mat<double> JG; // g(x) jacobian
        lin::Mat<double>* Xk; // state matrix
        lin::Mat<double>* Uk; // control matrix
        lin::Mat<double> Qk; // process noise covariance matrix
        lin::Mat<double> Rk; // sensor noise covariance matrix
        lin::Mat<double> Pk; // process covariance matrix
        lin::Mat<double> K; // kalman gain
        lin::Mat<double>* Zk; // sensor  matrix
        lin::Mat<double>(*jf)(lin::Mat<double>); // compute jacobian at x
        lin::Mat<double>(*jg)(lin::Mat<double>); // compute jacobian at x
        double delta = 0.01; // used to compute jacobians numericaly
        bool  jf_s = false;
        bool  jg_s = false;
        void computeJF(){ // compute jacobian of f(x,u)
                int i,j;
                lin::Mat<double> J(Xk->rows, Xk->rows);
                J = 0;
                JF = J;
                lin::Mat<double> X0, X1, x, Xf;

                for(i=0;i<Xk->rows;i++){
                    if(!jf_s){
                        x = (*Xk);
                        x(i,0) = x(i,0) - delta;
                        X0 = F(x, (*Uk));
                        x(i,0) = x(i,0) + 2*delta;
                        X1 = F(x, (*Uk));
                        Xf = X1 - X0;
                        for(j=0;j<Xk->rows;j++){
                            JF(j,i) = Xf(j,0)/(2*delta);
                        }
                    }else{
                        JF = jf((*Xk));
                        break;
                    }
                }
        }
        void computeJG(){ // compute jacobian of g(x)
                int i,j;
                lin::Mat<double> J(Zk->rows, Xk->rows);
                JG = J;
                lin::Mat<double> X0, X1, x, Xf;

                for(i=0;i<Xk->rows;i++){
                    if(!jg_s){
                        x = (*Xk);
                        x(i,0) = x(i,0) - delta;
                        X0 = G(x);
                        x(i,0) = x(i,0) + 2*delta;
                        X1 = G(x);
                        Xf = X1 - X0;
                        for(j=0;j<Zk->rows;j++){
                            JG(j,i) = Xf(j,0)/(2*delta);
                        }
                    }else{
                        JG = jg((*Xk));
                        break;
                    }
                }
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
        void setFJacobian(lin::Mat<double> (*f)(lin::Mat<double>)){jf = f; jf_s = true;}
        void setGJacobian(lin::Mat<double>(*g)(lin::Mat<double>)){jg = g;jg_s = true;}
        void init(double multiplier){
            lin::Mat<double> pkC(Xk->rows, Xk->rows);
            Pk = pkC.I()*multiplier;
        }
        void predict(){
            computeJF();
            std::cout<<"Xk_:"<<(*Xk);
            (*Xk) =  F((*Xk), (*Uk)); // diff filterpy
            std::cout<<"Xk:"<<(*Xk);
            Pk = JF*Pk*JF.T() + Qk; // diff
            std::cout<<"Pk:"<<Pk;
        }
        void update(){
            lin::Mat<double> Y, S, SI,KH, I_KH,  PHT;
            computeJG();
            std::cout<<"JG:"<<JG;
            PHT = Pk*JG.T();
            std::cout<<"PHT:"<<PHT;
            S = JG*PHT + Rk;
            std::cout<<"S:"<<S;
            bool s;
            SI = S.inverse(s);
            std::cout<<"SI:"<<SI;
            if(s)
                K = PHT*SI;
            else
                std::cout << "Cant compute Kalman gain!"<< std::endl;

            
            std::cout<<"K:"<<K;

            Y = (*Zk) - G((*Xk));
            std::cout<<"Y:"<<Y;
            (*Xk) = (*Xk) + K*Y;
            std::cout<<"Xk__:"<<(*Xk);

            KH = K*JG;
            I_KH = KH.I() - KH; 
            std::cout<<"I_KH:"<<I_KH;
            Pk = I_KH*Pk*I_KH.T() + K*Rk*K.T();
            std::cout<<"Pk_:"<<Pk;

        }
        void predict_update(){
            // predict
            lin::Mat<double> x;
            computeJF();
            std::cout<<"Xk_:"<<(*Xk);
            x =  F((*Xk), (*Uk)); // diff filterpy
            std::cout<<"Xk:"<<x;
            Pk = JF*Pk*JF.T() + Qk; // diff
            std::cout<<"Pk:"<<Pk;

            // update
            
            lin::Mat<double> Y, S, SI,KH, I_KH,  PHT;
            computeJG();
            std::cout<<"JG:"<<JG;
            PHT = Pk*JG.T();
            std::cout<<"PHT:"<<PHT;
            S = JG*PHT + Rk;
            std::cout<<"S:"<<S;
            bool s;
            SI = S.inverse(s);
            std::cout<<"SI:"<<SI;
            if(s)
                K = PHT*SI;
            else
                std::cout << "Cant compute Kalman gain!"<< std::endl;

            
            std::cout<<"K:"<<K;

            Y = (*Zk) - G(x);
            std::cout<<"Y:"<<Y;
            (*Xk) = x + K*Y;
            std::cout<<"Xk__:"<<(*Xk);

            KH = K*JG;
            I_KH = KH.I() - KH; 
            std::cout<<"I_KH:"<<I_KH;
            Pk = I_KH*Pk*I_KH.T() + K*Rk*K.T();
            std::cout<<"Pk_:"<<Pk;

        }
};



}
