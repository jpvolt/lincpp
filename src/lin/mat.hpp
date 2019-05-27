#pragma once 
#include <assert.h> 
#include <vector>
#include <iostream>

namespace lin{

template <typename t>
    class Mat{
        private:
            std::vector<std::vector<t>> data;
            t determinant(Mat<t> mat){ // determinant calculator algorithm (Laplace expansion)
                size_t i, j, c, v;
                int sinal = 1;
                t d = 0;
                for(c=0;c<mat.cols;c++){
                    Mat<t> inside(mat.rows-1, mat.cols-1);
                    for(i=1;i<mat.rows;i++){
                        v = 0;
                        for(j=0;j<mat.cols;j++){
                            if(j == c)
                                continue;
                            
                            inside(i-1, v) = data[i][j];
                            v++;
                        }
                    }
                    d += data[0][c] * sinal * inside.determinant();
                    sinal = sinal * -1;
                }
                return d;
            }
        public:
            unsigned int rows;
            unsigned int cols;
            Mat(unsigned int rows, unsigned int cols):rows(rows), cols(cols){ // constructor 
                data.resize(rows);
                for(size_t i = 0;i<data.size();i++){
                    data[i].resize(cols);
                }
            }
            // operators overloading
            void operator=(Mat<t> B){ 
                rows = B.rows;
                cols = B.cols;
                data.clear();
                data.resize(rows);
                for(size_t i = 0;i<data.size();i++){
                    data[i].resize(cols);
                    for(size_t j = 0;j<cols;j++)
                        data[i][j] = B(i,j);
                }
            }
            t& operator ()(int idx, int idy){
                return data[idx][idy];
            }
            t& operator ()(size_t idx, size_t idy){
                return data[idx][idy];
            }
            template<typename y>
            Mat<t> operator+(Mat<y> &B){//sum
                assert(this->rows == B.rows && this->cols == B.cols); // diferent size matrices
                Mat<t> sum(this->rows, this->cols);
                for(size_t i=0; i<this->rows;i++){
                    for(size_t j=0;j<this->cols;j++)
                        sum(i,j) = data[i][j] + B(i,j);
                }
                return sum;
            }
            template<typename y>
            Mat<t> operator-(Mat<y> &B){//sub
                assert(this->rows == B.rows && this->cols == B.cols); // diferent size matrices
                Mat<t> sum(this->rows, this->cols);
                for(size_t i=0; i<this->rows;i++){
                    for(size_t j=0;j<this->cols;j++)
                        sum(i,j) = data[i][j] - B(i,j);
                }
                return sum;
            }
            template<typename y>
            Mat<t> operator*(Mat<y> &B){//mult
                assert(this->cols == B.rows); // check if multiplication is possible
                Mat<t> mult(this->rows, B.cols);
                for(size_t i =0;i<this->rows;i++){
                    for(size_t j=0;j<B.cols;j++){
                        t sum = 0;
                        for(size_t k = 0;k<this->cols;k++)
                            sum += data[i][k]*B(k,j);
                        mult(i,j) = sum;
                    }
                }
                return mult;
            }
            template<typename y>
            Mat<t> operator*(y &B){// scalar mult
                Mat<t> mult(rows, cols);
                for(size_t i =0;i<rows;i++){
                    for(size_t j=0;j<cols;j++){
                        mult(i,j) = data[i][j]*B;
                    }
                }
                return mult;
            }
            template<typename y>
            Mat<t> operator/(y &B){// scalar division
                Mat<t> div(rows, cols);
                for(size_t i =0;i<rows;i++){
                    for(size_t j=0;j<cols;j++){
                        div(i,j) = data[i][j]/B;
                    }
                }
                return div;
            }
            friend std::ostream& operator<< (std::ostream& os, Mat<t> &mat){ // NEED FIX! - dont compile with << std::endl;
                int i,j;
                os << std::endl;
                for(i=0;i<mat.rows;i++){
                    os << "|";
                    for(j=0;j<mat.cols;j++){
                        os << mat(i,j)<<" ";
                    }
                    
                    os << "|"<<std::endl;
                }
                return os;
            }
            Mat<t> T(){
                Mat<t> transposed(rows, cols);
                for(size_t i=0;i<rows;i++){
                    for(size_t j=0;j<cols;j++)
                        transposed(i,j) = data[j][i];
                }
                return transposed;
            }
            // end of operators overloading
            t determinant(){ //determinant calculator function (Laplace expansion)
                if(rows==1)
                    return data[0][0];
                else if(rows==2){
                   t d = (data[0][0]*data[1][1]) - (data[0][1]*data[1][0]);
                   return d;
                }
                else{
                    return determinant(*this);
                }
            }
            Mat<t> getMinorMat(bool cofactor = true){ // compute minors/cofactors mat - used for inverse calculation
                size_t i, j, c, v, b, n;
                int sinal = 1;
                Mat<t> minors(rows, cols);
                for(i=0;i<rows;i++){ // for each row
                    for(j=0;j<cols;j++){ // for each col
                        // fill minor matrix
                        Mat<t> minor(rows-1, cols-1);
                        b = 0;
                        for(c=0;c<rows;c++){
                            if(c == i)
                                continue;
                            n = 0;
                            for(v=0;v<rows;v++){
                                if(v == j)
                                    continue;

                                minor(b,n) = data[c][v];
                                n++;

                            }
                            b++;
                        }
                        std::cout<<"minor:"<<minor;
                        std::cout<<" "<<std::endl;
                        minors(i,j) = minor.determinant();
                        if(cofactor){
                            minors(i,j) = minors(i,j)*sinal;
                        }
                        sinal = sinal * -1;
                    }
                }
                std::cout<<"minors:"<<minors;
                std::cout<<" "<<std::endl;
                return minors;
            }
            Mat<t> inverse(bool &success){ // using Minors, Cofactors and Adjugate
                t d = this->determinant();
                Mat<t> inverse(rows, cols);
                if(d == 0){ // inverse dont exist
                    success = false;
                    return inverse;
                }
                std::cout<<"determinant:"<<d<<std::endl;
                
                Mat<t> cofactors(rows, cols);
                cofactors = getMinorMat();
                cofactors = cofactors.T(); // transpose cofactors matrix to get adjoint matrix
                std::cout<<"cofactors:"<<cofactors;
                std::cout<<" "<<std::endl;

                inverse = cofactors / d;
                return inverse;

            }
    };

}
template<typename t>
bool operator==( lin::Mat<t> A, lin::Mat<t> B){ 
    bool ret = false;
    assert(A.rows == B.rows && A.cols == B.cols); // diferent size matrices
    for(int i = 0;i<A.rows;i++){
        for(int j = 0;j<A.cols;j++){
            if(A(i,j) == B(i,j))
                ret = true;
            else
                ret = false;
        }
    }
    return ret;
}
