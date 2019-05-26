#pragma once 
#include <assert.h> 
#include <vector>
#include <iostream>

namespace lin{

template <typename t>
    class Mat{
        private:
            std::vector<std::vector<t>> data;
        public:
            unsigned int rows;
            unsigned int cols;
            Mat(unsigned int rows, unsigned int cols):rows(rows), cols(cols){
                data.resize(rows);
                for(size_t i = 0;i<data.size();i++){
                    data[i].resize(cols);
                }
            }
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
            friend std::ostream& operator<< (std::ostream& os, Mat<t> &mat){
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
            t determinant(Mat<t> mat){
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
            t determinant(){
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
    };

}
