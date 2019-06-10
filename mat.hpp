#pragma once 
#include <assert.h> 
#include <vector>
#include <iostream>
#include <initializer_list>


namespace lin{

template <typename t>
    class Mat{
        private:
            std::vector<std::vector<t>> data;
            t determinant(const Mat<t> mat) const{ // determinant calculator algorithm (Laplace expansion)
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
            Mat(){}
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
            template<typename y>
            void operator=(y B){ 
                for(size_t i = 0;i<rows;i++){
                    for(size_t j = 0;j<cols;j++)
                        data[i][j] = B;
                }
            }
            template<typename y>
            void operator=(y d[]){
                for(int i=0;i<rows;i++){
                    for(int j=0;j<cols;j++)
                        data[i][j] = d[cols*i + j];
                }
            }
            void operator=(std::initializer_list<t> l){
                std::vector<t> temp;
                temp.insert(temp.end(), l.begin(), l.end());
                assert(l.size() == rows*cols);
                for(int i=0;i<rows;i++){
                    for(int j=0;j<cols;j++){
                        data[i][j] = temp[cols*i +j]; 
                    }
                }
            }
            t& operator ()(int idx, int idy){
                assert(idx<rows && idy<cols);
                return data[idx][idy];
            }
            t operator ()(int idx, int idy) const{
                assert(idx<rows && idy<cols);
                return data[idx][idy];
            }
            template<typename y>
            Mat<t> operator+(const Mat<y> &B) const{//sum
                assert(this->rows == B.rows && this->cols == B.cols); // diferent size matrices
                Mat<t> sum(this->rows, this->cols);
                for(size_t i=0; i<this->rows;i++){
                    for(size_t j=0;j<this->cols;j++)
                        sum(i,j) = data[i][j] + B(i,j);
                }
                return sum;
                }
            template<typename y>
            Mat<t> operator-(const Mat<y> &B) const{//sub
                assert(this->rows == B.rows && this->cols == B.cols); // diferent size matrices
                Mat<t> sum(this->rows, this->cols);
                for(size_t i=0; i<this->rows;i++){
                    for(size_t j=0;j<this->cols;j++)
                        sum(i,j) = data[i][j] - B(i,j);
                }
                return sum;
            }
            template<typename y>
            Mat<t> operator*(const Mat<y> &B) const{//mult
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
            Mat<t> operator*(const y &B) const{// scalar mult
                Mat<t> mult(rows, cols);
                for(size_t i =0;i<rows;i++){
                    for(size_t j=0;j<cols;j++){
                        mult(i,j) = data[i][j]*B;
                    }
                }
                return mult;
            }
            template<typename y>
            Mat<t> operator/(const y &B) const{// scalar division
                Mat<t> div(rows, cols);
                for(size_t i =0;i<rows;i++){
                    for(size_t j=0;j<cols;j++){
                        div(i,j) = data[i][j]/B;
                    }
                }
                return div;
            }
            friend std::ostream& operator<< (std::ostream& os, const Mat<t> &mat){ // NEED FIX! - dont compile with << std::endl;
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
            Mat<t> T() const{
                Mat<t> transposed(cols, rows);
                for(size_t i=0;i<rows;i++){
                    for(size_t j=0;j<cols;j++){
                        transposed(j,i) = data[i][j];
                }
            }
                return transposed;
            }
            Mat<t> I() const{
                assert(rows == cols); // square matrix
                Mat<t> indentity(rows, rows);
                indentity = 0; // fill matrix with 0s
                for(int i =0;i<rows;i++)
                    indentity(i,i) = 1;

                return indentity;
            }
            // end of operators overloading
            void pSize(){
                std::cout<<"[ "<<rows<<" x "<<cols<< " ]"<<std::endl;
            }
            t determinant() const{ //determinant calculator function (Laplace expansion)
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
            Mat<t> getMinorMat(bool cofactor = true) const{ // compute minors/cofactors mat - used for inverse calculation
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
                        minors(i,j) = minor.determinant();
                        if(cofactor){
                            minors(i,j) = minors(i,j)*sinal;
                        }
                        sinal = sinal * -1;
                    }
                }
                return minors;
            }
            Mat<t> inverse(bool &success) const{ // using Minors, Cofactors and Adjugate
                success = true;
                t d = this->determinant();
                if(rows == 1 && cols == 1){ // scalars
                    Mat<t> r(1,1);
                    r = 1/data[0][0];
                    return r;
                }
                Mat<t> inverse(rows, cols);
                if(d == 0){ // inverse dont exist
                    std::cout<<"invert failed!"<<std::endl;
                    success = false;
                    return inverse;
                }
                
                Mat<t> cofactors(rows, cols);
                cofactors = getMinorMat();
                cofactors = cofactors.T(); // transpose cofactors matrix to get adjoint matrix

                inverse = cofactors / d;
                return inverse;

            }
    };

}
template<typename t>
bool operator==(const lin::Mat<t> A, const lin::Mat<t> B){ 
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
