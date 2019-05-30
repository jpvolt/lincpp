# lincpp
simple linear algebra c++ header only library - educational purpose 


# Usage

```C++
#include "lin/mat.hpp"

using namespace lin;

int main(int argc, char* argv[]){
  
  //definition Mat<type> name(rows, cols);
  Mat<double> n(2,2);
  
  // fill matrix with value
  n = 0;
  // | 0.0 0.0 |
  // | 0.0 0.0 |
  
  n = {0.0, 1.0, 2.0, 3.0};
  // | 0.0 1.0 |
  // | 2.0 3.0 |
  
  n = n.T(); // matrix transpose
  // | 0.0 2.0 |
  // | 1.0 3.0 |
  
  n = n.I() // matrix identity
  // | 1.0 0.0 |
  // | 0.0 1.0 |
  
  // multiplication by scalar
  n = n*2.5;
  // | 2.5 0.0 |
  // | 0.0 2.5 |
  
  // division by scalar
  n = n/5;
  // | 0.5 0.0 |
  // | 0.0 0.5 |
  
  // matrix multiplication
  n = n*n.I();
  // | 0.5 0.0 |
  // | 0.0 0.5 |
  
  // Matrix sum and subtraction
  n = n - n*0.1;
  // | 0.4 0.0 |
  // | 0.0 0.4 |
  
  bool success;
  n = n.inverse(success); // matrix inverse
  
  // element access n(row, col)
  n(1,0) = 20;
  // | 0.4 0.0 |
  // | 20.0 0.4 |
  
  
}

```
