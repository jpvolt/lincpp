#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "mat.hpp"

TEST_CASE("Matrix multiplication", "[Multiplication]"){
    lin::Mat<double> a(3,3);
    a = 3;
    REQUIRE(a*a.I() == a);
    lin::Mat<double> b(3,2);
    lin::Mat<double> c(3,2);
    b = 1;
    b(0,0) = -1;
    b(1,0) = 3;
    b(2,0) = 0;
    b(0,1) = 2;
    c = 12;
    c(0,0) = 6;
    c(1,0) = 6;
    c(2,0) = 6;
    REQUIRE(a*b == c);
}
TEST_CASE("Matrix mutiplication scalar", "[Multiplication by scalar]"){
    lin::Mat<double> a(3,3);
    a = 3;
    lin::Mat<double> b(3,3);
    b = 9;
    REQUIRE(a*3 == b);
}
TEST_CASE("Matrix division scalar", "[Division by scalar]"){
    lin::Mat<double> a(3,3);
    a = 3;
    lin::Mat<double> b(3,3);
    b = 1;
    REQUIRE(a/3 == b);
}
TEST_CASE("Matrix transpose", "[Transpose]"){
    lin::Mat<double> a(3,2);
    double A[] = {1,2,3,4,5,6};
    a = A;
    lin::Mat<double> b(2,3);
    double B[] = {1,3,5,2,4,6};
    b = B;
    REQUIRE(a.T() == b);
}
TEST_CASE("Matrix indentity", "[indentity]"){
    lin::Mat<double> a(3,3);
    a = 12.2;
    lin::Mat<double> i(3,3);
    i = 0;
    i(0,0) = 1;
    i(1,1) = 1;
    i(2,2) = 1;
    REQUIRE(a.I() == i);
}
TEST_CASE("Determinants are computed", "[Determinant]"){
    lin::Mat<int> test(4,4);
    test(0,0) = 1;
    test(0,1) = 3;
    test(0,2) = 5;
    test(0,3) = 9;
    test(1,0) = 1;
    test(1,1) = 3;
    test(1,2) = 1;
    test(1,3) = 7;
    test(2,0) = 4;
    test(2,1) = 3;
    test(2,2) = 9;
    test(2,3) = 7;
    test(3,0) = 5;
    test(3,1) = 2;
    test(3,2) = 0;
    test(3,3) = 9;
    REQUIRE(test.determinant() == -376);
}
TEST_CASE("Inverse matrix are computed", "[Inverse]"){
    lin::Mat<double> test(3,3);
    lin::Mat<double> testinv(3,3);
    test(0,0) = 3;
    test(0,1) = 0;
    test(0,2) = 2;
    test(1,0) = 2;
    test(1,1) = 0;
    test(1,2) = -2;
    test(2,0) = 0;
    test(2,1) = 1;
    test(2,2) = 1;
    testinv(0,0) = 0.2;
    testinv(0,1) = 0.2;
    testinv(0,2) = -0;
    testinv(1,0) = -0.2;
    testinv(1,1) = 0.3;
    testinv(1,2) = 1;
    testinv(2,0) = 0.2;
    testinv(2,1) = -0.3;
    testinv(2,2) = 0;
    bool success;
    bool equal = test.inverse(success) == testinv;
    REQUIRE(equal);
}
