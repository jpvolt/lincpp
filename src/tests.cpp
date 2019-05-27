#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "lin/mat.hpp"
#include "lin/kalman.hpp"

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
