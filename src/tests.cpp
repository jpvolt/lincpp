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
    REQUIRE(test.determinant() == -476);
}
