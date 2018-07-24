// C++11 source code.
#include "catch.hpp"
#include "sundials_cxx.hpp" // realtype
#include <vector>

using sundials_cxx::nvector_serial::Vector;
using sundials_cxx::nvector_serial::VectorView;


TEST_CASE( "methods" "[sundials::nvector_serial::Vector]" ) {
    auto vec1 = Vector(3);
    REQUIRE(vec1.size() == 3);
    vec1.set_all(42);
    REQUIRE(vec1[2] == 42);
    auto vec2 = Vector(3);
    vec2 = vec1;
    REQUIRE(vec2[2] == 42);
    vec1[2] = 7;
    REQUIRE(vec1[2] == 7);
    REQUIRE(vec2[2] == 42);
    std::vector<realtype> a(3, 9);
    vec2.load(&a[0]);
    REQUIRE(vec2[0] == 9);
}

TEST_CASE( "methods" "[sundials::nvector_serial::VectorView]" ) {
    auto vec1 = Vector(3);
    REQUIRE(vec1.size() == 3);
    vec1.set_all(42);
    REQUIRE(vec1[2] == 42);
    auto vec2 = VectorView(vec1);
    REQUIRE(vec2[2] == 42);
    vec1[2] = 7;
    REQUIRE(vec1[2] == 7);
    REQUIRE(vec2[2] == 7);
    std::vector<realtype> a(3, 9);
    vec2.load(&a[0]);
    REQUIRE(vec1[0] == 9);
}
