// C++11 source code.
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "catch.hpp"
#include "cvodes_anyode.hpp"
#include "testing_utils.hpp"

TEST_CASE( "decay_adaptive", "[simple_adaptive]" ) {
    Decay odesys(1.0);
    double y0 = 1.0;
    std::vector<int> root_indices;
    auto tout_yout = cvodes_anyode::simple_adaptive(&odesys, {1e-10}, 1e-10, cvodes_cxx::LMM::Adams, &y0, 0.0, 1.0, root_indices);
    auto& tout = tout_yout.first;
    auto& yout = tout_yout.second;
    REQUIRE( tout.size() == yout.size() );
    for (uint i = 0; i < tout.size(); ++i){
        REQUIRE( std::abs(std::exp(-tout[i]) - yout[i]) < 1e-8 );
    }
    REQUIRE( odesys.last_integration_info["n_steps"] > 1 );
    REQUIRE( odesys.last_integration_info["n_steps"] < 997 );
}


TEST_CASE( "decay_adaptive_dx_max", "[simple_adaptive]" ) {
    Decay odesys(1.0);
    double y0 = 1.0;
    std::vector<int> root_indices;
    auto tout_yout = cvodes_anyode::simple_adaptive(&odesys, {1e-10}, 1e-10, cvodes_cxx::LMM::Adams, &y0, 0.0, 1.0, root_indices, 1100, 0.0, 0.0, 1e-3);
    auto& tout = tout_yout.first;
    auto& yout = tout_yout.second;
    REQUIRE( tout.size() == yout.size() );
    for (uint i = 0; i < tout.size(); ++i){
        REQUIRE( std::abs(std::exp(-tout[i]) - yout[i]) < 1e-8 );
    }
    REQUIRE( odesys.last_integration_info["n_steps"] > 998 );
}
