// C++11 source code.
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "catch.hpp"
#include "sundials_cxx.hpp"
#include "testing_utils.hpp"

using SVector = sundials_cxx::nvector_serial::Vector;


TEST_CASE( "methods", "[CVodeIntegrator]" ) {
    auto intgr = cvodes_cxx::CVodeIntegrator(cvodes_cxx::LMM::Adams, cvodes_cxx::IterType::Functional);
    std::vector<double> y(1, 1.0);
    double t, yref;
    SVector yout(1);
    intgr.init(rhs_cb, 0.0, &y[0], 1);
    intgr.set_tol(1e-10, 1e-10);
    intgr.step(1.0, yout, &t, cvodes_cxx::Task::Normal);
    yref = std::exp(-t);
    REQUIRE( t > 0 );
    REQUIRE( t <= 1 );
    REQUIRE( std::abs(yout[0] - yref) < 1e-8 );
}


TEST_CASE( "decay_adaptive", "[simple_adaptive]" ) {
    Decay odesys(1.0);
    double y0 = 1.0;
    std::vector<int> root_indices;
    auto tout_yout = cvodes_cxx::simple_adaptive(&odesys, {1e-10}, 1e-10, CV_ADAMS, &y0, 0.0, 1.0, root_indices);
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
    auto tout_yout = cvodes_cxx::simple_adaptive(&odesys, {1e-10}, 1e-10, CV_ADAMS, &y0, 0.0, 1.0, root_indices, 0.0, 0.0, 1e-3, 1100);
    auto& tout = tout_yout.first;
    auto& yout = tout_yout.second;
    REQUIRE( tout.size() == yout.size() );
    for (uint i = 0; i < tout.size(); ++i){
        REQUIRE( std::abs(std::exp(-tout[i]) - yout[i]) < 1e-8 );
    }
    REQUIRE( odesys.last_integration_info["n_steps"] > 998 );
}
