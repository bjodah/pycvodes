// C++11 source code.
#include "catch.hpp"
#include "cvodes_anyode.hpp"
#include "testing_utils.hpp"

TEST_CASE( "decay_adaptive", "[simple_adaptive]" ) {
    Decay<double> odesys(1.0);
    std::vector<int> root_indices;
    int td = 1;
    double * xyout = (double *)malloc(td*(odesys.get_ny()+1)*sizeof(double));
#define xout(ti) xyout[2*ti]
#define yout(ti) xyout[2*ti + 1]
    xout(0) = 0.0;
    yout(0) = 1.0;
    auto nout = cvodes_anyode::simple_adaptive(&xyout, &td, &odesys, {1e-10}, 1e-10, cvodes_cxx::LMM::Adams, 1.0, root_indices);
    for (int i = 0; i < nout; ++i){
        REQUIRE( std::abs(std::exp(-xout(i)) - yout(i)) < 1e-8 );
    }
    REQUIRE( odesys.last_integration_info["n_steps"] > 1 );
    REQUIRE( odesys.last_integration_info["n_steps"] >= nout );
    REQUIRE( odesys.last_integration_info["n_steps"] < 997 );
#undef xout
#undef yout
    free(xyout);
}

TEST_CASE( "decay_adaptive_get_dx_max", "[simple_adaptive]" ) {
    Decay<double> odesys(1.0);
    double y0 = 1.0;
    std::vector<int> root_indices;
    odesys.use_get_dx_max = true;
    int td = 1;
    double * xyout = (double *)malloc(td*(odesys.get_ny()+1)*sizeof(double));
#define xout(ti) xyout[2*ti]
#define yout(ti) xyout[2*ti + 1]
    xout(0) = 0.0;
    yout(0) = 1.0;
    auto nout = cvodes_anyode::simple_adaptive(&xyout, &td, &odesys, {1e-10}, 1e-10, cvodes_cxx::LMM::Adams, 1.0, root_indices, 1005);
    for (int i = 0; i < nout; ++i){
        REQUIRE( std::abs(std::exp(-xout(i)) - yout(i)) < 1e-8 );
    }
    REQUIRE( odesys.last_integration_info["n_steps"] > 1000 );  // dx_max == 1e-3
#undef xout
#undef yout
    free(xyout);
}


TEST_CASE( "decay_adaptive_dx_max", "[simple_adaptive]" ) {
    Decay<double> odesys(1.0);
    double y0 = 1.0;
    std::vector<int> root_indices;
    int td = 1;
    double * xyout = (double *)malloc(td*(odesys.get_ny()+1)*sizeof(double));
#define xout(ti) xyout[2*ti]
#define yout(ti) xyout[2*ti + 1]
    xout(0) = 0.0;
    yout(0) = 1.0;
    auto nout = cvodes_anyode::simple_adaptive(&xyout, &td, &odesys, {1e-10}, 1e-10, cvodes_cxx::LMM::Adams, 1.0, root_indices, 1100, 0.0, 0.0, 1e-3);
    for (int i = 0; i < nout; ++i){
        REQUIRE( std::abs(std::exp(-xout(i)) - yout(i)) < 1e-8 );
    }
    REQUIRE( odesys.last_integration_info["n_steps"] > 998 );
    REQUIRE( odesys.last_integration_info["n_steps"] >= nout );
    REQUIRE( nout >= 998 );
    free(xyout);
}
