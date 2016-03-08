// C++11 source code.
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "catch.hpp"
#include "cvodes_cxx.hpp"
#include "sundials_cxx.hpp"

using SVector = sundials_cxx::nvector_serial::Vector;

int rhs_cb(double t, N_Vector y, N_Vector f, void * user_data){
    cvodes_cxx::ignore(t); cvodes_cxx::ignore(user_data);
    NV_DATA_S(f)[0] = -NV_DATA_S(y)[0];
    return 0;
}

TEST_CASE( "methods" "[CVodeIntegrator]" ) {
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
