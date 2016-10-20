#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "catch.hpp"
#include <math.h>
#include <vector>
#include "anyode/anyode.hpp"
#include "cvodes_anyode.hpp"
#include "cetsa_case.hpp"


TEST_CASE( "adaptive_autorestart", "[simple_adaptive]" ) {
    std::vector<double> p = {{298.15, 39390, -135.3, 18010, 44960, 48.2, 65919.5, -93.8304, 1780, 3790, 57.44, 19700, -157.4}};
    std::vector<double> y0 = {{8.99937e-07, 0.000693731, 0.000264211, 0.000340312, 4.11575e-05}};
    double t0=0, tend=60;
    OdeSys odesys(&p[0]);
    std::vector<int> root_indices;

    const long int mxsteps=0;
    const realtype dx0=0.0;
    const realtype dx_min=0.0;
    const realtype dx_max=0.0;
    const bool with_jacobian=true;
    cvodes_cxx::IterType iter_type=cvodes_cxx::IterType::Undecided;
    int linear_solver=0;
    const int maxl=0;
    const realtype eps_lin=0.0;
    const unsigned nderiv=0;
    bool return_on_root=false;
    int autorestart=2;

    auto tout_yout = cvodes_anyode::simple_adaptive(&odesys, {1e-8}, 1e-8, cvodes_cxx::LMM::BDF, &y0[0], t0, tend, root_indices,
                                                    mxsteps, dx0, dx_min, dx_max, with_jacobian, iter_type, linear_solver,
                                                    maxl, eps_lin, nderiv, return_on_root, autorestart);
    auto& tout = tout_yout.first;
    auto& yout = tout_yout.second;
    const int ref = tout.size() * odesys.get_ny();
    REQUIRE( ref == yout.size() );
    REQUIRE( odesys.last_integration_info["n_steps"] > 1 );
    REQUIRE( odesys.last_integration_info["n_steps"] < 997 );
}
