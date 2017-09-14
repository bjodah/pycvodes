#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "catch.hpp"
#include <math.h>
#include <vector>
#include "anyode/anyode.hpp"
#include "cvodes_anyode.hpp"
#include "cetsa_case_iterative.hpp"

TEST_CASE( "adaptive_tricky_svd", "[simple_adaptive]" ) {
    std::vector<double> p = {{321.14999999999998, 39390, -135.30000000000001, 18010, 44960, 48.200000000000003,
                              49320, -114.59999999999999, 1780, -34400.547966379738, -2.865040967667511,
                              93065.338440593958, 5.7581184659305222}};
    std::vector<double> y0 = {{0.00043976285661326595, 0.00010031950134333891, 7.8448068109052759e-05,
                               8.4073190114798924e-11, 0.0003814695739342757}};
    double t0=0, tend=180;
    OdeSys odesys(&p[0]);
    std::vector<int> root_indices;

    const long int mxsteps=45000;
    const realtype dx0=0.0;
    const realtype dx_min=0.0;
    const realtype dx_max=0.0;
    const bool with_jacobian=true;
    cvodes_cxx::IterType iter_type=cvodes_cxx::IterType::Undecided;
    const int maxl=0;
    const realtype eps_lin=0.0;
    const unsigned nderiv=0;
    bool return_on_root=false;

    double atol=1e-8, rtol=1e-8;

    int linear_solver=10; // 10 => GMRES
    int autorestart=0;
    bool return_on_error = false;
    bool with_jtimes = false;
    auto tout_yout = cvodes_anyode::simple_adaptive(&odesys, {atol}, rtol, cvodes_cxx::LMM::BDF, &y0[0], t0, tend, root_indices,
                                                    mxsteps, dx0, dx_min, dx_max, with_jacobian, iter_type, linear_solver,
                                                    maxl, eps_lin, nderiv, return_on_root, autorestart, return_on_error, with_jtimes);
    auto& tout = tout_yout.first;
    auto& yout = tout_yout.second;
    const int ref = tout.size() * odesys.get_ny();
    REQUIRE( ref == yout.size() );
    REQUIRE( odesys.last_integration_info["n_steps"] > 1 );
    REQUIRE( odesys.last_integration_info["n_steps"] < 997 );
}
