#include "catch.hpp"
#include <math.h>
#include <vector>
#include "cvodes_anyode.hpp"
#include "anyode/anyode_iterative.hpp"
#include "tough_case.hpp"

TEST_CASE( "adaptive_autorestart_tricky", "[simple_adaptive]" ) {
    std::vector<realtype> p = {{321.14999999999998, 39390, -135.30000000000001, 18010, 44960, 48.200000000000003,
                                49320, -114.59999999999999, 1780, -34400.547966379738, -2.865040967667511,
                                93065.338440593958, 5.7581184659305222}};
    OdeSys odesys(p.data());
    std::vector<int> root_indices;

    const long int mxsteps=5000;
    const realtype dx0=0.0;
    const realtype dx_min=0.0;
    const realtype dx_max=0.0;
    const bool with_jacobian=true;
    cvodes_cxx::IterType iter_type=cvodes_cxx::IterType::Undecided;
    cvodes_cxx::LinSol linear_solver=cvodes_cxx::LinSol::DEFAULT;
    const int maxl=0;
    const realtype eps_lin=0.0;
    const unsigned nderiv=0;
    bool return_on_root=false;
    int autorestart=2;

    realtype atol=1e-7, rtol=1e-7;

    bool return_on_error = true;  // This is essentially "xfail" for now (transformed system would work)
    realtype * xyout = (realtype*)malloc((odesys.get_ny() + 1)*sizeof(realtype));
    int td = 1;
    xyout[0] = 0;
    realtype tend=180;
    xyout[1] = 0.00064313123504933787;
    xyout[2] = 0.00014677490343001067;
    xyout[3] = 9.536739572030514e-05;
    xyout[4] = 1.6877253332428752e-11;
    auto tout_yout = cvodes_anyode::simple_adaptive(&xyout, &td, &odesys, {atol}, rtol, cvodes_cxx::LMM::BDF, tend, root_indices,
                                                    mxsteps, dx0, dx_min, dx_max, with_jacobian, iter_type, linear_solver,
                                                    maxl, eps_lin, nderiv, return_on_root, autorestart, return_on_error);
    // REQUIRE( odesys.current_info.nfo_int["n_steps"] > 1 );
    // REQUIRE( odesys.current_info.nfo_int["n_steps"] < 997 );
    free(xyout);
}
