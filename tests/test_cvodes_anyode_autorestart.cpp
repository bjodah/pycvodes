#include "catch.hpp"
#include <math.h>
#include <vector>
#include "anyode/anyode.hpp"
#include "cvodes_anyode.hpp"
#include "cvodes_cxx.hpp"
#include "cetsa_case.hpp"


TEST_CASE( "adaptive_autorestart", "[simple_adaptive]" ) {
    std::vector<realtype> p = {{298.15, 39390, -135.3, 18010, 44960, 48.2, 65919.5, -93.8304, 1780, 3790, 57.44, 19700, -157.4}};
    OdeSys odesys(&p[0]);
    int td = 1;
    realtype * xyout = (realtype*)malloc(td*(odesys.get_ny()+1)*sizeof(realtype));
    xyout[0] = 0; // t0
    xyout[1] = 8.99937e-07;
    xyout[2] = 0.000693731;
    xyout[3] = 0.000264211;
    xyout[4] = 0.000340312;
    xyout[5] = 4.11575e-05;
    realtype t0=0, tend=60;
    std::vector<int> root_indices;

    const long int mxsteps=0;
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

    auto nout = cvodes_anyode::simple_adaptive(&xyout, &td, &odesys, {1e-8}, 1e-8, cvodes_cxx::LMM::BDF, tend, root_indices,
                                               mxsteps, dx0, dx_min, dx_max, with_jacobian, iter_type, linear_solver,
                                               maxl, eps_lin, nderiv, return_on_root, autorestart);
    REQUIRE( odesys.current_info.nfo_int["n_steps"] > 1 );
    REQUIRE( odesys.current_info.nfo_int["n_steps"] < 997 );
    REQUIRE( nout > 1 );
    REQUIRE( nout < 997 );
    free(xyout);
}
