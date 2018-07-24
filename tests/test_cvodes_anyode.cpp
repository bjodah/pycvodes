// C++11 source code.
#include "catch.hpp"
#include "cvodes_anyode.hpp"
#include "testing_utils.hpp"

TEST_CASE( "decay_adaptive", "[simple_adaptive]" ) {
    Decay odesys(1.0);
    std::vector<int> root_indices;
    int td = 1;
    realtype * xyout = (realtype *)malloc(td*(odesys.get_ny()+1)*sizeof(realtype));
#define xout(ti) xyout[2*ti]
#define yout(ti) xyout[2*ti + 1]
    xout(0) = 0.0;
    yout(0) = 1.0;
    auto nout = cvodes_anyode::simple_adaptive(&xyout, &td, &odesys, {1e-10}, 1e-10, cvodes_cxx::LMM::Adams, 1.0, root_indices);
    for (int i = 0; i < nout; ++i){
        REQUIRE( std::abs(std::exp(-xout(i)) - yout(i)) < 1e-8 );
    }
    REQUIRE( odesys.current_info.nfo_int["n_steps"] > 1 );
    REQUIRE( odesys.current_info.nfo_int["n_steps"] >= nout );
    REQUIRE( odesys.current_info.nfo_int["n_steps"] < 997 );
#undef xout
#undef yout
    free(xyout);
}

TEST_CASE( "decay_adaptive_ew_ele", "[simple_adaptive]" ) {
    Decay odesys(1.0);
    std::vector<int> root_indices;
    int td = 1;
    realtype * xyout = (realtype *)malloc(td*(odesys.get_ny()+1)*sizeof(realtype));
    realtype * ew_ele = (realtype *)malloc(2*td*odesys.get_ny()*sizeof(realtype));
#define xout(ti) xyout[2*ti]
#define yout(ti) xyout[2*ti + 1]
    xout(0) = 0.0;
    yout(0) = 1.0;
    const long int mxsteps=0;
    realtype dx0=0.0;
    const realtype dx_min=0.0;
    const realtype dx_max=0.0;
    const bool with_jacobian=false;
    cvodes_cxx::IterType iter_type=cvodes_cxx::IterType::Undecided;
    cvodes_cxx::LinSol linear_solver=cvodes_cxx::LinSol::DEFAULT;
    const int maxl=0;
    const realtype eps_lin=0.0;
    const unsigned nderiv=0;
    bool return_on_root=false;
    int autorestart=0;
    bool return_on_error=false;
    bool with_jtimes=false;
    int tidx=0;
    auto nout = cvodes_anyode::simple_adaptive(
        &xyout, &td, &odesys, {1e-10}, 1e-10, cvodes_cxx::LMM::Adams, 1.0, root_indices,
        mxsteps, dx0, dx_min, dx_max, with_jacobian, iter_type, linear_solver, maxl, eps_lin,
        nderiv, return_on_root, autorestart, return_on_error, with_jtimes, tidx, &ew_ele);
    for (int i = 0; i < nout; ++i){
        REQUIRE( std::abs(std::exp(-xout(i)) - yout(i)) < 1e-8 );
    }
    REQUIRE( odesys.current_info.nfo_int["n_steps"] > 1 );
    REQUIRE( odesys.current_info.nfo_int["n_steps"] >= nout );
    REQUIRE( odesys.current_info.nfo_int["n_steps"] < 997 );
#undef xout
#undef yout
    free(xyout);
    free(ew_ele);
}

TEST_CASE( "decay_predefined_ew_ele", "[simple_predefined]" ) {
    Decay odesys(1.0);
    int nt = 37;
    realtype t0=0, tend=4.0;
    std::vector<realtype> tout(nt);
    std::vector<realtype> yqout(nt*(odesys.get_ny()+odesys.get_nquads()));
    std::vector<realtype> ew_ele(2*nt*odesys.get_ny());

    for (int i=0; i<nt; ++i){
        tout[i] = t0 + i*(tend - t0)/(nt-1);
    }
    yqout[0] = 1.0;
    std::vector<int> root_indices;
    std::vector<realtype> root_out;
    const long int mxsteps=0;
    realtype dx0=0.0;
    const realtype dx_min=0.0;
    const realtype dx_max=0.0;
    const bool with_jacobian=false;
    cvodes_cxx::IterType iter_type=cvodes_cxx::IterType::Undecided;
    cvodes_cxx::LinSol linear_solver=cvodes_cxx::LinSol::DEFAULT;
    const int maxl=0;
    const realtype eps_lin=0.0;
    const unsigned nderiv=0;
    int autorestart=0;
    bool return_on_error=false;
    bool with_jtimes=false;
    auto nout = cvodes_anyode::simple_predefined(
        &odesys, {1e-10}, 1e-10, cvodes_cxx::LMM::BDF, yqout.data(), nt, tout.data(),
        yqout.data(), root_indices, root_out, mxsteps, dx0, dx_min, dx_max, with_jacobian, iter_type,
        linear_solver, maxl, eps_lin, nderiv, autorestart, return_on_error, with_jtimes, ew_ele.data()
        );
    REQUIRE(nout == nt);

    for (int i=0; i < nt; ++i){
        realtype t = tout[i];
        REQUIRE( std::abs(yqout[i] - 1.0*exp(-t)) < 1e-6 );
        REQUIRE( std::abs(ew_ele[2*i] * ew_ele[2*i + 1]) < 1.0 );
    }
    REQUIRE( odesys.current_info.nfo_int["n_steps"] > 1 );
    REQUIRE( odesys.current_info.nfo_int["n_steps"] < 997 );
}

TEST_CASE( "decay_adaptive_get_dx_max", "[simple_adaptive]" ) {
    Decay odesys(1.0);
    realtype y0 = 1.0;
    std::vector<int> root_indices;
    odesys.use_get_dx_max = true;
    int td = 1;
    realtype * xyout = (realtype *)malloc(td*(odesys.get_ny()+1)*sizeof(realtype));
#define xout(ti) xyout[2*ti]
#define yout(ti) xyout[2*ti + 1]
    xout(0) = 0.0;
    yout(0) = 1.0;
    auto nout = cvodes_anyode::simple_adaptive(&xyout, &td, &odesys, {1e-10}, 1e-10, cvodes_cxx::LMM::Adams, 1.0, root_indices, 1005);
    for (int i = 0; i < nout; ++i){
        REQUIRE( std::abs(std::exp(-xout(i)) - yout(i)) < 1e-8 );
    }
    REQUIRE( odesys.current_info.nfo_int["n_steps"] > 1000 );  // dx_max == 1e-3
#undef xout
#undef yout
    free(xyout);
}


TEST_CASE( "decay_adaptive_dx_max", "[simple_adaptive]" ) {
    Decay odesys(1.0);
    realtype y0 = 1.0;
    std::vector<int> root_indices;
    int td = 1;
    realtype * xyout = (realtype *)malloc(td*(odesys.get_ny()+1)*sizeof(realtype));
#define xout(ti) xyout[2*ti]
#define yout(ti) xyout[2*ti + 1]
    xout(0) = 0.0;
    yout(0) = 1.0;
    auto nout = cvodes_anyode::simple_adaptive(&xyout, &td, &odesys, {1e-10}, 1e-10, cvodes_cxx::LMM::Adams, 1.0, root_indices, 1100, 0.0, 0.0, 1e-3);
    for (int i = 0; i < nout; ++i){
        REQUIRE( std::abs(std::exp(-xout(i)) - yout(i)) < 1e-8 );
    }
    REQUIRE( odesys.current_info.nfo_int["n_steps"] > 998 );
    REQUIRE( odesys.current_info.nfo_int["n_steps"] >= nout );
    REQUIRE( nout >= 998 );
    free(xyout);
}
