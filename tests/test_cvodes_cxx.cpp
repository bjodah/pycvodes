// C++11 source code.
#include "catch.hpp"
#include "cvodes_cxx.hpp" // realtype
#include "testing_utils.hpp"

using SVector = sundials_cxx::nvector_serial::Vector;

int rhs_cb(realtype /* t */, N_Vector y, N_Vector f, void * /* user_data */){
    NV_DATA_S(f)[0] = -NV_DATA_S(y)[0];
    return 0;
}

int rhs_cb3(realtype /* t */, N_Vector y, N_Vector f, void * /* user_data */){
    NV_DATA_S(f)[0] = -NV_DATA_S(y)[0];
    NV_DATA_S(f)[1] = NV_DATA_S(y)[0]-NV_DATA_S(y)[1];
    NV_DATA_S(f)[2] = NV_DATA_S(y)[1]-NV_DATA_S(y)[2];
    return 0;
}


TEST_CASE( "methods", "[Integrator]" ) {
    auto intgr = cvodes_cxx::Integrator(cvodes_cxx::LMM::Adams, cvodes_cxx::IterType::Functional);
    std::vector<realtype> y(1, 1.0);
    realtype t, yref;
    SVector yout(1);
    intgr.init(rhs_cb, 0.0, &y[0], 1);
    intgr.set_tol(1e-10, 1e-10);
    intgr.step(1.0, yout, &t, cvodes_cxx::Task::Normal);
    yref = std::exp(-t);
    REQUIRE( t > 0 );
    REQUIRE( t <= 1 );
    REQUIRE( std::abs(yout[0] - yref) < 1e-8 );
}

realtype get_dx_max(realtype /* x */, const realtype * const /* y */){
    return 1e-3;
}

TEST_CASE( "adaptive", "[Integrator]" ) {
    const int ny = 1;
    auto intgr = cvodes_cxx::Integrator(cvodes_cxx::LMM::Adams, cvodes_cxx::IterType::Functional);
    std::vector<int> root_indices;
    int td = 1;  // trailing dimension
    realtype * xyout = (realtype*)malloc(td*(ny+1)*sizeof(realtype));
    xyout[0] = 0.0;
    xyout[1] = 1.0;
#define xout(idx) xyout[idx*2]
#define yout(idx) xyout[idx*2+1]
    REQUIRE( xout(0) == 0.0 );
    REQUIRE( yout(0) == 1.0 );
    bool return_on_root=false, return_on_error=false;
    int autorestart=0;
    intgr.init(rhs_cb, 0.0, xyout+1, 1);
    intgr.set_tol(1e-9, 1e-9);
    intgr.set_max_num_steps(1005);
    const realtype xend = 1.0;
    const int nderiv=0;
    auto nt = intgr.adaptive(&xyout, &td, xend, nderiv, root_indices, return_on_root, autorestart, return_on_error, get_dx_max);
    REQUIRE( xout(0) == 0.0 );
    REQUIRE( yout(0) == 1.0 );
    for (int idx=1; idx<nt; ++idx){
        const realtype yref = std::exp(-xout(idx));
        REQUIRE( xout(idx) > 0 );
        REQUIRE( xout(idx) <= 1 );
        REQUIRE( std::abs(yout(idx) - yref) < 1e-8 );
    }
    free(xyout);
#undef xout
#undef yout
}

realtype get_dx_max2(realtype /* x */, const realtype * const /* y */){
    return 1e-3;
}


TEST_CASE( "predefined", "[Integrator]" ) {
    const int ny = 1;
    auto intgr = cvodes_cxx::Integrator(cvodes_cxx::LMM::Adams, cvodes_cxx::IterType::Functional);
    std::vector<realtype> y(1, 1.0);
    std::vector<int> root_indices;
    bool return_on_error=false;
    int autorestart=0;
    int nderiv = 0;
    intgr.init(rhs_cb, 0.0, &y[0], 1);
    intgr.set_tol(1e-9, 1e-9);
    intgr.set_max_num_steps(1005);
    int nt = 702;
    realtype tend = 2.0;
    std::vector<realtype> tout(702);
    for (int idx=0; idx<nt; ++idx)
        tout[idx] = idx*tend/(nt - 1);
    std::vector<realtype> yout(nt*ny);
    std::vector<realtype> root_out;
    auto result = intgr.predefined(nt, &tout[0], &y[0], &yout[0], nderiv,
                                   root_indices, root_out, autorestart, return_on_error, get_dx_max2);
    REQUIRE( result == 702 );
    REQUIRE( intgr.get_n_steps() > 1000 );
    for (int idx=0; idx<nt; ++idx){
        const realtype yref = std::exp(-tout[idx]);
        REQUIRE( std::abs(yout[idx] - yref) < 1e-8 );
    }
}

TEST_CASE( "predefined_autorestart", "[Integrator]" ) {
    const int ny = 3;
    auto intgr = cvodes_cxx::Integrator(cvodes_cxx::LMM::Adams, cvodes_cxx::IterType::Functional);
    std::vector<realtype> y(ny, 1.0);
    std::vector<int> root_indices;
    bool return_on_error=false;
    int autorestart=4;
    int nderiv = 0;
    intgr.init(rhs_cb3, 0.0, &y[0], ny);
    intgr.set_tol(1e-9, 1e-9);
    intgr.set_max_num_steps(4);
    int nt = 174;
    realtype tend = 2.0;
    std::vector<realtype> tout(nt);
    for (int idx=0; idx<nt; ++idx)
        tout[idx] = idx*tend/(nt - 1);
    std::vector<realtype> yout(nt*ny);
    std::vector<realtype> root_out;
    auto result2 = intgr.predefined(nt, &tout[0], &y[0], &yout[0], nderiv,
                                    root_indices, root_out, autorestart, return_on_error, get_dx_max2);
    REQUIRE( result2 == nt );
    REQUIRE( intgr.get_n_steps() > 1000 );
    for (int idx=0; idx<nt; ++idx){
        const realtype yref = std::exp(-tout[idx]);
        //std::cout << yout[idx*ny] << " " << yref << " " << yout[idx*ny] - yref << "\n";
        REQUIRE( std::abs(yout[idx*ny] - yref) < 1e-7 );
    }
}
