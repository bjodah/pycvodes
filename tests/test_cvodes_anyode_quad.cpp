#include "catch.hpp"
#include <math.h>
#include <vector>
#include "anyode/anyode.hpp"
#include "cvodes_anyode.hpp"

struct OdeSys : public AnyODE::CvodesOdeSysBase {
    indextype get_ny() const override {
        return 1;
    }
    int get_nquads() const override {
        return 2;
    }
    AnyODE::Status rhs(realtype /* t */,
                       const realtype * const ANYODE_RESTRICT y,
                       realtype * const ANYODE_RESTRICT f) override
    {
        f[0] = -0.7*y[0]; // y(t) = y(0)*exp(-0.7*t)
        return AnyODE::Status::success;
    }
    AnyODE::Status dense_jac_cmaj(realtype /* t */,
                                  const realtype * const ANYODE_RESTRICT /* y */,
                                  const realtype * const ANYODE_RESTRICT /* fy */,
                                  realtype * const ANYODE_RESTRICT jac,
                                  long int /* ldim */,
                                  realtype * const ANYODE_RESTRICT /*dfdt*/=nullptr) override
    {
        jac[0] = -0.7;
        return AnyODE::Status::success;
    }

    AnyODE::Status quads(realtype xval, const realtype * const y, realtype * const out) override {
        out[0] = xval*y[0]; //
        out[1] = y[0]*y[0];
        return AnyODE::Status::success;
    }

};

realtype integral_A_t_exp_minus_k_t(realtype A, realtype k, realtype t){
    return A/pow(k, 2) + (-A*pow(k, 2)*t - A*k)*exp(-k*t)/pow(k, 3);
}

realtype integral_A_exp_minus_k_t__squared(realtype A, realtype k, realtype t){
    return (1.0L/2.0L)*pow(A, 2)/k - 1.0L/2.0L*pow(A, 2)*exp(-2*k*t)/k;
}

TEST_CASE( "quadrature_adaptive", "[simple_adaptive]" ) {
    auto odesys = OdeSys();
    int td = 1;
    realtype * xyqout = (realtype*)malloc(td*(odesys.get_ny()+odesys.get_nquads()+1)*sizeof(realtype));
    realtype A = 42.0;
    realtype k = 0.7;
    xyqout[0] = 0; // t0
    xyqout[1] = A; // y0
    xyqout[2] = 2.0; // q0
    xyqout[3] = 3.0; // q1
    realtype t0=0, tend=4.0;
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

    auto nout = cvodes_anyode::simple_adaptive(
        &xyqout, &td, &odesys, {1e-10, 1e-11, 1e-10}, 1e-10, cvodes_cxx::LMM::BDF, tend, root_indices,
        mxsteps, dx0, dx_min, dx_max, with_jacobian, iter_type, linear_solver,
        maxl, eps_lin, nderiv, return_on_root, autorestart);
    REQUIRE((nout + 1) == td);

    for (int i=0; i < nout; ++i){
        realtype t = xyqout[i*4];
        REQUIRE( std::abs(xyqout[i*4+1] - 42*exp(-k*t)) < 1e-6 );
        realtype q0 = 2 + integral_A_t_exp_minus_k_t(A, k, t);
        realtype q1 = 3 + integral_A_exp_minus_k_t__squared(A, k, t);
        REQUIRE( std::abs(xyqout[i*4+2] - q0) < 1e-4 );
        REQUIRE( std::abs(xyqout[i*4+3] - q1) < 1e-4 );
    }
    REQUIRE( odesys.current_info.nfo_int["n_steps"] > 1 );
    REQUIRE( odesys.current_info.nfo_int["n_steps"] < 997 );
    REQUIRE( nout > 1 );
    REQUIRE( nout < 997 );
    free(xyqout);
}

TEST_CASE( "quadrature_predefined", "[simple_predefined]" ) {
    auto odesys = OdeSys();
    int nt = 37;
    realtype t0=0, tend=4.0;
    std::vector<realtype> tout(nt);
    std::vector<realtype> yqout(nt*(odesys.get_ny()+odesys.get_nquads()));
    for (int i=0; i<nt; ++i){
        tout[i] = t0 + i*(tend - t0)/(nt-1);
    }
    realtype A = 42.0;
    realtype k = 0.7;
    yqout[0] = A; // y0
    yqout[1] = 2.0; // q0
    yqout[2] = 3.0; // q1
    std::vector<int> root_indices;
    std::vector<realtype> root_out;

    auto nout = cvodes_anyode::simple_predefined(
        &odesys, {1e-10, 1e-11, 1e-10}, 1e-10, cvodes_cxx::LMM::BDF, yqout.data(), nt, tout.data(),
        yqout.data(), root_indices, root_out);
    REQUIRE(nout == nt);

    for (int i=0; i < nt; ++i){
        realtype t = tout[i];
        REQUIRE( std::abs(yqout[i*3] - 42*exp(-k*t)) < 1e-6 );
        realtype q0 = 2 + integral_A_t_exp_minus_k_t(A, k, t);
        realtype q1 = 3 + integral_A_exp_minus_k_t__squared(A, k, t);
        REQUIRE( std::abs(yqout[i*3+1] - q0) < 1e-4 );
        REQUIRE( std::abs(yqout[i*3+2] - q1) < 1e-4 );
    }
    REQUIRE( odesys.current_info.nfo_int["n_steps"] > 1 );
    REQUIRE( odesys.current_info.nfo_int["n_steps"] < 997 );
}
