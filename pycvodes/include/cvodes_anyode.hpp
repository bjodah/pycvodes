#pragma once

#include <chrono>
#include <ctime>
#include <functional>
#include <memory>

#include "anyode/anyode.hpp"
#include "anyode/anyode_buffer.hpp"
#include "cvodes_cxx.hpp"

BEGIN_NAMESPACE(cvodes_anyode)

using cvodes_cxx::Integrator;
using cvodes_cxx::LMM;
using cvodes_cxx::IterType;
using cvodes_cxx::IterLinSolEnum;
using cvodes_cxx::PrecType;
using cvodes_cxx::GramSchmidtType;


inline int handle_status_(AnyODE::Status status){
    switch (status){
    case AnyODE::Status::success:
        return 0;
    case AnyODE::Status::recoverable_error:
        return 1;
    case AnyODE::Status::unrecoverable_error:
        return -1;
    default:
        throw std::runtime_error(cvodes_cxx::StreamFmt() <<
                                 "Got an unhandled status: " << static_cast<int>(status));
    }
}

template<class OdeSys>
int rhs_cb(realtype t, N_Vector y, N_Vector ydot, void *user_data){
    auto t_start = std::chrono::high_resolution_clock::now();
    auto& odesys = *static_cast<OdeSys*>(user_data);
    if (odesys.record_rhs_xvals)
        odesys.current_info.nfo_vecdbl["rhs_xvals"].push_back(t);
    AnyODE::Status status = odesys.rhs(t, NV_DATA_S(y), NV_DATA_S(ydot));
    static_cast<Integrator*>(odesys.integrator)->time_rhs += std::chrono::duration<double>(
        std::chrono::high_resolution_clock::now() - t_start).count();
    return handle_status_(status);
}

template<class OdeSys>
int roots_cb(realtype t, N_Vector y, realtype *gout, void *user_data){
    auto t_start = std::chrono::high_resolution_clock::now();
    auto& odesys = *static_cast<OdeSys*>(user_data);
    AnyODE::Status status = odesys.roots(t, NV_DATA_S(y), gout);
    if (status == AnyODE::Status::recoverable_error)
        throw std::runtime_error("There are only unrecoverable errors for roots().");
    static_cast<Integrator*>(odesys.integrator)->time_roots += std::chrono::duration<double>(
        std::chrono::high_resolution_clock::now() - t_start).count();
    return handle_status_(status);
}

template<class OdeSys>
int quads_cb(realtype t, N_Vector y, N_Vector yQdot, void *user_data){
    auto t_start = std::chrono::high_resolution_clock::now();
    auto& odesys = *static_cast<OdeSys*>(user_data);
    AnyODE::Status status = odesys.quads(t, NV_DATA_S(y), NV_DATA_S(yQdot));
    if (status == AnyODE::Status::recoverable_error)
        throw std::runtime_error("There are only unrecoverable errors for quads().");
    static_cast<Integrator*>(odesys.integrator)->time_quads += std::chrono::duration<double>(
        std::chrono::high_resolution_clock::now() - t_start).count();
    return handle_status_(status);
}

template <class OdeSys>
int jac_dense_cb(
#if SUNDIALS_VERSION_MAJOR < 3
    long int N,
#endif
    realtype t,
    N_Vector y, N_Vector fy,
#if SUNDIALS_VERSION_MAJOR < 3
    DlsMat Jac,
#else
    SUNMatrix Jac,
#endif
    void *user_data,
    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3
    ){
    // callback of req. signature wrapping OdeSys method.
#if SUNDIALS_VERSION_MAJOR < 3
    AnyODE::ignore(N);
#endif
    AnyODE::ignore(tmp1); AnyODE::ignore(tmp2); AnyODE::ignore(tmp3);
    auto t_start = std::chrono::high_resolution_clock::now();
    auto& odesys = *static_cast<OdeSys*>(user_data);
    if (odesys.record_jac_xvals)
        odesys.current_info.nfo_vecdbl["jac_xvals"].push_back(t);
    AnyODE::Status status = odesys.dense_jac_cmaj(t, NV_DATA_S(y), NV_DATA_S(fy),
#if SUNDIALS_VERSION_MAJOR < 3
                                                  DENSE_COL(Jac, 0), Jac->ldim
#else
                                                  SM_DATA_D(Jac), odesys.get_ny()
#endif
        );

    static_cast<Integrator*>(odesys.integrator)->time_jac += std::chrono::duration<double>(
        std::chrono::high_resolution_clock::now() - t_start).count();
    return handle_status_(status);
}

template <typename OdeSys>
int jac_band_cb(
#if SUNDIALS_VERSION_MAJOR < 3
    long int N, long int mupper, long int mlower,
#endif
    realtype t,
    N_Vector y, N_Vector fy,
#if SUNDIALS_VERSION_MAJOR < 3
    DlsMat Jac,
#else
    SUNMatrix Jac,
#endif
    void *user_data,
    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
    AnyODE::ignore(tmp1); AnyODE::ignore(tmp2); AnyODE::ignore(tmp3);
    auto t_start = std::chrono::high_resolution_clock::now();
    auto& odesys = *static_cast<OdeSys*>(user_data);
#if SUNDIALS_VERSION_MAJOR < 3
    if (odesys.get_mupper() != mupper)
        throw std::runtime_error("mupper mismatch");
    if (odesys.get_mlower() != mlower)
        throw std::runtime_error("mlower mismatch");
    auto Jac_ = Jac;
    AnyODE::ignore(N);
#else
    auto Jac_ = SM_CONTENT_B(Jac);
#endif
    if (odesys.record_jac_xvals)
        odesys.current_info.nfo_vecdbl["jac_xvals"].push_back(t);
    AnyODE::Status status = odesys.banded_jac_cmaj(t, NV_DATA_S(y), NV_DATA_S(fy), Jac_->data + Jac_->s_mu - Jac_->mu, Jac_->ldim);
    static_cast<Integrator*>(odesys.integrator)->time_jac += std::chrono::duration<double>(
        std::chrono::high_resolution_clock::now() - t_start).count();
    return handle_status_(status);
}


template <typename OdeSys>
int jac_times_vec_cb(N_Vector v, N_Vector Jv, realtype t, N_Vector y,
                     N_Vector fy, void *user_data, N_Vector tmp){
    // callback of req. signature wrapping OdeSys method.
    AnyODE::ignore(tmp);
    auto t_start = std::chrono::high_resolution_clock::now();
    auto& odesys = *static_cast<OdeSys*>(user_data);
    AnyODE::Status status = odesys.jac_times_vec(NV_DATA_S(v), NV_DATA_S(Jv), t, NV_DATA_S(y), NV_DATA_S(fy));
    if (status == AnyODE::Status::recoverable_error)
        throw std::runtime_error("There are only unrecoverable errors for JacTimesVec().");
    static_cast<Integrator*>(odesys.integrator)->time_jtimes += std::chrono::duration<double>(
        std::chrono::high_resolution_clock::now() - t_start).count();
    return handle_status_(status);
}

template <typename OdeSys> // Section 4.6.9 Preconditioning in cvs_guide.pdf
int prec_solve_cb(realtype t, N_Vector y, N_Vector fy, N_Vector r,
                  N_Vector z, realtype gamma, realtype delta, int lr,
                  void *user_data
#if SUNDIALS_VERSION_MAJOR < 3
                  , N_Vector tmp
#endif
    ){
    // callback of req. signature wrapping OdeSys method.
#if SUNDIALS_VERSION_MAJOR < 3
    AnyODE::ignore(tmp);  // delta used for iterative methods
#endif
    auto t_start = std::chrono::high_resolution_clock::now();
    double * ewt {nullptr};
    auto& odesys = *static_cast<OdeSys*>(user_data);
    if (lr != 1)
        throw std::runtime_error("Only left preconditioning implemented.");
    AnyODE::Status status =  odesys.prec_solve_left(t, NV_DATA_S(y), NV_DATA_S(fy), NV_DATA_S(r),
                                                    NV_DATA_S(z), gamma, delta, ewt);
    static_cast<Integrator*>(odesys.integrator)->time_prec += std::chrono::duration<double>(
        std::chrono::high_resolution_clock::now() - t_start).count();
    return handle_status_(status);
}

template <typename OdeSys>
int prec_setup_cb(realtype t, N_Vector y, N_Vector fy, booleantype jok,
                  booleantype *jcurPtr, realtype gamma, void *user_data
#if SUNDIALS_VERSION_MAJOR < 3
                  ,N_Vector tmp1, N_Vector tmp2, N_Vector tmp3
#endif
    ){
    // callback of req. signature wrapping OdeSys method.
#if SUNDIALS_VERSION_MAJOR < 3
    AnyODE::ignore(tmp1); AnyODE::ignore(tmp2); AnyODE::ignore(tmp3);
#endif
    auto t_start = std::chrono::high_resolution_clock::now();
    auto& odesys = *static_cast<OdeSys*>(user_data);
    bool jac_recomputed = false;
    AnyODE::Status status = odesys.prec_setup(t, NV_DATA_S(y), NV_DATA_S(fy), jok, jac_recomputed, gamma);
    (*jcurPtr) = (jac_recomputed) ? SUNTRUE : SUNFALSE;
    static_cast<Integrator*>(odesys.integrator)->time_prec += std::chrono::duration<double>(
        std::chrono::high_resolution_clock::now() - t_start).count();
    return handle_status_(status);
}

template <class OdeSys>
std::unique_ptr<Integrator> get_integrator(
    OdeSys * odesys,
    std::vector<realtype> &atol,
    const realtype rtol,
    const LMM lmm,
    const realtype * const yq0,
    const realtype t0,
    const long int mxsteps=0,
    const realtype dx0=0.0,
    const realtype dx_min=0.0,
    const realtype dx_max=0.0,
    const bool with_jacobian=false,
    const IterType iter_type=IterType::Newton,
    const int linear_solver=0,
    const int maxl=0,
    const realtype eps_lin=0.0,
    const bool with_jtimes=false,
    const std::vector<realtype> &constraints={})
{
    const int ny = odesys->get_ny();
    const int nroots = odesys->get_nroots();
    const int nq = odesys->get_nquads();
    auto integr_ptr = AnyODE::make_unique<Integrator>(lmm, iter_type);
    auto &integr = *integr_ptr;
    integr.autonomous_exprs = odesys->autonomous_exprs;
    integr.record_order = odesys->record_order;
    integr.record_fpe = odesys->record_fpe;
    integr.record_steps = odesys->record_steps;
    integr.set_user_data(static_cast<void *>(odesys));
    integr.init(rhs_cb<OdeSys>, t0, yq0, ny);
    if (nroots > 0)
        integr.root_init(nroots, roots_cb<OdeSys>);
    if (nq > 0){
        sundials_cxx::nvector_serial::VectorView q0(nq, const_cast<realtype*>(yq0)+ny);
        integr.quad_init(quads_cb<OdeSys>, q0.n_vec);
        if (atol.size() == (size_t)(ny + nq)){
            sundials_cxx::nvector_serial::VectorView quad_atol(nq, atol.data()+ny);
            integr.set_quad_err_con(true);
            integr.set_quad_tol(rtol, quad_atol.n_vec);
        } else if (atol.size() == 1) {
            integr.set_quad_err_con(true);
            integr.set_quad_tol(rtol, atol[0]);
        }
    }
    if (atol.size() == 1){
        integr.set_tol(rtol, atol[0]);
    } else if (atol.size() == (size_t)ny) {
        integr.set_tol(rtol, atol);
    } else if (atol.size() == (size_t)(ny+nq)) {
        sundials_cxx::nvector_serial::VectorView atol_(ny, atol.data());
        integr.set_tol(rtol, atol_.n_vec);
    } else {
        throw std::runtime_error("atol of incorrect length");
    }
    if (constraints.size())
        integr.set_constraints(constraints);
    integr.set_init_step(dx0);
    if (dx_min != 0.0)
        integr.set_min_step(dx_min);
    if (dx_max != 0.0)
        integr.set_max_step(dx_max);
    if (mxsteps)
        integr.set_max_num_steps(mxsteps);

    char * anyode_verbosity = std::getenv("ANYODE_VERBOSITY");
    integr.verbosity = (anyode_verbosity) ? std::atoi(anyode_verbosity) : 50;
    if (integr.verbosity == 0)
        integr.set_err_file_path("/dev/null"); // supposedly "NUL:" on windows

    if (iter_type == IterType::Newton){
        // Newton iteration --> we need a linear solver:
        switch(linear_solver){
        case 1:
            integr.set_linear_solver_to_dense(ny);
            if (with_jacobian)
                integr.set_dense_jac_fn(jac_dense_cb<OdeSys >);
            break;
        case 2:
            integr.set_linear_solver_to_banded(ny, odesys->get_mupper(), odesys->get_mlower());
            if (with_jacobian)
                integr.set_band_jac_fn(jac_band_cb<OdeSys >);
            break;
        case 10:
        case 11:
            integr.set_linear_solver_to_iterative(IterLinSolEnum::GMRES, maxl); break;
        case 20:
            integr.set_linear_solver_to_iterative(IterLinSolEnum::BICGSTAB, maxl); break;
        case 30:
            integr.set_linear_solver_to_iterative(IterLinSolEnum::TFQMR, maxl); break;
        default:
            throw std::runtime_error("Invalid linear_solver");
        }
        if (linear_solver >= 10){
            integr.set_prec_type(PrecType::Left);
            integr.set_iter_eps_lin(eps_lin);
            if (with_jtimes)
                integr.set_jac_times_vec_fn(jac_times_vec_cb<OdeSys>);
            integr.set_preconditioner(prec_setup_cb<OdeSys>, prec_solve_cb<OdeSys>);
            if (linear_solver == 10 || linear_solver == 11) // GMRES
                integr.set_gram_schmidt_type((linear_solver == 10) ? GramSchmidtType::Modified : GramSchmidtType::Classical);
            else if (linear_solver == 20 || linear_solver == 30) // BiCGStab, TFQMR
                ;
            else
                throw std::runtime_error("Unknown linear_solver.");
        }
    }
    return integr_ptr;
}

template <class OdeSys>
int
simple_adaptive(realtype ** xyqout,
                int * td,  // trailing dimension of xyqout ( == len(x) )
                OdeSys * const odesys,
                std::vector<realtype> atol,
                const realtype rtol,
                const LMM lmm,
                const realtype xend,
                std::vector<int>& root_indices,
                const long int mxsteps=0,
                realtype dx0=0.0,
                const realtype dx_min=0.0,
                const realtype dx_max=0.0,
                const bool with_jacobian=false,
                IterType iter_type=IterType::Undecided,
                int linear_solver=0,
                const int maxl=0,
                const realtype eps_lin=0.0,
                const unsigned nderiv=0,
                bool return_on_root=false,
                int autorestart=0,
                bool return_on_error=false,
                bool with_jtimes=false,
                int tidx=0,
                realtype ** ew_ele=nullptr,
                const std::vector<realtype> &constraints={}
    ){
    // iter_type == Undecided => Functional if lmm == Adams else Newton

    // linear_solver ==  0 => 1 if get_mlower() == -1 else 2
    // linear_solver ==  1 => Direct (dense LU-factorization)
    // linear_solver ==  2 => Direct (banded LU-factorization)
    // linear_solver == 10 => Iterative, GMRES, Modified Gram-Schmidt
    // linear_solver == 11 => Iterative, GMRES, Classical Gram-Schmidt
    // linear_solver == 20 => Iterative, Bi-CGStab (maxl => maximum dimension of Krylov subspace)
    // linear_solver == 30 => Iterative, TFQMR (maxl => maximum dimension of Krylov subspace)
    if (iter_type == IterType::Undecided)
        iter_type = (lmm == LMM::Adams) ? IterType::Functional : IterType::Newton;
    if (linear_solver == 0)
        linear_solver = (odesys->get_mlower() == -1) ? 1 : 2;
    realtype x0 = (*xyqout)[0];
    realtype * y0 = (*xyqout) + 1;
    if (dx0 == 0.0)
        dx0 = odesys->get_dx0(x0, y0);
    auto integr = get_integrator<OdeSys>(
        odesys, atol, rtol, lmm, y0, x0, mxsteps, dx0, dx_min, dx_max,
        with_jacobian, iter_type, linear_solver, maxl, eps_lin, with_jtimes, constraints);

    odesys->integrator = static_cast<void*>(integr.get());

    odesys->current_info.clear();
    if (odesys->record_rhs_xvals)
        odesys->current_info.nfo_vecdbl["rhs_xvals"] = {};
    if (odesys->record_jac_xvals)
        odesys->current_info.nfo_vecdbl["jac_xvals"] = {};

    std::time_t cput0 = std::clock();
    auto t_start = std::chrono::high_resolution_clock::now();

    int result = integr->adaptive(
        xyqout, td, xend, nderiv, root_indices, return_on_root, autorestart,
        return_on_error, (
            (odesys->use_get_dx_max) ?
            static_cast<cvodes_cxx::get_dx_max_fn>(std::bind(
                                                       &OdeSys::get_dx_max, odesys,
                                                       std::placeholders::_1,
                                                       std::placeholders::_2))
            :
            cvodes_cxx::get_dx_max_fn()
            ), tidx, ew_ele);

    odesys->current_info.nfo_dbl["time_cpu"] = (std::clock() - cput0) / (double)CLOCKS_PER_SEC;
    odesys->current_info.nfo_dbl["time_wall"] = std::chrono::duration<double>(
        std::chrono::high_resolution_clock::now() - t_start).count();
    cvodes_cxx::update_integration_info(
        odesys->current_info.nfo_int,
        odesys->current_info.nfo_dbl,
        odesys->current_info.nfo_vecdbl,
        odesys->current_info.nfo_vecint,
        *integr, iter_type, linear_solver);
    odesys->current_info.nfo_int["nfev"] = odesys->nfev;
    odesys->current_info.nfo_int["njev"] = odesys->njev;
    return result;
}

template <class OdeSys>
int simple_predefined(OdeSys * const odesys,
                      std::vector<realtype> atol,
                      const realtype rtol,
                      const LMM lmm,
                      const realtype * const yq0,
                      const std::size_t nout,
                      const realtype * const xout,
                      realtype * const yqout,
                      std::vector<int>& root_indices,
                      std::vector<double>& root_out,
                      const long int mxsteps=0,
                      realtype dx0=0.0,
                      const realtype dx_min=0.0,
                      const realtype dx_max=0.0,
                      const bool with_jacobian=false,
                      IterType iter_type=IterType::Undecided,
                      int linear_solver=0,
                      const int maxl=0,
                      const realtype eps_lin=0.0,
                      const unsigned nderiv=0,
                      int autorestart=0,
                      bool return_on_error=false,
                      bool with_jtimes=false,
                      realtype * ew_ele=nullptr,
                      const std::vector<realtype> &constraints={}
    ){
    // iter_type == Undecided => Functional if lmm == Adams else Newton

    // linear_solver ==  0 => 1 if get_mlower() == -1 else 2
    // linear_solver ==  1 => Direct (dense LU-factorization)
    // linear_solver ==  2 => Direct (banded LU-factorization)
    // linear_solver == 10 => Iterative, GMRES, Modified Gram-Schmidt
    // linear_solver == 11 => Iterative, GMRES, Classical Gram-Schmidt
    // linear_solver == 20 => Iterative, Bi-CGStab (maxl => maximum dimension of Krylov subspace)
    // linear_solver == 30 => Iterative, TFQMR (maxl => maximum dimension of Krylov subspace)
    if (iter_type == IterType::Undecided)
        iter_type = (lmm == LMM::Adams) ? IterType::Functional : IterType::Newton;
    if (linear_solver == 0)
        linear_solver = (odesys->get_mlower() == -1) ? 1 : 2;
    if (dx0 == 0.0)
        dx0 = odesys->get_dx0(xout[0], yq0);
    auto integr = get_integrator<OdeSys>(odesys, atol, rtol, lmm, yq0, xout[0], mxsteps, dx0, dx_min, dx_max,
                                         with_jacobian, iter_type, linear_solver, maxl, eps_lin, with_jtimes, constraints);
    odesys->integrator = static_cast<void*>(integr.get());

    odesys->current_info.clear();
    if (odesys->record_rhs_xvals)
        odesys->current_info.nfo_vecdbl["rhs_xvals"] = {};
    if (odesys->record_jac_xvals)
        odesys->current_info.nfo_vecdbl["jac_xvals"] = {};

    std::time_t cput0 = std::clock();
    auto t_start = std::chrono::high_resolution_clock::now();

    auto nreached = integr->predefined(
        nout, xout, yq0, yqout, nderiv, root_indices, root_out, autorestart, return_on_error,
        ((odesys->use_get_dx_max) ? static_cast<cvodes_cxx::get_dx_max_fn>(std::bind(&OdeSys::get_dx_max, odesys, std::placeholders::_1 , std::placeholders::_2))
         : cvodes_cxx::get_dx_max_fn()), ew_ele);

    odesys->current_info.nfo_dbl["time_cpu"] = (std::clock() - cput0) / (double)CLOCKS_PER_SEC;
    odesys->current_info.nfo_dbl["time_wall"] = std::chrono::duration<double>(
        std::chrono::high_resolution_clock::now() - t_start).count();
    if (odesys->record_order)
        odesys->current_info.nfo_vecint["orders"] = integr->orders_seen;
    if (odesys->record_fpe)
        odesys->current_info.nfo_vecint["fpes"] = integr->fpes_seen;

    cvodes_cxx::update_integration_info(
        odesys->current_info.nfo_int,
        odesys->current_info.nfo_dbl,
        odesys->current_info.nfo_vecdbl,
        odesys->current_info.nfo_vecint,
        *integr, iter_type, linear_solver);
    odesys->current_info.nfo_int["nfev"] = odesys->nfev;
    odesys->current_info.nfo_int["njev"] = odesys->njev;
    return nreached;
}


struct SolverSettings{
    double rtol {1e-8};
    std::vector<double> atol {1e-8};
    int mxsteps {0};
    double dx0{0}, dx_min{0}, dx_max{0};
    std::string method {"bdf"}, iter_type {"undecided"};
    int linear_solver {0};
    int maxl {0};
    realtype eps_lin {0.0};
    unsigned nderiv {0};
    bool return_on_root {false};
    int autorestart {0};
    bool return_on_error {true};
    bool with_jacobian {true};
    bool with_jtimes {true};
    std::vector<realtype> constraints={};
};

template <class OdeSys>
void check_atol(const OdeSys * odesys, const SolverSettings& settings){
    if (settings.atol.size() != 1 && static_cast<int>(settings.atol.size()) != (odesys->ny + odesys->nquads))
        throw std::logic_error("atol of incorrect length");
}

template <class OdeSys>
std::unique_ptr<AnyODE::Result> chained_predefined(
    OdeSys * const odesys,
    const std::vector<double> durations,
    const std::vector<double> &yq0,
    const std::vector<double> &varied_values,
    const std::vector<int> &varied_indices,
    int npoints,
    const SolverSettings &settings)
{
    LMM lmm = cvodes_cxx::lmm_from_name(settings.method);
    IterType iter_type = cvodes_cxx::iter_type_from_name(settings.iter_type);
    if (iter_type == IterType::Undecided)
        iter_type = (lmm == LMM::Adams) ? IterType::Functional : IterType::Newton;
    int linear_solver = settings.linear_solver;
    if (linear_solver == 0)
        linear_solver = (odesys->get_mlower() == -1) ? 1 : 2;
    double dx0 = settings.dx0;
    double x0 = 0.0;
    if (dx0 == 0.0)
        dx0 = odesys->get_dx0(x0, yq0);
    auto integr = get_integrator<OdeSys>(
        odesys, settings.atol, settings.rtol, lmm, yq0, x0, settings.mxsteps, dx0,
        settings.dx_min, settings.dx_max,
        settings.with_jacobian, iter_type, linear_solver, settings.maxl, settings.eps_lin,
        settings.with_jtimes, settings.constraints);
    odesys->integrator = static_cast<void*>(integr.get());

    std::time_t cput0 = std::clock();
    auto t_start = std::chrono::high_resolution_clock::now();

    auto ny = odesys->ny;
    auto nq = odesys->nquads;
    if (static_cast<int>(yq0.size()) != (ny + nq)) {
        throw std::logic_error("Incorrect length");
    }
    check_atol(odesys, settings);
    if (durations.size() != varied_values.size())
        throw std::logic_error("durations and varied_values have different lengths");
    int ndur = durations.size();
    double * xyqout = (double *)std::malloc((ndur*npoints+1)*(1+ny+nq)*sizeof(double));
    xyqout[0] = 0;
    for (int i=0; i < ny + nq; ++i){
        xyqout[i+1] = yq0[i];
    }
    for (int i=0; i < ndur; ++i){
        for (int j=0; j < npoints; ++j)
            xyqout[i*(1+ny+nq)*npoints + (j+1)*(1+ny+nq)] = xyqout[i*(1+ny+nq)*npoints] + (j+1)*durations[i]/npoints;
    }
    std::vector<double> tbuffer(npoints+1, 0.0), yqbuffer((npoints+1)*(ny+nq));
    std::vector<int> root_indices;
    std::vector<double> roots_out;
    bool success = true;
    AnyODE::Info info;
    int ntot = 1;
    for (int i=0; i<ndur; ++i){
#if defined(PYCVODES_VERBOSE)
        printf("\rProgress in %s: %.3f %%", __FUNCTION__, i*100.0/ndur);
        fflush(stdout);
#endif
        for (int j=0; j<npoints; ++j){
            tbuffer[j+1] = (j+1)*durations[i]/npoints;
        }
        for (int j=0; j<ny+nq; ++j){
            yqbuffer[j] = xyqout[i*npoints*(ny+nq+1) + 1 + j];
        }
        for (int j=0; j < static_cast<int>(varied_indices.size()); ++j){
            odesys->set_param_value(varied_indices[j], varied_values[i*varied_indices.size() + j]);
        }
        int nreached = cvodes_anyode::simple_predefined(
            odesys, settings.atol, settings.rtol, cvodes_cxx::lmm_from_name(settings.method),
            yqbuffer.data(), npoints+1, tbuffer.data(), yqbuffer.data(), root_indices, roots_out,
            settings.mxsteps, settings.dx0, settings.dx_min, settings.dx_max, settings.with_jacobian,
            cvodes_cxx::iter_type_from_name(settings.iter_type), settings.linear_solver, settings.maxl,
            settings.eps_lin, settings.nderiv, settings.autorestart, settings.return_on_error, settings.with_jtimes);
        if (nreached != static_cast<int>(npoints+1)){
            success = false;
            break;
        } else {
            ntot += nreached - 1;
        }
        info.update(odesys->current_info.nfo_int,
                    odesys->current_info.nfo_dbl,
                    odesys->current_info.nfo_vecdbl,
                    odesys->current_info.nfo_vecint);
        for (int j=0; j<npoints; ++j){
            xyqout[i*npoints*(1+ny+nq) + (j+1)*(1+ny+nq)] = xyqout[i*npoints*(1+ny+nq)] + tbuffer[j+1];
            for (int k=0; k<ny+nq; ++k)
                xyqout[(i*npoints+1)*(1+ny+nq) + j*(1+ny+nq) + k + 1] = yqbuffer[(j+1)*(ny+nq) + k];
        }
    }

    odesys->last_integration_info_dbl["time_cpu"] = (std::clock() - cput0) / (double)CLOCKS_PER_SEC;
    odesys->last_integration_info_dbl["time_wall"] = std::chrono::duration<double>(
        std::chrono::high_resolution_clock::now() - t_start).count();
    if (odesys->record_order)
        odesys->last_integration_info_vecint["orders"] = integr->orders_seen;
    if (odesys->record_fpe)
        odesys->last_integration_info_vecint["fpes"] = integr->fpes_seen;

    cvodes_cxx::update_integration_info(
        odesys->last_integration_info,
        odesys->last_integration_info_dbl,
        odesys->last_integration_info_vecdbl,
        odesys->last_integration_info_vecint,
        *integr, iter_type, linear_solver);
    odesys->current_info.nfo_int["nfev"] = odesys->nfev;
    odesys->current_info.nfo_int["njev"] = odesys->njev;

    auto result = AnyODE::make_unique<AnyODE::Result>(
        ntot, odesys->ny, odesys->nquads, odesys->nroots, xyqout);
    info.nfo_int["success"] = success;
    result->info = info;
    return result;
}

END_NAMESPACE(cvodes_anyode)
