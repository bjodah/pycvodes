#pragma once

#include <chrono>
#include "anyode/anyode.hpp"
#include "cvodes_cxx.hpp"

namespace cvodes_anyode {

    using cvodes_cxx::CVodeIntegrator;
    using cvodes_cxx::LMM;
    using cvodes_cxx::IterType;
    using cvodes_cxx::IterLinSolEnum;
    using cvodes_cxx::PrecType;
    using cvodes_cxx::GramSchmidtType;

    int handle_status_(AnyODE::Status status){
        switch (status){
        case AnyODE::Status::success:
            return 0;
        case AnyODE::Status::recoverable_error:
            return 1;
        case AnyODE::Status::unrecoverable_error:
            return -1;
        default:
            throw std::runtime_error("impossible (this is for silencing -Wreturn-type)");
        }
    }

    template<class OdeSys>
    int rhs_cb(realtype t, N_Vector y, N_Vector ydot, void *user_data){
        auto& odesys = *static_cast<OdeSys*>(user_data);
        AnyODE::Status status = odesys.rhs(t, NV_DATA_S(y), NV_DATA_S(ydot));
        return handle_status_(status);
    }

    template<class OdeSys>
    int roots_cb(realtype t, N_Vector y, realtype *gout, void *user_data){
        auto& odesys = *static_cast<OdeSys*>(user_data);
        AnyODE::Status status = odesys.roots(t, NV_DATA_S(y), gout);
        if (status == AnyODE::Status::recoverable_error)
            throw std::runtime_error("There are only unrecoverable errors for roots().");
        return handle_status_(status);
    }

    template <class OdeSys>
    int jac_dense_cb(long int N, realtype t,
                     N_Vector y, N_Vector fy, DlsMat Jac, void *user_data,
                     N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
        // callback of req. signature wrapping OdeSys method.
        AnyODE::ignore(N); AnyODE::ignore(tmp1); AnyODE::ignore(tmp2); AnyODE::ignore(tmp3);
        auto& odesys = *static_cast<OdeSys*>(user_data);
        AnyODE::Status status = odesys.dense_jac_cmaj(t, NV_DATA_S(y), NV_DATA_S(fy), DENSE_COL(Jac, 0), Jac->ldim);
        return handle_status_(status);
    }

    template <typename OdeSys>
    int jac_band_cb(long int N, long int mupper, long int mlower, realtype t,
                    N_Vector y, N_Vector fy, DlsMat Jac, void *user_data,
                    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
        AnyODE::ignore(N); AnyODE::ignore(tmp1); AnyODE::ignore(tmp2); AnyODE::ignore(tmp3);
        auto& odesys = *static_cast<OdeSys*>(user_data);
        if (odesys.get_mupper() != mupper)
            throw std::runtime_error("mupper mismatch");
        if (odesys.get_mlower() != mlower)
            throw std::runtime_error("mlower mismatch");
        AnyODE::Status status = odesys.banded_jac_cmaj(t, NV_DATA_S(y), NV_DATA_S(fy), Jac->data + Jac->s_mu - Jac->mu, Jac->ldim);
        return handle_status_(status);
    }


    template <typename OdeSys>
    int jac_times_vec_cb(N_Vector v, N_Vector Jv, realtype t, N_Vector y,
                         N_Vector fy, void *user_data, N_Vector tmp){
        // callback of req. signature wrapping OdeSys method.
        AnyODE::ignore(tmp);
        auto& odesys = *static_cast<OdeSys*>(user_data);
        AnyODE::Status status = odesys.jac_times_vec(NV_DATA_S(v), NV_DATA_S(Jv), t, NV_DATA_S(y), NV_DATA_S(fy));
        if (status == AnyODE::Status::recoverable_error)
            throw std::runtime_error("There are only unrecoverable errors for JacTimesVec().");
        return handle_status_(status);
    }

    template <typename OdeSys> // Section 4.6.9 Preconditioning in cvs_guide.pdf
    int jac_prec_solve_cb(realtype t, N_Vector y, N_Vector fy, N_Vector r,
                          N_Vector z, realtype gamma, realtype delta, int lr,
                          void *user_data, N_Vector tmp){
        // callback of req. signature wrapping OdeSys method.
        AnyODE::ignore(tmp);  // delta used for iterative methods
        double * ewt {nullptr};
        auto& odesys = *static_cast<OdeSys*>(user_data);
        if (lr != 1)
            throw std::runtime_error("Only left preconditioning implemented.");
        AnyODE::Status status =  odesys.prec_solve_left(t, NV_DATA_S(y), NV_DATA_S(fy), NV_DATA_S(r),
                                                        NV_DATA_S(z), gamma, delta, ewt);
        return handle_status_(status);
    }

    template <typename OdeSys>
    int prec_setup_cb(realtype t, N_Vector y, N_Vector fy, booleantype jok,
                      booleantype *jcurPtr, realtype gamma, void *user_data,
                      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
        // callback of req. signature wrapping OdeSys method.
        AnyODE::ignore(tmp1); AnyODE::ignore(tmp2); AnyODE::ignore(tmp3);
        auto& odesys = *static_cast<OdeSys*>(user_data);
        bool jac_recomputed = false;
        bool compute_jac = (jok == TRUE) ? false : true; // TRUE, FALSE sundials macros
        AnyODE::Status status = odesys.prec_setup(t, NV_DATA_S(y), NV_DATA_S(fy), compute_jac, jac_recomputed, gamma);
        (*jcurPtr) = (jac_recomputed) ? TRUE : FALSE;
        return handle_status_(status);
    }

    template <class OdeSys>
    CVodeIntegrator get_integrator(OdeSys * odesys,
                                   const std::vector<realtype> atol,
                                   const realtype rtol,
                                   const LMM lmm,  // Cython does not supported enum class...
                                   const realtype * const y0,
                                   const realtype t0,
                                   const long int mxsteps=0,
                                   const realtype dx0=0.0,
                                   const realtype dx_min=0.0,
                                   const realtype dx_max=0.0,
                                   const bool with_jacobian=false,
                                   const IterType iter_type=IterType::Newton,
                                   const int linear_solver=0,
                                   const int maxl=0,
                                   const realtype eps_lin=0.0
                                   )
    {
        const int ny = odesys->get_ny();
        CVodeIntegrator integr {lmm, iter_type};
        integr.set_user_data(static_cast<void *>(odesys));
        integr.init(rhs_cb<OdeSys>, t0, y0, ny);
        integr.root_init(odesys->get_nroots(), roots_cb<OdeSys>);
        if (atol.size() == 1){
            integr.set_tol(rtol, atol[0]);
        }else{
            integr.set_tol(rtol, atol);
        }
        integr.set_init_step(dx0);
        if (dx_min != 0.0)
            integr.set_min_step(dx_min);
        if (dx_max != 0.0)
            integr.set_max_step(dx_max);
        if (mxsteps)
            integr.set_max_num_steps(mxsteps);


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
                if (!with_jacobian)
                    throw std::runtime_error("Iterative method requires an (approximate) jacobian");
                integr.set_prec_type(PrecType::Left);
                integr.set_iter_eps_lin(eps_lin);
                integr.set_jac_times_vec_fn(jac_times_vec_cb<OdeSys>);
                integr.set_preconditioner(prec_setup_cb<OdeSys>, jac_prec_solve_cb<OdeSys>);
                if (linear_solver == 10 || linear_solver == 11) // GMRES
                    integr.set_gram_schmidt_type((linear_solver == 10) ? GramSchmidtType::Modified : GramSchmidtType::Classical);
                else if (linear_solver == 20 or linear_solver == 30) // BiCGStab, TFQMR
                    ;
                else
                    throw std::runtime_error("Unknown linear_solver.");
            }
        }
        return integr;
    }

    template <class OdeSys>
    std::pair<std::vector<realtype>, std::vector<realtype> >
    simple_adaptive(OdeSys * const odesys,
                    const std::vector<realtype> atol,
                    const realtype rtol,
                    const LMM lmm,
                    const realtype * const y0,
                    const realtype t0,
                    const realtype tend,
                    std::vector<int>& root_indices,
                    const long int mxsteps=0,
                    const realtype dx0=0.0,
                    const realtype dx_min=0.0,
                    const realtype dx_max=0.0,
                    const bool with_jacobian=false,
                    IterType iter_type=IterType::Undecided,
                    int linear_solver=0,
                    const int maxl=0,
                    const realtype eps_lin=0.0,
                    const unsigned nderiv=0,
                    bool return_on_root=false
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

        auto integr = get_integrator<OdeSys>(odesys, atol, rtol, lmm, y0, t0, mxsteps, dx0, dx_min, dx_max,
                                             with_jacobian, iter_type, linear_solver, maxl, eps_lin);
        odesys->integrator = static_cast<void*>(&integr);
        std::time_t cput0 = std::clock();
        auto t_start = std::chrono::high_resolution_clock::now();

        auto result = integr.adaptive(t0, tend, y0, nderiv, root_indices, return_on_root);

        odesys->last_integration_info_dbl["time_cpu"] = (std::clock() - cput0) / (double)CLOCKS_PER_SEC;
        odesys->last_integration_info_dbl["time_wall"] = std::chrono::duration<double>(
                std::chrono::high_resolution_clock::now() - t_start).count();

        odesys->last_integration_info.clear();
        cvodes_cxx::set_integration_info(odesys->last_integration_info, integr,
                                         iter_type, linear_solver);
        odesys->last_integration_info["nfev"] = odesys->nfev;
        odesys->last_integration_info["njev"] = odesys->njev;
        return result;
    }

    template <class OdeSys>
    void simple_predefined(OdeSys * const odesys,
                           const std::vector<realtype> atol,
                           const realtype rtol,
                           const LMM lmm,
                           const realtype * const y0,
                           const std::size_t nout,
                           const realtype * const tout,
                           realtype * const yout,
                           std::vector<int>& root_indices,
                           std::vector<double>& root_out,
                           const long int mxsteps=0,
                           const realtype dx0=0.0,
                           const realtype dx_min=0.0,
                           const realtype dx_max=0.0,
                           const bool with_jacobian=false,
                           IterType iter_type=IterType::Undecided,
                           int linear_solver=0,
                           const int maxl=0,
                           const realtype eps_lin=0.0,
                           const unsigned nderiv=0
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
        auto integr = get_integrator<OdeSys>(odesys, atol, rtol, lmm, y0, tout[0], mxsteps, dx0, dx_min, dx_max,
                                             with_jacobian, iter_type, linear_solver, maxl, eps_lin);
        odesys->integrator = static_cast<void*>(&integr);

        std::time_t cput0 = std::clock();
        auto t_start = std::chrono::high_resolution_clock::now();

        integr.predefined(nout, tout, y0, yout, nderiv, root_indices, root_out);

        odesys->last_integration_info_dbl["time_cpu"] = (std::clock() - cput0) / (double)CLOCKS_PER_SEC;
        odesys->last_integration_info_dbl["time_wall"] = std::chrono::duration<double>(
                std::chrono::high_resolution_clock::now() - t_start).count();

        odesys->last_integration_info.clear();
        cvodes_cxx::set_integration_info(odesys->last_integration_info, integr,
                                         iter_type, linear_solver);
        odesys->last_integration_info["nfev"] = odesys->nfev;
        odesys->last_integration_info["njev"] = odesys->njev;
    }
}
