#pragma once
// Thin C++11 wrapper around CVODES v2.8.2 (SUNDIALS v2.6.2)
// far from all functionality is available yet.
// sundials-2.6.2.tar.gz (MD5: 3deeb0ede9f514184c6bd83ecab77d95)

#include <assert.h>
#include <cmath>
#include <cstring>
#include <memory>
#include <utility>
#include <vector>

#include "sundials_cxx.hpp" // sundials_cxx::nvector_serial::Vector
#include <cvodes/cvodes_spils.h>
#include <cvodes/cvodes_spgmr.h>
#include <cvodes/cvodes_spbcgs.h>
#include <cvodes/cvodes_sptfqmr.h>
#include <cvodes/cvodes.h> /* CVODE fcts., CV_BDF, CV_ADAMS */
#include <cvodes/cvodes_impl.h> /* CVodeMem */
#include <cvodes/cvodes_lapack.h>       /* prototype for CVDense */

// Overflows easily, check your architecture
int factorial(int n) {
    int res = 1;
    while (n > 1) {
        res *= n--;
    }
    return res;
}

namespace cvodes_cxx {

    //using sundials_cxx::nvector_serial::N_Vector; // native sundials vector
    using SVector = sundials_cxx::nvector_serial::Vector; // serial vector

    // Wrapper for Revision 4306 of cvodes.c

    enum class LMM : int {ADAMS=CV_ADAMS, BDF=CV_BDF}; // Linear multistep method

    enum class Task : int {NORMAL=CV_NORMAL, ONE_STEP=CV_ONE_STEP};
    enum class IterType : int {FUNCTIONAL=CV_FUNCTIONAL, NEWTON=CV_NEWTON};
    enum class IterLinSolEnum : int {GMRES=1, BICGSTAB=2, TFQMR=3};
    enum class PrecType : int {NONE=PREC_NONE, LEFT=PREC_LEFT,
            RIGHT=PREC_RIGHT, BOTH=PREC_BOTH};
    enum class GramSchmidtType : int {CLASSICAL=CLASSICAL_GS, MODIFIED=MODIFIED_GS};

    void check_flag(int flag) {
        switch (flag){
        case CV_SUCCESS:
            break;
        case CV_MEM_NULL:
            throw std::runtime_error("cvode_mem is NULL");
        default:
            throw std::runtime_error("Unhandled flag");
        }
    }

    class Integrator{ // Thin wrapper class of CVode in CVODES
    public:
        void *mem {nullptr};
        Integrator(const LMM lmm, const IterType iter) {
            this->mem = CVodeCreate(static_cast<int>(lmm), static_cast<int>(iter));
        }
        // init
        void init(CVRhsFn cb, realtype t0, N_Vector y) {
            int status = CVodeInit(this->mem, cb, t0, y);
            if (status == CV_ILL_INPUT)
                throw std::runtime_error("CVodeInit failed (CV_ILL_INPUT).");
            else if (status == CV_MEM_FAIL)
                throw std::runtime_error("CVodeInit failed (allocation failed).");
            else
                check_flag(status);
        }
        void init(CVRhsFn cb, realtype t0, SVector &y) {
            this->init(cb, t0, y.n_vec);
        }
        void init(CVRhsFn cb, realtype t0, const realtype * const y, int ny) {
            SVector yvec (ny, const_cast<realtype*>(y));
            this->init(cb, t0, yvec.n_vec);
            // it is ok that yvec is destructed here
            // (see CVodeInit in cvodes.c which at L848 calls cvAllocVectors (L3790))
        }

        // reinit
        void reinit(realtype t0, N_Vector y){
            CVodeReInit(this->mem, t0, y);
        }
        void reinit(realtype t0, SVector &y){
            reinit(t0, y.n_vec);
        }
        void reinit(realtype t0, const realtype * const y, int ny){
            SVector yvec (ny, const_cast<realtype*>(y));
            reinit(t0, yvec.n_vec);
        }
        // Root finding
        void root_init(const int nrtfn, CVRootFn g){
            if (!nrtfn)
                return;
            int status = CVodeRootInit(this->mem, nrtfn, g);
            if (status == CV_ILL_INPUT)
                throw std::runtime_error("CVodeRootInit failed (CV_ILL_INPUT).");
            else if (status == CV_MEM_FAIL)
                throw std::runtime_error("CVodeRootInit failed (allocation failed).");
            else
                check_flag(status);
        }


        // Step specifications
        void set_init_step(realtype h0){
            int status = CVodeSetInitStep(this->mem, h0);
            check_flag(status);
        }
        void set_min_step(realtype hmin){
            int status = CVodeSetMinStep(this->mem, hmin);
            if (status == CV_ILL_INPUT)
                throw std::runtime_error("hmin non-positive or exceeding maximum allowable step size");
            else
                check_flag(status);
        }
        void set_max_step(realtype hmax){
            int status = CVodeSetMaxStep(this->mem, hmax);
            if (status == CV_ILL_INPUT)
                throw std::runtime_error("hmax non-positive or smaller than minimumem allowable step size");
            else
                check_flag(status);
        }
        void set_max_num_steps(long int mxsteps){
            int status = CVodeSetMaxNumSteps(this->mem, mxsteps);
            check_flag(status);
        }

        // set_tol
        void set_tol(realtype rtol, realtype atol){
            int status = CVodeSStolerances(this->mem, rtol, atol);
            if (status < 0)
                throw std::runtime_error("CVodeSStolerances failed.");
        }
        void set_tol(realtype rtol, N_Vector atol){
            int status = CVodeSVtolerances(this->mem, rtol, atol);
            if (status < 0)
                throw std::runtime_error("CVodeSVtolerances failed.");
        }
        void set_tol(realtype rtol, std::vector<realtype> atol){
            SVector atol_(atol.size(), &atol[0]);
            set_tol(rtol, atol_.n_vec);
        }

        // set stop time
        void set_stop_time(realtype tend){
            int status = CVodeSetStopTime(this->mem, tend);
            if (status == CV_ILL_INPUT)
                throw std::runtime_error("tend not beyond current t");
            else
                check_flag(status);
        }

        // user data
        void set_user_data(void *user_data){
            int status = CVodeSetUserData(this->mem, user_data);
            if (status < 0)
                throw std::runtime_error("CVodeSetUserData failed.");
        }

        // dense jacobian
        void set_linear_solver_to_dense(int ny){
            int status = CVLapackDense(this->mem, ny);
            if (status != CVDLS_SUCCESS)
                throw std::runtime_error("CVLapackDense failed");
        }
        void set_dense_jac_fn(CVDlsDenseJacFn djac){
            int status = CVDlsSetDenseJacFn(this->mem, djac);
            if (status < 0)
                throw std::runtime_error("CVDlsSetDenseJacFn failed.");
        }

        // banded jacobian
        void set_linear_solver_to_banded(int N, int mupper, int mlower){
            int status = CVLapackBand(this->mem, N, mupper, mlower);
            if (status != CVDLS_SUCCESS)
                throw std::runtime_error("CVLapackBand failed");
        }
        void set_band_jac_fn(CVDlsBandJacFn djac){
            int status = CVDlsSetBandJacFn(this->mem, djac);
            if (status < 0)
                throw std::runtime_error("CVDlsSetBandJacFn failed.");
        }

        // Iterative Linear solvers
        void set_linear_solver_to_iterative(IterLinSolEnum solver, int maxl=0){
            int flag;
            switch (solver) {
            case IterLinSolEnum::GMRES:
                flag = CVSpgmr(this->mem, static_cast<int>(PrecType::LEFT), maxl);
                break;
            case IterLinSolEnum::BICGSTAB:
                flag = CVSpbcg(this->mem, static_cast<int>(PrecType::LEFT), maxl);
                break;
            case IterLinSolEnum::TFQMR:
                flag = CVSptfqmr(this->mem, static_cast<int>(PrecType::LEFT), maxl);
                break;
            }
            switch (flag){
            case CVSPILS_SUCCESS:
                break;
            case CVSPILS_MEM_NULL:
                throw std::runtime_error("set_linear_solver_to_iterative failed (cvode_mem is NULL)");
            case CVSPILS_ILL_INPUT:
                throw std::runtime_error("PREC_LEFT invalid.");
            case CVSPILS_MEM_FAIL:
                throw std::runtime_error("Memory allocation request failed.");
            }
        }
        void cvspils_check_flag(int flag, bool check_ill_input=false) {
            switch (flag){
            case CVSPILS_SUCCESS:
                break;
            case CVSPILS_MEM_NULL:
                throw std::runtime_error("cvode_mem is NULL");
            case CVSPILS_LMEM_NULL:
                throw std::runtime_error("CVSPILS linear solver has not been initialized)");
            }
            if ((check_ill_input) && (flag == CVSPILS_ILL_INPUT))
                throw std::runtime_error("Bad input.");
        }
        void set_jac_times_vec_fn(CVSpilsJacTimesVecFn jac_times_vec_fn){
            int flag = CVSpilsSetJacTimesVecFn(this->mem, jac_times_vec_fn);
            this->cvspils_check_flag(flag);
        }
        void set_preconditioner(CVSpilsPrecSetupFn setup_fn, CVSpilsPrecSolveFn solve_fn){
            int flag = CVSpilsSetPreconditioner(this->mem, setup_fn, solve_fn);
            this->cvspils_check_flag(flag);
        }
        void set_iter_eps_lin(realtype delta){
            int flag = CVSpilsSetEpsLin(this->mem, delta);
            this->cvspils_check_flag(flag, true);
        }
        void set_prec_type(PrecType pretyp){
            int flag = CVSpilsSetPrecType(this->mem, (int)pretyp);
            this->cvspils_check_flag(flag, true);
        }
        void set_gram_schmidt_type(GramSchmidtType gs_type){
            int flag = CVSpilsSetGSType(this->mem, (int)gs_type);
            this->cvspils_check_flag(flag, true);
        }
        void set_krylov_max_len(int maxl){
            int flag = CVSpilsSetMaxl(this->mem, maxl);
            this->cvspils_check_flag(flag);
        }

        // get info
        long int get_n_lin_iters(){
            long int res=0;
            int flag;
            flag = CVSpilsGetNumLinIters(this->mem, &res);
            this->cvspils_check_flag(flag);
            return res;
        }

        long int get_n_prec_evals(){
            long int res=0;
            int flag;
            flag = CVSpilsGetNumPrecEvals(this->mem, &res);
            this->cvspils_check_flag(flag);
            return res;
        }

        long int get_n_prec_solves(){
            long int res=0;
            int flag;
            flag = CVSpilsGetNumPrecSolves(this->mem, &res);
            this->cvspils_check_flag(flag);
            return res;
        }

        long int get_n_conv_fails(){
            long int res=0;
            int flag;
            flag = CVSpilsGetNumConvFails(this->mem, &res);
            this->cvspils_check_flag(flag);
            return res;
        }

        long int get_n_jac_times_evals(){
            long int res=0;
            int flag;
            flag = CVSpilsGetNumJtimesEvals(this->mem, &res);
            this->cvspils_check_flag(flag);
            return res;
        }

        long int get_n_iter_rhs(){
            long int res=0;
            int flag;
            flag = CVSpilsGetNumRhsEvals(this->mem, &res);
            this->cvspils_check_flag(flag);
            return res;
        }

        long int get_n_steps(){
            long int res=0;
            int flag = CVodeGetNumSteps(this->mem, &res);
            check_flag(flag);
            return res;
        }

        long int get_n_rhs_evals(){
            long int res=0;
            int flag = CVodeGetNumRhsEvals(this->mem, &res);
            check_flag(flag);
            return res;
        }

        long int get_n_lin_solv_setups(){
            long int res=0;
            int flag = CVodeGetNumLinSolvSetups(this->mem, &res);
            check_flag(flag);
            return res;
        }

        long int get_n_err_test_fails(){
            long int res=0;
            int flag = CVodeGetNumErrTestFails(this->mem, &res);
            check_flag(flag);
            return res;
        }

        long int get_n_nonlin_solv_iters(){
            long int res=0;
            int flag = CVodeGetNumNonlinSolvIters(this->mem, &res);
            check_flag(flag);
            return res;
        }

        long int get_n_nonlin_solv_conv_fails(){
            long int res=0;
            int flag = CVodeGetNumNonlinSolvConvFails(this->mem, &res);
            check_flag(flag);
            return res;
        }

        void cvdls_check_flag(int flag) {
            switch (flag){
            case CVDLS_SUCCESS:
                break;
            case CVDLS_MEM_NULL:
                throw std::runtime_error("cvode_mem is NULL");
            case CVDLS_LMEM_NULL:
                throw std::runtime_error("CVDLS linear solver has not been initialized)");
            }
        }
        long int get_n_dls_jac_evals(){
            long int res=0;
            int flag = CVDlsGetNumJacEvals(this->mem, &res);
            cvdls_check_flag(flag);
            return res;
        }

        long int get_n_dls_rhs_evals(){
            long int res=0;
            int flag = CVDlsGetNumRhsEvals(this->mem, &res);
            cvdls_check_flag(flag);
            return res;
        }

        void get_dky(realtype t, int k, SVector &dky) {
            int flag = CVodeGetDky(this->mem, t, k, dky.n_vec);
            switch(flag){
            case CV_SUCCESS:
                // CVodeGetDky succeeded.
                break;
            case CV_BAD_K:
                throw std::runtime_error("CVodeGetDky failed with (invalid k)");
            case CV_BAD_T:
                throw std::runtime_error("CVodeGetDky failed with (invalid t)");
            case CV_BAD_DKY:
                throw std::runtime_error("CVodeGetDky failed with (dky.n_vec was NULL)");
            case CV_MEM_NULL:
                throw std::runtime_error("CVodeGetDky failed with (cvode_mem was NULL)");
            default:
                throw std::runtime_error("Undhandled flag in get_dky()");
            }
        }
        int step(realtype tout, SVector &yout, realtype *tret, Task task){
            return CVode(this->mem, tout, yout.n_vec, tret, static_cast<int>(task));
        }
        void call_rhs(realtype t, SVector y, SVector &ydot){
            CVodeMem cv_mem = (CVodeMem) this->mem;
            int status = cv_mem->cv_f(t, y.n_vec, ydot.n_vec, cv_mem->cv_user_data);
            if (status)
                throw std::runtime_error("call_rhs failed.");
        }

        std::pair<std::vector<realtype>, std::vector<realtype> >
        adaptive(long int ny, const realtype x0, const realtype xend,
                 const realtype * const y0, int nderiv, std::vector<int>& root_indices,
                 bool return_on_root=false){
            std::vector<realtype> xout;
            std::vector<realtype> yout;
            realtype cur_t;
            int status;
            int idx = 0;
            SVector y {ny};
            SVector work {ny};
            xout.push_back(x0);
            for (int i=0; i<ny; ++i){
                y[i] = y0[i];
                yout.push_back(y0[i]);
            }
            this->reinit(x0, y);
            if (nderiv >= 1){
                this->call_rhs(x0, y, work);
                for (int i=0; i<ny; ++i)
                    yout.push_back(work[i]);
            }
            for (int di=1; di<nderiv; ++di){
                for (int i=0; i<ny; ++i)  // higher order too expensive
                    yout.push_back(0);
            }
            this->set_stop_time(xend);
            do {
                idx++;
                status = this->step(xend, y, &cur_t, Task::ONE_STEP);
                if(status != CV_SUCCESS && status != CV_TSTOP_RETURN){
                    if (status == CV_ROOT_RETURN){
                        root_indices.push_back(idx);
                    }else{
                        throw std::runtime_error("Unsuccessful CVode step.");
                    }
                }
                xout.push_back(cur_t);
                for (int i=0; i<ny; ++i)
                    yout.push_back(y[i]);
                // Derivatives for interpolation
                for (int di=0; di<nderiv; ++di){
                    if (this->get_n_steps() < 2*(nderiv+1))
                        // Too few points collected
                        work.zero_out();
                    else
                        this->get_dky(cur_t, di+1, work);
                    for (int i=0; i<ny; ++i)
                        yout.push_back(work[i]);
                }
                if (return_on_root && status == CV_ROOT_RETURN)
                    break;
            } while (status != CV_TSTOP_RETURN);
            return std::pair<std::vector<realtype>, std::vector<realtype>>(xout, yout);
        }

        void predefined(int nt, int ny, const realtype * const tout, const realtype * const y0,
                        realtype * const yout, int nderiv, std::vector<int>& root_indices){
            realtype cur_t;
            int status;
            SVector y {ny};
            SVector work {ny};

            for (int i=0; i<ny; ++i)
                y[i] = y0[i];
            this->reinit(tout[0], y);
            y.dump(yout);
            if (nderiv >= 1){
                this->call_rhs(tout[0], y, work);
                work.dump(&yout[ny]);
            }
            for (int di=1; di<nderiv; ++di){
                for (int i=0; i<ny; ++i)  // too expensive
                    yout[ny*(di+1) + i] = 0;
            }

            for(int iout=1; (iout < nt); iout++) {
                status = this->step(tout[iout], y, &cur_t, Task::NORMAL);
                if(status != CV_SUCCESS){
                    if (status == CV_ROOT_RETURN){
                        root_indices.push_back(iout);
                        iout--;
                        continue;
                    }else{
                        throw std::runtime_error("Unsuccessful CVode step.");
                    }
                }
                y.dump(&yout[ny*(iout*(nderiv+1))]);
                // Derivatives for interpolation
                for (int di=0; di<nderiv; ++di){
                    if (this->get_n_steps() < 2*(nderiv+1))
                        // Too few points collected
                        work.zero_out();
                    else
                        this->get_dky(tout[iout], di+1, work);
                    work.dump(&yout[ny*(di+1+(iout*(nderiv+1)))]);
                }
            }
        }

        ~Integrator(){
            if (this->mem)
                CVodeFree(&(this->mem));
        }
    };

    template<class OdeSys>
    int f_cb(realtype t, N_Vector y, N_Vector ydot, void *user_data){
        auto& odesys = *static_cast<OdeSys*>(user_data);
        odesys.rhs(t, NV_DATA_S(y), NV_DATA_S(ydot));
        return 0;
    }

    template<class OdeSys>
    int roots_cb(realtype t, N_Vector y, realtype *gout, void *user_data){
        auto& odesys = *static_cast<OdeSys*>(user_data);
        odesys.roots(t, NV_DATA_S(y), gout);
        return 0;
    }

    template <class OdeSys>
    int jac_dense_cb(long int N, realtype t,
                     N_Vector y, N_Vector fy, DlsMat Jac, void *user_data,
                     N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
        // callback of req. signature wrapping OdeSys method.
        auto& odesys = *static_cast<OdeSys*>(user_data);
        odesys.dense_jac_cmaj(t, NV_DATA_S(y), NV_DATA_S(fy), DENSE_COL(Jac, 0),
                               Jac->ldim);
        return 0;
    }

    template <typename OdeSys>
    int jac_band_cb(long int N, long int mupper, long int mlower, realtype t,
                    N_Vector y, N_Vector fy, DlsMat Jac, void *user_data,
                    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
        // callback of req. signature wrapping OdeSys method.
        auto& odesys = *static_cast<OdeSys*>(user_data);
        if (odesys.get_mupper() != mupper)
            throw std::runtime_error("mupper mismatch");
        if (odesys.get_mlower() != mlower)
            throw std::runtime_error("mlower mismatch");
        odesys.banded_padded_jac_cmaj(t, NV_DATA_S(y), NV_DATA_S(fy), Jac->data, Jac->ldim);
        return 0;
    }


    template <typename OdeSys>
    int jac_times_vec_cb(N_Vector v, N_Vector Jv, realtype t, N_Vector y,
                         N_Vector fy, void *user_data, N_Vector tmp){
        // callback of req. signature wrapping OdeSys method.
        auto& odesys = *static_cast<OdeSys*>(user_data);
        odesys.jac_times_vec(NV_DATA_S(v), NV_DATA_S(Jv), t, NV_DATA_S(y), NV_DATA_S(fy));
        return 0;
    }

    template <typename OdeSys>
    int jac_prec_solve_cb(realtype t, N_Vector y, N_Vector fy, N_Vector r,
                          N_Vector z, realtype gamma, realtype delta, int lr,
                          void *user_data, N_Vector tmp){
        // callback of req. signature wrapping OdeSys method.
        auto& odesys = *static_cast<OdeSys*>(user_data);
        if (lr != 1)
            throw std::runtime_error("Only left preconditioning implemented.");
        odesys.prec_solve_left(t, NV_DATA_S(y), NV_DATA_S(fy), NV_DATA_S(r),
                                NV_DATA_S(z), gamma);
        return 0; // Direct solver give no hint on success, hence report success.
    }

    template <typename OdeSys>
    int prec_setup_cb(realtype t, N_Vector y, N_Vector fy, booleantype jok,
                      booleantype *jcurPtr, realtype gamma, void *user_data,
                      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
        // callback of req. signature wrapping OdeSys method.
        auto& odesys = *static_cast<OdeSys*>(user_data);
        bool jac_recomputed = false;
        bool compute_jac = (jok == TRUE) ? false : true; // TRUE, FALSE sundials macros
        odesys.prec_setup(t, NV_DATA_S(y), NV_DATA_S(fy), compute_jac, jac_recomputed, gamma);
        (*jcurPtr) = (jac_recomputed) ? TRUE : FALSE;
        return 0;
    }

    template <class OdeSys>
    Integrator get_integrator(OdeSys * odesys,
                              const std::vector<realtype> atol,
                              const realtype rtol, const int lmm,
                              const realtype * const y0,
                              const realtype t0,
                              const realtype dx0=0.0,
                              const realtype dx_min=0.0,
                              const realtype dx_max=0.0,
                              const int mxsteps=0,
                              const bool with_jacobian=false,
                              const int iter_type=0,
                              const int linear_solver=0,
                              const int maxl=0,
                              const realtype eps_lin=0.0
                              )
    {
        const int ny = odesys->ny;
        Integrator integr {(lmm == CV_BDF) ? LMM::BDF : LMM::ADAMS,
                (iter_type == 1) ? IterType::FUNCTIONAL : IterType::NEWTON};
        integr.set_user_data(static_cast<void *>(odesys));
        integr.init(f_cb<OdeSys>, t0, y0, ny);
        integr.root_init(odesys->nroots, roots_cb<OdeSys>);
        if (atol.size() == 1){
            integr.set_tol(rtol, atol[0]);
        }else{
            integr.set_tol(rtol, atol);
        }
        integr.set_init_step(dx0);
        if (dx_min != 0.0)
            integr.set_min_step(dx0);
        if (dx_max != 0.0)
            integr.set_max_step(dx0);
        if (mxsteps)
            integr.set_max_num_steps(mxsteps);


        if (iter_type == 2){
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
                integr.set_prec_type(PrecType::LEFT);
                integr.set_iter_eps_lin(eps_lin);
                integr.set_jac_times_vec_fn(jac_times_vec_cb<OdeSys>);
                integr.set_preconditioner(prec_setup_cb<OdeSys>, jac_prec_solve_cb<OdeSys>);
                if (linear_solver == 10 || linear_solver == 11) // GMRES
                    integr.set_gram_schmidt_type((linear_solver == 10) ? GramSchmidtType::MODIFIED : GramSchmidtType::CLASSICAL);
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
    simple_adaptive(OdeSys * odesys,
                    const std::vector<realtype> atol,
                    const realtype rtol,
                    const int lmm,
                    const realtype * const y0,
                    const realtype t0,
                    const realtype tend,
                    std::vector<int>& root_indices,
                    const realtype dx0=0.0,
                    const realtype dx_min=0.0,
                    const realtype dx_max=0.0,
                    const long int mxsteps=0,
                    const bool with_jacobian=false,
                    int iter_type=0,
                    int linear_solver=0,
                    const int maxl=0,
                    const realtype eps_lin=0.0,
                    const int nderiv=0,
                    bool return_on_root=false
                    ){
        // iter_type == 0 => 1 if lmm == CV_ADAMS else 2
        // iter_type == 1 => Functional (ignore linear_solver)
        // iter_type == 2 => Newton

        // linear_solver ==  0 => 1 if get_mlower() == -1 else 2
        // linear_solver ==  1 => Direct (dense LU-factorization)
        // linear_solver ==  2 => Direct (banded LU-factorization)
        // linear_solver == 10 => Iterative, GMRES, Modified Gram-Schmidt
        // linear_solver == 11 => Iterative, GMRES, Classical Gram-Schmidt
        // linear_solver == 20 => Iterative, Bi-CGStab (maxl => maximum dimension of Krylov subspace)
        // linear_solver == 30 => Iterative, TFQMR (maxl => maximum dimension of Krylov subspace)
        if (iter_type == 0)
            iter_type = (lmm == CV_ADAMS) ? 1 : 2;
        if (iter_type != 1 and iter_type != 2)
            throw std::runtime_error("Invalid iter_type");
        if (linear_solver == 0)
            linear_solver = (odesys->get_mlower() == -1) ? 1 : 2;

        auto integr = get_integrator<OdeSys>(odesys, atol, rtol, lmm, y0, t0, dx0, dx_min, dx_max, mxsteps,
                                             with_jacobian, iter_type, linear_solver, maxl, eps_lin);
        return integr.adaptive(odesys->ny, t0, tend, y0, nderiv, root_indices, return_on_root);
    }

    template <class OdeSys>
    void simple_predefined(OdeSys * odesys,
                           const std::vector<realtype> atol,
                           const realtype rtol,
                           const int lmm,
                           const realtype * const y0,
                           const std::size_t nout,
                           const realtype * const tout,
                           realtype * const yout,
                           std::vector<int>& root_indices,
                           const realtype dx0=0.0,
                           const realtype dx_min=0.0,
                           const realtype dx_max=0.0,
                           const long int mxsteps=0,
                           const bool with_jacobian=false,
                           int iter_type=0,
                           int linear_solver=0,
                           const int maxl=0,
                           const realtype eps_lin=0.0,
                           const int nderiv=0
                           ){
        // iter_type == 0 => 1 if lmm == CV_ADAMS else 2
        // iter_type == 1 => Functional (ignore linear_solver)
        // iter_type == 2 => Newton

        // linear_solver ==  0 => 1 if get_mlower() == -1 else 2
        // linear_solver ==  1 => Direct (dense LU-factorization)
        // linear_solver ==  2 => Direct (banded LU-factorization)
        // linear_solver == 10 => Iterative, GMRES, Modified Gram-Schmidt
        // linear_solver == 11 => Iterative, GMRES, Classical Gram-Schmidt
        // linear_solver == 20 => Iterative, Bi-CGStab (maxl => maximum dimension of Krylov subspace)
        // linear_solver == 30 => Iterative, TFQMR (maxl => maximum dimension of Krylov subspace)
        if (iter_type == 0)
            iter_type = (lmm == CV_ADAMS) ? 1 : 2;
        if (iter_type != 1 and iter_type != 2)
            throw std::runtime_error("Invalid iter_type");
        if (linear_solver == 0)
            linear_solver = (odesys->get_mlower() == -1) ? 1 : 2;
        auto integr = get_integrator<OdeSys>(odesys, atol, rtol, lmm, y0, tout[0], dx0, dx_min, dx_max, mxsteps,
                                             with_jacobian, iter_type, linear_solver, maxl, eps_lin);
        integr.predefined(nout, odesys->ny, tout, y0, yout, nderiv, root_indices);
        odesys->last_integration_info.clear();
        odesys->last_integration_info["n_steps"] = integr.get_n_steps();
        odesys->last_integration_info["n_rhs_evals"] = integr.get_n_rhs_evals();
        odesys->last_integration_info["n_lin_solv_setups"] = integr.get_n_lin_solv_setups();
        odesys->last_integration_info["n_err_test_fails"] = integr.get_n_err_test_fails();
        odesys->last_integration_info["n_nonlin_solv_iters"] = integr.get_n_nonlin_solv_iters();
        odesys->last_integration_info["n_nonlin_solv_conv_fails"] = integr.get_n_nonlin_solv_conv_fails();
        if (iter_type == 2){
            if (linear_solver >= 10) {  // iterative linear solver
                odesys->last_integration_info["krylov_n_lin_iters"] = integr.get_n_lin_iters();
                odesys->last_integration_info["krylov_n_prec_evals"] = integr.get_n_prec_evals();
                odesys->last_integration_info["krylov_n_prec_solves"] = integr.get_n_prec_solves();
                odesys->last_integration_info["krylov_n_conv_fails"] = integr.get_n_conv_fails();
                odesys->last_integration_info["krylov_n_jac_times_evals"] = integr.get_n_jac_times_evals();
                odesys->last_integration_info["krylov_n_iter_rhs"] = integr.get_n_iter_rhs();
            } else { // direct linear solver
                odesys->last_integration_info["dense_n_dls_jac_evals"] = integr.get_n_dls_jac_evals();
                odesys->last_integration_info["dense_n_dls_rhs_evals"] = integr.get_n_dls_rhs_evals();
            }
        }
    }

}
