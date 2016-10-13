#pragma once
// Thin C++11 wrapper around CVODES v2.8.2 (SUNDIALS v2.6.2)
// far from all functionality is available yet.
// sundials-2.6.2.tar.gz (MD5: 3deeb0ede9f514184c6bd83ecab77d95)

#include <assert.h>
#include <cmath>
#include <cstring>
#include <memory>
#include <new> // bad_alloc
#include <utility>
#include <vector>
#include <unordered_map> // std::unordered_map
#include <sstream>

#include "sundials_cxx.hpp" // sundials_cxx::nvector_serial::Vector
#include <cvodes/cvodes_spils.h>
#include <cvodes/cvodes_spgmr.h>
#include <cvodes/cvodes_spbcgs.h>
#include <cvodes/cvodes_sptfqmr.h>
#include <cvodes/cvodes.h> /* CVODE fcts., CV_BDF, CV_ADAMS */
#include <cvodes/cvodes_impl.h> /* CVodeMem */
#include <cvodes/cvodes_lapack.h>       /* prototype for CVDense */


namespace {
    class StreamFmt
    {
        std::stringstream m_s;
    public:
        StreamFmt() {}
        ~StreamFmt() {}

        template <typename T>
        StreamFmt& operator << (const T& v) {
            this->m_s << v;
            return *this;
        }

        std::string str() const {
            return this->m_s.str();
        }
        operator std::string() const {
            return this->m_s.str();
        }

    };
}

namespace cvodes_cxx {

    //using sundials_cxx::nvector_serial::N_Vector; // native sundials vector
    using SVector = sundials_cxx::nvector_serial::Vector; // serial vector

    // Wrapper for Revision 4306 of cvodes.c

    enum class LMM : int {Adams=CV_ADAMS, BDF=CV_BDF}; // Linear multistep method
    LMM lmm_from_name(std::string name){
        if (name == "adams")
            return LMM::Adams;
        else if (name == "bdf")
            return LMM::BDF;
        else
            throw std::runtime_error(StreamFmt() << "Unknown linear multistep method: " << name);
    }

    enum class Task : int {Normal=CV_NORMAL, One_Step=CV_ONE_STEP};
    enum class IterType : int {Functional=CV_FUNCTIONAL, Newton=CV_NEWTON, Undecided=-1};
    IterType iter_type_from_name(std::string name){
        if (name == "functional")
            return IterType::Functional;
        else if (name == "newton")
            return IterType::Newton;
        else if (name == "undecided")
            return IterType::Undecided;
        else
            throw std::runtime_error(StreamFmt() << "Unknown iter_type: " << name);
    }
    enum class IterLinSolEnum : int {GMRES=1, BICGSTAB=2, TFQMR=3};
    enum class PrecType : int {None=PREC_NONE, Left=PREC_LEFT,
            Right=PREC_RIGHT, Both=PREC_BOTH};
    enum class GramSchmidtType : int {Classical=CLASSICAL_GS, Modified=MODIFIED_GS};

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

    class CVodeIntegrator{ // Thin wrapper class of CVode in CVODES
    public:
        void *mem {nullptr};
        long int ny {0};  // cvodes uses a signed data type for this...
        CVodeIntegrator(const LMM lmm, const IterType iter) {
            this->mem = CVodeCreate(static_cast<int>(lmm), static_cast<int>(iter));
            if (!this->mem)
                throw std::bad_alloc(); // "CVodeCreate failed (allocation failed)."
        }
        ~CVodeIntegrator(){
            if (this->mem)
                CVodeFree(&(this->mem));
        }
        // init
        void init(CVRhsFn cb, realtype t0, N_Vector y) {
            int status = CVodeInit(this->mem, cb, t0, y);
            if (status == CV_ILL_INPUT)
                throw std::runtime_error("CVodeInit failed (CV_ILL_INPUT).");
            else if (status == CV_MEM_FAIL)
                throw std::bad_alloc(); // "CVodeInit failed (allocation failed).";
            else
                check_flag(status);
            this->ny = NV_LENGTH_S(y);
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
                throw std::bad_alloc(); // "CVodeRootInit failed (allocation failed)."
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
        long int get_max_num_steps(){
            return static_cast<CVodeMem>(this->mem)->cv_mxstep;
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
                flag = CVSpgmr(this->mem, static_cast<int>(PrecType::Left), maxl);
                break;
            case IterLinSolEnum::BICGSTAB:
                flag = CVSpbcg(this->mem, static_cast<int>(PrecType::Left), maxl);
                break;
            case IterLinSolEnum::TFQMR:
                flag = CVSptfqmr(this->mem, static_cast<int>(PrecType::Left), maxl);
                break;
            default:
                throw std::runtime_error("unknown solver kind.");
            }
            switch (flag){
            case CVSPILS_SUCCESS:
                break;
            case CVSPILS_MEM_NULL:
                throw std::runtime_error("set_linear_solver_to_iterative failed (cvode_mem is NULL)");
            case CVSPILS_ILL_INPUT:
                throw std::runtime_error("PREC_LEFT invalid.");
            case CVSPILS_MEM_FAIL:
                throw std::bad_alloc(); // "Memory allocation for sparse iterative solver request failed."
            }
        }
        void cvspils_check_flag(int flag, bool check_ill_input=false) const {
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
        int get_current_order() {
            int qcur;
            CVodeGetCurrentOrder(this->mem, &qcur);
            return qcur;
        }
        realtype get_current_step() {
            realtype hcur;
            CVodeGetCurrentStep(this->mem, &hcur);
            return hcur;
        }
        realtype get_current_time() {
            realtype tcur;
            CVodeGetCurrentTime(this->mem, &tcur);
            return tcur;
        }
        void get_err_weights(N_Vector ew) const {
            check_flag(CVodeGetErrWeights(this->mem, ew));
        }
        void get_err_weights(SVector &ew) const {
            this->get_err_weights(ew.n_vec);
        }
        void get_err_weights(realtype * ew) const {
            SVector ew_(ny, ew);
            this->get_err_weights(ew_);
        }

        // get info
        long int get_n_lin_iters() const {
            long int res=0;
            int flag;
            flag = CVSpilsGetNumLinIters(this->mem, &res);
            this->cvspils_check_flag(flag);
            return res;
        }

        long int get_n_prec_evals() const {
            long int res=0;
            int flag;
            flag = CVSpilsGetNumPrecEvals(this->mem, &res);
            this->cvspils_check_flag(flag);
            return res;
        }

        long int get_n_prec_solves() const {
            long int res=0;
            int flag;
            flag = CVSpilsGetNumPrecSolves(this->mem, &res);
            this->cvspils_check_flag(flag);
            return res;
        }

        long int get_n_conv_fails() const {
            long int res=0;
            int flag;
            flag = CVSpilsGetNumConvFails(this->mem, &res);
            this->cvspils_check_flag(flag);
            return res;
        }

        long int get_n_jac_times_evals() const {
            long int res=0;
            int flag;
            flag = CVSpilsGetNumJtimesEvals(this->mem, &res);
            this->cvspils_check_flag(flag);
            return res;
        }

        long int get_n_iter_rhs() const {
            long int res=0;
            int flag;
            flag = CVSpilsGetNumRhsEvals(this->mem, &res);
            this->cvspils_check_flag(flag);
            return res;
        }

        long int get_n_steps() const {
            long int res=0;
            int flag = CVodeGetNumSteps(this->mem, &res);
            check_flag(flag);
            return res;
        }

        long int get_n_rhs_evals() const {
            long int res=0;
            int flag = CVodeGetNumRhsEvals(this->mem, &res);
            check_flag(flag);
            return res;
        }

        long int get_n_lin_solv_setups() const {
            long int res=0;
            int flag = CVodeGetNumLinSolvSetups(this->mem, &res);
            check_flag(flag);
            return res;
        }

        long int get_n_err_test_fails() const {
            long int res=0;
            int flag = CVodeGetNumErrTestFails(this->mem, &res);
            check_flag(flag);
            return res;
        }

        long int get_n_nonlin_solv_iters() const {
            long int res=0;
            int flag = CVodeGetNumNonlinSolvIters(this->mem, &res);
            check_flag(flag);
            return res;
        }

        long int get_n_nonlin_solv_conv_fails() const {
            long int res=0;
            int flag = CVodeGetNumNonlinSolvConvFails(this->mem, &res);
            check_flag(flag);
            return res;
        }

        void cvdls_check_flag(int flag) const {
            switch (flag){
            case CVDLS_SUCCESS:
                break;
            case CVDLS_MEM_NULL:
                throw std::runtime_error("cvode_mem is NULL");
            case CVDLS_LMEM_NULL:
                throw std::runtime_error("CVDLS linear solver has not been initialized)");
            }
        }
        long int get_n_dls_jac_evals() const {
            long int res=0;
            int flag = CVDlsGetNumJacEvals(this->mem, &res);
            cvdls_check_flag(flag);
            return res;
        }

        long int get_n_dls_rhs_evals() const {
            long int res=0;
            int flag = CVDlsGetNumRhsEvals(this->mem, &res);
            cvdls_check_flag(flag);
            return res;
        }

        void get_dky(realtype t, int k, SVector &dky) const {
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
        void unsuccessful_step_(int flag){
            throw std::runtime_error(StreamFmt() << std::scientific << "Unsuccessful step (t="
                                     << get_current_time() << ", h=" << get_current_step() << "): " <<
                                     CVodeGetReturnFlagName(flag));
        }
        std::pair<std::vector<realtype>, std::vector<realtype> >
        adaptive(const realtype x0,
                 const realtype xend,
                 const realtype * const y0,
                 const unsigned nderiv,
                 std::vector<int>& root_indices,
                 bool return_on_root=false){
            std::vector<realtype> xout;
            std::vector<realtype> yout;
            realtype cur_t;
            int status;
            int idx = 0;
            SVector y {ny};
            SVector work {ny};
            long int mxsteps = get_max_num_steps();
            if (mxsteps == 0) { mxsteps = 500; } // cvodes default (MXSTEP_DEFAULT)
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
            for (unsigned di=1; di<nderiv; ++di){
                for (int i=0; i<ny; ++i)  // higher order too expensive
                    yout.push_back(0);
            }
            this->set_stop_time(xend);
            do {
                idx++;
                if (idx > 0 and idx > mxsteps)
                    throw std::runtime_error(StreamFmt() << std::scientific << "Maximum number of steps reached (at t="
                                             << cur_t <<"): " << mxsteps);
                status = this->step(xend, y, &cur_t, Task::One_Step);
                if(status != CV_SUCCESS && status != CV_TSTOP_RETURN){
                    if (status == CV_ROOT_RETURN){
                        root_indices.push_back(idx);
                    }else{
                        unsuccessful_step_(status);
                    }
                }
                xout.push_back(cur_t);
                for (int i=0; i<ny; ++i)
                    yout.push_back(y[i]);
                // Derivatives for interpolation
                for (unsigned di=0; di<nderiv; ++di){
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

        void predefined(const long int nt,  // sundials does not use unsigned types here...
                        const realtype * const tout,
                        const realtype * const y0,
                        realtype * const yout,
                        const unsigned nderiv,
                        std::vector<int>& root_indices,
                        std::vector<realtype>& root_out){
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
            for (unsigned di=1; di<nderiv; ++di){
                for (int i=0; i<ny; ++i)  // too expensive
                    yout[ny*(di+1) + i] = 0;
            }

            for(int iout=1; (iout < nt); iout++) {
                status = this->step(tout[iout], y, &cur_t, Task::Normal);
                if(status != CV_SUCCESS){
                    if (status == CV_ROOT_RETURN){
                        root_out.push_back(cur_t);
                        for (int i=0; i<ny; ++i)
                            root_out.push_back(y[i]);
                        root_indices.push_back(iout);
                        iout--;
                        continue;
                    }else{
                        unsuccessful_step_(status);
                    }
                }
                y.dump(&yout[ny*(iout*(nderiv+1))]);
                // Derivatives for interpolation
                for (unsigned di=0; di<nderiv; ++di){
                    if (this->get_n_steps() < 2*(nderiv+1))
                        // Too few points collected
                        work.zero_out();
                    else
                        this->get_dky(tout[iout], di+1, work);
                    work.dump(&yout[ny*(di+1+(iout*(nderiv+1)))]);
                }
            }
        }

    };

    void set_integration_info(std::unordered_map<std::string, int>& info,
                              const CVodeIntegrator& integrator,
                              IterType iter_type, int linear_solver){
        info["n_steps"] = integrator.get_n_steps();
        info["n_rhs_evals"] = integrator.get_n_rhs_evals();
        info["n_lin_solv_setups"] = integrator.get_n_lin_solv_setups();
        info["n_err_test_fails"] = integrator.get_n_err_test_fails();
        info["n_nonlin_solv_iters"] = integrator.get_n_nonlin_solv_iters();
        info["n_nonlin_solv_conv_fails"] = integrator.get_n_nonlin_solv_conv_fails();
        if (iter_type == IterType::Newton){
            if (linear_solver >= 10) {  // iterative linear solver
                info["krylov_n_lin_iters"] = integrator.get_n_lin_iters();
                info["krylov_n_prec_evals"] = integrator.get_n_prec_evals();
                info["krylov_n_prec_solves"] = integrator.get_n_prec_solves();
                info["krylov_n_conv_fails"] = integrator.get_n_conv_fails();
                info["krylov_n_jac_times_evals"] = integrator.get_n_jac_times_evals();
                info["krylov_n_iter_rhs"] = integrator.get_n_iter_rhs();
            } else { // direct linear solver
                info["dense_n_dls_jac_evals"] = integrator.get_n_dls_jac_evals();
                info["dense_n_dls_rhs_evals"] = integrator.get_n_dls_rhs_evals();
            }
        }
    }

}
