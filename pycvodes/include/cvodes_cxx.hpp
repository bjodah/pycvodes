#pragma once
// Thin C++11 wrapper around CVODES from (SUNDIALS v2.7.0 and v3.1.x)
// far from all functionality has been wrapped yet.

#include <assert.h>
#include <cfenv>
#include <cmath>
#include <cstring>
#include <functional>
#include <memory>
#include <new> // bad_alloc
#include <utility>
#include <vector>
#include <unordered_map> // std::unordered_map
#include <sstream>
#include <stdlib.h>
#include <iostream>

#include <sundials/sundials_config.h>
#include "sundials_cxx.hpp" // sundials_cxx::nvector_serial::Vector
#include <cvodes/cvodes_spils.h>
#if SUNDIALS_VERSION_MAJOR >= 3
#  include <cvodes/cvodes_direct.h> /* CVODE fcts., CV_BDF, CV_ADAMS */
#  include <sunmatrix/sunmatrix_dense.h>
#  include <sunmatrix/sunmatrix_band.h>
#  include <sunmatrix/sunmatrix_sparse.h>
#  if !defined(USE_LAPACK)
#    if defined(SUNDIALS_BLAS_LAPACK)
#      define USE_LAPACK 1
#    else
#      define USE_LAPACK 0
#    endif
#  endif
#  if USE_LAPACK == 1
#    include <sunlinsol/sunlinsol_lapackdense.h>
#    include <sunlinsol/sunlinsol_lapackband.h>
#  else
#    include <sunlinsol/sunlinsol_dense.h>
#    include <sunlinsol/sunlinsol_band.h>
#  endif
#  include <sunlinsol/sunlinsol_spgmr.h>
#  include <sunlinsol/sunlinsol_spbcgs.h>
#  include <sunlinsol/sunlinsol_sptfqmr.h>
#else
#  if defined(SUNDIALS_PACKAGE_VERSION)   /* == 2.7.0 */
#    include <cvodes/cvodes_spgmr.h>
#    include <cvodes/cvodes_spbcgs.h>
#    include <cvodes/cvodes_sptfqmr.h>
#    if !defined(USE_LAPACK)
#      if defined(SUNDIALS_BLAS_LAPACK)
#        define USE_LAPACK 1
#      else
#        define USE_LAPACK 0
#      endif
#    endif
#    if USE_LAPACK == 1
#      include <cvodes/cvodes_lapack.h>       /* prototype for CVDense */
#    else
#      include <cvodes/cvodes_dense.h>
#      include <cvodes/cvodes_band.h>
#    endif
#    define SUNTRUE TRUE
#    define SUNFALSE FALSE
#  else
#    error "Unkown sundials version"
#  endif
#endif
#include <cvodes/cvodes.h> /* CVODE fcts., CV_BDF, CV_ADAMS */
#include <cvodes/cvodes_impl.h> /* CVodeMem */
#include <cvodes/cvodes_diag.h>       /* prototype for CVDiag */

#if !defined(BEGIN_NAMESPACE)
#define BEGIN_NAMESPACE(s) namespace s{
#endif
#if !defined(END_NAMESPACE)
#define END_NAMESPACE(s) }
#endif

//pragma STDC FENV_ACCESS on  // GCC 5.4 does not seem to support the pragma

BEGIN_NAMESPACE(cvodes_cxx)
class StreamFmt {
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

template<class T> void ignore( const T& ) { }

static const std::unordered_map<std::string, int> fpes {{
        {"FE_INEXACT", FE_INEXACT},
        {"FE_UNDERFLOW", FE_UNDERFLOW},
        {"FE_OVERFLOW", FE_OVERFLOW},
        {"FE_INVALID", FE_INVALID},
        {"FE_DIVBYZERO", FE_DIVBYZERO}
    }};
//using sundials_cxx::nvector_serial::N_Vector; // native sundials vector
using SVector = sundials_cxx::nvector_serial::Vector; // serial vector
using SVectorV = sundials_cxx::nvector_serial::VectorView; // view of serial vector
//using get_dx_max_fn = double(double, const double * const) *;
using get_dx_max_fn = std::function<double(double, const double * const)>;
// Wrapper for Revision 4306 of cvodes.c

enum class LMM : int {Adams=CV_ADAMS, BDF=CV_BDF}; // Linear multistep method
inline LMM lmm_from_name(std::string name){
    if (name == "adams")
        return LMM::Adams;
    else if (name == "bdf")
        return LMM::BDF;
    else
        throw std::runtime_error(StreamFmt() << "Unknown linear multistep method: " << name);
}

enum class Task : int {Normal=CV_NORMAL, One_Step=CV_ONE_STEP};
enum class IterType : int {Functional=CV_FUNCTIONAL, Newton=CV_NEWTON, Undecided=-1};
inline IterType iter_type_from_name(std::string name){
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

inline void check_flag(int flag) {
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
    FILE *errfp = nullptr;
#if SUNDIALS_VERSION_MAJOR >= 3
    SUNMatrix A_ = nullptr;
    SUNLinearSolver LS_ = nullptr;
    N_Vector y_ = nullptr;
    IterLinSolEnum solver_;
    std::vector<realtype> constraints_;
#endif
public:
    void *mem {nullptr};
    long int ny {0};
    int nq {0};
    int verbosity = 50;  // "50%" -- plenty of room for future tuning.
    bool autonomous_exprs = false;
    bool record_order = false, record_fpe = false, record_steps = false;
    double time_rhs {0}, time_jac {0}, time_roots {0}, time_quads {0}, time_prec {0}, time_jtimes {0};
    std::vector<int> orders_seen, fpes_seen;
    std::vector<double> steps_seen;  // Conversion from float / long double not a problem.
    Integrator(const LMM lmm, const IterType iter) {
        this->mem = CVodeCreate(static_cast<int>(lmm), static_cast<int>(iter));
        if (!this->mem)
            throw std::bad_alloc(); // "CVodeCreate failed (allocation failed)."
    }
    ~Integrator(){
        if (this->mem)
            CVodeFree(&(this->mem));
        if (this->errfp)
            fclose(errfp);
#if SUNDIALS_VERSION_MAJOR >= 3
        if (this->y_)
            N_VDestroy(this->y_);
        if (this->LS_)
            SUNLinSolFree(this->LS_);
        if (this->A_)
            SUNMatDestroy(this->A_);
#endif
    }
    // init
    void init(CVRhsFn cb, realtype t0, N_Vector y) {
        this->ny = NV_LENGTH_S(y);
#if SUNDIALS_VERSION_MAJOR >= 3
        if (y_)
            throw std::runtime_error("y_ already allocated");
        y_ = N_VNew_Serial(NV_LENGTH_S(y));
        std::memcpy(NV_DATA_S(y_), NV_DATA_S(y), ny*sizeof(realtype));
#else
        auto y_ = y;
#endif
        int status = CVodeInit(this->mem, cb, t0, y_);
        if (status == CV_ILL_INPUT)
            throw std::runtime_error("CVodeInit failed (CV_ILL_INPUT).");
        else if (status == CV_MEM_FAIL)
            throw std::bad_alloc(); // "CVodeInit failed (allocation failed).";
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
#if SUNDIALS_VERSION_MAJOR >= 3
        std::memcpy(NV_DATA_S(y_), NV_DATA_S(y), NV_LENGTH_S(y)*sizeof(realtype));
#else
        auto y_ = y;
#endif
        CVodeReInit(this->mem, t0, y_);
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
    // Quadrature integartion
    void quad_init(CVQuadRhsFn fQ, N_Vector yQ0){
        int flag = CVodeQuadInit(this->mem, fQ, yQ0);
        this->nq = NV_LENGTH_S(yQ0);
        if (flag == CV_MEM_FAIL)
            throw std::bad_alloc(); // "CVodeQuadInit failed (allocation failed)."
        else
            check_flag(flag);
    }
    void quad_init(CVQuadRhsFn fQ, std::vector<double>& yQ0){
        SVector yQ0_(yQ0.size(), yQ0.data());
        quad_init(fQ, yQ0_.n_vec);
    }
    void quad_reinit(const N_Vector yQ0){
        int flag = CVodeQuadReInit(this->mem, yQ0);
        if (flag == CV_MEM_FAIL)
            throw std::bad_alloc(); // "CVodeQuadInit failed (allocation failed)."
        else if (flag == CV_NO_QUAD)
            throw std::runtime_error("CVodeQuadReInit failed (CV_NO_QUAD).");
        else
            check_flag(flag);
    }
    void quad_reinit(SVector& yQ0){
        quad_reinit(yQ0.n_vec);
    };
    void quad_reinit(std::vector<double>& yQ0){
        SVector yQ0_(yQ0.size(), yQ0.data());
        quad_reinit(yQ0_.n_vec);
    }
    // Main solver optional input functions
    void set_err_file(FILE * errfp){
        int status = CVodeSetErrFile(this->mem, errfp);
        check_flag(status);
    }

    void set_err_file_path(const char * filename) {
        this->errfp = fopen(filename, "o");
        set_err_file(this->errfp);
    }

    // Step specifications
    void set_init_step(realtype h0){
        int status = CVodeSetInitStep(this->mem, h0);
        check_flag(status);
    }
    void set_min_step(realtype hmin){
        int status = CVodeSetMinStep(this->mem, hmin);
        if (status == CV_ILL_INPUT)
            throw std::runtime_error(StreamFmt() << "hmin=" << hmin << " non-positive or exceeding maximum allowable step size.");
        else
            check_flag(status);
    }
    void set_max_step(realtype hmax){
        int status = CVodeSetMaxStep(this->mem, hmax);
        if (status == CV_ILL_INPUT)
            throw std::runtime_error(StreamFmt() << "hmax=" << hmax << " non-positive or smaller than minimumem allowable step size");
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
        if (NV_LENGTH_S(atol) != ny)
            throw std::runtime_error("atol of incorrect length");
        int status = CVodeSVtolerances(this->mem, rtol, atol);
        if (status < 0)
            throw std::runtime_error("CVodeSVtolerances failed.");
    }
    void set_tol(realtype rtol, std::vector<realtype> &atol){
        SVector atol_(atol.size(), atol.data());
        set_tol(rtol, atol_.n_vec);
    }
    void set_quad_err_con(bool errconQ){ // quad_init must have been called
        int flag = CVodeSetQuadErrCon(this->mem, errconQ ? SUNTRUE : SUNFALSE);
        if (flag == CV_NO_QUAD)
            throw std::runtime_error("Quadrature integation has not been initialized");
        check_flag(flag);
    }
    void set_quad_tol(realtype reltolQ, realtype abstolQ){
        int flag = CVodeQuadSStolerances(this->mem, reltolQ, abstolQ);
        switch(flag){
        case CV_ILL_INPUT:
            throw std::runtime_error("One of the input tolerances was negative");
        case CV_NO_QUAD:
            throw std::runtime_error("Quadrature integation has not been initialized");
        default:
            check_flag(flag);
        }
    }
    void set_quad_tol(realtype reltolQ, const N_Vector abstolQ){
        if (NV_LENGTH_S(abstolQ) != nq)
            throw std::runtime_error("abstolQ of incorrect length");
        int flag = CVodeQuadSVtolerances(this->mem, reltolQ, abstolQ);
        switch(flag){
        case CV_ILL_INPUT:
            throw std::runtime_error("One of the input tolerances was negative");
        case CV_NO_QUAD:
            throw std::runtime_error("Quadrature integation has not been initialized");
        default:
            check_flag(flag);
        }
    }
    // void set_quad_tol(realtype reltolQ, SVectorV &abstolQ){
    //     set_quad_tol(reltolQ, abstolQ.n_vec);
    // }
    void set_quad_tol(realtype reltolQ, std::vector<realtype> &abstolQ){
        SVector atol_(abstolQ.size(), abstolQ.data());
        set_quad_tol(reltolQ, atol_.n_vec);
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
    // diagonal  solver
    void set_linear_solver_to_diag(){
        int status;
        status = CVDiag(this->mem);
        if (status < 0)
            throw std::runtime_error("CVDiag failed.");
    }

    // dense jacobian
    void set_linear_solver_to_dense(int ny){
        int status;
#if SUNDIALS_VERSION_MAJOR >= 3
        if (A_ == nullptr){
            if (A_)
                throw std::runtime_error("matrix already set");
            A_ = SUNDenseMatrix(ny, ny);
            if (!A_)
                throw std::runtime_error("SUNDenseMatrix failed.");
        }
        if (LS_ == nullptr){
            if (LS_)
                throw std::runtime_error("linear solver already set");
#if USE_LAPACK == 1
            LS_ = SUNLapackDense(y_, A_);
#else
            LS_ = SUNDenseLinearSolver(y_, A_);
#endif
            if (!LS_)
                throw std::runtime_error("SUNDenseLinearSolver failed.");
        }
        status = CVDlsSetLinearSolver(this->mem, LS_, A_);
        if (status < 0)
            throw std::runtime_error("CVDlsSetLinearSolver failed.");
#else
#if USE_LAPACK == 1
        status = CVLapackDense(this->mem, ny);
#else
        status = CVDense(this->mem, ny);
#endif
        if (status != CVDLS_SUCCESS) {
            throw std::runtime_error(
#if USE_LAPACK == 1
                "CVLapackDense failed"
#else
                "CVDense failed"
#endif
                );
        }
#endif
    }
    void set_dense_jac_fn(
#if SUNDIALS_VERSION_MAJOR >= 3
        CVDlsJacFn
#else
        CVDlsDenseJacFn
#endif
        djac
        ){
        int status;
#if SUNDIALS_VERSION_MAJOR >= 3
        status = CVDlsSetJacFn(this->mem, djac);
        if (status < 0)
            throw std::runtime_error("CVDlsSetJacFn failed.");
#else
        status = CVDlsSetDenseJacFn(this->mem, djac);
        if (status < 0)
            throw std::runtime_error("CVDlsSetDenseJacFn failed.");
#endif
    }

    // banded jacobian
    void set_linear_solver_to_banded(int N, int mupper, int mlower){
        int status;
#if SUNDIALS_VERSION_MAJOR >= 3
        ignore(N);
        if (A_ == nullptr){
            if (A_)
                throw std::runtime_error("matrix already set");
            A_ = SUNBandMatrix(ny, mupper, mlower, mlower+mupper);
            if (!A_)
                throw std::runtime_error("SUNDenseMatrix failed.");
        }
        if (LS_ == nullptr){
            if (LS_)
                throw std::runtime_error("linear solver already set");
#if USE_LAPACK == 1
            LS_ = SUNLapackBand(y_, A_);
#else
            LS_ = SUNBandLinearSolver(y_, A_);
#endif
            if (!LS_)
                throw std::runtime_error("SUNDenseLinearSolver failed.");
        }
        status = CVDlsSetLinearSolver(this->mem, LS_, A_);
        if (status < 0)
            throw std::runtime_error("CVDlsSetLinearSolver failed.");
#else
        status =
#if USE_LAPACK == 1
        CVLapackBand(this->mem, N, mupper, mlower);
#else
        CVBand(this->mem, N, mupper, mlower);
#endif
        if (status != CVDLS_SUCCESS)
            throw std::runtime_error(
#if USE_LAPACK == 1
                "CVLapackBand failed"
#else
                "CVBand failed"
#endif
                );
#endif
    }
    void set_band_jac_fn(
#if SUNDIALS_VERSION_MAJOR >= 3
        CVDlsJacFn
#else
        CVDlsBandJacFn
#endif
        djac){
        int status;
#if SUNDIALS_VERSION_MAJOR >= 3
        if (A_ == nullptr || LS_ == nullptr)
            throw std::runtime_error("set_linear_solver_to_banded not called?");
        status = CVDlsSetJacFn(this->mem, djac);
        if (status < 0)
            throw std::runtime_error("CVDlsSetJacFn failed.");
#else
        status = CVDlsSetBandJacFn(this->mem, djac);
        if (status < 0)
            throw std::runtime_error("CVDlsSetBandJacFn failed.");
#endif
    }

    // Iterative Linear solvers
    void set_linear_solver_to_iterative(IterLinSolEnum solver, int maxl=0, PrecType ptyp=PrecType::Left){
        int flag;
#if SUNDIALS_VERSION_MAJOR >= 3
        if (LS_)
            throw std::runtime_error("linear solver already set.");
#endif
        switch (solver) {
        case IterLinSolEnum::GMRES:
#if SUNDIALS_VERSION_MAJOR >= 3
            LS_ = SUNSPGMR(y_, static_cast<int>(ptyp), maxl);
#else
            flag = CVSpgmr(this->mem, static_cast<int>(ptyp), maxl);
#endif
            break;
        case IterLinSolEnum::BICGSTAB:
#if SUNDIALS_VERSION_MAJOR >= 3
            LS_ = SUNSPBCGS(y_, static_cast<int>(ptyp), maxl);
#else
            flag = CVSpbcg(this->mem, static_cast<int>(ptyp), maxl);
#endif
            break;
        case IterLinSolEnum::TFQMR:
#if SUNDIALS_VERSION_MAJOR >= 3
            LS_ = SUNSPTFQMR(y_, static_cast<int>(ptyp), maxl);
#else
            flag = CVSptfqmr(this->mem, static_cast<int>(ptyp), maxl);
#endif
            break;
        default:
            throw std::runtime_error(StreamFmt() << "Unknown solver kind: "
                                     << static_cast<int>(solver));
        }
#if SUNDIALS_VERSION_MAJOR >= 3
        flag = CVSpilsSetLinearSolver(this->mem, LS_);
        solver_ = solver;
#endif
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
    void set_jac_times_vec_fn(CVSpilsJacTimesVecFn jtimes){
#if SUNDIALS_VERSION_MAJOR >= 3
        int flag = CVSpilsSetJacTimes(this->mem, NULL, jtimes);
#else
        int flag = CVSpilsSetJacTimesVecFn(this->mem, jtimes);
#endif
        this->cvspils_check_flag(flag);
    }
#if SUNDIALS_VERSION_MAJOR >=3
    void set_jac_times_vec_fn(CVSpilsJacTimesSetupFn jtsetup,
                              CVSpilsJacTimesVecFn jtimes){
        int flag = CVSpilsSetJacTimes(this->mem, jtsetup, jtimes);
        this->cvspils_check_flag(flag);
    }
#endif
    void set_preconditioner(CVSpilsPrecSetupFn setup_fn, CVSpilsPrecSolveFn solve_fn){
        int flag = CVSpilsSetPreconditioner(this->mem, setup_fn, solve_fn);
        this->cvspils_check_flag(flag);
    }
    void set_iter_eps_lin(realtype delta){
        int flag = CVSpilsSetEpsLin(this->mem, delta);
        this->cvspils_check_flag(flag, true);
    }
    void set_prec_type(PrecType pretyp){
        int flag;
#if SUNDIALS_VERSION_MAJOR >= 3
        switch(solver_){
        case IterLinSolEnum::GMRES:
            flag = SUNSPGMRSetPrecType(LS_, (int)pretyp); break;
        case IterLinSolEnum::BICGSTAB:
            flag = SUNSPBCGSSetPrecType(LS_, (int)pretyp); break;
        case IterLinSolEnum::TFQMR:
            flag = SUNSPTFQMRSetPrecType(LS_, (int)pretyp); break;
        default:
            throw std::runtime_error("unknown solver kind.");
        }
#else
        flag = CVSpilsSetPrecType(this->mem, (int)pretyp);
#endif
        this->cvspils_check_flag(flag, true);
    }
    void set_gram_schmidt_type(GramSchmidtType gs_type){
        int flag;
#if SUNDIALS_VERSION_MAJOR >= 3
        if (solver_ == IterLinSolEnum::GMRES)
            flag = SUNSPGMRSetGSType(LS_, (int)gs_type);
        else
            throw std::runtime_error("Setting Gram-Schmidt type only makes sense for GMRES");
#else
        flag = CVSpilsSetGSType(this->mem, (int)gs_type);
#endif
        this->cvspils_check_flag(flag, true);
    }
    void set_krylov_max_len(int maxl){
        int flag;
#if SUNDIALS_VERSION_MAJOR >= 3
        switch(solver_){
        case IterLinSolEnum::GMRES:
            throw std::runtime_error("GMRES has no max length option");
        case IterLinSolEnum::BICGSTAB:
            flag = SUNSPBCGSSetMaxl(LS_, maxl); break;
        case IterLinSolEnum::TFQMR:
            flag = SUNSPTFQMRSetMaxl(LS_, maxl); break;
        default:
            throw std::runtime_error("unknown solver kind.");
        }
#else
        flag = CVSpilsSetMaxl(this->mem, maxl);
#endif
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
        SVectorV ew_(ny, ew);
        this->get_err_weights(ew_.n_vec);
    }
    void get_est_local_errors(N_Vector ele) const {
        check_flag(CVodeGetEstLocalErrors(this->mem, ele));
    }
    void get_est_local_errors(SVector &ele) const {
        this->get_est_local_errors(ele.n_vec);
    }
    void get_est_local_errors(realtype * ele) const {
        SVectorV ele_(ny, ele);
        this->get_est_local_errors(ele_.n_vec);
    }
    void set_constraints(N_Vector constraints) {
#if SUNDIALS_VERSION_MAJOR >= 3 && SUNDIALS_VERSION_MINOR >= 2
        if (NV_LENGTH_S(constraints) != ny)
            throw std::runtime_error("constraints of incorrect length");
        int status = CVodeSetConstraints(this->mem, constraints);
        if (status == CV_ILL_INPUT) {
            throw std::runtime_error("constraints vector contains illegal values");
        } else {
            check_flag(status);
        }
        constraints_.insert(constraints_.begin(),
                            NV_DATA_S(constraints),
                            NV_DATA_S(constraints) + NV_LENGTH_S(constraints));
#else
        ignore(constraints);
        throw std::runtime_error("setting constraints requires sundials >=3.2.0");
#endif
    }
    void set_constraints(const std::vector<realtype> &constraints) {
        SVector constraints_(constraints.size(), constraints.data());
        set_constraints(constraints_.n_vec);
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
    long int get_n_root_evals() const {
        long int res=0;
        int flag = CVodeGetNumGEvals(this->mem, &res);
        check_flag(flag);
        return res;
    }
    long int get_quad_num_rhs_evals() const {
        long int nfQevals = 0;
        int flag = CVodeGetQuadNumRhsEvals(this->mem, &nfQevals);
        if (flag == CV_NO_QUAD)
            throw std::runtime_error("Quadrature integation has not been initialized");
        check_flag(flag);
        return nfQevals;
    }
    long int get_quad_num_err_test_fails() const {
        long int nQetfails = 0;
        int flag = CVodeGetQuadNumErrTestFails(this->mem, &nQetfails);
        if (flag == CV_NO_QUAD)
            throw std::runtime_error("Quadrature integation has not been initialized");
        check_flag(flag);
        return nQetfails;
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
    int step(realtype tout, SVector &yout, realtype *tret, Task task){
        return CVode(this->mem, tout, yout.n_vec, tret, static_cast<int>(task));
    }
    void get_dky(realtype t, int k, SVector &dky) const {
        int flag = CVodeGetDky(this->mem, t, k, dky.n_vec);
        switch(flag){
        case CV_BAD_K:
            throw std::runtime_error("CVodeGetDky failed with (invalid k)");
        case CV_BAD_T:
            throw std::runtime_error("CVodeGetDky failed with (invalid t)");
        case CV_BAD_DKY:
            throw std::runtime_error("CVodeGetDky failed with (dky.n_vec was NULL)");
        default:
            check_flag(flag);
        }
    }
    void get_quad(realtype * const tret, N_Vector yQ){
        int flag = CVodeGetQuad(this->mem, tret, yQ);
        if (flag == CV_NO_QUAD)
            throw std::runtime_error("CVodeGetQuad failed (quadrature integration not initialized)");
        else if (flag == CV_BAD_DKY)
            throw std::runtime_error("CVodeGetQuad failed (yQ is NULL)");
        else
            check_flag(flag);
    }
    void get_quad_dky(realtype t, int k, N_Vector dkyQ){
        int flag = CVodeGetQuadDky(this->mem, t, k, dkyQ);
        switch(flag){
        case CV_NO_QUAD:
            throw std::runtime_error("CVodeGetQuadDky failed (CV_NO_QUAD)");
        case CV_BAD_DKY:
            throw std::runtime_error("CVodeGetQuadDky failed (BAD_DKY)");
        case CV_BAD_K:
            throw std::runtime_error("CVodeGetQuadDky failed (BAD_K)");
        case CV_BAD_T:
            throw std::runtime_error("CVodeGetQuadDky failed (BAD_T)");
        default:
            check_flag(flag);
        }
    }
    void call_rhs(realtype t, SVector y, SVector &ydot){
        CVodeMem cv_mem = (CVodeMem) this->mem;
        int status = cv_mem->cv_f(t, y.n_vec, ydot.n_vec, cv_mem->cv_user_data);
        if (status)
            throw std::runtime_error("call_rhs failed.");
    }
    void unsuccessful_step_throw_(int flag){
        throw std::runtime_error(StreamFmt() << std::scientific << "Unsuccessful step (t="
                                 << get_current_time() << ", h=" << get_current_step() << "): " <<
                                 CVodeGetReturnFlagName(flag));
    }
#define row(ti) (1 + ny*(nderiv+1) + nq)*(ti)
#define datalen(nt, nd, ny, nq) ((1 + (ny)*((nd)+1) + nq)*(nt)*sizeof(realtype))
#define ew_ele_len(nt, ny) (2*ny*nt*sizeof(realtype))
#define xout(ti) (*xyqout)[row(ti)]
#define yout(ti, di, yi) (*xyqout)[row(ti) + 1 + ny*(di) + yi]
#define qout(ti, qi) (*xyqout)[row(ti) + 1 + ny*(nderiv+1) + qi]
    int adaptive(realtype ** xyqout, // allocated using malloc (may be realloc:ed)
                 int * const td, // trailing dimension of xyqout
                 const realtype xend,
                 const unsigned nderiv,
                 std::vector<int>& root_indices,
                 bool return_on_root=false,
                 int autorestart=0,
                 bool return_on_error=false,
                 get_dx_max_fn get_dx_max = get_dx_max_fn(),
                 int tidx=0,
                 realtype ** ew_ele=nullptr // length(ew_ele) must be == 2*td*ny if not nullptr
        ){
        int status;
        SVector y {ny};
        SVector yQ {nq};
        SVector work {ny};
        long int mxsteps = get_max_num_steps();
        if (*td < 1)
            throw std::logic_error("Expected td >= 1");
        if (xyqout == nullptr)
            throw std::logic_error("xyqout cannot be a nullptr");
        if (*xyqout == nullptr)
            throw std::logic_error("xyqout cannot point to a nullptr");
        if (ew_ele){
            if (*ew_ele == nullptr) {
                throw std::logic_error("ew_ele cannot point to a nullptr");
            }
            for (int i=0; i<2*ny; ++i){
                (*ew_ele)[tidx*2*ny+i] = 0.0;
            }
        }
        realtype cur_t = xout(tidx);
        if (mxsteps == 0) { mxsteps = 500; } // cvodes default (MXSTEP_DEFAULT)
        if (record_steps)
            steps_seen.push_back(get_current_step());
        if (record_order)
            orders_seen.push_back(get_current_order()); // len(orders_seen) == len(xout)
        if (record_fpe){
            std::feclearexcept(FE_ALL_EXCEPT);
            fpes_seen.push_back(std::fetestexcept(FE_ALL_EXCEPT));  // gives equal length of output
        }
        for (int i=0; i<ny; ++i)
            y[i] = yout(tidx, 0, i);
        this->reinit(xout(tidx), y);
        if (nq){
            auto yQ0 = SVector(nq, &qout(tidx, 0));
            this->quad_reinit(yQ0);
        }
        if (nderiv >= 1){
            this->call_rhs(xout(tidx), y, work);
            for (int i=0; i<ny; ++i)
                yout(tidx, 1, i) = work[i];
            for (unsigned di=2; di<=nderiv; ++di){
                for (int i=0; i<ny; ++i)  // higher order too expensive
                    yout(tidx, di, i) = 0;
            }
        }
        if (xout(tidx) >= xend)  // negative step-sizes NOT supported (trade-off wrt. rounding errors)
            return tidx;
        this->set_stop_time(xend);
        do {
            tidx++;
            if (tidx >= *td){
                (*td) *= 2;
                {
                    void * new_xyqout = realloc(*xyqout, datalen(*td, nderiv, ny, nq));
                    if (new_xyqout == nullptr){
                        throw std::bad_alloc();
                    } else {
                        *xyqout = (realtype *)new_xyqout;
                    }
                }
                if (ew_ele)
                {
                    void * new_ew_ele = realloc(*ew_ele, ew_ele_len(*td, 2*ny));
                    if (new_ew_ele == nullptr){
                        throw std::bad_alloc();
                    } else {
                        *ew_ele = (realtype *)new_ew_ele;
                    }
                }
            }
            if (get_dx_max)
                this->set_max_step(get_dx_max(cur_t, y.get_data_ptr()));
            status = this->step(xend, y, &cur_t, Task::One_Step);
            if((status != CV_SUCCESS && status != CV_TSTOP_RETURN) || (tidx > mxsteps)){
                if (status == CV_ROOT_RETURN){
                    root_indices.push_back(tidx);
                }else{
                    if (autorestart == 0) {
                        if (return_on_error) {
                            --tidx;
                            break;
                        } else if (tidx > mxsteps) {
                            throw std::runtime_error(StreamFmt() << std::scientific << "Maximum number of steps reached (at t="
                                                     << cur_t <<"): " << mxsteps);
                        } else {
                            unsuccessful_step_throw_(status);
                        }
                    } else {
                        if (this->verbosity > 0){
                            std::cerr << "cvodes_cxx.hpp:" << __LINE__ << ": Autorestart (" << autorestart << ") t=" << cur_t << "\n";
                            if (status >= 0) {
                                N_Vector ele, ew;
                                ele = N_VNew_Serial(ny);
                                ew = N_VNew_Serial(ny);
                                get_est_local_errors(ele);
                                get_err_weights(ew);
                                realtype mx = 0.0;
                                int mxi = -1;
                                for (unsigned i=0; i < ny; ++i){
                                    const double cur = NV_DATA_S(ele)[i]*NV_DATA_S(ew)[i];
                                    if (cur > mx){
                                        mxi = i;
                                        mx = cur;
                                    }
                                }
                                N_VDestroy_Serial(ele);
                                N_VDestroy_Serial(ew);
                                std::cerr << "cvodes_cxx.hpp:" << __LINE__ << ":     max(ew[i]*ele[i]) = " << mx << ", i=" << mxi << "\n";
                            }
                        }
                        if (status == CV_CONV_FAILURE && autorestart == 1) { // Most likely close to singular matrix
                            if (this->verbosity > 0)
                                std::cerr << "    Singular Jacobian?";
                            this->set_tol(1e-3, 1e-3); if (this->verbosity > 0) std::cerr << " - using atol=1e-3, rtol=1e-3";
                            this->set_dense_jac_fn(nullptr); if (this->verbosity > 0) std::cerr << " - using finite differences.\n"; // Hail Mary
                        }
                        if (this->verbosity > 0) std::cerr << '\n';
                        if (tidx > mxsteps){
                            this->set_max_num_steps(this->get_max_num_steps() + 500);
                        }
                        const int step_back = (tidx > 1) ? 2 : 1;
                        const double last_x = xout(tidx - step_back);
                        if (autonomous_exprs)
                            xout(tidx - step_back) = 0; // allows for smaller step sizes
                        auto inner = this->adaptive(xyqout, td, xend - last_x, nderiv, root_indices, return_on_root,
                                                    autorestart-1, return_on_error, get_dx_max, tidx - step_back, ew_ele);
                        if (autonomous_exprs){
                            for (int i=tidx - step_back; i<=inner; ++i)
                                xout(i) += last_x;
                        }
                        tidx = inner;
                        break;
                    }
                }
            }
            xout(tidx) = cur_t;
            if (record_steps)
                steps_seen.push_back(get_current_step());
            if (record_order)
                orders_seen.push_back(get_current_order());
            if (record_fpe){
                fpes_seen.push_back(std::fetestexcept(FE_ALL_EXCEPT));
                std::feclearexcept(FE_ALL_EXCEPT);
            }

            for (int i=0; i<ny; ++i) {
#if SUNDIALS_VERSION_MAJOR >= 3 && SUNDIALS_VERSION_MINOR >= 2
                if (constraints_.size() and constraints_[i] == 1.0 and y[i] < 0) {
                    if (this->verbosity > 60) std::clog << "clipping y[" << i << "] to zero.\n";
                    yout(tidx, 0, i) = 0.0;
                } else {
                    yout(tidx, 0, i) = y[i];
                }
#else
                yout(tidx, 0, i) = y[i];
#endif
            }
            // Derivatives for interpolation
            for (unsigned di=1; di<=nderiv; ++di){
                if (this->get_n_steps() < 2*(nderiv+1))
                    // Too few points collected
                    work.zero_out();
                else
                    this->get_dky(cur_t, di, work);
                for (int i=0; i<ny; ++i)
                    yout(tidx, di, i) = work[i];
            }
            if (nq)
                get_quad_dky(cur_t, 0, yQ.n_vec);
            for (int i=0; i<nq; ++i)
                qout(tidx, i) = yQ[i];
            if (ew_ele) {
                this->get_err_weights(work);
                work.dump((*ew_ele) + 2*ny*tidx);
                this->get_est_local_errors(work);
                work.dump((*ew_ele) + 2*ny*tidx + ny);
            }
            if (return_on_root && status == CV_ROOT_RETURN)
                break;
        } while (status != CV_TSTOP_RETURN);

        if (*td > tidx + 1) { // Shrink xyqout:
            *td = tidx + 1;
            void * new_xyqout = realloc(*xyqout, datalen(*td, nderiv, ny, nq));
            if (new_xyqout == nullptr){
                throw std::bad_alloc();
            } else {
                *xyqout = (realtype *)new_xyqout;
            }
        }
        return tidx;
    }
#undef datalen
#undef xout
#undef yout
#undef row
    int predefined(const long int nt,
                   const realtype * const tout,
                   const realtype * const yq0,
                   realtype * const yqout,
                   const unsigned nderiv,
                   std::vector<int>& root_indices,
                   std::vector<realtype>& root_out,
                   int autorestart=0, // must be autonomous if >0b
                   bool return_on_error=false,
                   get_dx_max_fn get_dx_max=get_dx_max_fn(),
                   realtype * ew_ele=nullptr
        ){
        int iout = 0;
        realtype cur_t=tout[0];
        int status;
        SVector y {ny};
        SVector yQ {nq};
        SVector work {ny};
        long int mxsteps = get_max_num_steps();
        for (int i=0; i<ny; ++i)
            y[i] = yq0[i];
        this->reinit(tout[0], y);
        if (nq){
            auto yQ0 = SVector(nq, const_cast<realtype*>(yq0) + (nderiv+1)*ny);
            this->quad_reinit(yQ0);
        }
        if (record_order)
            orders_seen.push_back(get_current_order()); // len(orders_seen) == len(xout)
        if (record_fpe){
            std::feclearexcept(FE_ALL_EXCEPT);
            fpes_seen.push_back(std::fetestexcept(FE_ALL_EXCEPT));
        }
        y.dump(yqout);
        if (nderiv >= 1){
            this->call_rhs(tout[0], y, work);
            work.dump(&yqout[ny]);
        }
        for (unsigned di=1; di<nderiv; ++di){
            for (int i=0; i<ny; ++i)  // too expensive
                yqout[ny*(di+1) + i] = 0;
        }
        for (int i=0; i<nq; ++i){
            yqout[ny*(nderiv+1) + i] = yq0[ny+i];
        }
        if (ew_ele) {
            for (int i=0; i<2*ny; ++i){
                ew_ele[i] = 0.0;
            }
        }
        for(iout=1; (iout < nt); iout++) {
            if (get_dx_max)
                this->set_max_step(get_dx_max(cur_t, y.get_data_ptr()));
            status = this->step(tout[iout], y, &cur_t, Task::Normal);
            if(status == CV_SUCCESS){
                if (record_order)
                    orders_seen.push_back(get_current_order());
                if (record_fpe){
                    fpes_seen.push_back(std::fetestexcept(FE_ALL_EXCEPT));
                    std::feclearexcept(FE_ALL_EXCEPT);
                }
                y.dump(&yqout[iout*(ny*(nderiv+1) + nq)]);
                // Derivatives for interpolation
                for (unsigned di=0; di<nderiv; ++di){
                    if (this->get_n_steps() < 2*(nderiv+1))
                        // Too few points collected
                        work.zero_out();
                    else
                        this->get_dky(tout[iout], di+1, work);
                    work.dump(&yqout[iout*(ny*(nderiv+1) + nq) + (di+1)*ny]);
                }
                if (nq)
                    get_quad_dky(cur_t, 0, yQ.n_vec);
                for (int i=0; i<nq; ++i)
                    yqout[iout*(ny*(nderiv+1) + nq) + (nderiv+1)*ny + i] = yQ[i];
                if (ew_ele){
                    this->get_err_weights(work);
                    work.dump(ew_ele + 2*ny*iout);
                    this->get_est_local_errors(work);
                    work.dump(ew_ele + 2*ny*iout + ny);
                }
            } else if (status == CV_ROOT_RETURN) {
                root_out.push_back(cur_t);
                for (int i=0; i<ny; ++i)
                    root_out.push_back(y[i]);
                root_indices.push_back(iout);
                iout--;
            } else {  // unsuccessful step
                if (autorestart == 0){
                    if (return_on_error){
                        break;
                    } else {
                        unsuccessful_step_throw_(status);
                    }
                } else {
                    if (this->verbosity > 0){
                        std::cerr << "cvodes_cxx.hpp:" << __LINE__ << ": Autorestart (" << autorestart << ") t=" << cur_t << "\n";
                        if (status >= 0) {
                            N_Vector ele, ew;
                            ele = N_VNew_Serial(ny);
                            ew = N_VNew_Serial(ny);
                            get_est_local_errors(ele);
                            get_err_weights(ew);
                            realtype mx = 0.0;
                            int mxi = -1;
                            for (unsigned i=0; i < ny; ++i){
                                const double cur = NV_DATA_S(ele)[i]*NV_DATA_S(ew)[i];
                                if (cur > mx){
                                    mxi = i;
                                    mx = cur;
                                }
                            }
                            N_VDestroy_Serial(ele);
                            N_VDestroy_Serial(ew);
                            std::cerr << "cvodes_cxx.hpp:" << __LINE__ << ":     max(ew[i]*ele[i]) = " << mx << ", i=" << mxi << "\n";
                        }
                    }
                    this->set_max_num_steps(mxsteps + this->get_max_num_steps());
                    std::vector<double> tout_;
                    long int nleft = nt - iout + 1;
                    for (int i=0; i<nleft; ++i)
                        tout_.push_back(tout[iout + i - 1] - tout[iout - 1]);
                    std::vector<int> root_indices_;
                    std::vector<realtype> root_out_;
#if SUNDIALS_VERSION_MAJOR >= 3 && SUNDIALS_VERSION_MINOR >= 2
                    if (constraints_.size()) {
                        for (int i=0; i<ny; ++i){
                            if (constraints_[i] == 1.0 and yqout[i + (iout-1)*((nderiv+1)*ny+nq)] < 0) {
                                if (this->verbosity > 60) std::clog << "clipping to y[" << i << "] to zero.\n";
                                yqout[i + (iout-1)*((nderiv+1)*ny+nq)] = 0;
                            }
                        }
                    }
#endif
                    int n_reached = this->predefined(nleft, &tout_[0],
                                                     yqout + (iout-1)*((nderiv+1)*ny+nq),
                                                     yqout + (iout-1)*((nderiv+1)*ny+nq),
                                                     nderiv, root_indices_, root_out_,
                                                     autorestart-1,
                                                     return_on_error,
                                                     get_dx_max);
                    for (const auto &ri : root_indices_)
                        root_indices.push_back(iout+ri-1);
                    root_out.insert(root_out.end(), root_out_.begin(), root_out_.end());
                    iout += n_reached - 1;
                    break;
                }
            }
        }
        return iout;
    }

};

template<typename T>
void extend_vector(std::vector<T> &dest, const std::vector<T> &source)
{
    dest.reserve(dest.size() + std::distance(source.begin(), source.end()));
    dest.insert(dest.end(), source.begin(), source.end());
}

inline void update_integration_info(std::unordered_map<std::string, int> &info_int,
                                    std::unordered_map<std::string, double> &info_dbl,
                                    std::unordered_map<std::string, std::vector<double>> &info_vecdbl,
                                    std::unordered_map<std::string, std::vector<int>> &info_vecint,
                                    const Integrator& integrator,
                                    const IterType iter_type, const int linear_solver)
{
    info_int["n_steps"] += integrator.get_n_steps();
    info_int["n_root_evals"] += integrator.get_n_root_evals();
    info_int["n_rhs_evals"] += integrator.get_n_rhs_evals();
    info_int["n_lin_solv_setups"] += integrator.get_n_lin_solv_setups();
    info_int["n_err_test_fails"] += integrator.get_n_err_test_fails();
    info_int["n_nonlin_solv_iters"] += integrator.get_n_nonlin_solv_iters();
    info_int["n_nonlin_solv_conv_fails"] += integrator.get_n_nonlin_solv_conv_fails();
    if (iter_type == IterType::Newton){
        if (linear_solver >= 10) {  // iterative linear solver
            info_int["krylov_n_lin_iters"] += integrator.get_n_lin_iters();
            info_int["krylov_n_prec_evals"] += integrator.get_n_prec_evals();
            info_int["krylov_n_prec_solves"] += integrator.get_n_prec_solves();
            info_int["krylov_n_conv_fails"] += integrator.get_n_conv_fails();
            info_int["krylov_n_jac_times_evals"] += integrator.get_n_jac_times_evals();
            info_int["krylov_n_iter_rhs"] += integrator.get_n_iter_rhs();
        } else { // direct linear solver
            info_int["dense_n_dls_jac_evals"] += integrator.get_n_dls_jac_evals();
            info_int["dense_n_dls_rhs_evals"] += integrator.get_n_dls_rhs_evals();
        }
    }
    info_dbl["time_rhs"] += integrator.time_rhs;
    info_dbl["time_quads"] += integrator.time_quads;
    info_dbl["time_roots"] += integrator.time_roots;
    info_dbl["time_jac"] += integrator.time_jac;
    info_dbl["time_jtimes"] += integrator.time_jtimes;
    info_dbl["time_prec"] += integrator.time_prec;
    extend_vector(info_vecint["orders"], integrator.orders_seen);
    extend_vector(info_vecint["fpes"], integrator.fpes_seen);
    extend_vector(info_vecdbl["steps"], integrator.steps_seen);
}
END_NAMESPACE(cvodes_cxx)
