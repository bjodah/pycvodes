#pragma once

#include <cstdlib>
#include "anyode/anyode.hpp"
#include "anyode/anyode_parallel.hpp"
#include "cvodes_anyode.hpp"

namespace cvodes_anyode_parallel {

    using cvodes_cxx::LMM;
    using cvodes_cxx::IterType;
    using cvodes_cxx::LinSol;
    using cvodes_anyode::simple_adaptive;
    using cvodes_anyode::simple_predefined;

    template <class OdeSys>
    std::vector<std::pair<int, std::vector<int>>>
    multi_adaptive(realtype ** xyout_arr, // vectorized
                   int * td_arr, // vectorized
                   std::vector<OdeSys *> odesys, // vectorized
                   std::vector<realtype> atol,
                   const realtype rtol,
                   const LMM lmm,
                   const realtype * tend,  // vectorized
                   const long int mxsteps,
                   const realtype * dx0,  // vectorized
                   const realtype * dx_min,  // vectorized
                   const realtype * dx_max,  // vectorized
                   const bool with_jacobian=false,
                   IterType iter_type=IterType::Undecided,
                   LinSol linear_solver=LinSol::DEFAULT,
                   const int maxl=0,
                   const realtype eps_lin=0.0,
                   const unsigned nderiv=0,
                   const bool return_on_root=false,
                   const int autorestart=0, // must be autonomous!
                   const bool return_on_error=false,
                   const bool with_jtimes=false,
                   int tidx=0,
                   realtype *** ew_ele=nullptr,
                   const std::vector<realtype> &constraints={}
                   ){
        const indextype ny = odesys[0]->get_ny();
        AnyODE::ignore(ny);
        const int nsys = odesys.size();
        auto results = std::vector<std::pair<int, std::vector<int>>>(nsys);
        anyode_parallel::ThreadException te;
        char * num_threads_var = std::getenv("ANYODE_NUM_THREADS");
        int nt = (num_threads_var) ? std::atoi(num_threads_var) : 1;
        if (nt < 0)
            nt = 1;
        #pragma omp parallel for num_threads(nt) // OMP_NUM_THREADS should be 1 for openblas LU (small matrices)
        for (int idx=0; idx<nsys; ++idx){
            te.run([&]{
                results[idx].first = simple_adaptive<OdeSys>(
                    xyout_arr + idx, td_arr + idx,
                    odesys[idx], atol, rtol, lmm, tend[idx],
                    results[idx].second, mxsteps, dx0[idx], dx_min[idx], dx_max[idx],
                    with_jacobian, iter_type, linear_solver, maxl, eps_lin, nderiv,
                    return_on_root, autorestart, return_on_error, with_jtimes, tidx, (ew_ele) ? ew_ele[idx] : nullptr, constraints);
            });
        }
        te.rethrow();

        return results;
    }

    template <class OdeSys>
    std::vector<std::pair<int, std::pair<std::vector<int>, std::vector<realtype>>>>
    multi_predefined(std::vector<OdeSys *> odesys,  // vectorized
                     std::vector<realtype> atol,
                     const realtype rtol,
                     const LMM lmm,
                     realtype * y0, // vectorized
                     const std::size_t nout,
                     realtype * tout, // vectorized
                     realtype * yout,  // vectorized
                     const long int mxsteps,
                     const realtype * dx0,  // vectorized
                     const realtype * dx_min,  // vectorized
                     const realtype * dx_max,  // vectorized
                     const bool with_jacobian=false,
                     IterType iter_type=IterType::Undecided,
                     LinSol linear_solver=LinSol::DEFAULT,
                     const int maxl=0,
                     const realtype eps_lin=0.0,
                     const unsigned nderiv=0,
                     const int autorestart=0, // must be autonomous!
                     const bool return_on_error=false,
                     const bool with_jtimes=false,
                     realtype ** ew_ele=nullptr,
                     const std::vector<realtype> &constraints={}
                     ){
        const indextype ny = odesys[0]->get_ny();
        const int nsys = odesys.size();

        auto nreached_roots = std::vector<std::pair<int, std::pair<std::vector<int>, std::vector<realtype>>>>(nsys);

        anyode_parallel::ThreadException te;
        char * num_threads_var = std::getenv("ANYODE_NUM_THREADS");
        int nt = (num_threads_var) ? std::atoi(num_threads_var) : 1;
        if (nt < 0)
            nt = 1;

        #pragma omp parallel for num_threads(nt) // OMP_NUM_THREADS should be 1 for openblas LU (small matrices)
        for (int idx=0; idx<nsys; ++idx){
            te.run([&]{
                    nreached_roots[idx].first = simple_predefined<OdeSys>(
                        odesys[idx], atol, rtol, lmm, y0 + idx*ny,
                        nout, tout + idx*nout, yout + idx*ny*nout*(nderiv+1),
                        nreached_roots[idx].second.first, nreached_roots[idx].second.second,
                        mxsteps, dx0[idx], dx_min[idx], dx_max[idx], with_jacobian,
                        iter_type, linear_solver, maxl, eps_lin, nderiv,
                        autorestart, return_on_error, with_jtimes, (ew_ele) ? ew_ele[idx] : nullptr, constraints);
            });
        }
        if (!return_on_error)
            te.rethrow();

        return nreached_roots;
    }

}
