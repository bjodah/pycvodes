#pragma once

#include "anyode/anyode_parallel.hpp"
#include "cvodes_anyode.hpp"

namespace cvodes_anyode_parallel {

    using cvodes_cxx::LMM;
    using cvodes_cxx::IterType;
    using cvodes_anyode::simple_adaptive;
    using cvodes_anyode::simple_predefined;

    using sa_t = std::pair<std::vector<realtype>, std::vector<realtype> >;

    template <class OdeSys>
    std::vector<std::pair<sa_t, std::vector<int>>>
    multi_adaptive(std::vector<OdeSys *> odesys, // vectorized
                   const std::vector<realtype> atol,
                   const realtype rtol,
                   const LMM lmm,
                   const realtype * const y0,  // vectorized
                   const realtype * t0,  // vectorized
                   const realtype * tend,  // vectorized
                   const long int mxsteps,
                   const realtype * dx0,  // vectorized
                   const realtype * dx_min,  // vectorized
                   const realtype * dx_max,  // vectorized
                   const bool with_jacobian=false,
                   IterType iter_type=IterType::Undecided,
                   int linear_solver=0,
                   const int maxl=0,
                   const realtype eps_lin=0.0,
                   const unsigned nderiv=0,
                   bool return_on_root=false,
                   int autorestart=0, // must be autonomous!
                   bool return_on_error=false
                   ){
        const int ny = odesys[0]->get_ny();
        const int nsys = odesys.size();
        auto results = std::vector<std::pair<sa_t, std::vector<int>>>(nsys);

        anyode_parallel::ThreadException te;
        #pragma omp parallel for
        for (int idx=0; idx<nsys; ++idx){
            std::pair<sa_t, std::vector<int>> local_result;
            te.run([&]{
                local_result.first = simple_adaptive<OdeSys>(
                    odesys[idx], atol, rtol, lmm, y0 + idx*ny, t0[idx], tend[idx],
                    local_result.second, mxsteps, dx0[idx], dx_min[idx], dx_max[idx], with_jacobian,
                    iter_type, linear_solver, maxl, eps_lin, nderiv, return_on_root,
                    autorestart, return_on_error);
            });
            results[idx] = local_result;
        }
        te.rethrow();

        return results;
    }

    template <class OdeSys>
    std::vector<std::pair<std::vector<int>, std::vector<double>>>
    multi_predefined(std::vector<OdeSys *> odesys,  // vectorized
                     const std::vector<realtype> atol,
                     const realtype rtol,
                     const LMM lmm,
                     const realtype * const y0, // vectorized
                     const std::size_t nout,
                     const realtype * const tout, // vectorized
                     realtype * const yout,  // vectorized
                     const long int mxsteps,
                     const realtype * dx0,  // vectorized
                     const realtype * dx_min,  // vectorized
                     const realtype * dx_max,  // vectorized
                     const bool with_jacobian=false,
                     IterType iter_type=IterType::Undecided,
                     int linear_solver=0,
                     const int maxl=0,
                     const realtype eps_lin=0.0,
                     const unsigned nderiv=0,
                     int autorestart=0 // must be autonomous!
                     ){
        const int ny = odesys[0]->get_ny();
        const int nsys = odesys.size();

        auto roots = std::vector<std::pair<std::vector<int>, std::vector<double>>>(nsys);

        anyode_parallel::ThreadException te;
        #pragma omp parallel for
        for (int idx=0; idx<nsys; ++idx){
            te.run([&]{
               simple_predefined<OdeSys>(odesys[idx], atol, rtol, lmm, y0 + idx*ny,
                                         nout, tout + idx*nout, yout + idx*ny*nout*(nderiv+1),
                                         roots[idx].first, roots[idx].second,
                                         mxsteps, dx0[idx], dx_min[idx], dx_max[idx], with_jacobian,
                                         iter_type, linear_solver, maxl, eps_lin, nderiv,
                                         autorestart);
            });
        }
        te.rethrow();

        return roots;
    }

}
