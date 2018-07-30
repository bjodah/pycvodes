#include "catch.hpp"
#include <math.h>
#include <vector>
#include <cstdlib>
#include <stdio.h>
#include "anyode/anyode.hpp"
#include "cvodes_anyode.hpp"
#include "cetsa_case.hpp"

TEST_CASE( "adaptive_sparse_jac", "[simple_adaptive]" ) {
    std::vector<realtype> p = {{298.15, 39390, -135.3, 18010, 44960, 48.2, 65919.5, -93.8304, 1780, 3790, 57.44, 19700, -157.4}};
    OdeSys odesys(&p[0]);
    std::vector<int> root_indices;

    const long int mxsteps=45000;
    const realtype dx0=0.0;
    const realtype dx_min=0.0;
    const realtype dx_max=0.0;
    const bool with_jacobian=true;
    cvodes_cxx::IterType iter_type=cvodes_cxx::IterType::Newton;
    const int maxl=0;
    const realtype eps_lin=0.0;
    const unsigned nderiv=0;
    bool return_on_root=false;

    realtype atol=1e-8, rtol=1e-8;

    cvodes_cxx::LinSol linear_solver=cvodes_cxx::LinSol::KLU;
    int autorestart=0;
    bool return_on_error = false;
    bool with_jtimes = false;
    realtype * xyout = (realtype*)malloc((odesys.get_ny() + 1)*sizeof(realtype));
    int td = 1;
    realtype tend=5.0;
    xyout[0] = 0; // t0
    xyout[1] = 8.99937e-07;
    xyout[2] = 0.000693731;
    xyout[3] = 0.000264211;
    xyout[4] = 0.000340312;
    xyout[5] = 4.11575e-05;

    // Test equivalence of dense and sparse jacs
    const indextype ny = odesys.get_ny();
    const int ldim = ny;
    const indextype nnz = odesys.get_nnz();
    indextype * colptrs = (indextype *) malloc((ny + 1) * sizeof(indextype));
    indextype * rowvals = (indextype *) malloc(nnz * sizeof(indextype));
    realtype * data = (realtype *) malloc(nnz * sizeof(realtype));
    realtype * J = (realtype *) malloc(ny * ny * sizeof(realtype));

    odesys.dense_jac_cmaj(0.0, &xyout[1], nullptr, J, ldim);
    odesys.sparse_jac_csc(0.0, &xyout[1], nullptr, data, colptrs, rowvals);

    for (int i=0; i < ny; i++) {
        for (int j = colptrs[i]; j < colptrs[i+1]; j++) {
            int k = rowvals[j];
            // both dense jac and sparse CSC matrix are column-major
            REQUIRE( J[ldim*i + k] == data[j]);
        }
    }

    auto nout = cvodes_anyode::simple_adaptive(&xyout, &td, &odesys, {atol}, rtol, cvodes_cxx::LMM::BDF, tend, root_indices,
                                               mxsteps, dx0, dx_min, dx_max, with_jacobian, iter_type, linear_solver,
                                               maxl, eps_lin, nderiv, return_on_root, autorestart, return_on_error,
                                               with_jtimes);
    REQUIRE( odesys.current_info.nfo_int["njev"] > 0 );
    free(colptrs);
    free(rowvals);
    free(data);
    free(J);
    free(xyout);
}
