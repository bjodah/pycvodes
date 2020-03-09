// C++11 source code.
#include "catch.hpp"
#include "cvodes_anyode_parallel.hpp"
#include "testing_utils.hpp"

TEST_CASE( "decay_adaptive", "[multi_adaptive]" ) {
    std::vector<realtype> k {{ 2.0, 3.0}};
    Decay odesys1(k[0]);
    Decay odesys2(k[1]);
    std::vector<Decay *> systems {{ &odesys1, &odesys2 }};
    const int nsys = systems.size();
    std::vector<realtype> y0 {{ 5.0, 7.0 }};
    std::vector<realtype> t0 {{ 1.0, 3.0 }};  // delta = 2
    std::vector<realtype> tend {{ 2.0, 5.0 }};  // delta = 3
    const int mxsteps = 0;
    std::vector<realtype> dx0 {{ 0, 0 }};
    std::vector<realtype> dx_min {{ 0, 0 }};
    std::vector<realtype> dx_max {{ 0, 0 }};
    realtype ** xyout_arr = (realtype **)malloc(nsys*sizeof(realtype*));
    int * td_arr = (int*)malloc(nsys*sizeof(int));
    for (int i=0; i<nsys; ++i){
        xyout_arr[i] = (realtype*)malloc((systems[i]->get_ny()+1)*sizeof(realtype));
        xyout_arr[i][0] = t0[i];
        xyout_arr[i][1] = y0[i];
        td_arr[i] = 1;
    }
    auto results = cvodes_anyode_parallel::multi_adaptive(
        xyout_arr, td_arr, systems, {1e-10}, 1e-10, cvodes_cxx::LMM::Adams, &tend[0],
        mxsteps, &dx0[0], &dx_min[0], &dx_max[0]);
    for (int idx=0; idx<nsys; ++idx){
        const auto& roots = results[idx].second;
        for (int j=0; j <= results[idx].first; ++j){
            REQUIRE( std::abs(y0[idx]*std::exp(t0[idx]-xyout_arr[idx][2*j])
                              - xyout_arr[idx][2*j+1]) < 1e-8 );
        }
        REQUIRE( systems[idx]->current_info.nfo_int["n_steps"] > 1 );
        REQUIRE( systems[idx]->current_info.nfo_int["n_steps"] < 997 );
        REQUIRE( systems[idx]->current_info.nfo_dbl["time_wall"] > 0 );
        REQUIRE( systems[idx]->current_info.nfo_dbl["time_wall"] < 100 );
        free(xyout_arr[idx]);
    }
    free(td_arr);
    free(xyout_arr);
}

TEST_CASE( "decay_predefined", "[multi_predefined]" ) {
    std::vector<realtype> k {{ 2.0, 3.0}};
    Decay odesys1(k[0]);
    Decay odesys2(k[1]);
    std::vector<Decay *> systems {{ &odesys1, &odesys2 }};
    std::vector<realtype> y0 {{ 5.0, 7.0 }};
    std::vector<realtype> yout(2*10);
    std::vector<realtype> tout(2*10);

    yout[0] = 1.0;
    yout[10] = 1.0;
    for (int i=0; i<2; ++i)
        for (int j=0; j<10; ++j)
            tout[10*i + j] = j/10.0;
    const int mxsteps = 0;
    std::vector<realtype> dx0 {{ 0, 0 }};
    std::vector<realtype> dx_min {{ 0, 0 }};
    std::vector<realtype> dx_max {{ 0, 0 }};

    auto result = cvodes_anyode_parallel::multi_predefined(
        systems, {1e-10}, 1e-10, cvodes_cxx::LMM::Adams, &y0[0], 10, &tout[0], &yout[0],
        mxsteps, &dx0[0], &dx_min[0], &dx_max[0]);
    REQUIRE(result.size() == 2);
    REQUIRE(result[0].first == 10);
    REQUIRE(result[1].first == 10);

    for (int idx=0; idx<2; ++idx){
        for (unsigned j=0; j<10; ++j){
            REQUIRE( std::abs(y0[idx]*std::exp(tout[10*idx + 0]-tout[10*idx + j]) - yout[10*idx + j]) < 1e-8 );
        }
        REQUIRE( systems[idx]->current_info.nfo_int["n_steps"] > 1 );
        REQUIRE( systems[idx]->current_info.nfo_int["n_steps"] < 997 );
        REQUIRE( systems[idx]->current_info.nfo_dbl["time_wall"] > 0 );
        REQUIRE( systems[idx]->current_info.nfo_dbl["time_wall"] < 100 );
    }
}
