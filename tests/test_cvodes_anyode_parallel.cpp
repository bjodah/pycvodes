// C++11 source code.
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "catch.hpp"
#include "cvodes_anyode_parallel.hpp"
#include "testing_utils.hpp"

TEST_CASE( "decay_adaptive", "[multi_adaptive]" ) {
    std::vector<double> k {{ 2.0, 3.0}};
    Decay<double> odesys1(k[0]);
    Decay<double> odesys2(k[1]);
    std::vector<Decay<double> *> systems {{ &odesys1, &odesys2 }};
    std::vector<double> y0 {{ 5.0, 7.0 }};
    std::vector<double> t0 {{ 1.0, 3.0 }};  // delta = 2
    std::vector<double> tend {{ 2.0, 5.0 }};  // delta = 3
    const int mxsteps = 0;
    std::vector<double> dx0 {{ 0, 0 }};
    std::vector<double> dx_min {{ 0, 0 }};
    std::vector<double> dx_max {{ 0, 0 }};

    auto result = cvodes_anyode_parallel::multi_adaptive(
        systems, {1e-10}, 1e-10, cvodes_cxx::LMM::Adams, &y0[0], &t0[0], &tend[0],
        mxsteps, &dx0[0], &dx_min[0], &dx_max[0]);
    for (int idx=0; idx<2; ++idx){
        const auto& tout_yout_roots = result[idx];
        const auto& tout = tout_yout_roots.first.first;
        const auto& yout = tout_yout_roots.first.second;
        const auto& roots = tout_yout_roots.second;
        for (unsigned j=0; j<tout.size(); ++j){
            REQUIRE( std::abs(y0[idx]*std::exp(tout[0]-tout[j]) - yout[j]) < 1e-8 );
        }
        REQUIRE( systems[idx]->last_integration_info["n_steps"] > 1 );
        REQUIRE( systems[idx]->last_integration_info["n_steps"] < 997 );
        REQUIRE( systems[idx]->last_integration_info_dbl["time_wall"] > 0 );
        REQUIRE( systems[idx]->last_integration_info_dbl["time_wall"] < 100 );
    }
}

TEST_CASE( "decay_predefined", "[multi_predefined]" ) {
    std::vector<double> k {{ 2.0, 3.0}};
    Decay<double> odesys1(k[0]);
    Decay<double> odesys2(k[1]);
    std::vector<Decay<double> *> systems {{ &odesys1, &odesys2 }};
    std::vector<double> y0 {{ 5.0, 7.0 }};
    std::vector<double> yout(2*10);
    std::vector<double> tout(2*10);

    yout[0] = 1.0;
    yout[10] = 1.0;
    for (int i=0; i<2; ++i)
        for (int j=0; j<10; ++j)
            tout[10*i + j] = j/10.0;
    const int mxsteps = 0;
    std::vector<double> dx0 {{ 0, 0 }};
    std::vector<double> dx_min {{ 0, 0 }};
    std::vector<double> dx_max {{ 0, 0 }};

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
        REQUIRE( systems[idx]->last_integration_info["n_steps"] > 1 );
        REQUIRE( systems[idx]->last_integration_info["n_steps"] < 997 );
        REQUIRE( systems[idx]->last_integration_info_dbl["time_wall"] > 0 );
        REQUIRE( systems[idx]->last_integration_info_dbl["time_wall"] < 100 );
    }
}
