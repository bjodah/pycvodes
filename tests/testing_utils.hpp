#include "cvodes_cxx.hpp"

int rhs_cb(double t, N_Vector y, N_Vector f, void * user_data){
    cvodes_cxx::ignore(t); cvodes_cxx::ignore(user_data);
    NV_DATA_S(f)[0] = -NV_DATA_S(y)[0];
    return 0;
}

struct Decay : public cvodes_cxx::OdeSysBase {
    double m_k;

    Decay(double k) : m_k(k) {}
    int get_ny() { return 1; }
    void rhs(double t, const double * const y, double * const f) override {
        cvodes_cxx::ignore(t);
        f[0] = -y[0];
    }
};
