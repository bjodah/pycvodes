//#include <unordered_map>
#include "cvodes_cxx.hpp"
#include "anyode/anyode.hpp"

struct Decay : public AnyODE::OdeSysBase {
    double m_k;

    Decay(double k) : m_k(k) {}
    int get_ny() const override { return 1; }
    AnyODE::Status rhs(double t, const double * const __restrict__ y, double * const __restrict__ f) override {
        AnyODE::ignore(t);
        f[0] = -y[0];
        return AnyODE::Status::success;
    }
    double get_dx_max(double /* t */, const double * const /* y */){
        return 1e-3;
    }
};
