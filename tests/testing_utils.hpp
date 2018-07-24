//#include <unordered_map>
#include "cvodes_anyode.hpp" // CvodesOdeSysBase, realtype, indextype
#include "anyode/anyode.hpp"

struct Decay : public AnyODE::CvodesOdeSysBase {
    realtype m_k;

    Decay(realtype k) : m_k(k) {}
    indextype get_ny() const override { return 1; }
    AnyODE::Status rhs(realtype t, const realtype * const ANYODE_RESTRICT y, realtype * const ANYODE_RESTRICT f) override {
        AnyODE::ignore(t);
        f[0] = -y[0];
        return AnyODE::Status::success;
    }
    realtype get_dx_max(realtype /* t */, const realtype * const /* y */) override{
        return 1e-3;
    }
};
