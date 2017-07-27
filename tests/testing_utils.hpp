//#include <unordered_map>
#include "cvodes_cxx.hpp"
#include "anyode/anyode.hpp"

template <typename T>
struct Decay : public AnyODE::OdeSysBase<T> {
    T m_k;

    Decay(T k) : m_k(k) {}
    int get_ny() const override { return 1; }
    AnyODE::Status rhs(T t, const T * const __restrict__ y, T * const __restrict__ f) override {
        AnyODE::ignore(t);
        f[0] = -y[0];
        return AnyODE::Status::success;
    }
    T get_dx_max(T /* t */, const T * const /* y */) override{
        return 1e-3;
    }
};
