// -*- coding: utf-8 -*-
// -*- mode: read-only -*-
// This is a rendered source file (from template).

// User provided system description: CETSA (no Aggregate)
// Names of dependent variables: ['N', 'U', 'NL', 'L']

#pragma once
#include <math.h>
#include <vector>
#include "cvodes_cxx.hpp" // realtype, indextype

struct OdeSys : public AnyODE::CvodesOdeSysBase {
    std::vector<realtype> m_p;
    OdeSys(const realtype * const params);
    indextype get_ny() const override;
    AnyODE::Status rhs(realtype t,
                       const realtype * const ANYODE_RESTRICT y,
                       realtype * const ANYODE_RESTRICT f) override;
    AnyODE::Status dense_jac_cmaj(realtype t,
                                  const realtype * const ANYODE_RESTRICT y,
                                  const realtype * const ANYODE_RESTRICT fy,
                                  realtype * const ANYODE_RESTRICT jac,
                                  long int ldim,
                                  realtype * const ANYODE_RESTRICT dfdt=nullptr) override;
    AnyODE::Status dense_jac_rmaj(realtype t,
                                  const realtype * const ANYODE_RESTRICT y,
                                  const realtype * const ANYODE_RESTRICT fy,
                                  realtype * const ANYODE_RESTRICT jac,
                                  long int ldim,
                                  realtype * const ANYODE_RESTRICT dfdt=nullptr) override;
};

OdeSys::OdeSys(const realtype * const params) {
    m_p.assign(params, params + 13);
}
indextype OdeSys::get_ny() const {
    return 4;
}
AnyODE::Status OdeSys::rhs(realtype t,
                           const realtype * const ANYODE_RESTRICT y,
                           realtype * const ANYODE_RESTRICT f) {
    AnyODE::ignore(t);
    const realtype x0 = m_p[0]*y[0];
    const realtype x1 = 0.120272395808565/m_p[0];
    const realtype x2 = -m_p[9];
    const realtype x3 = 20836.6122251252*x0*y[3]*exp(x1*(m_p[0]*m_p[10] + x2));
    const realtype x4 = m_p[0] - 298.15;
    const realtype x5 = log(0.00335401643468053*m_p[0]);
    const realtype x6 = 20836612225.1252*m_p[0]*y[2]*exp(x1*(m_p[0]*(m_p[10] + m_p[12] + m_p[8]*x5) - m_p[11] - m_p[8]*x4 + x2));
    const realtype x7 = -x3 + x6;
    const realtype x8 = 20836612225.1252*m_p[0]*y[1];
    const realtype x9 = -m_p[1];
    const realtype x10 = x8*exp(x1*(m_p[0]*m_p[2] + x9));
    const realtype x11 = 20836612225.1252*x0*exp(x1*(m_p[0]*(m_p[2] + m_p[3]*x5 + m_p[4]/(m_p[5] + 273.15)) - m_p[3]*x4 - m_p[4] + x9));

    f[0] = x10 - x11 + x7;
    f[1] = -x10 + x11 - x8*exp(x1*(m_p[0]*m_p[7] - m_p[6]));
    f[2] = x3 - x6;
    f[3] = x7;
    this->nfev++;
    return AnyODE::Status::success;
}


AnyODE::Status OdeSys::dense_jac_cmaj(realtype t,
                                      const realtype * const ANYODE_RESTRICT y,
                                      const realtype * const ANYODE_RESTRICT fy,
                                      realtype * const ANYODE_RESTRICT jac,
                                      long int ldim,
                                      realtype * const ANYODE_RESTRICT dfdt) {
    // The AnyODE::ignore(...) calls below are used to generate code free from compiler warnings.
    AnyODE::ignore(fy);  // Currently we are not using fy (could be done through extensive pattern matching)
    AnyODE::ignore(t);


    const realtype x0 = 0.120272395808565/m_p[0];
    const realtype x1 = -m_p[9];
    const realtype x2 = 20836.6122251252*m_p[0]*exp(x0*(m_p[0]*m_p[10] + x1));
    const realtype x3 = x2*y[3];
    const realtype x4 = -x3;
    const realtype x5 = 20836612225.1252*m_p[0];
    const realtype x6 = -m_p[1];
    const realtype x7 = m_p[0] - 298.15;
    const realtype x8 = log(0.00335401643468053*m_p[0]);
    const realtype x9 = x5*exp(x0*(m_p[0]*(m_p[2] + m_p[3]*x8 + m_p[4]/(m_p[5] + 273.15)) - m_p[3]*x7 - m_p[4] + x6));
    const realtype x10 = x5*exp(x0*(m_p[0]*m_p[2] + x6));
    const realtype x11 = x5*exp(x0*(m_p[0]*(m_p[10] + m_p[12] + m_p[8]*x8) - m_p[11] - m_p[8]*x7 + x1));
    const realtype x12 = x2*y[0];
    const realtype x13 = -x12;

    jac[ldim*0 + 0] = x4 - x9;
    jac[ldim*0 + 1] = x9;
    jac[ldim*0 + 2] = x3;
    jac[ldim*0 + 3] = x4;

    jac[ldim*1 + 0] = x10;
    jac[ldim*1 + 1] = -x10 - x5*exp(x0*(m_p[0]*m_p[7] - m_p[6]));
    jac[ldim*1 + 2] = 0;
    jac[ldim*1 + 3] = 0;

    jac[ldim*2 + 0] = x11;
    jac[ldim*2 + 1] = 0;
    jac[ldim*2 + 2] = -x11;
    jac[ldim*2 + 3] = x11;

    jac[ldim*3 + 0] = x13;
    jac[ldim*3 + 1] = 0;
    jac[ldim*3 + 2] = x12;
    jac[ldim*3 + 3] = x13;

    if (dfdt){
        dfdt[0] = 0;
        dfdt[1] = 0;
        dfdt[2] = 0;
        dfdt[3] = 0;
    }
    this->njev++;
    return AnyODE::Status::success;
}

AnyODE::Status OdeSys::dense_jac_rmaj(realtype t,
                                      const realtype * const ANYODE_RESTRICT y,
                                      const realtype * const ANYODE_RESTRICT fy,
                                      realtype * const ANYODE_RESTRICT jac,
                                      long int ldim,
                                      realtype * const ANYODE_RESTRICT dfdt) {
    // The AnyODE::ignore(...) calls below are used to generate code free from compiler warnings.
    AnyODE::ignore(fy);  // Currently we are not using fy (could be done through extensive pattern matching)
    AnyODE::ignore(t);


    const realtype x0 = 0.120272395808565/m_p[0];
    const realtype x1 = -m_p[9];
    const realtype x2 = 20836.6122251252*m_p[0]*exp(x0*(m_p[0]*m_p[10] + x1));
    const realtype x3 = x2*y[3];
    const realtype x4 = -x3;
    const realtype x5 = 20836612225.1252*m_p[0];
    const realtype x6 = -m_p[1];
    const realtype x7 = m_p[0] - 298.15;
    const realtype x8 = log(0.00335401643468053*m_p[0]);
    const realtype x9 = x5*exp(x0*(m_p[0]*(m_p[2] + m_p[3]*x8 + m_p[4]/(m_p[5] + 273.15)) - m_p[3]*x7 - m_p[4] + x6));
    const realtype x10 = x5*exp(x0*(m_p[0]*m_p[2] + x6));
    const realtype x11 = x5*exp(x0*(m_p[0]*(m_p[10] + m_p[12] + m_p[8]*x8) - m_p[11] - m_p[8]*x7 + x1));
    const realtype x12 = x2*y[0];
    const realtype x13 = -x12;

    jac[ldim*0 + 0] = x4 - x9;
    jac[ldim*0 + 1] = x10;
    jac[ldim*0 + 2] = x11;
    jac[ldim*0 + 3] = x13;

    jac[ldim*1 + 0] = x9;
    jac[ldim*1 + 1] = -x10 - x5*exp(x0*(m_p[0]*m_p[7] - m_p[6]));
    jac[ldim*1 + 2] = 0;
    jac[ldim*1 + 3] = 0;

    jac[ldim*2 + 0] = x3;
    jac[ldim*2 + 1] = 0;
    jac[ldim*2 + 2] = -x11;
    jac[ldim*2 + 3] = x12;

    jac[ldim*3 + 0] = x4;
    jac[ldim*3 + 1] = 0;
    jac[ldim*3 + 2] = x11;
    jac[ldim*3 + 3] = x13;

    if (dfdt){
        dfdt[0] = 0;
        dfdt[1] = 0;
        dfdt[2] = 0;
        dfdt[3] = 0;
    }
    this->njev++;
    return AnyODE::Status::success;
}
