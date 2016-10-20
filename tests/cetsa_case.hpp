#pragma once
// This is a real-world based test example

struct OdeSys : public AnyODE::OdeSysBase {
    std::vector<double> m_p;
    OdeSys(const double * const params);
    int get_ny() const override;
    AnyODE::Status rhs(double t,
                       const double * const __restrict__ y,
                       double * const __restrict__ f) override;
    AnyODE::Status dense_jac_cmaj(double t,
                                  const double * const __restrict__ y,
                                  const double * const __restrict__ fy,
                                  double * const __restrict__ jac,
                                  long int ldim,
                                  double * const __restrict__ dfdt=nullptr) override;
    AnyODE::Status dense_jac_rmaj(double t,
                                  const double * const __restrict__ y,
                                  const double * const __restrict__ fy,
                                  double * const __restrict__ jac,
                                  long int ldim,
                                  double * const __restrict__ dfdt=nullptr) override;
};


OdeSys::OdeSys(const double * const params) {
    m_p.assign(params, params + 13);
}
int OdeSys::get_ny() const {
    return 5;
}
AnyODE::Status OdeSys::rhs(double t,
                           const double * const __restrict__ y,
                           double * const __restrict__ f) {
    AnyODE::ignore(t);
    const double x0 = m_p[0]*y[0];
    const double x1 = 0.120272395808565/m_p[0];
    const double x2 = -m_p[9];
    const double x3 = 20836.6122251252*x0*y[3]*exp(x1*(m_p[0]*m_p[10] + x2));
    const double x4 = m_p[0] - 298.15;
    const double x5 = log(0.00335401643468053*m_p[0]);
    const double x6 = 20836612225.1252*m_p[0]*y[4]*exp(x1*(m_p[0]*(m_p[10] + m_p[12] + m_p[8]*x5) - m_p[11] - m_p[8]*x4 + x2));
    const double x7 = -x3 + x6;
    const double x8 = 20836612225.1252*m_p[0]*y[1];
    const double x9 = -m_p[1];
    const double x10 = x8*exp(x1*(m_p[0]*m_p[2] + x9));
    const double x11 = 20836612225.1252*x0*exp(x1*(m_p[0]*(m_p[2] + m_p[3]*x5 + m_p[4]/(m_p[5] + 273.15)) - m_p[3]*x4 - m_p[4] + x9));
    const double x12 = x8*exp(x1*(m_p[0]*m_p[7] - m_p[6]));

    f[0] = x10 - x11 + x7;
    f[1] = -x10 + x11 - x12;
    f[2] = x12;
    f[3] = x7;
    f[4] = x3 - x6;
    this->nfev++;
    return AnyODE::Status::success;
}


AnyODE::Status OdeSys::dense_jac_cmaj(double t,
                                      const double * const __restrict__ y,
                                      const double * const __restrict__ fy,
                                      double * const __restrict__ jac,
                                      long int ldim,
                                      double * const __restrict__ dfdt) {
    // The AnyODE::ignore(...) calls below are used to generate code free from compiler warnings.
    AnyODE::ignore(fy);  // Currently we are not using fy (could be done through extensive pattern matching)
    AnyODE::ignore(t);


    const double x0 = 0.120272395808565/m_p[0];
    const double x1 = -m_p[9];
    const double x2 = 20836.6122251252*m_p[0]*exp(x0*(m_p[0]*m_p[10] + x1));
    const double x3 = x2*y[3];
    const double x4 = -x3;
    const double x5 = 20836612225.1252*m_p[0];
    const double x6 = -m_p[1];
    const double x7 = m_p[0] - 298.15;
    const double x8 = log(0.00335401643468053*m_p[0]);
    const double x9 = x5*exp(x0*(m_p[0]*(m_p[2] + m_p[3]*x8 + m_p[4]/(m_p[5] + 273.15)) - m_p[3]*x7 - m_p[4] + x6));
    const double x10 = x5*exp(x0*(m_p[0]*m_p[2] + x6));
    const double x11 = x2*y[0];
    const double x12 = -x11;
    const double x13 = x5*exp(x0*(m_p[0]*(m_p[10] + m_p[12] + m_p[8]*x8) - m_p[11] - m_p[8]*x7 + x1));
    const double x14 = x5*exp(x0*(m_p[0]*m_p[7] - m_p[6]));

    jac[ldim*0 + 0] = x4 - x9;
    jac[ldim*0 + 1] = x9;
    jac[ldim*0 + 2] = 0;
    jac[ldim*0 + 3] = x4;
    jac[ldim*0 + 4] = x3;

    jac[ldim*1 + 0] = x10;
    jac[ldim*1 + 1] = -x10 - x14;
    jac[ldim*1 + 2] = x14;
    jac[ldim*1 + 3] = 0;
    jac[ldim*1 + 4] = 0;

    jac[ldim*2 + 0] = 0;
    jac[ldim*2 + 1] = 0;
    jac[ldim*2 + 2] = 0;
    jac[ldim*2 + 3] = 0;
    jac[ldim*2 + 4] = 0;

    jac[ldim*3 + 0] = x12;
    jac[ldim*3 + 1] = 0;
    jac[ldim*3 + 2] = 0;
    jac[ldim*3 + 3] = x12;
    jac[ldim*3 + 4] = x11;

    jac[ldim*4 + 0] = x13;
    jac[ldim*4 + 1] = 0;
    jac[ldim*4 + 2] = 0;
    jac[ldim*4 + 3] = x13;
    jac[ldim*4 + 4] = -x13;

    if (dfdt){
        dfdt[0] = 0;
        dfdt[1] = 0;
        dfdt[2] = 0;
        dfdt[3] = 0;
        dfdt[4] = 0;
    }
    this->njev++;
    return AnyODE::Status::success;
}

AnyODE::Status OdeSys::dense_jac_rmaj(double t,
                                      const double * const __restrict__ y,
                                      const double * const __restrict__ fy,
                                      double * const __restrict__ jac,
                                      long int ldim,
                                      double * const __restrict__ dfdt) {
    // The AnyODE::ignore(...) calls below are used to generate code free from compiler warnings.
    AnyODE::ignore(fy);  // Currently we are not using fy (could be done through extensive pattern matching)
    AnyODE::ignore(t);


    const double x0 = 0.120272395808565/m_p[0];
    const double x1 = -m_p[9];
    const double x2 = 20836.6122251252*m_p[0]*exp(x0*(m_p[0]*m_p[10] + x1));
    const double x3 = x2*y[3];
    const double x4 = -x3;
    const double x5 = 20836612225.1252*m_p[0];
    const double x6 = -m_p[1];
    const double x7 = m_p[0] - 298.15;
    const double x8 = log(0.00335401643468053*m_p[0]);
    const double x9 = x5*exp(x0*(m_p[0]*(m_p[2] + m_p[3]*x8 + m_p[4]/(m_p[5] + 273.15)) - m_p[3]*x7 - m_p[4] + x6));
    const double x10 = x5*exp(x0*(m_p[0]*m_p[2] + x6));
    const double x11 = x2*y[0];
    const double x12 = -x11;
    const double x13 = x5*exp(x0*(m_p[0]*(m_p[10] + m_p[12] + m_p[8]*x8) - m_p[11] - m_p[8]*x7 + x1));
    const double x14 = x5*exp(x0*(m_p[0]*m_p[7] - m_p[6]));

    jac[ldim*0 + 0] = x4 - x9;
    jac[ldim*0 + 1] = x10;
    jac[ldim*0 + 2] = 0;
    jac[ldim*0 + 3] = x12;
    jac[ldim*0 + 4] = x13;

    jac[ldim*1 + 0] = x9;
    jac[ldim*1 + 1] = -x10 - x14;
    jac[ldim*1 + 2] = 0;
    jac[ldim*1 + 3] = 0;
    jac[ldim*1 + 4] = 0;

    jac[ldim*2 + 0] = 0;
    jac[ldim*2 + 1] = x14;
    jac[ldim*2 + 2] = 0;
    jac[ldim*2 + 3] = 0;
    jac[ldim*2 + 4] = 0;

    jac[ldim*3 + 0] = x4;
    jac[ldim*3 + 1] = 0;
    jac[ldim*3 + 2] = 0;
    jac[ldim*3 + 3] = x12;
    jac[ldim*3 + 4] = x13;

    jac[ldim*4 + 0] = x3;
    jac[ldim*4 + 1] = 0;
    jac[ldim*4 + 2] = 0;
    jac[ldim*4 + 3] = x11;
    jac[ldim*4 + 4] = -x13;

    if (dfdt){
        dfdt[0] = 0;
        dfdt[1] = 0;
        dfdt[2] = 0;
        dfdt[3] = 0;
        dfdt[4] = 0;
    }
    this->njev++;
    return AnyODE::Status::success;
}
