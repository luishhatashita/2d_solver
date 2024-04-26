#ifndef CFD_PARAMETERS_PARAMETERS_H_
#define CFD_PARAMETERS_PARAMETERS_H_

struct RunTime {
    int nrest, nit, nlog;
    double CFL, dt;
};

struct Refs {
    double pref, uref, Tref, lref, Mref;
};

struct BCs {
    int s, e, n, w;
};

struct Thermo {
    double cp, R, gamma, c;
};

struct MUSCL {
    double kappa, epsilon;
};

struct AUSM {
    double beta, kp, ku, sigma;
};

struct Parameters {
    RunTime rt; 
    Refs    ref;
    BCs     bcs;
    Thermo  td;
    MUSCL   muscl;
    AUSM    ausm;
};

//void sampleFunction();

#endif // CFD_PARAMETERS_PARAMETERS_H_
