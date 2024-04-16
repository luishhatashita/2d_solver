#ifndef CFD_PARAMETERS_PARAMETERS_H_
#define CFD_PARAMETERS_PARAMETERS_H_

struct RunTime {
    int nrest, nit;
    double CFL, dt;
};

struct Refs {
    double pref, uref, Tref, lref;
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

struct Parameters {
    RunTime rt; 
    Refs    ref;
    BCs     bcs;
    Thermo  td;
    MUSCL   muscl;
};

//void sampleFunction();

#endif // CFD_PARAMETERS_PARAMETERS_H_
