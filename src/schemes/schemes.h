#ifndef CFD_SCHEMES_SCHEMES_H_
#define CFD_SCHEMES_SCHEMES_H_

#include "grid.h"
#include "parameters.h"
#include "solution.h"

// MUSCL interpolation
void interpolateMUSCL(
    double kappa, double epsilon,
    int nx, int ny, int nhc, 
    double***& Q,
    double***& QuL, double***& QuR,
    double***& QvL, double***& QvR
);
double fluxLimiter(double r);

// Shock capturing scheme
void computeFirstOrderUpwindFluxes(
    const Parameters& par,
    int nx, int ny, int nhc,
    double***& Su, double***& Sv, double**& V,
    double***& QuL, double***& QuR, double***& Qu, 
    double***& EhL, double***& EhR, double***& Eh,
    double***& QvL, double***& QvR, double***& Qv, 
    double***& FhL, double***& FhR, double***& Fh
);
void computeRoeFluxes(
    const Parameters& par,
    int nx, int ny, int nhc,
    double***& Su, double***& Sv, double**& V,
    double***& QuL, double***& QuR, double***& Qu, 
    double***& EhL, double***& EhR, double***& Eh,
    double***& QvL, double***& QvR, double***& Qv, 
    double***& FhL, double***& FhR, double***& Fh
);

// Time integration
void computeTimeSteps(Solution& solution);
void integrateTime(Solution& solution);

#endif // CFD_SCHEMES_SCHEMES_H_
