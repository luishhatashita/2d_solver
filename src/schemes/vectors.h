#ifndef CFD_SCHEMES_VECTORS_H
#define CFD_SCHEMES_VECTORS_H

#include "parameters.h"
#include "grid.h"

enum BoundaryConditions {
    INFLOW_SUBSONIC    = 1,
    INFLOW_SUPERSONIC  = 2,
    OUTFLOW_SUBSONIC   = 3,
    OUTFLOW_SUPERSONIC = 4,
    SLIP_ADIABATIC     = 5,
    NOSLIP_ADIABATIC   = 6,
    SLIP_ISOTHERMAL    = 7,
    NOSLIP_ISOTHERMAL  = 8,
    SLIP_GRADIENT      = 9,
    NOSLIP_GRADIENT    = 10,
};

void initializeField(
    const Parameters& par, const Grid& grid,
    int nx, int ny, int nhc,
    double***& Qn, double***& Qvn
);
void computeBoundaryConditions(
    const Parameters& par, const Grid& grid,
    int nx, int ny, int nhc,
    double***& Qvn
);
void primitiveToConservative(
    const Parameters& par,
    int nx, int ny, int nhc,
    double**& V,
    double***& Qn, double***& Qvn
);
void conservativeToPrimitive(
    const Parameters& par,
    int nx, int ny, int nhc,
    double**& V,
    double***& Qn, double***& Qvn
);
void computeXiFaceQuantities(
    int nx, int ny, int nhc,
    double***& QuL, double***& QuR, double***& Qu 
);
void computeXiFlux(
    const Parameters& par,
    int nx, int ny, int nhc,
    double***& Su, double**& V,
    double***& Qu, double***& Eh
);
void computeXiDissipation(
    const Parameters& par,
    int nx, int ny, int nhc,
    double***& Su, double**& V,
    double***& QuL, double***& QuR, double***& Qu,
    double***& AQu
);
void computeEtaFaceQuantities(
    int nx, int ny, int nhc,
    double***& QvL, double***& QvR, double***& Qv 
);
void computeEtaFlux(
    const Parameters& par,
    int nx, int ny, int nhc,
    double***& Sv, double**& V,
    double***& Qv, double***& Fh
);
void computeEtaDissipation(
    const Parameters& par,
    int nx, int ny, int nhc,
    double***& Sv, double**& V,
    double***& QvL, double***& QvR, double***& Qv,
    double***& BQv
);

#endif
