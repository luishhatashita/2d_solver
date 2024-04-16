#ifndef CFD_NORMS_ERROR_H_
#define CFD_NORMS_ERROR_H_

#include "parameters.h"

void computeL2(
    const Parameters& par,
    int nx, int ny, int nhc, 
    double***& Qn, double***& Qnp1, 
    double*& L2
);
void computeLinfinity(
    const Parameters& par,
    int nx, int ny, int nhc, 
    double***& Qn, double***& Qnp1, 
    double*& Linfty
);

#endif // CFD_NORMS_ERROR_H_ 
