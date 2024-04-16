#ifndef CFD_NORMS_ERROR_H_
#define CFD_NORMS_ERROR_H_

void computeL2(
    int nx, int ny, int nhc, 
    double***& Qn, double***& Qnp1, 
    double*& L2
);
void computeLinfinity(
    int nx, int ny, int nhc, 
    double***& Qn, double***& Qnp1, 
    double*& Linfty
);

#endif // CFD_NORMS_ERROR_H_ 
