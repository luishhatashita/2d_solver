#include "error.h"

#include <cmath>

void computeL2(
    int nx, int ny, int nhc, 
    double ***&Qn, double ***&Qnp1, 
    double *&L2
)
{
    double L2l;
    for (int l=0; l<4; l++) {
        L2l = 0.0;
        for (int i=nhc; i<(nx+nhc-1); i++) {
            for (int j=nhc; j<(ny+nhc-1); j++) {
                L2l += (Qnp1[i][j][l]-Qn[i][j][l])*(Qnp1[i][j][l]-Qn[i][j][l]);
            }
        }
        L2[l] = sqrt(L2l);
    }
}

void computeLinfinity(
    int nx, int ny, int nhc, 
    double ***&Qn, double ***&Qnp1, 
    double *&Linfty
)
{
    double Linftyl;
    for (int l=0; l<4; l++) {
        Linftyl = 0.0;
        for (int i=nhc; i<(nx+nhc-1); i++) {
            for (int j=nhc; j<(ny+nhc-1); j++) {
                Linftyl = fmax(Linftyl, Qnp1[i][j][l]-Qn[i][j][l]);
            }
        }
        Linfty[l] = Linftyl;
    }
}
