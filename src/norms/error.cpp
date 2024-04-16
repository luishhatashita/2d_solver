#include "error.h"

#include <cmath>

#include "parameters.h"

void computeL2(
    const Parameters& par, 
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

    // Normalize error values:
    double rho_ref;
    rho_ref = par.ref.pref/(par.td.R*par.ref.Tref);
    L2[0] = L2[0]/rho_ref;
    L2[1] = L2[1]/(rho_ref*par.ref.uref);
    L2[2] = L2[2]/(rho_ref*par.ref.uref);
    L2[3] = L2[3]/(rho_ref*par.td.c*par.td.c);
}

void computeLinfinity(
    const Parameters& par, 
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

    // Normalize error values:
    double rho_ref;
    rho_ref = par.ref.pref/(par.td.R*par.ref.Tref);
    Linfty[0] = Linfty[0]/rho_ref;
    Linfty[1] = Linfty[1]/(rho_ref*par.ref.uref);
    Linfty[2] = Linfty[2]/(rho_ref*par.ref.uref);
    Linfty[3] = Linfty[3]/(rho_ref*par.td.c*par.td.c);
}
