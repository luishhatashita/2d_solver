#include "schemes.h"

#include <iostream>
#include <cmath>

#include "parameters.h"
#include "solution.h"
#include "vectors.h"

double fluxLimiter(double r);

/* Method to left and right states at cell faces based on MUSCL formulation.
 *
 * Parameters:
 * -----------
 *  double     kappa   : MUSCL parameter 1;
 *  double     epsilon : MUSCL parameter 2;
 *  int        nx      : size of xi dimension;
 *  int        ny      : size of eta dimension;
 *  int        nhc     : number of halo cells;
 *  double***& Q       : conservative variables at cell centers at n-th time step;
 *  double***& QuL     : left state interpolated value for face in the xi direction;
 *  double***& QuR     : right state interpolated value for face in the xi direction;
 *  double***& QvL     : left state interpolated value for face in the eta direction;
 *  double***& QvR     : right state interpolated value for face in the eta direction;
 */
void interpolateMUSCL(
    double kappa, double epsilon,
    int nx, int ny, int nhc, 
    double***& Q,
    double***& QuL, double***& QuR,
    double***& QvL, double***& QvR
)
{
    // xi direction:
    // a) high order stencils at the inner faces
    double r_L, r_R;
    for (int i=nhc+1; i<(nx+nhc-1)  ; i++) {
        for (int j=nhc; j<(ny+nhc-1); j++) {
            for (int l=0; l<4; l++) {
                r_L = (Q[i][j][l]-Q[i-1][j][l])/(Q[i-1][j][l]-Q[i-2][j][l]);
                QuL[i][j][l] = Q[i-1][j][l] \
                    + 0.25*epsilon*(
                        (1.0-kappa)*(Q[i-1][j][l]-Q[i-2][j][l])*fluxLimiter(r_L)
                       +(1.0+kappa)*(Q[i][j][l]  -Q[i-1][j][l])*fluxLimiter(1.0/r_L));
                r_R = (Q[i][j][l]-Q[i-1][j][l])/(Q[i+1][j][l]-Q[i][j][l]);
                QuR[i][j][l] = Q[i][j][l]; \
                    - 0.25*epsilon*(
                        (1.0+kappa)*(Q[i][j][l]  -Q[i-1][j][l])*fluxLimiter(1.0/r_R)
                       +(1.0-kappa)*(Q[i+1][j][l]-Q[i][j][l]  )*fluxLimiter(r_R));
            }
        }
    }
    // b) first order for boundaries:
    // Do averages here already?
    for (int j=nhc; j<(ny+nhc-1); j++) {
        for (int l=0; l<4; l++) {
            //// West boundary:
            //QuL[nhc][j][l]      = Q[nhc-1][j][l];
            //QuR[nhc][j][l]      = Q[nhc][j][l];
            //// East boundary:
            //QuL[nx+nhc-1][j][l] = Q[nx+nhc-2][j][l];
            //QuR[nx+nhc-1][j][l] = Q[nx+nhc-1][j][l];
            // West boundary:
            QuL[nhc][j][l]      = 0.5*(Q[nhc-1][j][l]+Q[nhc][j][l]);
            QuR[nhc][j][l]      = 0.5*(Q[nhc-1][j][l]+Q[nhc][j][l]);
            // East boundary:
            QuL[nx+nhc-1][j][l] = 0.5*(Q[nx+nhc-2][j][l]+Q[nx+nhc-1][j][l]);
            QuR[nx+nhc-1][j][l] = 0.5*(Q[nx+nhc-2][j][l]+Q[nx+nhc-1][j][l]);
        }
    }

    // eta direction:
    // a) high order stencils at the inner faces
    for (int i=nhc; i<(nx+nhc-1); i++) {
        for (int j=nhc+1; j<(ny+nhc-1)  ; j++) {
            for (int l=0; l<4; l++) {
                r_L = (Q[i][j][l]-Q[i][j-1][l])/(Q[i][j-1][l]-Q[i][j-2][l]);
                QvL[i][j][l] = Q[i][j-1][l] \
                    + 0.25*epsilon*(
                        (1.0-kappa)*(Q[i][j-1][l]-Q[i][j-2][l])*fluxLimiter(r_L)
                       +(1.0+kappa)*(Q[i][j][l]  -Q[i][j-1][l])*fluxLimiter(1.0/r_L));
                r_R = (Q[i][j][l]-Q[i][j-1][l])/(Q[i][j+1][l]-Q[i][j][l]);
                QvR[i][j][l] = Q[i][j][l] \
                    - 0.25*epsilon*(
                        (1.0+kappa)*(Q[i][j][l]  -Q[i][j-1][l])*fluxLimiter(1.0/r_R)
                       +(1.0-kappa)*(Q[i][j+1][l]-Q[i][j][l]  )*fluxLimiter(r_R));
            }
        }
    }
    // b) first order for boundaries:
    // Do averages here already?
    for (int i=nhc; i<(nx+nhc-1); i++) {
        for (int l=0; l<4; l++) {
            //// South boundary:
            //QvL[i][nhc][l]      = Q[i][nhc-1][l];
            //QvR[i][nhc][l]      = Q[i][nhc][l];
            //// East boundary:
            //QvL[i][ny+nhc-1][l] = Q[i][ny+nhc-2][l];
            //QvR[i][ny+nhc-1][l] = Q[i][ny+nhc-1][l];
            // South boundary:
            QvL[i][nhc][l]      = 0.5*(Q[i][nhc-1][l]+Q[i][nhc][l]);
            QvR[i][nhc][l]      = 0.5*(Q[i][nhc-1][l]+Q[i][nhc][l]);
            // East boundary:
            QvL[i][ny+nhc-1][l] = 0.5*(Q[i][ny+nhc-2][l]+Q[i][ny+nhc-1][l]);
            QvR[i][ny+nhc-1][l] = 0.5*(Q[i][ny+nhc-2][l]+Q[i][ny+nhc-1][l]);
        }
    }
}

double fluxLimiter(double r)
{
    double max, min;
    min = fmin(1.0, r);
    max = fmax(0.0, min);
    return max;
    //return 0.0; 
}

/* Method to manually upwind the fluxes, given that the flow comes from the left
 * and before the airfoil upwards and after downwards.
 *
 * Parameters:
 * -----------
 *  struct Parameters& par :
 *  int                nx      : size of xi dimension;
 *  int                ny      : size of eta dimension;
 *  int                nhc     : number of halo cells;
 *  double***&         Su      : projected cell face areas in the xi direction;
 *  double***&         Sv      : projected cell face areas in the eta direction;
 *  double***&         QuL     : left state interpolated value for face in the 
 *                               xi direction;
 *  double***&         QuR     : right state interpolated value for face in the 
 *                               xi direction;
 *  double***&         EhL     : xi direction fluxes computed based on the left 
 *                               state;
 *  double***&         EhR     : xi direction fluxes computed based on the right 
 *                               state;
 *  double***&         Eh      : xi direction fluxes upwinded;
 *  double***&         QvL     : left state interpolated value for face in the 
 *                               eta direction;
 *  double***&         QvR     : right state interpolated value for face in the 
 *                               eta direction;
 *  double***&         FhL     : eta direction fluxes computed based on the left 
 *                               state;
 *  double***&         FhR     : eta direction fluxes computed based on the right 
 *                               state;
 *  double***&         Fh      : eta direction fluxes upwinded;
 */
void computeFirstOrderUpwindFluxes(
    const Parameters& par,
    int nx, int ny, int nhc,
    double***& Su, double***& Sv, double**& V,
    double***& QuL, double***& QuR, double***& Qu, 
    double***& EhL, double***& EhR, double***& Eh,
    double***& QvL, double***& QvR, double***& Qv, 
    double***& FhL, double***& FhR, double***& Fh
)
{

    // Already constructing Left and Right states such that the boundaries are
    // the averages of the adjecent cells;
    // TODO:
    // -----
    //  - can remove some of the loop;
    // Manually upwind:
    // 1. xi fluxes:
    // a) interior faces:
    for (int i=nhc+1; i<(nx+nhc-1); i++) {
        for (int j=nhc; j<(ny+nhc-1); j++) {
            for (int l=0; l<4; l++) {
                //Eh[i][j][l] = EhL[i][j][l];
                Qu[i][j][l] = QuL[i][j][l];
            }
        }
    }
    // b) boundaries - central:
    for (int j=nhc; j<(ny+nhc-1); j++) {
        for (int l=0; l<4; l++) {
            //Eh[nhc][j][l]      = 0.5*(EhL[nhc][j][l]     +EhR[nhc][j][l]);
            //Eh[nx+nhc-1][j][l] = 0.5*(EhL[nx+nhc-1][j][l]+EhR[nx+nhc-1][j][l]);
            Qu[nhc][j][l]      = 0.5*(QuL[nhc][j][l]     +QuR[nhc][j][l]     );
            Qu[nx+nhc-1][j][l] = 0.5*(QuL[nx+nhc-1][j][l]+QuR[nx+nhc-1][j][l]);
        }
    }
    // 2. eta fluxes:
    // a) interior faces before airfoil: 
    /*
    int mid = (nx-1)/2;// + 5;
    for (int i=nhc; i<(mid+nhc); i++) {
        for (int j=nhc+1; j<(ny+nhc-1); j++) {
            for (int l=0; l<4; l++) {
                //Fh[i][j][l] = FhL[i][j][l];
                Qv[i][j][l] = QvL[i][j][l];
                //Qv[i][j][l] = QvR[i][j][l];
            }
        }
    }
    // a) interior faces after airfoil:
    for (int i=(mid+nhc); i<(nx+nhc-1); i++) {
        for (int j=nhc+1; j<(ny+nhc-1); j++) {
            for (int l=0; l<4; l++) {
                //Fh[i][j][l] = FhR[i][j][l];
                //Qv[i][j][l] = QvL[i][j][l];
                Qv[i][j][l] = QvR[i][j][l];
            }
        }
    }
    */
    ///* upwind from down to upwards
    //double u, v, S_et, v_et, b_v_et;
    for (int i=nhc; i<(nx+nhc-1); i++) {
        for (int j=nhc+1; j<(ny+nhc-1); j++) {
            for (int l=0; l<4; l++) {
                //Fh[i][j][l] = FhL[i][j][l];
                //Qv[i][j][l] = QvL[i][j][l];
                //Qv[i][j][l] = QvR[i][j][l];
                Qv[i][j][l] = 0.5*(QvL[i][j][l]+QvR[i][j][l]);
            }
            //u    = Qv[i][j][1]/Qv[i][j][0];
            //v    = Qv[i][j][2]/Qv[i][j][0];
            //S_et = sqrt(Sv[i][j][0]*Sv[i][j][0] + Sv[i][j][1]*Sv[i][j][1]);
            //v_et = u*Sv[i][j][0]/S_et + v*Sv[i][j][1]/S_et;
            //b_v_et = (v_et + fabs(v_et))/(2*fabs(v_et));
            //for (int l=0; l<4; l++) {
            //    //Qv[i][j][l] = b_v_et*QvL[i][j][l] + (1.0-b_v_et)*QvR[i][j][l];
            //    //Qv[i][j][l] = b_v_et*QvR[i][j][l] + (1.0-b_v_et)*QvL[i][j][l];
            //}
        }
    }
    //*/
    // b) boundaries:
    for (int i=nhc; i<(nx+nhc-1); i++) {
        for (int l=0; l<4; l++) {
            //Fh[i][nhc][l]      = 0.5*(FhL[i][nhc][l]     +FhR[i][nhc][l]);
            //Fh[i][ny+nhc-1][l] = 0.5*(FhL[i][ny+nhc-1][l]+FhR[i][ny+nhc-1][l]);
            Qv[i][nhc][l]      = 0.5*(QvL[i][nhc][l]     +QvR[i][nhc][l]     );
            Qv[i][ny+nhc-1][l] = 0.5*(QvL[i][ny+nhc-1][l]+QvR[i][ny+nhc-1][l]);
        }
    }
    
    computeXiFlux (par, nx, ny, nhc, Su, V, Qu, Eh);
    computeEtaFlux(par, nx, ny, nhc, Sv, V, Qv, Fh);
}

void computeRoeFluxes(
    const Parameters& par,
    int nx, int ny, int nhc,
    double***& Su, double***& Sv, double**& V,
    double***& QuL, double***& QuR, double***& Qu, 
    double***& EhL, double***& EhR, double***& Eh,
    double***& QvL, double***& QvR, double***& Qv, 
    double***& FhL, double***& FhR, double***& Fh,
    double***& AQu, double***& BQv
)
{
    // Compute left and right state fluxes:
    computeXiFlux (par, nx, ny, nhc, Su, V, QuL, EhL);
    computeXiFlux (par, nx, ny, nhc, Su, V, QuR, EhR);
    computeEtaFlux(par, nx, ny, nhc, Sv, V, QvL, FhL);
    computeEtaFlux(par, nx, ny, nhc, Sv, V, QvR, FhR);
    
    // Compute dissipation:
    computeXiDissipation(par, nx, ny, nhc, Su, V, QuL, QuR, Qu, AQu);
    computeEtaDissipation(par, nx, ny, nhc, Sv, V, QvL, QvR, Qv, BQv);

    // Compute fluxes:
    for (int i=nhc; i<(nx+nhc); i++) {
        for (int j=nhc; j<(ny+nhc-1); j++) {
            for (int l=0; l<4; l++) {
                Eh[i][j][l] = 0.5*(EhL[i][j][l]+EhR[i][j][l]) \
                            - 0.5*AQu[i][j][l];
            }
        }
    }
    for (int i=nhc; i<(nx+nhc-1); i++) {
        for (int j=nhc; j<(ny+nhc); j++) {
            for (int l=0; l<4; l++) {
                Fh[i][j][l] = 0.5*(FhL[i][j][l]+FhR[i][j][l]) \
                            - 0.5*BQv[i][j][l];
            }
        }
    }
}

void computeTimeSteps(Solution &solution)
{
}

void integrateTime(Solution &solution)
{
}
