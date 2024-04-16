#include "vectors.h"

#include <iostream>
#include <string>
#include <cmath>

#include "parameters.h"
#include "allocate.h"
#include "grid.h"
#include "writers.h"

/* Method to inialize the flow field, based on .
 *
 * Parameters:
 * -----------
 *  struct Parameters& par   : const reference to struct of parameters;
 *  class  Grid&       grid  : const reference to struct of parameters;
 *  int                nx    : number of nodes in xi direction;
 *  int                ny    : number of nodes in eta direction;
 *  int                nhc   : number of halo cells;
 *  double***&         Qn    : reference to conservative variable vector at the
 *                             n-th time step;
 *  double***&         Qvn   : reference to primitive  variable vector at the
 *                             n-th time step;
 *  double***&         Qvnp1 : reference to primitive  variable vector at the
 *                             n+1-th time step;
 */
void initializeField(
    const Parameters& par, const Grid& grid,
    int nx, int ny, int nhc, 
    double***& Qn, double***& Qvn, double***& Qvnp1
)
{
    std::cout << "--Flow field initialization:" << std::endl;
    // Loop inside grid:
    for (int i=nhc; i<(nx+nhc-1); i++) {
        for (int j=nhc; j<(ny+nhc-1); j++) {
            Qvn[i][j][0] = par.ref.pref;
            Qvn[i][j][1] = par.ref.uref;
            Qvn[i][j][2] = 0;
            Qvn[i][j][3] = par.ref.Tref;
        }
    }
    std::cout << "----Inner cells computed." << std::endl;
    // Copy to Qvnp1 to use the updated computeBCs;
    for (int i=nhc; i<(nx+nhc-1); i++) {
        for (int j=nhc; j<(ny+nhc-1); j++) {
            for (int l=0; l<4; l++) {
                Qvnp1[i][j][l] = Qvn[i][j][l];
            }
        }
    }
    computeBoundaryConditions(par, grid, nx, ny, nhc, Qvn, Qvnp1);
    std::cout << "----Boundary conditions computed." << std::endl;
    //std::string wfpath = "./out/rest/Qv00000_original.bin";
    //writeBinary3DArray(wfpath, nx+2*nhc-1, ny+2*nhc-1, 4, Qvn);

    double** V;
    V = grid.getCellVolumes();
    primitiveToConservative(par, nx, ny, nhc, V, Qn, Qvnp1);
    //std::cout << "----Conservative vector converted." << std::endl;
    //wfpath = "./out/rest/Q00000_original.bin";
    //writeBinary3DArray(wfpath, nx+2*nhc-1, ny+2*nhc-1, 4, Qn);

    //conservativeToPrimitive(par, nx, ny, nhc, V, Qn, Qvn);
    //wfpath = "./out/rest/Qv00000_converted.bin";
    //writeBinary3DArray(wfpath, nx+2*nhc-1, ny+2*nhc-1, 4, Qvn);
}

/* Method to compute boundary conditions in primitive variable space based on
 * the predefined parameters.
 *
 * Parameters:
 * -----------
 *  struct Parameters& par  : const reference to struct of parameters;
 *  class  Grid&       grid : const reference to grid class with metrics;
 *  int                nx   : number of nodes in xi direction;
 *  int                ny   : number of nodes in eta direction;
 *  int                nhc  : number of halo cells;
 *  double***&         Qvn  : reference to primitive variable vector at the
 *                            n-th time step;
 */
void computeBoundaryConditions(
    const Parameters &par, const Grid& grid,
    int nx, int ny, int nhc, 
    double***& Qvn, double***& Qvnp1
)
{
    double ***su, ***sv;
    su = grid.getProjectedFaceAreasXi(); 
    sv = grid.getProjectedFaceAreasEta(); 

    // South BC:
    double svx, svy, svx2, svy2;
    switch((BoundaryConditions)par.bcs.s) {
        case INFLOW_SUBSONIC: break;
        case INFLOW_SUPERSONIC: break;
        case OUTFLOW_SUBSONIC: break;
        case OUTFLOW_SUPERSONIC: break;
        case SLIP_ADIABATIC:
            for (int j=0; j<nhc; j++) {
                for (int i=nhc; i<(nx+nhc-1); i++) {
                    Qvnp1[i][(nhc-1)-j][0] = Qvnp1[i][nhc+j][0];
                    svx = sv[i][nhc][0];
                    svy = sv[i][nhc][1];
                    svx2 = svx*svx;
                    svy2 = svy*svy;
                    Qvnp1[i][(nhc-1)-j][1] = ((svy2-svx2)*Qvnp1[i][nhc+j][1] \
                                            -2*svx*svy*Qvnp1[i][nhc+j][2]) \
                                            / (svx2 + svy2);
                    Qvnp1[i][(nhc-1)-j][2] = (-2*svx*svy*Qvnp1[i][nhc+j][1]    \
                                            +(svx2-svy2)*Qvnp1[i][nhc+j][2]) \
                                            / (svx2 + svy2);
                    Qvnp1[i][(nhc-1)-j][3] = Qvnp1[i][nhc+j][3];
                }
            }
            break;
        case NOSLIP_ADIABATIC:
            for (int j=0; j<nhc; j++) {
                for (int i=nhc; i<(nx+nhc-1); i++) {
                    Qvnp1[i][(nhc-1)-j][0] =  Qvnp1[i][nhc+j][0];
                    Qvnp1[i][(nhc-1)-j][1] = -Qvnp1[i][nhc+j][1];
                    Qvnp1[i][(nhc-1)-j][2] = -Qvnp1[i][nhc+j][2];
                    Qvnp1[i][(nhc-1)-j][3] =  Qvnp1[i][nhc+j][3];
                }
            }
            break;
        case SLIP_ISOTHERMAL: break;
        case NOSLIP_ISOTHERMAL: break;
        case SLIP_GRADIENT: break;
        case NOSLIP_GRADIENT: break;
        default:
            std::cout << "Invalid South Boundary Condition." << std::endl;
    }
    //std::cout << "----South boundary condition computed." << std::endl;

    // East BC:
    switch((BoundaryConditions)par.bcs.e) {
        case INFLOW_SUBSONIC: break;
        case INFLOW_SUPERSONIC: break;
        case OUTFLOW_SUBSONIC: break;
        case OUTFLOW_SUPERSONIC:
            for (int i=0; i<nhc; i++) {
                for (int j=nhc; j<(ny+nhc-1); j++) {
                    Qvnp1[(nx+nhc-1)+i][j][0] = Qvn[(nx+nhc-2)-i][j][0];
                    Qvnp1[(nx+nhc-1)+i][j][1] = Qvn[(nx+nhc-2)-i][j][1];
                    Qvnp1[(nx+nhc-1)+i][j][2] = Qvn[(nx+nhc-2)-i][j][2];
                    Qvnp1[(nx+nhc-1)+i][j][3] = Qvn[(nx+nhc-2)-i][j][3];
                }
            }
            break;
        case SLIP_ADIABATIC:
            //for (int j=0; j<nhc; j++) {
            //    for (int i=nhc; i<(nx+nhc-1); i++) {
            //        Qvn[i][(nhc-1)-j][0] = Qvn[i][nhc+j][0];
            //        svx = sv[i][nhc][0];
            //        svy = sv[i][nhc][1];
            //        svx2 = svx*svx;
            //        svy2 = svy*svy;
            //        Qvn[i][(nhc-1)-j][1] = ((svy2-svx2)*Qvn[i][nhc+j][1] \
            //                                -2*svx*svy*Qvn[i][nhc+j][2]) \
            //                                / (svx2 + svy2);
            //        Qvn[i][(nhc-1)-j][2] = (-2*svx*svy*Qvn[i][nhc+j][1]    \
            //                                +(svx2-svy2)*Qvn[i][nhc+j][2]) \
            //                                / (svx2 + svy2);
            //        Qvn[i][(nhc-1)-j][3] = Qvn[i][nhc+j][3];
            //    }
            //}
            break;
        case NOSLIP_ADIABATIC:
            //for (int j=0; j<nhc; j++) {
            //    for (int i=nhc; i<(nx+nhc-1); i++) {
            //        Qvn[i][(nhc-1)-j][0] =  Qvn[i][nhc+j][0];
            //        Qvn[i][(nhc-1)-j][1] = -Qvn[i][nhc+j][1];
            //        Qvn[i][(nhc-1)-j][2] = -Qvn[i][nhc+j][2];
            //        Qvn[i][(nhc-1)-j][3] =  Qvn[i][nhc+j][3];
            //    }
            //}
            break;
        case SLIP_ISOTHERMAL: break;
        case NOSLIP_ISOTHERMAL: break;
        case SLIP_GRADIENT: break;
        case NOSLIP_GRADIENT: break;
        default:
            std::cout << "Invalid East Boundary Condition." << std::endl;
    }
    //std::cout << "----East boundary condition computed." << std::endl;
    
    // North BC:
    switch((BoundaryConditions)par.bcs.n) {
        case INFLOW_SUBSONIC: break;
        case INFLOW_SUPERSONIC: break;
        case OUTFLOW_SUBSONIC: break;
        case OUTFLOW_SUPERSONIC:
            //for (int i=0; i<nhc; i++) {
            //    for (int j=nhc; j<(ny+nhc-1); j++) {
            //        Qvn[(nx+nhc)+i][j][0] = Qvn[(nx+nhc-1)-i][j][0];
            //        Qvn[(nx+nhc)+i][j][1] = Qvn[(nx+nhc-1)-i][j][1];
            //        Qvn[(nx+nhc)+i][j][2] = Qvn[(nx+nhc-1)-i][j][2];
            //        Qvn[(nx+nhc)+i][j][3] = Qvn[(nx+nhc-1)-i][j][3];
            //    }
            //}
            break;
        case SLIP_ADIABATIC:
            for (int j=0; j<nhc; j++) {
                for (int i=nhc; i<(nx+nhc-1); i++) {
                    Qvnp1[i][(ny+nhc-1)+j][0] = Qvnp1[i][(ny+nhc-2)-j][0];
                    svx = sv[i][ny+nhc][0];
                    svy = sv[i][ny+nhc][1];
                    svx2 = svx*svx;
                    svy2 = svy*svy;
                    Qvnp1[i][(ny+nhc-1)+j][1] = ((svy2-svx2)*Qvnp1[i][(ny+nhc-2)-j][1] \
                                            -2*svx*svy*Qvnp1[i][(ny+nhc-2)-j][2]) \
                                            / (svx2 + svy2);
                    Qvnp1[i][(ny+nhc-1)+j][2] = (-2*svx*svy*Qvnp1[i][(ny+nhc-2)-j][1]    \
                                            +(svx2-svy2)*Qvnp1[i][(ny+nhc-2)-j][2]) \
                                            / (svx2 + svy2);
                    Qvnp1[i][(ny+nhc-1)+j][3] = Qvnp1[i][(ny+nhc-2)-j][3];
                }
            }
            break;
        case NOSLIP_ADIABATIC:
            for (int j=0; j<nhc; j++) {
                for (int i=nhc; i<(nx+nhc-1); i++) {
                    Qvnp1[i][(ny+nhc)+j][0] =  Qvnp1[i][(ny+nhc-1)-j][0];
                    Qvnp1[i][(ny+nhc)+j][1] = -Qvnp1[i][(ny+nhc-1)-j][1];
                    Qvnp1[i][(ny+nhc)+j][2] = -Qvnp1[i][(ny+nhc-1)-j][2];
                    Qvnp1[i][(ny+nhc)+j][3] =  Qvnp1[i][(ny+nhc-1)-j][3];
                }
            }
            break;
        case SLIP_ISOTHERMAL: break;
        case NOSLIP_ISOTHERMAL: break;
        case SLIP_GRADIENT: break;
        case NOSLIP_GRADIENT: break;
        default:
            std::cout << "Invalid North Boundary Condition." << std::endl;
    }
    //std::cout << "----North boundary condition computed." << std::endl;
    
    // West BC:
    switch((BoundaryConditions)par.bcs.w) {
        case INFLOW_SUBSONIC: break;
        case INFLOW_SUPERSONIC:
            for (int i=0; i<nhc; i++) {
                for (int j=nhc; j<(ny+nhc-1); j++) {
                    Qvnp1[(nhc-1)-i][j][0] = par.ref.pref;
                    Qvnp1[(nhc-1)-i][j][1] = par.ref.uref;
                    Qvnp1[(nhc-1)-i][j][2] = 0.0;
                    Qvnp1[(nhc-1)-i][j][3] = par.ref.Tref;
                }
            }
            break;
        case OUTFLOW_SUBSONIC: break;
        case OUTFLOW_SUPERSONIC:
            //for (int i=0; i<nhc; i++) {
            //    for (int j=nhc; j<(ny+nhc-1); j++) {
            //        Qvn[(nx+nhc)+i][j][0] = Qvn[(nx+nhc-1)-i][j][0];
            //        Qvn[(nx+nhc)+i][j][1] = Qvn[(nx+nhc-1)-i][j][1];
            //        Qvn[(nx+nhc)+i][j][2] = Qvn[(nx+nhc-1)-i][j][2];
            //        Qvn[(nx+nhc)+i][j][3] = Qvn[(nx+nhc-1)-i][j][3];
            //    }
            //}
            break;
        case SLIP_ADIABATIC:
            //for (int j=0; j<nhc; j++) {
            //    for (int i=nhc; i<(nx+nhc-1); i++) {
            //        Qvn[i][(ny+nhc)+j][0] = Qvn[i][(ny+nhc-1)-j][0];
            //        svx = sv[i][ny+nhc][0];
            //        svy = sv[i][ny+nhc][1];
            //        svx2 = svx*svx;
            //        svy2 = svy*svy;
            //        Qvn[i][(ny+nhc)+j][1] = ((svy2-svx2)*Qvn[i][(ny+nhc-1)-j][1] \
            //                                -2*svx*svy*Qvn[i][(ny+nhc-1)-j][2]) \
            //                                / (svx2 + svy2);
            //        Qvn[i][(ny+nhc)+j][2] = (-2*svx*svy*Qvn[i][(ny+nhc-1)-j][1]    \
            //                                +(svx2-svy2)*Qvn[i][(ny+nhc-1)-j][2]) \
            //                                / (svx2 + svy2);
            //        Qvn[i][(ny+nhc)+j][3] = Qvn[i][(ny+nhc-1)-j][3];
            //    }
            //}
            break;
        case NOSLIP_ADIABATIC:
            //for (int j=0; j<nhc; j++) {
            //    for (int i=nhc; i<(nx+nhc-1); i++) {
            //        Qvn[i][(ny+nhc)+j][0] =  Qvn[i][(ny+nhc-1)-j][0];
            //        Qvn[i][(ny+nhc)+j][1] = -Qvn[i][(ny+nhc-1)-j][1];
            //        Qvn[i][(ny+nhc)+j][2] = -Qvn[i][(ny+nhc-1)-j][2];
            //        Qvn[i][(ny+nhc)+j][3] =  Qvn[i][(ny+nhc-1)-j][3];
            //    }
            //}
            break;
        case SLIP_ISOTHERMAL: break;
        case NOSLIP_ISOTHERMAL: break;
        case SLIP_GRADIENT: break;
        case NOSLIP_GRADIENT: break;
        default:
            std::cout << "Invalid West Boundary Condition." << std::endl;
    }
    //std::cout << "----West boundary condition computed." << std::endl;
}

/* Method to convert primitive variable vector into conservative variables.
 *
 * Parameters:
 * -----------
 *  struct Parameters& par  : const reference to struct of parameters;
 *  int                nx   : number of nodes in xi direction;
 *  int                ny   : number of nodes in eta direction;
 *  int                nhc  : number of halo cells;
 *  double***&         Qn   : reference to conservative variable vector at the
 *                            n-th time step;
 *  double***&         Qvn  : reference to primitive variable vector at the
 *                            n-th time step;
 */
void primitiveToConservative(
    const Parameters& par, 
    int nx, int ny, int nhc, 
    double **&V,
    double ***&Qn, double ***&Qvn
)
{
    double rho, v;
    for (int i=0; i<(nx+2*nhc-1); i++) {
        for (int j=0; j<(ny+2*nhc-1); j++) {
            //v   = V[i][j];
            rho = Qvn[i][j][0]/(par.td.R*Qvn[i][j][3]); // p/RT
            Qn[i][j][0] = rho;
            Qn[i][j][1] = rho*Qvn[i][j][1];
            Qn[i][j][2] = rho*Qvn[i][j][2];
            Qn[i][j][3] = Qvn[i][j][0]/(par.td.gamma - 1.0) \
                        + 0.5*rho*(Qvn[i][j][1]*Qvn[i][j][1] + Qvn[i][j][2]*Qvn[i][j][2]);
            //Qhn[i][j][0] = v*rho;
            //Qhn[i][j][1] = v*rho*Qvn[i][j][1];
            //Qhn[i][j][2] = v*rho*Qvn[i][j][2];
            //Qhn[i][j][3] = v*(Qvn[i][j][0]/(par.td.gamma - 1.0) \
            //            + 0.5*rho*(Qvn[i][j][1]*Qvn[i][j][1] + Qvn[i][j][2]*Qvn[i][j][2]));
        }
    }
}

/* Method to convert conservative variable vector into primitive variables.
 *
 * Parameters:
 * -----------
 *  struct Parameters& par  : const reference to struct of parameters;
 *  int                nx   : number of nodes in xi direction;
 *  int                ny   : number of nodes in eta direction;
 *  int                nhc  : number of halo cells;
 *  double***&         Qn   : reference to conservative variable vector at the
 *                            n-th time step;
 *  double***&         Qvn  : reference to primitive variable vector at the
 *                            n-th time step;
 */
void conservativeToPrimitive(
    const Parameters& par, 
    int nx, int ny, int nhc, 
    double**& V,
    double***& Qn, double***& Qvn
)
{
    double p, v;
    for (int i=0; i<(nx+2*nhc-1); i++) {
        for (int j=0; j<(ny+2*nhc-1); j++) {
            //v            = V[i][j];
            p            = (Qn[i][j][3] \
                           - 0.5*(Qn[i][j][1]*Qn[i][j][1]/Qn[i][j][0] \
                                 +Qn[i][j][2]*Qn[i][j][2]/Qn[i][j][0]))*(par.td.gamma-1.0);
            Qvn[i][j][0] = p;
            Qvn[i][j][1] = Qn[i][j][1]/Qn[i][j][0];
            Qvn[i][j][2] = Qn[i][j][2]/Qn[i][j][0];
            Qvn[i][j][3] = p/(Qn[i][j][0]*par.td.R);
            //p            = (Qhn[i][j][3]/v \
            //               - 0.5*(Qhn[i][j][1]*Qhn[i][j][1]/Qhn[i][j][0]/v \
            //                     +Qhn[i][j][2]*Qhn[i][j][2]/Qhn[i][j][0]/v))*(par.td.gamma-1.0);
            //Qvn[i][j][0] = p;
            //Qvn[i][j][1] = Qhn[i][j][1]/Qhn[i][j][0];
            //Qvn[i][j][2] = Qhn[i][j][2]/Qhn[i][j][0];
            //Qvn[i][j][3] = p/(Qhn[i][j][0]*par.td.R/v);
        }
    }
}

/* Compute quantities at xi faces within and include primary grid faces.
 *
 * Parameters:
 * -----------
 *  int        nx  : number of nodes in xi direction;
 *  int        ny  : number of nodes in eta direction;
 *  int        nhc : number of halo cells;
 *  double***& QuL : conservative variables in the xi cell faces from the left 
 *                   state;
 *  double***& QuR : conservative variables in the xi cell faces from the right
 *                   state;
 *  double***& Qu  : conservative variables in the xi cell faces;
 */
void computeXiFaceQuantities(
    int nx, int ny, int nhc, 
    double ***&QuL, double ***&QuR, double ***&Qu
)
{
    for (int i=nhc; i<(nx+nhc); i++) {
        for (int j=nhc; j<(ny+nhc-1); j++) {
            for (int l=0; l<4; l++) {
                Qu[i][j][l] = 0.5*(QuL[i][j][l] + QuR[i][j][l]);
            }
        }
    }
}

/* Method to compute xi flux at cell faces.
 *
 * Parameters:
 * -----------
 *  struct Parameters& par  : const reference to struct of parameters;
 *  int                nx   : number of nodes in xi direction;
 *  int                ny   : number of nodes in eta direction;
 *  int                nhc  : number of halo cells;
 *  double***&         Su   : projected cell face areas in the xi direction;
 *  double**&          V    : cell volumes;
 *  double***&         Qu   : conservative variables in the xi cell faces;
 *  double***&         Eh   : xi flux at cell faces;
 */
void computeXiFlux(
    const Parameters& par,
    int nx, int ny, int nhc, 
    double***& Su, double**& V,
    double***& Qu, double***&Eh
)
{
    double u, v, u_xi, S_xi, p, vc;
    for (int i=nhc; i<(nx+nhc); i++) {
        for (int j=nhc; j<(ny+nhc-1); j++) {
            //vc   = 0.5*(V[i-1][j]+V[i][j]);
            u    = Qu[i][j][1]/Qu[i][j][0];
            v    = Qu[i][j][2]/Qu[i][j][0];
            S_xi = sqrt(Su[i][j][0]*Su[i][j][0] + Su[i][j][1]*Su[i][j][1]);
            u_xi = u*Su[i][j][0]/S_xi + v*Su[i][j][1]/S_xi;
            p    = (Qu[i][j][3] \
                   - 0.5*(Qu[i][j][1]*u + Qu[i][j][2]*v))*(par.td.gamma-1.0);
            Eh[i][j][0] = Qu[i][j][0]*u_xi;
            Eh[i][j][1] = Qu[i][j][1]*u_xi + p*Su[i][j][0]/S_xi;
            Eh[i][j][2] = Qu[i][j][2]*u_xi + p*Su[i][j][1]/S_xi;
            Eh[i][j][3] = (Qu[i][j][3]+p)*u_xi;
            //u    = Qhu[i][j][1]/Qhu[i][j][0];
            //v    = Qhu[i][j][2]/Qhu[i][j][0];
            //S_xi = sqrt(Su[i][j][0]*Su[i][j][0] + Su[i][j][1]*Su[i][j][1]);
            //u_xi = u*Su[i][j][0]/S_xi + v*Su[i][j][1]/S_xi;
            //p    = (Qhu[i][j][3]/vc \
            //       - 0.5*(Qhu[i][j][1]*u/vc + Qhu[i][j][2]*v/vc))*(par.td.gamma-1.0);
            //Eh[i][j][0] = Qhu[i][j][0]*u_xi/vc;
            //Eh[i][j][1] = Qhu[i][j][1]*u_xi/vc + p*Su[i][j][0]/S_xi;
            //Eh[i][j][2] = Qhu[i][j][2]*u_xi/vc + p*Su[i][j][1]/S_xi;
            //Eh[i][j][3] = (Qhu[i][j][3]/vc+p)*u_xi;
        }
    }
}

void computeXiDissipation(
    const Parameters& par,
    int nx, int ny, int nhc, 
    double***& Su, double**& V,
    double***& QuL, double***& QuR, double***& Qu,
    double***& AQu
)
{
}

/* Compute quantities at eta faces within and include primary grid faces.
 *
 * Parameters:
 * -----------
 *  int        nx  : number of nodes in xi direction;
 *  int        ny  : number of nodes in eta direction;
 *  int        nhc : number of halo cells;
 *  double***& QvL : conservative variables in the eta cell faces from the left 
 *                   state;
 *  double***& QvR : conservative variables in the eta cell faces from the right
 *                   state;
 *  double***& Qv  : conservative variables in the eta cell faces;
 */
void computeEtaFaceQuantities(
    int nx, int ny, int nhc, 
    double ***&QvL, double ***&QvR, double ***&Qv
)
{
    for (int i=nhc; i<(nx+nhc-1); i++) {
        for (int j=nhc; j<(ny+nhc); j++) {
            for (int l=0; l<4; l++) {
                Qv[i][j][l] = 0.5*(QvL[i][j][l] + QvR[i][j][l]);
            }
        }
    }
}

/* Method to compute eta flux at cell faces.
 *
 * Parameters:
 * -----------
 *  struct Parameters& par  : const reference to struct of parameters;
 *  int                nx   : number of nodes in xi direction;
 *  int                ny   : number of nodes in eta direction;
 *  int                nhc  : number of halo cells;
 *  double***&         Sv   : projected cell face areas in the eta direction;
 *  double***&         Qv   : conservative variables in the eta cell faces;
 *  double***&         Fh   : eta flux at cell faces;
 */
void computeEtaFlux(
    const Parameters& par,
    int nx, int ny, int nhc, 
    double***& Sv, double**& V,
    double***& Qv, double***&Fh
)
{
    double u, v, v_et, S_et, p, vc;
    double b_ep, epsilon = 1.0e1;
    for (int i=nhc; i<(nx+nhc-1); i++) {
        for (int j=nhc; j<(ny+nhc); j++) {
            //vc   = V[i][j];
            //vc   = 0.5*(V[i][j]+V[i][j-1]);
            u    = Qv[i][j][1]/Qv[i][j][0];
            v    = Qv[i][j][2]/Qv[i][j][0];
            S_et = sqrt(Sv[i][j][0]*Sv[i][j][0] + Sv[i][j][1]*Sv[i][j][1]);
            v_et = u*Sv[i][j][0]/S_et + v*Sv[i][j][1]/S_et;
            // to comment below afterwards;
            //b_ep = -((fabs(v_et)-epsilon)-fabs(fabs(v_et)-epsilon))/(2.0*fabs(fabs(v_et)-epsilon));
            //v_et = b_ep*copysign(epsilon, v_et) + (1.0-b_ep)*v_et;
            //
            p    = (Qv[i][j][3] \
                   - 0.5*(Qv[i][j][1]*u + Qv[i][j][2]*v))*(par.td.gamma-1.0);
            Fh[i][j][0] = Qv[i][j][0]*v_et;
            Fh[i][j][1] = Qv[i][j][1]*v_et + p*Sv[i][j][0]/S_et;
            Fh[i][j][2] = Qv[i][j][2]*v_et + p*Sv[i][j][1]/S_et;
            Fh[i][j][3] = (Qv[i][j][3]+p)*v_et;
            //u    = Qhv[i][j][1]/Qhv[i][j][0];
            //v    = Qhv[i][j][2]/Qhv[i][j][0];
            //S_et = sqrt(Sv[i][j][0]*Sv[i][j][0] + Sv[i][j][1]*Sv[i][j][1]);
            //v_et = u*Sv[i][j][0]/S_et + v*Sv[i][j][1]/S_et;
            //p    = (Qhv[i][j][3]/vc \
            //       - 0.5*(Qhv[i][j][1]*u/vc + Qhv[i][j][2]*v/vc))*(par.td.gamma-1.0);
            //Fh[i][j][0] = Qhv[i][j][0]*v_et/vc;
            //Fh[i][j][1] = Qhv[i][j][1]*v_et/vc + p*Sv[i][j][0]/S_et;
            //Fh[i][j][2] = Qhv[i][j][2]*v_et/vc + p*Sv[i][j][1]/S_et;
            //Fh[i][j][3] = (Qhv[i][j][3]/vc+p)*v_et;
        }
    }
}

void computeEtaDissipation(
    const Parameters& par,
    int nx, int ny, int nhc, 
    double***& Sv, double**& V,
    double***& QvL, double***& QvR, double***& Qv,
    double***& BQv
)
{
}
