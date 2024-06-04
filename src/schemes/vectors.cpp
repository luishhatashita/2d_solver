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

/* Method to compute xi flux at cell face.
 *
 * Parameters:
 * -----------
 *  struct Parameters& par  : const reference to struct of parameters;
 *  double*&           Su   : projected cell face areas in the xi direction;
 *  double*&           Qu   : conservative variables in the xi cell face;
 *  double*&           Eh   : xi flux at cell face;
 */
void computeXiFlux(
    const Parameters& par,
    double*& Su,
    double* Qu, double*Eh
)
{
    double u, v, u_xi, S_xi, p;
    u    = Qu[1]/Qu[0];
    v    = Qu[2]/Qu[0];
    S_xi = sqrt(Su[0]*Su[0] + Su[1]*Su[1]);
    u_xi = u*Su[0]/S_xi + v*Su[1]/S_xi;
    p    = (Qu[3] - 0.5*(Qu[1]*u + Qu[2]*v))*(par.td.gamma-1.0);
    Eh[0] = Qu[0]*u_xi;
    Eh[1] = Qu[1]*u_xi + p*Su[0]/S_xi;
    Eh[2] = Qu[2]*u_xi + p*Su[1]/S_xi;
    Eh[3] = (Qu[3]+p)*u_xi;
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
void computeXiFluxes(
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
    double*& Su,
    double* QuL, double* QuR,
    double* AdQu
)
{
    // Left state quantities:
    double uL, vL, rhoL, srrhoL, pL, htL, TL;
    uL     =  QuL[1]/QuL[0];
    vL     =  QuL[2]/QuL[0];
    rhoL   =  QuL[0];
    srrhoL = sqrt(rhoL);
    pL     = (QuL[3] - 0.5*(QuL[1]*uL + QuL[2]*vL))*(par.td.gamma-1.0);
    TL     = pL/(rhoL*par.td.R);
    htL    = (QuL[3] + pL)/rhoL;
    
    // Right state quantities:
    double uR, vR, rhoR, srrhoR, pR, htR, TR;
    uR     =  QuR[1]/QuR[0];
    vR     =  QuR[2]/QuR[0];
    rhoR   =  QuR[0];
    srrhoR = sqrt(rhoR);
    pR     = (QuR[3] - 0.5*(QuR[1]*uR + QuR[2]*vR))*(par.td.gamma-1.0);
    TR     = pR/(rhoR*par.td.R);
    htR    = (QuR[3] + pR)/rhoR;

    // Averaged quanties:
    double uh, vh, ch, rhoh, hth, Hth, Th;
    rhoh = sqrt(rhoL*rhoR); 
    uh   = (srrhoL*uL + srrhoR*uR)/(srrhoL + srrhoR);
    vh   = (srrhoL*vL + srrhoR*vR)/(srrhoL + srrhoR);
    hth  = (srrhoL*htL + srrhoR*htR)/(srrhoL + srrhoR);
    Hth  = hth;
    //Hth  = rhoh*hth;
    Th   = (srrhoL*TL + srrhoR*TR)/(srrhoL + srrhoR);
    //ch   = sqrt(par.td.gamma*par.td.R*Th);
    ch   = sqrt((par.td.gamma - 1.0)*(hth - 0.5*(uh*uh + vh*vh)));

    // Dissipation matrix as in Swanson and Turkel, JCP (1992):
    double q, a1, a2, phi;
    double lambda1, lambda2, lambda3, rhoA, Vn, Vl;
    a1  = Su[0]; // S_xi_x = J^-1 \xi_x
    a2  = Su[1]; // S_xi_y = J^-1 \xi_y
    q   = a1*uh + a2*vh;
    phi = 0.5*(uh*uh + vh*vh);
    Vn  = 0.2;
    Vl  = 0.2;
    // Can rearrange for \tilde{lambda_i} later
    lambda1 = q + sqrt(a1*a1 + a2*a2)*ch;
    lambda2 = q - sqrt(a1*a1 + a2*a2)*ch;
    lambda3 = q;
    lambda1 = fabs(lambda1);
    lambda2 = fabs(lambda2);
    lambda3 = fabs(lambda3);
    rhoA    = fabs(q) + ch*sqrt(a1*a1 + a2*a2);
    lambda1 = fmax(lambda1, Vn*rhoA);
    lambda2 = fmax(lambda2, Vn*rhoA);
    lambda3 = fmax(lambda3, Vl*rhoA);

    // Compute each term of the matrix |A|;
    double A11, A12, A13, A14,
           A21, A22, A23, A24,
           A31, A32, A33, A34,
           A41, A42, A43, A44;
    // premultipiers:
    double CE1, CE2, CE12, CE4, CE34;
    CE12 = 0.5*(lambda1 + lambda2) - lambda3; 
    CE1  = (par.td.gamma - 1.0)/(ch*ch);
    CE2  = 1.0/(a1*a1 + a2*a2);
    CE34 = 0.5*(lambda1 - lambda2)/(sqrt(a1*a1 + a2*a2)*ch);
    CE4  = par.td.gamma - 1.0;
    A11  = lambda3 + CE12*CE1*phi                      - CE34*q;
    A12  =         - CE12*CE1*uh                       + CE34*a1;
    A13  =         - CE12*CE1*vh                       + CE34*a2;
    A14  =         + CE12*CE1;
    A21  =         + CE12*CE1*uh*phi  - CE12*CE2*a1*q  - CE34*uh*q   + CE34*CE4*a1*phi;
    A22  = lambda3 - CE12*CE1*uh*uh   + CE12*CE2*a1*a1 + CE34*uh*a1  - CE34*CE4*a1*uh;
    A23  =         - CE12*CE1*uh*vh   + CE12*CE2*a1*a2 + CE34*uh*a2  - CE34*CE4*a1*vh;
    A24  =         + CE12*CE1*uh                                     + CE34*CE4*a1;
    A31  =         + CE12*CE1*vh*phi  - CE12*CE2*a2*q  - CE34*vh*q   + CE34*CE4*a2*phi;
    A32  =         - CE12*CE1*uh*vh   + CE12*CE2*a2*a1 + CE34*vh*a1  - CE34*CE4*a2*uh;
    A33  = lambda3 - CE12*CE1*vh*vh   + CE12*CE2*a2*a2 + CE34*vh*a2  - CE34*CE4*a2*vh;
    A34  =         + CE12*CE1*vh                                     + CE34*CE4*a2;
    A41  =         + CE12*CE1*Hth*phi - CE12*CE2*q*q   - CE34*Hth*q  + CE34*CE4*q*phi;
    A42  =         - CE12*CE1*uh*Hth  + CE12*CE2*q*a1  + CE34*Hth*a1 - CE34*CE4*q*uh;
    A43  =         - CE12*CE1*vh*Hth  + CE12*CE2*q*a2  + CE34*Hth*a2 - CE34*CE4*q*vh;
    A44  = lambda3 + CE12*CE1*Hth                                    + CE34*CE4*q;
    
    // Manually multiply |A|(QR - QL):
    double dQ[4];
    for (int l=0; l<4; l++){
        dQ[l] = QuR[l] - QuL[l];
    }
    AdQu[0] = A11*dQ[0] + A12*dQ[1] + A13*dQ[2] + A14*dQ[3];
    AdQu[1] = A21*dQ[0] + A22*dQ[1] + A23*dQ[2] + A24*dQ[3];
    AdQu[2] = A31*dQ[0] + A32*dQ[1] + A33*dQ[2] + A34*dQ[3];
    AdQu[3] = A41*dQ[0] + A42*dQ[1] + A43*dQ[2] + A44*dQ[3];
}

void computeXiDissipations(
    const Parameters& par,
    int nx, int ny, int nhc, 
    double***& Su, double**& V,
    double***& QuL, double***& QuR, double***& Qu,
    double***& AQu
)
{
}

double computeM4p(double M, double beta)
{
    double M1p, M2p, M2m, M4p;
    M1p =  0.5*(M + fabs(M)); 
    M2p =  0.25*(M+1.0)*(M+1.0);
    M2m = -0.25*(M-1.0)*(M-1.0);
    M4p =  ((fabs(M)-1.0)+fabs(fabs(M)-1.0))/(2.0*fabs(fabs(M)-1.0))*M1p \
          -((fabs(M)-1.0)-fabs(fabs(M)-1.0))/(2.0*fabs(fabs(M)-1.0)) \
                *(M2p*(1.0-16.0*beta*M2m));
    return M4p;
}

double computeM4m(double M, double beta)
{
    double M1m, M2p, M2m, M4m;
    M1m =  0.5*(M - fabs(M)); 
    M2p =  0.25*(M+1.0)*(M+1.0);
    M2m = -0.25*(M-1.0)*(M-1.0);
    M4m =  ((fabs(M)-1.0)+fabs(fabs(M)-1.0))/(2.0*fabs(fabs(M)-1.0))*M1m \
          -((fabs(M)-1.0)-fabs(fabs(M)-1.0))/(2.0*fabs(fabs(M)-1.0)) \
                *(M2m*(1.0+16.0*beta*M2p));
    return M4m;
}

double computep5p(double M, double alpha)
{
    double M1p, M2p, M2m, p5p;
    double epsilon = 1.0e-5;
    M  += epsilon;

    M1p =  0.5*(M + fabs(M)); 
    M2p =  0.25*(M+1.0)*(M+1.0);
    M2m = -0.25*(M-1.0)*(M-1.0);
    //if (fabs(M) >= 1) {
    //    p5p = M1p/M;
    //} else {
    //    p5p = M2p*((2.0-M)-16.0*alpha*M*M2m);
    //}
    p5p =  ((fabs(M)-1.0)+fabs(fabs(M)-1.0))/(2.0*fabs(fabs(M)-1.0))*M1p/M \
          -((fabs(M)-1.0)-fabs(fabs(M)-1.0))/(2.0*fabs(fabs(M)-1.0)) \
                *(M2p*((2.0-M)-16.0*alpha*M*M2m));
    return p5p;
}

double computep5m(double M, double alpha)
{
    double M1m, M2p, M2m, p5m;
    double epsilon = 1.0e-5;
    M  += epsilon;

    M1m =  0.5*(M - fabs(M)); 
    M2p =  0.25*(M+1.0)*(M+1.0);
    M2m = -0.25*(M-1.0)*(M-1.0);
    //if (fabs(M) >= 1) {
    //    p5m = M1m/M;
    //} else {
    //    p5m = M2m*((-2.0-M)+16.0*alpha*M*M2p);
    //}
    p5m =  ((fabs(M)-1.0)+fabs(fabs(M)-1.0))/(2.0*fabs(fabs(M)-1.0))*M1m/M \
          -((fabs(M)-1.0)-fabs(fabs(M)-1.0))/(2.0*fabs(fabs(M)-1.0)) \
                *(M2m*((-2.0-M)+16.0*alpha*M*M2p));
    return p5m;
}

void computeXiAUSMFlux(
    const Parameters &par, 
    double*& Su, 
    double* QuL, double* QuR,
    double* Eh
)
{
    // AUSM Parameters:
    double beta  = par.ausm.beta; 
    double kp    = par.ausm.kp;
    double ku    = par.ausm.ku;
    double sigma = par.ausm.sigma;
    double epsilon = 1.0e-15;

    // Metrics
    double S_xi_x, S_xi_y, S_xi;
    S_xi_x = Su[0];
    S_xi_y = Su[1];
    S_xi   = sqrt(S_xi_x*S_xi_x + S_xi_y*S_xi_y);
    // Left state quantities:
    double uL, vL, u_xiL, rhoL, pL, htL, cs2L, chL, ML;
    double psiL[4];
    rhoL   =  QuL[0];
    uL     =  QuL[1]/QuL[0];
    vL     =  QuL[2]/QuL[0];
    u_xiL  =  uL*S_xi_x/S_xi + vL*S_xi_y/S_xi;
    pL     = (QuL[3] - 0.5*(QuL[1]*uL + QuL[2]*vL))*(par.td.gamma-1.0);
    htL    = (QuL[3] + pL)/rhoL;
    cs2L   =  2.0*(par.td.gamma - 1.0)*htL/(par.td.gamma + 1.0);
    chL    =  cs2L/fmax(sqrt(cs2L), u_xiL);
    //chL    =  cs2L/fmax(sqrt(cs2L), fabs(u_xiL));
    ML     =  u_xiL/chL;
    ML    += epsilon;
    psiL[0] = 1.0;
    psiL[1] = uL;
    psiL[2] = vL;
    psiL[3] = htL;
    
    // Right state quantities:
    double uR, vR, u_xiR, rhoR, pR, htR, cs2R, chR, MR;
    double psiR[4];
    rhoR   =  QuR[0];
    uR     =  QuR[1]/QuR[0];
    vR     =  QuR[2]/QuR[0];
    u_xiR  =  uR*S_xi_x/S_xi + vR*S_xi_y/S_xi;
    pR     = (QuR[3] - 0.5*(QuR[1]*uR + QuR[2]*vR))*(par.td.gamma-1.0);
    htR    = (QuR[3] + pR)/rhoR;
    cs2R   =  2.0*(par.td.gamma - 1.0)*htR/(par.td.gamma + 1.0);
    chR    =  cs2R/fmax(sqrt(cs2R), -u_xiR);
    //chR    =  cs2R/fmax(sqrt(cs2R), fabs(u_xiR));
    MR     =  u_xiR/chR;
    MR    += epsilon;
    psiR[0] = 1.0;
    psiR[1] = uR;
    psiR[2] = vR;
    psiR[3] = htR;

    // Interface quantities;
    double rho, c, Mb2, M0, M02, fa, alpha, Mp, M, pu, p;
    double pv[4];
    rho   = 0.5*(rhoL + rhoR);
    c     = fmin(chL, chR);
    Mb2   = 0.5*(u_xiL*u_xiL + u_xiR*u_xiR)/(c*c); 
    M02   = fmin(1.0, fmax(Mb2, par.ref.Mref*par.ref.Mref));
    M0    = sqrt(M02);
    fa    = M0*(2-M0);
    Mp    = -kp*fmax(1.0-sigma*Mb2, 0.0)*(pR-pL)/(fa*rho*c*c);
    M     = computeM4p(ML, beta) + computeM4m(MR, beta) \
            + Mp;
    M    += epsilon;
    alpha = 3.0*(-4.0+5.0*fa*fa)/16.0;
    pu    = -ku*computep5p(ML, alpha)*computep5m(MR, alpha) \
             *(rhoL+rhoR)*c*(u_xiR-u_xiL);
    p     = computep5p(ML, alpha)*pL + computep5m(MR, alpha)*pR \
            + pu;
    pv[0] = 0.0;
    pv[1] = S_xi_x*p/S_xi;
    pv[2] = S_xi_y*p/S_xi;
    pv[3] = 0.0;

    // Compute flux terms
    double mdot;//, b_mdotp;
    //if (M > 0.0) {
    //    mdot = c*M*rhoL;
    //} else {
    //    mdot = c*M*rhoR;
    //}
    mdot = (M + fabs(M))/(2*fabs(M))*c*M*rhoL \
          -(M - fabs(M))/(2*fabs(M))*c*M*rhoR;
    //mdot += epsilon;
    //b_mdotp = (mdot + fabs(mdot))/(2.0*fabs(mdot));
    //for (int l=0; l<4; l++) {
    //    Eh[l] = mdot*(b_mdotp*psiL[l] + (1.0-b_mdotp)*psiR[l]) \
    //          + pv[l];
    //}
    for (int l=0; l<4; l++) {
        Eh[l] = 0.5*mdot*(psiR[l]+psiL[l]) \
              - 0.5*fabs(mdot)*(psiR[l]-psiL[l]) \
              + pv[l];
    }
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

/* Method to compute eta flux at cell face.
 *
 * Parameters:
 * -----------
 *  struct Parameters& par  : const reference to struct of parameters;
 *  double***&         Sv   : projected cell face areas in the eta direction;
 *  double***&         Qv   : conservative variables in the eta cell faces;
 *  double***&         Fh   : eta flux at cell faces;
 */
void computeEtaFlux(
    const Parameters& par,
    double*& Sv,
    double* Qv, double*Fh
)
{
    double u, v, v_et, S_et, p;
    u    = Qv[1]/Qv[0];
    v    = Qv[2]/Qv[0];
    S_et = sqrt(Sv[0]*Sv[0] + Sv[1]*Sv[1]);
    v_et = u*Sv[0]/S_et + v*Sv[1]/S_et;
    p    = (Qv[3] - 0.5*(Qv[1]*u + Qv[2]*v))*(par.td.gamma-1.0);
    Fh[0] = Qv[0]*v_et;
    Fh[1] = Qv[1]*v_et + p*Sv[0]/S_et;
    Fh[2] = Qv[2]*v_et + p*Sv[1]/S_et;
    Fh[3] = (Qv[3]+p)*v_et;
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
void computeEtaFluxes(
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
    double*& Sv,
    double* QvL, double* QvR,
    double* BdQv
)
{
    // Left state quantities:
    double uL, vL, rhoL, srrhoL, pL, htL, TL;
    uL     =  QvL[1]/QvL[0];
    vL     =  QvL[2]/QvL[0];
    rhoL   =  QvL[0];
    srrhoL = sqrt(rhoL);
    pL     = (QvL[3] - 0.5*(QvL[1]*uL + QvL[2]*vL))*(par.td.gamma-1.0);
    TL     = pL/(rhoL*par.td.R);
    htL    = (QvL[3] + pL)/rhoL;
    
    // Right state quantities:
    double uR, vR, rhoR, srrhoR, pR, htR, TR;
    uR     =  QvR[1]/QvR[0];
    vR     =  QvR[2]/QvR[0];
    rhoR   =  QvR[0];
    srrhoR = sqrt(rhoR);
    pR     = (QvR[3] - 0.5*(QvR[1]*uR + QvR[2]*vR))*(par.td.gamma-1.0);
    TR     = pR/(rhoR*par.td.R);
    htR    = (QvR[3] + pR)/rhoR;

    // Averaged quanties:
    double uh, vh, ch, rhoh, hth, Hth, Th;
    rhoh = sqrt(rhoL*rhoR); 
    uh   = (srrhoL*uL + srrhoR*uR)/(srrhoL + srrhoR);
    vh   = (srrhoL*vL + srrhoR*vR)/(srrhoL + srrhoR);
    hth  = (srrhoL*htL + srrhoR*htR)/(srrhoL + srrhoR);
    Hth  = hth;
    //Hth  = rhoh*hth;
    Th   = (srrhoL*TL + srrhoR*TR)/(srrhoL + srrhoR);
    //ch   = sqrt(par.td.gamma*par.td.R*Th);
    ch   = sqrt((par.td.gamma - 1.0)*(hth - 0.5*(uh*uh + vh*vh)));

    // Dissipation matrix as in Swanson and Turkel, JCP (1992):
    double q, a1, a2, phi;
    double lambda1, lambda2, lambda3, rhoB, Vn, Vl;
    a1  = Sv[0]; // S_eta_x = J^-1 \eta_x
    a2  = Sv[1]; // S_eta_y = J^-1 \eta_y
    q   = a1*uh + a2*vh;
    phi = 0.5*(uh*uh + vh*vh);
    Vn  = 0.2;
    Vl  = 0.2;
    // Can rearrange for \tilde{lambda_i} later
    lambda1 = q + sqrt(a1*a1 + a2*a2)*ch;
    lambda2 = q - sqrt(a1*a1 + a2*a2)*ch;
    lambda3 = q;
    lambda1 = fabs(lambda1);
    lambda2 = fabs(lambda2);
    lambda3 = fabs(lambda3);
    rhoB    = fabs(q) + ch*sqrt(a1*a1 + a2*a2);
    lambda1 = fmax(lambda1, Vn*rhoB);
    lambda2 = fmax(lambda2, Vn*rhoB);
    lambda3 = fmax(lambda3, Vl*rhoB);

    // Compute each term of the matrix |A|;
    double B11, B12, B13, B14,
           B21, B22, B23, B24,
           B31, B32, B33, B34,
           B41, B42, B43, B44;
    // premultipiers:
    double CE1, CE2, CE12, CE4, CE34;
    CE12 = 0.5*(lambda1 + lambda2) - lambda3; 
    CE1  = (par.td.gamma - 1.0)/(ch*ch);
    CE2  = 1.0/(a1*a1 + a2*a2);
    CE34 = 0.5*(lambda1 - lambda2)/(sqrt(a1*a1 + a2*a2)*ch);
    CE4  = par.td.gamma - 1.0;
    B11  = lambda3 + CE12*CE1*phi                      - CE34*q;
    B12  =         - CE12*CE1*uh                       + CE34*a1;
    B13  =         - CE12*CE1*vh                       + CE34*a2;
    B14  =         + CE12*CE1;
    B21  =         + CE12*CE1*uh*phi  - CE12*CE2*a1*q  - CE34*uh*q   + CE34*CE4*a1*phi;
    B22  = lambda3 - CE12*CE1*uh*uh   + CE12*CE2*a1*a1 + CE34*uh*a1  - CE34*CE4*a1*uh;
    B23  =         - CE12*CE1*uh*vh   + CE12*CE2*a1*a2 + CE34*uh*a2  - CE34*CE4*a1*vh;
    B24  =         + CE12*CE1*uh                                     + CE34*CE4*a1;
    B31  =         + CE12*CE1*vh*phi  - CE12*CE2*a2*q  - CE34*vh*q   + CE34*CE4*a2*phi;
    B32  =         - CE12*CE1*uh*vh   + CE12*CE2*a2*a1 + CE34*vh*a1  - CE34*CE4*a2*uh;
    B33  = lambda3 - CE12*CE1*vh*vh   + CE12*CE2*a2*a2 + CE34*vh*a2  - CE34*CE4*a2*vh;
    B34  =         + CE12*CE1*vh                                     + CE34*CE4*a2;
    B41  =         + CE12*CE1*Hth*phi - CE12*CE2*q*q   - CE34*Hth*q  + CE34*CE4*q*phi;
    B42  =         - CE12*CE1*uh*Hth  + CE12*CE2*q*a1  + CE34*Hth*a1 - CE34*CE4*q*uh;
    B43  =         - CE12*CE1*vh*Hth  + CE12*CE2*q*a2  + CE34*Hth*a2 - CE34*CE4*q*vh;
    B44  = lambda3 + CE12*CE1*Hth                                    + CE34*CE4*q;
    
    // Manually multiply |A|(QR - QL):
    double dQ[4];
    for (int l=0; l<4; l++){
        dQ[l] = QvR[l] - QvL[l];
    }
    BdQv[0] = B11*dQ[0] + B12*dQ[1] + B13*dQ[2] + B14*dQ[3];
    BdQv[1] = B21*dQ[0] + B22*dQ[1] + B23*dQ[2] + B24*dQ[3];
    BdQv[2] = B31*dQ[0] + B32*dQ[1] + B33*dQ[2] + B34*dQ[3];
    BdQv[3] = B41*dQ[0] + B42*dQ[1] + B43*dQ[2] + B44*dQ[3];
}

void computeEtaDissipations(
    const Parameters& par,
    int nx, int ny, int nhc, 
    double***& Sv, double**& V,
    double***& QvL, double***& QvR, double***& Qv,
    double***& BQv
)
{
}

void computeEtaAUSMFlux(
    const Parameters &par, 
    double*& Sv, 
    double* QvL, double* QvR,
    double* Fh
)
{
    // AUSM Parameters:
    double beta  = par.ausm.beta; 
    double kp    = par.ausm.kp;
    double ku    = par.ausm.ku;
    double sigma = par.ausm.sigma;
    double epsilon = 1.0e-15;

    // Metrics
    double S_eta_x, S_eta_y, S_eta;
    S_eta_x = Sv[0];
    S_eta_y = Sv[1];
    S_eta   = sqrt(S_eta_x*S_eta_x + S_eta_y*S_eta_y);
    // Left state quantities:
    double uL, vL, u_etaL, rhoL, pL, htL, cs2L, chL, ML;
    double psiL[4];
    rhoL   =  QvL[0];
    uL     =  QvL[1]/QvL[0];
    vL     =  QvL[2]/QvL[0];
    u_etaL =  uL*S_eta_x/S_eta + vL*S_eta_y/S_eta;
    pL     = (QvL[3] - 0.5*(QvL[1]*uL + QvL[2]*vL))*(par.td.gamma-1.0);
    htL    = (QvL[3] + pL)/rhoL;
    cs2L   =  2.0*(par.td.gamma - 1.0)*htL/(par.td.gamma + 1.0);
    chL    =  cs2L/fmax(sqrt(cs2L), u_etaL);
    //chL    =  cs2L/fmax(sqrt(cs2L), fabs(u_etaL));
    ML     =  u_etaL/chL;
    ML    += epsilon;
    psiL[0] = 1.0;
    psiL[1] = uL;
    psiL[2] = vL;
    psiL[3] = htL;
    
    // Right state quantities:
    double uR, vR, u_etaR, rhoR, pR, htR, cs2R, chR, MR;
    double psiR[4];
    rhoR   =  QvR[0];
    uR     =  QvR[1]/QvR[0];
    vR     =  QvR[2]/QvR[0];
    u_etaR =  uR*S_eta_x/S_eta + vR*S_eta_y/S_eta;
    pR     = (QvR[3] - 0.5*(QvR[1]*uR + QvR[2]*vR))*(par.td.gamma-1.0);
    htR    = (QvR[3] + pR)/rhoR;
    cs2R   =  2.0*(par.td.gamma - 1.0)*htR/(par.td.gamma + 1.0);
    chR    =  cs2R/fmax(sqrt(cs2R), -u_etaR);
    //chR    =  cs2R/fmax(sqrt(cs2R), fabs(u_etaR));
    MR     =  u_etaR/chR;
    MR    += epsilon;
    psiR[0] = 1.0;
    psiR[1] = uR;
    psiR[2] = vR;
    psiR[3] = htR;

    // Interface quantities;
    double rho, c, Mb2, M0, M02, fa, alpha, Mp, M, pu, p;
    double pv[4];
    rho   = 0.5*(rhoL + rhoR);
    c     = fmin(chL, chR);
    Mb2   = 0.5*(u_etaL*u_etaL + u_etaR*u_etaR)/(c*c); 
    M02   = fmin(1.0, fmax(Mb2, par.ref.Mref*par.ref.Mref));
    M0    = sqrt(M02);
    fa    = M0*(2.0 - M0);
    Mp    = -kp*fmax(1.0 - sigma*Mb2, 0.0)*(pR - pL)/(fa*rho*c*c);
    M     = computeM4p(ML, beta) + computeM4m(MR, beta) \
            + Mp;
    M    += epsilon;
    alpha = 3.0*(-4.0 + 5.0*fa*fa)/16.0;
    pu    = -ku*computep5p(ML, alpha)*computep5m(MR, alpha) \
             *(rhoL+rhoR)*c*(u_etaR-u_etaL);
    p     = computep5p(ML, alpha)*pL + computep5m(MR, alpha)*pR \
            + pu;
    pv[0] = 0.0;
    pv[1] = S_eta_x*p/S_eta;
    pv[2] = S_eta_y*p/S_eta;
    pv[3] = 0.0;

    // Compute flux terms
    double mdot;//, b_mdotp;

    mdot = (M + fabs(M))/(2*fabs(M))*c*M*rhoL \
          -(M - fabs(M))/(2*fabs(M))*c*M*rhoR;
    //mdot += epsilon;
    //b_mdotp = (mdot + fabs(mdot))/(2.0*fabs(mdot));
    //for (int l=0; l<4; l++) {
    //    Fh[l] = mdot*(b_mdotp*psiL[l] + (1.0-b_mdotp)*psiR[l]) \
    //          + pv[l];
    //}
    for (int l=0; l<4; l++) {
        Fh[l] = 0.5*mdot*(psiR[l]+psiL[l]) \
              - 0.5*fabs(mdot)*(psiR[l]-psiL[l]) \
              + pv[l];
    }
}

