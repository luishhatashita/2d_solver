#include "solution.h"

#include <iostream>
#include <cmath>

#include "allocate.h"
#include "grid.h"
#include "parameters.h"
#include "vectors.h"
#include "schemes.h"
#include "error.h"
#include "writers.h"

/* Constructor method to initialize and allocate memory for all arrays.                                                                
 *                                                                             
 * Parameters:                                                                 
 * -----------                                                                 
 *  class  Grid&     grid  : const reference to the grid class with all array 
 *                           size information and metrics;
 *  struct Parameters& par : const reference to struct will all case related para-
 *                           meters;
 *  int              nrest : number of restart file to write or read;
 *  bool             write : boolean flag to write restart files;
 */
Solution::Solution(
    const class Grid& grid, const Parameters& par, 
    int nrest, 
    bool write
)
{
    // Get dimensions and number of halos
    int *dims = new int[3];
    grid.getArraySizes(dims);
    m_nx  = dims[0];
    m_ny  = dims[1];
    m_nhc = dims[2];

    s_nit = par.rt.nit;
    //s_grid = grid;
    //s_par  = par;

    // Allocate memory:
    allocate3D(m_nx+2*m_nhc-1, m_ny+2*m_nhc-1, 4, m_Qn);
    allocate3D(m_nx+2*m_nhc-1, m_ny+2*m_nhc-1, 4, m_Qvn);
    allocate3D(m_nx+2*m_nhc-1, m_ny+2*m_nhc-1, 4, m_Qvnp1);
    
    // Solution start either from file or scratch:
    m_write = write;
    if (nrest >= 0) {
        // read rest file;
    } else {
        // create initial restart;
        initializeField(par, grid, m_nx, m_ny, m_nhc, m_Qn, m_Qvn, m_Qvnp1);
        //for (int i=0; i<(m_nx+2*m_nhc-1); i++) {
        //    for (int j=0; j<(m_ny+2*m_nhc-1); j++) {
        //        for (int l=0; l<4; l++) {
        //            m_Qvn[i][j][l] = m_Qvnp1[i][j][l];
        //        }
        //    }
        //}
        if (m_write) {
            writePrimitives(par, grid, 0);
            writeConservatives(grid, 0);
        }
    }

    // Additional memory allocation needed:
    allocate3D(m_nx+2*m_nhc,   m_ny+2*m_nhc-1, 4, m_QnuL );
    allocate3D(m_nx+2*m_nhc,   m_ny+2*m_nhc-1, 4, m_QnuR );
    allocate3D(m_nx+2*m_nhc,   m_ny+2*m_nhc-1, 4, m_Qnu  );
    allocate3D(m_nx+2*m_nhc-1, m_ny+2*m_nhc,   4, m_QnvL );
    allocate3D(m_nx+2*m_nhc-1, m_ny+2*m_nhc,   4, m_QnvR );
    allocate3D(m_nx+2*m_nhc-1, m_ny+2*m_nhc,   4, m_Qnv  );
    // Fluxes
    allocate3D(m_nx+2*m_nhc,   m_ny+2*m_nhc-1, 4, m_EhnL );
    allocate3D(m_nx+2*m_nhc,   m_ny+2*m_nhc-1, 4, m_EhnR );
    allocate3D(m_nx+2*m_nhc,   m_ny+2*m_nhc-1, 4, m_Ehn  );
    allocate3D(m_nx+2*m_nhc,   m_ny+2*m_nhc-1, 4, m_AQu  );
    allocate3D(m_nx+2*m_nhc-1, m_ny+2*m_nhc,   4, m_FhnL );
    allocate3D(m_nx+2*m_nhc-1, m_ny+2*m_nhc,   4, m_FhnR );
    allocate3D(m_nx+2*m_nhc-1, m_ny+2*m_nhc,   4, m_Fhn  );
    allocate3D(m_nx+2*m_nhc-1, m_ny+2*m_nhc,   4, m_BQv  );
    // Time step and integration
    allocate2D(m_nx+2*m_nhc-1, m_ny+2*m_nhc-1,    m_dt   );
    allocate3D(m_nx+2*m_nhc-1, m_ny+2*m_nhc-1, 4, m_Qnp1 );
    // Error norms;
    allocate2D(s_nit-1, 4,m_L2);
    allocate2D(s_nit-1, 4,m_Linfty);
    
    delete[] dims;
}

// Destructor
Solution::~Solution()
{
    //std::cout << "Solution destructor" << std::endl;
    deallocate3D(m_nx+2*m_nhc-1, m_ny+2*m_nhc-1, m_Qn    );
    deallocate3D(m_nx+2*m_nhc-1, m_ny+2*m_nhc-1, m_Qvn   );
    deallocate3D(m_nx+2*m_nhc,   m_ny+2*m_nhc-1, m_QnuL  );
    deallocate3D(m_nx+2*m_nhc,   m_ny+2*m_nhc-1, m_QnuR  );
    deallocate3D(m_nx+2*m_nhc,   m_ny+2*m_nhc-1, m_Qnu   );
    deallocate3D(m_nx+2*m_nhc-1, m_ny+2*m_nhc,   m_QnvL  );
    deallocate3D(m_nx+2*m_nhc-1, m_ny+2*m_nhc,   m_QnvR  );
    deallocate3D(m_nx+2*m_nhc-1, m_ny+2*m_nhc,   m_Qnv   );
    deallocate3D(m_nx+2*m_nhc,   m_ny+2*m_nhc-1, m_EhnL  );
    deallocate3D(m_nx+2*m_nhc,   m_ny+2*m_nhc-1, m_EhnR  );
    deallocate3D(m_nx+2*m_nhc,   m_ny+2*m_nhc-1, m_Ehn   );
    deallocate3D(m_nx+2*m_nhc,   m_ny+2*m_nhc-1, m_AQu   );
    deallocate3D(m_nx+2*m_nhc-1, m_ny+2*m_nhc,   m_FhnL  );
    deallocate3D(m_nx+2*m_nhc-1, m_ny+2*m_nhc,   m_FhnR  );
    deallocate3D(m_nx+2*m_nhc-1, m_ny+2*m_nhc,   m_Fhn   );
    deallocate3D(m_nx+2*m_nhc-1, m_ny+2*m_nhc,   m_BQv   );
    deallocate2D(m_nx+2*m_nhc-1,                 m_dt    );
    deallocate2D(s_nit-1,                        m_L2    );
    deallocate2D(s_nit-1,                        m_Linfty);
}

/* Method that calls writer method for the primitive variables array for the 
 * specific iteration.
 *
 * Parameters:
 * -----------
 *  int it : iteration;
 */
void Solution::writePrimitives(const Parameters& par, const Grid& grid, int it)
{
    double** V;
    V = grid.getCellVolumes();
    conservativeToPrimitive(par, m_nx, m_ny, m_nhc, V, m_Qn, m_Qvn);
    char buffer[50];
    sprintf(buffer, "./out/rest/Qv%05d.bin", it);
    writeBinary3DArray(buffer, m_nx+2*m_nhc-1, m_ny+2*m_nhc-1, 4, m_Qvn);
}

/* Method that calls writer method for the conservative variables array for the 
 * specific iteration.
 *
 * Parameters:
 * -----------
 *  int it : iteration;
 */
void Solution::writeConservatives(const Grid& grid, int it)
{
    char buffer[50];
    sprintf(buffer, "./out/rest/Q%05d.bin", it);
    //double*** Qn, **V;
    //V = grid.getCellVolumes();
    //allocate3D(m_nx+2*m_nhc-1, m_ny+2*m_nhc-1, 4, Qn);
    //for (int i=0; i<(m_nx+2*m_nhc-1); i++) {
    //    for (int j=0; j<(m_ny+2*m_nhc-1); j++) {
    //        for (int l=0; l<4; l++) {
    //            Qn[i][j][l] = m_Qn[i][j][l]/V[i][j];
    //        }
    //    }
    //}
    writeBinary3DArray(buffer, m_nx+2*m_nhc-1, m_ny+2*m_nhc-1, 4, m_Qn);
    //writeBinary3DArray(buffer, m_nx+2*m_nhc-1, m_ny+2*m_nhc-1, 4, Qn);
    //deallocate3D(m_nx+2*m_nhc-1, m_ny+2*m_nhc-1, Qn);
}

// Method that calls MUSCL interpolation.
void Solution::constructLeftRightStates(double kappa, double epsilon)
{
    interpolateMUSCL(
        kappa, epsilon,
        m_nx, m_ny, m_nhc, 
        m_Qn, 
        m_QnuL, m_QnuR, 
        m_QnvL, m_QnvR
    );
    //char buffer[50];
    //sprintf(buffer, "./out/rest/QuL%05d.bin", 0);
    //writeBinary3DArray(buffer, m_nx+2*m_nhc, m_ny+2*m_nhc-1, 4, m_QnuL);
    //sprintf(buffer, "./out/rest/QuR%05d.bin", 0);
    //writeBinary3DArray(buffer, m_nx+2*m_nhc, m_ny+2*m_nhc-1, 4, m_QnuR);
    //sprintf(buffer, "./out/rest/QvL%05d.bin", 0);
    //writeBinary3DArray(buffer, m_nx+2*m_nhc-1, m_ny+2*m_nhc, 4, m_QnvL);
    //sprintf(buffer, "./out/rest/QvR%05d.bin", 0);
    //writeBinary3DArray(buffer, m_nx+2*m_nhc-1, m_ny+2*m_nhc, 4, m_QnvR);
}

/* Method that call the flux construction methods.
 *
 * Parameters:
 * -----------
 *  struct Parameters& par  : reference to simulation parameters;
 *  class  Grid&       grid : reference to grid with finite volume metrics;
 */
void Solution::computeFluxes(const Parameters& par, const Grid& grid)
{
    // Compute fluxes:
    double ***Su, ***Sv, **V;
    Su = grid.getProjectedFaceAreasXi();
    Sv = grid.getProjectedFaceAreasEta();
    V  = grid.getCellVolumes();
    computeFirstOrderUpwindFluxes(
        par,
        m_nx, m_ny, m_nhc,
        Su, Sv, V,
        m_QnuL, m_QnuR, m_Qnu, m_EhnL, m_EhnR, m_Ehn,
        m_QnvL, m_QnvR, m_Qnv, m_FhnL, m_FhnR, m_Fhn
    );
    //char buffer[50];
    //sprintf(buffer, "./out/rest/Eh%05d.bin", 0);
    //writeBinary3DArray(buffer, m_nx+2*m_nhc, m_ny+2*m_nhc-1, 4, m_Ehn);
    //sprintf(buffer, "./out/rest/Fh%05d.bin", 0);
    //writeBinary3DArray(buffer, m_nx+2*m_nhc-1, m_ny+2*m_nhc, 4, m_Fhn);
    //computeRoeFluxes(
    //    par,
    //    m_nx, m_ny, m_nhc,
    //    Su, Sv, V,
    //    m_QnuL, m_QnuR, m_Qnu, m_EhnL, m_EhnR, m_Ehn,
    //    m_QnvL, m_QnvR, m_Qnv, m_FhnL, m_FhnR, m_Fhn
    //);
}

/* Method to compute local time step based on maximum CFL.
 *
 * Parameters:
 * -----------
 *  struct Parameters& par  : reference to simulation parameters;
 *  class  Grid&       grid : reference to grid with finite volume metrics;
 */
void Solution::computeTimeSteps(
    const Parameters& par, const Grid& grid, 
    double& dt_max
)
{
    // Grid metrics:
    double ***Su, ***Sv, **V;
    Su = grid.getProjectedFaceAreasXi();
    Sv = grid.getProjectedFaceAreasEta();
    V  = grid.getCellVolumes();

    conservativeToPrimitive(par, m_nx, m_ny, m_nhc, V, m_Qn, m_Qvn);

    // Loop through interior cells:
    double u, v, u_xi, v_et,
           S_xi_x, S_xi_y, S_xi, S_et_x, S_et_y, S_et,
           rho_xi, rho_et, c,
           min, max;
    double b_ep, epsilon = 1.0e1;
    for (int i=m_nhc; i<(m_nx+m_nhc-1); i++) {
        for (int j=m_nhc; j<(m_ny+m_nhc-1); j++) {
            // Interpolate metrics at cell centers:
            S_xi_x = 0.5*(Su[i][j][0]+Su[i+1][j][0]);
            S_xi_y = 0.5*(Su[i][j][1]+Su[i+1][j][1]);
            S_et_x = 0.5*(Sv[i][j][0]+Sv[i][j+1][0]);
            S_et_y = 0.5*(Sv[i][j][1]+Sv[i][j+1][1]);
            // Compute contra-variant velocities:
            u      = m_Qn[i][j][1]/m_Qn[i][j][0];
            v      = m_Qn[i][j][2]/m_Qn[i][j][0];
            S_xi   = sqrt(S_xi_x*S_xi_x + S_xi_y*S_xi_y);
            u_xi   = u*S_xi_x/S_xi + v*S_xi_y/S_xi;
            S_et   = sqrt(S_et_x*S_et_x + S_et_y*S_et_y);
            v_et   = u*S_et_x/S_et + v*S_et_y/S_et;
            //b_ep = -((fabs(v_et)-epsilon)-fabs(fabs(v_et)-epsilon))/(2.0*fabs(fabs(v_et)-epsilon));
            //v_et = b_ep*copysign(epsilon, v_et) + (1.0-b_ep)*v_et;
            // Local time step within CFL constraint:
            // Const c or function of T as \sqrt{\gamma R T}?
            c      = sqrt(par.td.gamma*par.td.R*m_Qvn[i][j][3]);
            //c      = par.td.c;
            rho_xi = fabs(u_xi) + c;
            rho_et = fabs(v_et) + c;
            //rho_xi = fabs(u_xi) + par.td.c;
            //rho_et = fabs(v_et) + par.td.c;
            // This first option looks incorrect;
            //m_dt[i][j] = par.rt.CFL*fmin(1.0/rho_xi, 1.0/rho_et);
            // The second option provides convergence for CFL closer to 1.0
            m_dt[i][j] = par.rt.CFL \
                *fmin(V[i][j]/(rho_xi*S_xi), V[i][j]/(rho_et*S_et));
        }
    }
    min = m_dt[0][0];
    max = m_dt[0][0];
    for (int i=m_nhc; i<(m_nx+m_nhc-1); i++) {
        for (int j=m_nhc; j<(m_ny+m_nhc-1); j++) {
            min = fmin(min, m_dt[i][j]);
            max = fmax(max, m_dt[i][j]);
        }
    }
    dt_max = max;
    m_dt_min = min;
    m_dt_max = max;
    //char buffer[50];
    //sprintf(buffer, "./out/rest/dt%05d.bin", 0);
    //writeBinary2DArray(buffer, m_nx+2*m_nhc-1, m_ny+2*m_nhc-1, m_dt);
}

void Solution::checkTimeStep(
    const Parameters& par, const Grid& grid, 
    double dt
)
{
    // Grid metrics:
    double ***Su, ***Sv, **V;
    Su = grid.getProjectedFaceAreasXi();
    Sv = grid.getProjectedFaceAreasEta();
    V  = grid.getCellVolumes();

    conservativeToPrimitive(par, m_nx, m_ny, m_nhc, V, m_Qn, m_Qvn);

    // Loop through interior cells:
    double u, v, u_xi, v_et,
           S_xi_x, S_xi_y, S_xi, S_et_x, S_et_y, S_et,
           rho_xi, rho_et, c,
           min, max_CFL, local_CFL;
    max_CFL = 0.0;
    for (int i=m_nhc; i<(m_nx+m_nhc-1); i++) {
        for (int j=m_nhc; j<(m_ny+m_nhc-1); j++) {
            // Interpolate metrics at cell centers:
            S_xi_x = 0.5*(Su[i][j][0]+Su[i+1][j][0]);
            S_xi_y = 0.5*(Su[i][j][1]+Su[i+1][j][1]);
            S_et_x = 0.5*(Sv[i][j][0]+Sv[i][j+1][0]);
            S_et_y = 0.5*(Sv[i][j][1]+Sv[i][j+1][1]);
            // Compute contra-variant velocities:
            u      = m_Qn[i][j][1]/m_Qn[i][j][0];
            v      = m_Qn[i][j][2]/m_Qn[i][j][0];
            S_xi   = sqrt(S_xi_x*S_xi_x + S_xi_y*S_xi_y);
            u_xi   = u*S_xi_x/S_xi + v*S_xi_y/S_xi;
            S_et   = sqrt(S_et_x*S_et_x + S_et_y*S_et_y);
            v_et   = u*S_et_x/S_et + v*S_et_y/S_et;
            // Local time step within CFL constraint:
            // Const c or function of T as \sqrt{\gamma R T}?
            //c      = sqrt(par.td.gamma*par.td.R*m_Qvn[i][j][3]);
            c      = par.td.c;
            rho_xi = fabs(u_xi) + c;
            rho_et = fabs(v_et) + c;
            //rho_xi = fabs(u_xi) + par.td.c;
            //rho_et = fabs(v_et) + par.td.c;
            //local_CFL = dt*fmax(rho_xi, rho_et);
            local_CFL = dt*fmax(rho_xi*S_xi/V[i][j], rho_et*S_et/V[i][j]);
            max_CFL = fmax(max_CFL, local_CFL);
        }
    }
    m_CFL_max = max_CFL;
}

/* Method to integrate over time based on the time-marching approach, thus non-
 * uniform time steps.
 *
 * Parameters:
 * -----------
 *  struct Parameters& par  : reference to simulation parameters;
 *  class  Grid&       grid : reference to grid with finite volume metrics;
 */
void Solution::integrateLocalTime(const Parameters& par, const Grid& grid)
{
    // Grid metrics:
    double ***Su, ***Sv, **V;
    Su = grid.getProjectedFaceAreasXi();
    Sv = grid.getProjectedFaceAreasEta();
    V  = grid.getCellVolumes();
    
    // Loop through interior cells:
    // left face in xi 
    double S_xi_L, S_xi_R, S_et_L, S_et_R;
    for (int i=m_nhc; i<(m_nx+m_nhc-1); i++) {
        for (int j=m_nhc; j<(m_ny+m_nhc-1); j++) {
            S_xi_L = sqrt(Su[i][j][0]*Su[i][j][0]     + Su[i][j][1]*Su[i][j][1]);
            S_xi_R = sqrt(Su[i+1][j][0]*Su[i+1][j][0] + Su[i+1][j][1]*Su[i+1][j][1]);
            S_et_L = sqrt(Sv[i][j][0]*Sv[i][j][0]     + Sv[i][j][1]*Sv[i][j][1]);
            S_et_R = sqrt(Sv[i][j+1][0]*Sv[i][j+1][0] + Sv[i][j+1][1]*Sv[i][j+1][1]);
            for (int l=0; l<4; l++) {
                m_Qnp1[i][j][l] = m_Qn[i][j][l] \
                    - m_dt[i][j]*((m_Ehn[i+1][j][l]*S_xi_R-m_Ehn[i][j][l]*S_xi_L)
                                 +(m_Fhn[i][j+1][l]*S_et_R-m_Fhn[i][j][l]*S_et_L))/V[i][j];
            }
        }
    }

    //char buffer[50];
    //sprintf(buffer, "./out/rest/Qnp1%05d.bin", 0);
    //writeBinary3DArray(buffer, m_nx+2*m_nhc-1, m_ny+2*m_nhc-1, 4, m_Qnp1);

    // Update boundary conditions:
    conservativeToPrimitive(par, m_nx, m_ny, m_nhc, V, m_Qnp1, m_Qvnp1);
    computeBoundaryConditions(par, grid, m_nx, m_ny, m_nhc, m_Qvn, m_Qvnp1);
    primitiveToConservative(par, m_nx, m_ny, m_nhc, V, m_Qnp1, m_Qvnp1);
    //sprintf(buffer, "./out/rest/Qnp1%05d_bc.bin", 0);
    //writeBinary3DArray(buffer, m_nx+2*m_nhc-1, m_ny+2*m_nhc-1, 4, m_Qnp1);
}

/* Method to integrate over time based on the time-accurate approach, thus uni-
 * form time steps.
 *
 * Parameters:
 * -----------
 *  struct Parameters& par  : reference to simulation parameters;
 *  class  Grid&       grid : reference to grid with finite volume metrics;
 */
void Solution::integrateFixedTime(const Parameters& par, const Grid& grid)
{
    // Grid metrics:
    double ***Su, ***Sv, **V;
    Su = grid.getProjectedFaceAreasXi();
    Sv = grid.getProjectedFaceAreasEta();
    V  = grid.getCellVolumes();
    
    // Loop through interior cells:
    // left face in xi 
    double S_xi_L, S_xi_R, S_et_L, S_et_R;
    for (int i=m_nhc; i<(m_nx+m_nhc-1); i++) {
        for (int j=m_nhc; j<(m_ny+m_nhc-1); j++) {
            S_xi_L = sqrt(Su[i  ][j  ][0]*Su[i  ][j  ][0] + Su[i  ][j  ][1]*Su[i  ][j  ][1]);
            S_xi_R = sqrt(Su[i+1][j  ][0]*Su[i+1][j  ][0] + Su[i+1][j  ][1]*Su[i+1][j  ][1]);
            S_et_L = sqrt(Sv[i  ][j  ][0]*Sv[i  ][j  ][0] + Sv[i  ][j  ][1]*Sv[i  ][j  ][1]);
            S_et_R = sqrt(Sv[i  ][j+1][0]*Sv[i  ][j+1][0] + Sv[i  ][j+1][1]*Sv[i  ][j+1][1]);
            for (int l=0; l<4; l++) {
                m_Qnp1[i][j][l] = m_Qn[i][j][l] \
                    - par.rt.dt*((m_Ehn[i+1][j][l]*S_xi_R-m_Ehn[i][j][l]*S_xi_L)
                                +(m_Fhn[i][j+1][l]*S_et_R-m_Fhn[i][j][l]*S_et_L))/V[i][j];
            }
        }
    }

    //char buffer[50];
    //sprintf(buffer, "./out/rest/Qnp1%05d.bin", 0);
    //writeBinary3DArray(buffer, m_nx+2*m_nhc-1, m_ny+2*m_nhc-1, 4, m_Qnp1);

    // Update boundary conditions:
    conservativeToPrimitive  (par, m_nx, m_ny, m_nhc, V, m_Qnp1, m_Qvnp1);
    computeBoundaryConditions(par, grid, m_nx, m_ny, m_nhc, m_Qvn, m_Qvnp1);
    primitiveToConservative  (par, m_nx, m_ny, m_nhc, V, m_Qnp1, m_Qvnp1);
    //sprintf(buffer, "./out/rest/Qnp1%05d_bc.bin", 0);
    //writeBinary3DArray(buffer, m_nx+2*m_nhc-1, m_ny+2*m_nhc-1, 4, m_Qnp1);
}

// Copies values from Qn+1 to Qn.
void Solution::updateVectors()
{
    for (int i=0; i<(m_nx+2*m_nhc-1); i++) {
        for (int j=0; j<(m_ny+2*m_nhc-1); j++) {
            for (int l=0; l<4; l++) {
                m_Qn[i][j][l] = m_Qnp1[i][j][l];
            }
        }
    }
}

void Solution::computeErrorNorms(const Parameters& par)
{
    computeL2(par, m_nx, m_ny, m_nhc, m_Qn, m_Qnp1, m_L2[m_i-1]);
    computeLinfinity(par, m_nx, m_ny, m_nhc, m_Qn, m_Qnp1, m_Linfty[m_i-1]);

    // Find min/max of L2 and L_\infty
    double minL2, maxL2, minLinfty, maxLinfty;
    minL2     = m_L2[m_i-1][0];
    maxL2     = m_L2[m_i-1][0];
    minLinfty = m_Linfty[m_i-1][0];
    maxLinfty = m_Linfty[m_i-1][0];
    for (int l=1; l<4; l++) {
        minL2     = fmin(minL2, m_L2[m_i-1][l]);
        maxL2     = fmax(maxL2, m_L2[m_i-1][l]);
        minLinfty = fmin(minLinfty, m_Linfty[m_i-1][l]);
        maxLinfty = fmax(maxLinfty, m_Linfty[m_i-1][l]);
    }
    m_L2_min     = minL2;
    m_L2_max     = maxL2;
    m_Linfty_min = minLinfty;
    m_Linfty_max = maxLinfty;
}

void Solution::writeErrorNorms()
{
    char buffer[50];
    sprintf(buffer, "./out/rest/L2.bin");
    writeBinary2DArray(buffer, s_nit-1, 4, m_L2);
    sprintf(buffer, "./out/rest/Linfty.bin");
    writeBinary2DArray(buffer, s_nit-1, 4, m_Linfty);
}

void Solution::setIteration(int i)
{
    m_i = i;
}

void Solution::log()
{
    std::cout << "Iteration " << m_i << " :" << std::endl;
    // Time stepping:
    //std::cout << "CFL:     max = " << m_CFL_max << std::endl;
    std::cout << "dt:      min = " << m_dt_min 
              << "; max = " << m_dt_max << std::endl;
    // Error norms;
    std::cout << "L_2:     min = " << m_L2_min 
              << "; max = " << m_L2_max << std::endl;
    std::cout << "L_infty: min = " << m_Linfty_min 
              << "; max = " << m_Linfty_max << std::endl;
}

