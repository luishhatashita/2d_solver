#include "metrics_gc.h"

#include <iostream>
#include <cmath>

#include "allocate.h"

/* Method to compute cell areas which are analytically equal to cell volumes in
 * 2D for a unit length width.
 *
 * From geometry the area of the quadrilateral A-B-C-D, may be computed by the
 * cross product of its diagonals herein defined as A-C, B-D, following the con-
 * vention of counter-clockwise nodes orientation.
 *
 *     C-----B  The formula for the area is:
 *    /     /       
 *   /     /        Area = 0.5 | \delta x_{CA}\delta y_{db}
 *  D-----A                    - \delta x_{DB}\delta y_{ca} |,
 *
 *  where, \delta_{CA} = x_C - x_A, etc.
 *
 * Parameters:
 * -----------
 *  int       gm_nx       : size of xi dimension;
 *  int       gm_ny      : size of eta dimension;
 *  int       m_nhc        : number of halo cells;
 *  double*** x : 3d array of nodes with halo cells;
 *  double**& ga_whc     : 2d array of cells areas with halo cells;
 */
void computeCellAreas(
    int nx, int ny, int nhc,
    double***& x, double**& a
) 
{
    // Allocate memory;
    allocate2D(nx+2*nhc-1, ny+2*nhc-1, a);
    for (int j=0; j<(ny+2*nhc-1); j++) {
        for (int i=0; i<(nx+2*nhc-1); i++) {
            a[i][j] = 0.5*fabs(
                ((x[i][j+1][0] - x[i+1][j  ][0])  \
                *(x[i][j  ][1] - x[i+1][j+1][1])) \
               -((x[i][j  ][0] - x[i+1][j+1][0])  \
                *(x[i][j+1][1] - x[i+1][j  ][1]))
            );
        }
    }
    std::cout << "--2D array of grid cell areas computed." << std::endl;
}

/* Method to compute physical coordinate metrics, ie.
 * A = [[x_xi, x_eta], [y_xi, y_eta]];
 *
 * Parameters:
 * -----------
 *  int        nx  : size of xi dimension;
 *  int        ny  : size of eta dimension;
 *  int        nhc : number of halo cells;
 *  double***  xc  : 3d array of cell centers with halo cells;
 *  double***& mpc : 3d array of physical coordinate metrics considering 
 *                   halo cells;
 */
void computePhysicalCoordinateMetrics(
    int nx, int ny, int nhc,
    double***& xc, double***& mpc 
) 
{
    std::cout << "--Computing physical coordinate metrics:" << std::endl;
    // Compute physical coordinate metrics by second order central differences 
    // inside the domain and second order forward/backward at the boundaries.
    // 1. Inside - central in xi and central in eta:
    for (int i=1; i<(nx+2*nhc-2); i++) {
        for (int j=1; j<(ny+2*nhc-2); j++) {
            mpc[i][j][0] = 0.5*(xc[i+1][j][0] - xc[i-1][j][0]);
            mpc[i][j][1] = 0.5*(xc[i][j+1][0] - xc[i][j-1][0]);
            mpc[i][j][2] = 0.5*(xc[i+1][j][1] - xc[i-1][j][1]);
            mpc[i][j][3] = 0.5*(xc[i][j+1][1] - xc[i][j-1][1]);
        }
    }
    std::cout << "----INSIDE ok." << std::endl;
    // 2. Faces: 
    // the indexed gm_nx+2*m_nhc-2, gm_ny+2*m_nhc-2 exist.
    //std::cout << gxc_whc[gm_nx+2*m_nhc-2][gm_ny+2*m_nhc-2][0] << ","
    //     << gxc_whc[gm_nx+2*m_nhc-2][gm_ny+2*m_nhc-2][1] << std::endl;
    for (int j=1; j<(ny+2*nhc-2); j++) {
        // a) West - forward in xi and central in eta:
        mpc[0][j][0] = 0.5*(
            -3*xc[0][j][0] + 4*xc[1][j][0] - xc[2][j][0]
        );
        mpc[0][j][1] = 0.5*(xc[0][j+1][0] - xc[0][j-1][0]);
        mpc[0][j][2] = 0.5*(
            -3*xc[0][j][1] + 4*xc[1][j][1] - xc[2][j][1]
        );
        mpc[0][j][3] = 0.5*(xc[0][j+1][1] - xc[0][j-1][1]);
        // b) East - backward in xi and central in eta:
        mpc[nx+2*nhc-2][j][0] = 0.5*(
             3*xc[nx+2*nhc-2][j][0] \
           - 4*xc[nx+2*nhc-3][j][0] \
           +   xc[nx+2*nhc-4][j][0]
        );
        mpc[nx+2*nhc-2][j][1] = 0.5*(
            xc[nx+2*nhc-2][j+1][0] - xc[nx+2*nhc-2][j-1][0]
        );
        mpc[nx+2*nhc-2][j][2] = 0.5*(
             3*xc[nx+2*nhc-2][j][1] \
           - 4*xc[nx+2*nhc-3][j][1] \
           +   xc[nx+2*nhc-4][j][1]
        );
        mpc[nx+2*nhc-2][j][3] = 0.5*(
            xc[nx+2*nhc-2][j+1][1] - xc[nx+2*nhc-2][j-1][1]
        );
    }
    std::cout << "----WE FACES ok." << std::endl;
    for (int i=1; i<(nx+2*nhc-2); i++) {
        // c. South - central in xi and forward in eta:
        mpc[i][0][0] = 0.5*(xc[i+1][0][0] - xc[i-1][0][0]);
        mpc[i][0][1] = 0.5*(
            - 3*xc[i][0][0] + 4*xc[i][1][0] - xc[i][2][0]
        ); 
        mpc[i][0][2] = 0.5*(xc[i+1][0][1] - xc[i-1][0][1]); 
        mpc[i][0][3] = 0.5*(
            - 3*xc[i][0][1] + 4*xc[i][1][1] - xc[i][2][1]
        ); 
        // d. North - central in xi and backward in eta:
        mpc[i][ny+2*nhc-2][0] = 0.5*(
            xc[i+1][ny+2*nhc-2][0] - xc[i-1][ny+2*nhc-2][0]
        );
        mpc[i][ny+2*nhc-2][1] = 0.5*(
              3*xc[i][ny+2*nhc-2][0] \
            - 4*xc[i][ny+2*nhc-3][0] \
            +   xc[i][ny+2*nhc-4][0]
        ); 
        mpc[i][ny+2*nhc-2][2] = 0.5*(
            xc[i+1][ny+2*nhc-2][1] - xc[i-1][ny+2*nhc-2][1]
        ); 
        mpc[i][ny+2*nhc-2][3] = 0.5*(
              3*xc[i][ny+2*nhc-2][1] \
            - 4*xc[i][ny+2*nhc-3][1] \
            +   xc[i][ny+2*nhc-4][1]
        ); 
    }
    std::cout << "----SN FACES ok." << std::endl;
    // 3. Corners:
    // a) SW - forward in xi and forward in eta:
    mpc[0][0][0] = 0.5*(
        - 3*xc[0][0][0] \
        + 4*xc[1][0][0] \
        -   xc[2][0][0]
    );
    mpc[0][0][1] = 0.5*(
        - 3*xc[0][0][0] \
        + 4*xc[0][1][0] \
        -   xc[0][2][0]
    );
    mpc[0][0][2] = 0.5*(
        - 3*xc[0][0][1] \
        + 4*xc[1][0][1] \
        -   xc[2][0][1]
    );
    mpc[0][0][3] = 0.5*(
        - 3*xc[0][0][1] \
        + 4*xc[0][1][1] \
        -   xc[0][2][1]
    );
    // b) SE - backward in xi and forward in eta:
    mpc[nx+2*nhc-2][0][0] = 0.5*(
        + 3*xc[nx+2*nhc-2][0][0] \
        - 4*xc[nx+2*nhc-3][0][0] \
        +   xc[nx+2*nhc-4][0][0]
    );
    mpc[nx+2*nhc-2][0][1] = 0.5*(
        - 3*xc[nx+2*nhc-2][0][0] \
        + 4*xc[nx+2*nhc-2][1][0] \
        -   xc[nx+2*nhc-2][2][0]
    );
    mpc[nx+2*nhc-2][0][2] = 0.5*(
        + 3*xc[nx+2*nhc-2][0][1] \
        - 4*xc[nx+2*nhc-3][0][1] \
        +   xc[nx+2*nhc-4][0][1]
    );
    mpc[nx+2*nhc-2][0][3] = 0.5*(
        - 3*xc[nx+2*nhc-2][0][1] \
        + 4*xc[nx+2*nhc-2][1][1] \
        -   xc[nx+2*nhc-2][2][1]
    );
    // c) NW - forward in xi and backward in eta:
    mpc[0][ny+2*nhc-2][0] = 0.5*(
        - 3*xc[0][ny+2*nhc-2][0] \
        + 4*xc[1][ny+2*nhc-2][0] \
        -   xc[2][ny+2*nhc-2][0]
    );
    mpc[0][ny+2*nhc-2][1] = 0.5*(
        + 3*xc[0][ny+2*nhc-2][0] \
        - 4*xc[0][ny+2*nhc-3][0] \
        +   xc[0][ny+2*nhc-4][0]
    );
    mpc[0][ny+2*nhc-2][2] = 0.5*(
        - 3*xc[0][ny+2*nhc-2][1] \
        + 4*xc[1][ny+2*nhc-2][1] \
        -   xc[2][ny+2*nhc-2][1]
    );
    mpc[0][ny+2*nhc-2][3] = 0.5*(
        + 3*xc[0][ny+2*nhc-2][1] \
        - 4*xc[0][ny+2*nhc-3][1] \
        +   xc[0][ny+2*nhc-4][1]
    );
    // d) NE - backward in xi and backward in eta:
    mpc[nx+2*nhc-2][ny+2*nhc-2][0] = 0.5*(
        + 3*xc[nx+2*nhc-2][ny+2*nhc-2][0] \
        - 4*xc[nx+2*nhc-3][ny+2*nhc-2][0] \
        +   xc[nx+2*nhc-4][ny+2*nhc-2][0]
    );
    mpc[nx+2*nhc-2][ny+2*nhc-2][1] = 0.5*(
        + 3*xc[nx+2*nhc-2][ny+2*nhc-2][0] \
        - 4*xc[nx+2*nhc-2][ny+2*nhc-3][0] \
        +   xc[nx+2*nhc-2][ny+2*nhc-4][0]
    );
    mpc[nx+2*nhc-2][ny+2*nhc-2][2] = 0.5*(
        + 3*xc[nx+2*nhc-2][ny+2*nhc-2][1] \
        - 4*xc[nx+2*nhc-3][ny+2*nhc-2][1] \
        +   xc[nx+2*nhc-4][ny+2*nhc-2][1]
    );
    mpc[nx+2*nhc-2][ny+2*nhc-2][3] = 0.5*(
        + 3*xc[nx+2*nhc-2][ny+2*nhc-2][1] \
        - 4*xc[nx+2*nhc-2][ny+2*nhc-3][1] \
        +   xc[nx+2*nhc-2][ny+2*nhc-4][1]
    );
    std::cout << "----CORNERS ok." << std::endl;
}

/* Method to compute generalized coordinate metrics.
 * 1. calls method to compute A = [[x_xi, x_eta], [y_xi, y_eta]];
 * 2. compute 1/J  = det([[x_xi, x_eta],[y_xi, y_eta]])
 *                 = x_xi*y_eta - x_eta*y_xi;
 * 3. compute A^-1 = 1/J [[y_eta, -x_eta],[-y_xi, x_xi]]
 *
 * Parameters:
 * -----------
 *  int        nx   : size of xi dimension;
 *  int        ny   : size of eta dimension;
 *  int        nhc  : number of halo cells;
 *  double***  xc   : 3d array of cell centers with halo cells;
 *  double***& mgc  : 3d array of generalized coordinate metrics considering 
 *                    halo cells;
 *  double**&  invj : 2d array of the inverse of the jacobian. 
 */
void computeGeneralCoordinateMetrics(
    int nx, int ny, int nhc,
    double***& xc, double***& mgc, double**& invj 
)
{
    double*** mpc;
    // Allocate memory for gpcm_whc (physical coordinate metrics), ggcm_whc (ge-
    // neralized coordinate metrics), ginvj_whc (jacobian inverse).
    allocate3D(nx+2*nhc-1, ny+2*nhc-1, 4, mpc );
    allocate3D(nx+2*nhc-1, ny+2*nhc-1, 4, mgc );
    allocate2D(nx+2*nhc-1, ny+2*nhc-1   , invj);

    // Calls method to compute physical coordinate metrics, separate because 
    // needs proper treatments near the boundaries.
    computePhysicalCoordinateMetrics(nx, ny, nhc, xc, mpc);
    //print3Darray(gpcm_whc, gm_nx+2*m_nhc-1, gm_ny+2*m_nhc-1, 4);
    std::cout << "--3D array of physical coordinate metrics computed." << std::endl;

    // Compute the inverse of the Jacobian:
    for (int i=0; i<(nx+2*nhc-1); i++) {
        for (int j=0; j<(ny+2*nhc-1); j++) {
            invj[i][j] = mpc[i][j][0]*mpc[i][j][3] \
                            - mpc[i][j][1]*mpc[i][j][2];
        }
    }

    // Compute the generalized coordinate metrics:
    for (int i=0; i<(nx+2*nhc-1); i++) {
        for (int j=0; j<(ny+2*nhc-1); j++) {
            mgc[i][j][0] =   invj[i][j]*mpc[i][j][3];
            mgc[i][j][1] = - invj[i][j]*mpc[i][j][1];
            mgc[i][j][2] = - invj[i][j]*mpc[i][j][2];
            mgc[i][j][3] =   invj[i][j]*mpc[i][j][0];
        }
    }
    
    std::cout << "--2D array of generalized coordinate metrics computed." << std::endl;
    deallocate3D(nx+2*nhc-1, ny+2*nhc-1, mpc);
}
