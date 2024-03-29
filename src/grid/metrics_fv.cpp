#include "metrics_fv.h"

#include <iostream>

#include "allocate.h"

/* Method to compute cell center based on padded nodes.
 *
 * Parameters:
 * -----------
 *  int        nx  : size of xi dimension;
 *  int        ny  : size of eta dimension;
 *  int        nhc : number of halo cells;
 *  double***  x   : 3d array of nodes with halo cells;
 *  double***& xc  : 3d array of cells centers with halo cells;
 */
void computeCellCenters(
    int nx, int ny, int nhc, 
    double***& x, double***& xc
)
{
    // Allocate memory;
    allocate3D(nx+2*nhc-1, ny+2*nhc-1, 2, xc);
    for (int j=0; j<(ny+2*nhc-1); j++) {
        for (int i=0; i<(nx+2*nhc-1); i++) {
            xc[i][j][0] = 0.5*(x[i][j][0] + x[i+1][j][0]);
            xc[i][j][1] = 0.5*(x[i][j][1] + x[i][j+1][1]);
        }
    }
    std::cout << "--2D array of grid cell centers computed." << std::endl;
}
