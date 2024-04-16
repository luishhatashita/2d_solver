#include "metrics_fv.h"

#include <iostream>
#include <cmath>

#include "allocate.h"
#include "writers.h"

/* Method to compute cell center based on padded nodes.
 *
 * Parameters:
 * -----------
 *  int        nx    : size of xi dimension;
 *  int        ny    : size of eta dimension;
 *  int        nhc   : number of halo cells;
 *  double***& x     : 3d array of nodes with halo cells;
 *  double***& xc    : 3d array of cells centers with halo cells;
 *  string     fpath : path and file name;
 *  bool       write : boolean flag to write binary file;
 */
void computeCellCenters(
    int nx, int ny, int nhc, 
    double***& x, double***& xc,
    std::string fpath, bool write
)
{
    // Allocate memory;
    allocate3D(nx+2*nhc-1, ny+2*nhc-1, 2, xc);

    // Cell centroids based on average of 4 enclosing nodes:
    for (int i=0; i<(nx+2*nhc-1); i++) {
        for (int j=0; j<(ny+2*nhc-1); j++) {
            xc[i][j][0] = 0.25*(x[i  ][j  ][0] + x[i+1][j  ][0]
                               +x[i+1][j+1][0] + x[i  ][j+1][0]);
            xc[i][j][1] = 0.25*(x[i  ][j  ][1] + x[i+1][j  ][1]
                               +x[i+1][j+1][1] + x[i  ][j+1][1]);
        }
    }

    std::cout << "--2D array of grid cell centers computed." << std::endl;

    if (write) {
        std::string wfpath = "./grid/" + fpath + "_xc.bin";
        writeBinary3DArray(wfpath, nx+2*nhc-1, ny+2*nhc-1, 2, xc);
    }
}

/* Method to compute cell face centroids based on padded nodes.
 *
 * Parameters:
 * -----------
 *  int        nx  : size of xi dimension;
 *  int        ny  : size of eta dimension;
 *  int        nhc : number of halo cells;
 *  double***& x   : 3d array of nodes with halo cells;
 *  double***& xu  : 3d array of nodes with halo cells of face centroids, for 
 *                   which the normal is the xi direction;
 *  double***& xv  : 3d array of nodes with halo cells of face centroids, for
 *                   which the normal is the eta direction;
 *  string     fpath : path and file name;
 *  bool       write : boolean flag to write binary file;
 */
void computeFaceCenters(
    int nx, int ny, int nhc, 
    double***& x, double***& xu, double***& xv,
    std::string fpath, bool write
)
{
    // Allocate memory;
    allocate3D(nx+2*nhc  , ny+2*nhc-1, 2, xu);
    allocate3D(nx+2*nhc-1, ny+2*nhc  , 2, xv);

    // Face normal to xi direction:
    for (int i=0; i<(nx+2*nhc); i++) {
        for (int j=0; j<(ny+2*nhc-1); j++) {
            xu[i][j][0] = 0.5*(x[i][j][0] + x[i][j+1][0]);
            xu[i][j][1] = 0.5*(x[i][j][1] + x[i][j+1][1]);
        }
    }

    // Face normal to eta direction:
    for (int i=0; i<(nx+2*nhc-1); i++) {
        for (int j=0; j<(ny+2*nhc); j++) {
            xv[i][j][0] = 0.5*(x[i][j][0] + x[i+1][j][0]);
            xv[i][j][1] = 0.5*(x[i][j][1] + x[i+1][j][1]);
        }
    }

    std::cout << "--2D arrays of grid cell face centroids computed." << std::endl;

    if (write) {
        std::string wfpath = "./grid/" + fpath + "_xu.bin";
        writeBinary3DArray(wfpath, nx+2*nhc, ny+2*nhc-1, 2, xu);
        wfpath = "./grid/" + fpath + "_xv.bin";
        writeBinary3DArray(wfpath, nx+2*nhc-1, ny+2*nhc, 2, xv);
    }
}

/* Method to compute projected cell face areas based on padded nodes.
 * 
 *      S_{\xi, x} =  y,\eta
 *      S_{\xi, y} = -x,\eta
 *      S_{\eta,x} = -y,\xi
 *      S_{\eta,y} =  x,\xi
 *
 * Parameters:
 * -----------
 *  int        nx  : size of xi dimension;
 *  int        ny  : size of eta dimension;
 *  int        nhc : number of halo cells;
 *  double***& x   : 3d array of nodes with halo cells;
 *  double***& su  : 3d array of nodes with halo projected cell face areas, for 
 *                   which the normal is the xi direction;
 *  double***& sv  : 3d array of nodes with halo projected cell face areas, for
 *                   which the normal is the eta direction;
 *  string     fpath : path and file name;
 *  bool       write : boolean flag to write binary file;
 */
void computeProjectedFaceAreas(
    int nx, int ny, int nhc, 
    double***& x, double***& su, double***& sv,
    std::string fpath, bool write
)
{
    // Allocate memory;
    allocate3D(nx+2*nhc  , ny+2*nhc-1, 2, su);
    allocate3D(nx+2*nhc-1, ny+2*nhc  , 2, sv);

    // Face normal to xi direction:
    for (int i=0; i<(nx+2*nhc); i++) {
        for (int j=0; j<(ny+2*nhc-1); j++) {
            su[i][j][0] =   x[i][j+1][1] - x[i][j][1];
            su[i][j][1] = -(x[i][j+1][0] - x[i][j][0]);
        }
    }

    // Face normal to eta direction:
    for (int i=0; i<(nx+2*nhc-1); i++) {
        for (int j=0; j<(ny+2*nhc); j++) {
            sv[i][j][0] = -(x[i+1][j][1] - x[i][j][1]);
            sv[i][j][1] =   x[i+1][j][0] - x[i][j][0];
        }
    }

    std::cout << "--2D arrays of grid cell face centroids computed." << std::endl;

    if (write) {
        std::string wfpath = "./grid/" + fpath + "_su.bin";
        writeBinary3DArray(wfpath, nx+2*nhc, ny+2*nhc-1, 2, su);
        wfpath = "./grid/" + fpath + "_sv.bin";
        writeBinary3DArray(wfpath, nx+2*nhc-1, ny+2*nhc, 2, sv);
    }
}

/* Method to compute cell volumes (areas in 2d) based on padded nodes.
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
 *  int        nx  : size of xi dimension;
 *  int        ny  : size of eta dimension;
 *  int        nhc : number of halo cells;
 *  double***& x   : 3d array of nodes with halo cells;
 *  double**&  v   : 2d array of cells volumes with halo cells;
 *  string     fpath : path and file name;
 *  bool       write : boolean flag to write binary file;
 */
void computeCellVolumes(
    int nx, int ny, int nhc, 
    double***& x, double**& v,
    std::string fpath, bool write
)
{
    // Allocate memory;
    allocate2D(nx+2*nhc-1, ny+2*nhc-1, v);
    for (int i=0; i<(nx+2*nhc-1); i++) {
        for (int j=0; j<(ny+2*nhc-1); j++) {
            v[i][j] = 0.5*fabs(
                ((x[i][j+1][0] - x[i+1][j  ][0])  \
                *(x[i][j  ][1] - x[i+1][j+1][1])) \
               -((x[i][j  ][0] - x[i+1][j+1][0])  \
                *(x[i][j+1][1] - x[i+1][j  ][1]))
            );
        }
    }

    std::cout << "--2D array of grid cell volumes computed." << std::endl;

    if (write) {
        std::string wfpath = "./grid/" + fpath + "_v.bin";
        writeBinary2DArray(wfpath, nx+2*nhc-1, ny+2*nhc-1, v);
    }
}
