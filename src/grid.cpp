#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <cmath>
#include "readers.h"
#include "allocate.h"
#include "grid.h"

using namespace std;

/* Constructor method that reads the Oefelein formatted grid file based on spe-
 * cified path.                                                                
 *                                                                             
 * Parameters:                                                                 
 * -----------                                                                 
 *  string cfpath : relative path + filename of grid file.                     
 */
Grid::Grid(string cfpath) 
{
    // Implementation using the pointers properly;
    gfpath  = cfpath;
    // Reader function for Oefelein specific grid file;
    readGrid(cfpath, &gcoords, &gnxi, &gneta);
    
    // Allocate memory;
    allocate3D(gnxi, gneta, 2, gnodes);

    // Populate 2D array from list of coordinates;
    for (int j=0; j<gneta; j++) {
        for (int i=0; i<gnxi; i++) {
            gnodes[i][j][0] = gcoords[i+gnxi*j][0];
            gnodes[i][j][1] = gcoords[i+gnxi*j][1];
        }
    }
    cout << "--2D array of grid nodes arranged." << endl;
}

// Destructor of Grid
Grid::~Grid() 
{
    deallocate3D(gnxi        , gneta        , gnodes    );
    deallocate3D(gnxi+2*nhc  , gneta+2*nhc  , gnodes_whc);
    deallocate3D(gnxi+2*nhc-1, gneta+2*nhc-1, gxc_whc   );
    deallocate2D(gnxi+2*nhc-1,                ga_whc    );
    deallocate3D(gnxi+2*nhc-1, gneta+2*nhc-1, ggcm_whc  );
    deallocate2D(gnxi+2*nhc-1,                ginvj_whc );
}

// Method to get computational dimension and size.
array<int, 2> Grid::getCompDim()
{
    // Return array from private object variables;
    array<int, 2> acompdim;
    acompdim[0] = gnxi;
    acompdim[1] = gneta;
    return acompdim;
}

// Method to get nodes as 3d vector.
double*** Grid::getNodes()
{
    // Return array of nodes;
    return gnodes;
}

// Method to get generalized metrics matrix as 3d vector.
double*** Grid::getGeneralizedCoordinateMetrics()
{
    // Return array of nodes;
    return ggcm_whc;
}

// Method to get inverse jacobian as 2d vector.
double** Grid::getInverseJacobian()
{
    // Return array of nodes;
    return ginvj_whc;
}

/* Method to add n halo cells around the grid.
 *
 * Parameters:
 * -----------
 *  int    inhc     : number of halos cells as input;
 *  bool   writecsv : flag to write grid nodes into list of points in csv file;
 *  string wfpath   : file name for grid filw with hallo cells.
 */
void Grid::addHaloCells(int inhc, bool writecsv)
{
    // Set grid halo cells whilst initializing the padded domain;
    nhc = inhc;
    
    // Allocate memory;
    allocate3D(gnxi+2*nhc, gneta+2*nhc, 2, gnodes_whc);

    // Copy previous node values first;
    for (int j=0; j<gneta; j++) {
        for (int i=0; i<gnxi; i++) {
            gnodes_whc[nhc+i][nhc+j][0] = gnodes[i][j][0];
            gnodes_whc[nhc+i][nhc+j][1] = gnodes[i][j][1];
        }
    }

    // Update padded cells with information from boundary nodes;
    // 1. west and east;
    double dxw, dxe;
    for (int i=0; i<nhc; i++) {
        for (int j=0; j<gneta; j++) {
            dxw                              = gnodes[1][j][0] - gnodes[0][j][0];
            gnodes_whc[nhc-(i+1)][nhc+j][0]  = gnodes[0][j][0] - (i+1)*dxw;
            gnodes_whc[nhc-(i+1)][nhc+j][1]  = gnodes[0][j][1];
            dxe                              = gnodes[gnxi-1][j][0] - gnodes[gnxi-2][j][0];
            gnodes_whc[nhc+gnxi+i][nhc+j][0] = gnodes[gnxi-1][j][0] + (i+1)*dxe;
            gnodes_whc[nhc+gnxi+i][nhc+j][1] = gnodes[gnxi-1][j][1];
        }
    }

    // 2. south and north;
    double dys, dyn;
    for (int j=0; j<nhc; j++) {
        for (int i=0; i<gnxi; i++) {
            dys                               = gnodes[i][1][1] - gnodes[i][0][1];
            gnodes_whc[nhc+i][nhc-(j+1)][0]   = gnodes[i][0][0];
            gnodes_whc[nhc+i][nhc-(j+1)][1]   = gnodes[i][0][1] - (j+1)*dys;
            dyn                               = gnodes[i][gneta-1][1] - gnodes[i][gneta-2][1];
            gnodes_whc[nhc+i][nhc+gneta+j][0] = gnodes[i][gneta-1][0];
            gnodes_whc[nhc+i][nhc+gneta+j][1] = gnodes[i][gneta-1][1] + (j+1)*dyn;
        }
    }

    // 3. corners - SW, SE, NW, NE;
    for (int i=0; i<nhc; i++) {
        for (int j=0; j<nhc; j++) {
            // SW - march left and down, copying x from above and y from the right
            gnodes_whc[nhc-(i+1)][nhc-(j+1)][0]    = gnodes_whc[nhc-(i+1)][nhc][0];
            gnodes_whc[nhc-(i+1)][nhc-(j+1)][1]    = gnodes_whc[nhc][nhc-(j+1)][1];
            // SE - march right and down, copying x from above and y from the left
            gnodes_whc[nhc+gnxi+i][nhc-(j+1)][0]   = gnodes_whc[nhc+gnxi+i][nhc][0];
            gnodes_whc[nhc+gnxi+i][nhc-(j+1)][1]   = gnodes_whc[nhc+gnxi-1][nhc-(j+1)][1];
            // NW - march left and up, copying x from below and y from the right
            gnodes_whc[nhc-(i+1)][nhc+gneta+j][0]  = gnodes_whc[nhc-(i+1)][nhc+gneta-1][0];
            gnodes_whc[nhc-(i+1)][nhc+gneta+j][1]  = gnodes_whc[nhc][nhc+gneta+j][1];
            // NE - march right and up, copying x from below and y from the left
            gnodes_whc[nhc+gnxi+i][nhc+gneta+j][0] = gnodes_whc[nhc+gnxi+i][nhc+gneta-1][0];
            gnodes_whc[nhc+gnxi+i][nhc+gneta+j][1] = gnodes_whc[nhc+gnxi-1][nhc+gneta+j][1];
        }
    }

    // Write file of grid nodes with halo cells;
    if (writecsv) {
        string wfpath = gfpath + "_whc.dat";
        writeGridWithHalos(&wfpath, gnodes_whc, &nhc, &gnxi, &gneta);
    }
    cout << "--2D array of grid nodes with halo cells computed." << endl;
}

/* Method to compute cell center based on padded nodes.
 *
 * Parameters:
 * -----------
 *  int        gnxi       : size of xi dimension;
 *  int        gneta      : size of eta dimension;
 *  int        nhc        : number of halo cells;
 *  double***  gnodes_whc : 3d array of nodes with halo cells;
 *  double***& gxc_whc    : 3d array of cells centers with halo cells;
 */
void computeCellCenters(
    int gnxi, int gneta, int nhc, 
    double*** gnodes_whc, double***& gxc_whc
)
{
    // Allocate memory;
    allocate3D(gnxi+2*nhc-1, gneta+2*nhc-1, 2, gxc_whc);
    for (int j=0; j<(gneta+2*nhc-1); j++) {
        for (int i=0; i<(gnxi+2*nhc-1); i++) {
            gxc_whc[i][j][0] = 0.5*(gnodes_whc[i][j][0] + gnodes_whc[i+1][j][0]);
            gxc_whc[i][j][1] = 0.5*(gnodes_whc[i][j][1] + gnodes_whc[i][j+1][1]);
        }
    }
    cout << "--2D array of grid cell centers computed." << endl;
}

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
 *  int       gnxi       : size of xi dimension;
 *  int       gneta      : size of eta dimension;
 *  int       nhc        : number of halo cells;
 *  double*** gnodes_whc : 3d array of nodes with halo cells;
 *  double**& ga_whc     : 2d array of cells areas with halo cells;
 */
void computeCellAreas(
    int gnxi, int gneta, int nhc,
    double*** gnodes_whc, double**& ga_whc
) 
{
    // Allocate memory;
    allocate2D(gnxi+2*nhc-1, gneta+2*nhc-1, ga_whc);
    for (int j=0; j<(gneta+2*nhc-1); j++) {
        for (int i=0; i<(gnxi+2*nhc-1); i++) {
            ga_whc[i][j] = 0.5*fabs(
                ((gnodes_whc[i][j+1][0] - gnodes_whc[i+1][j  ][0])  \
                *(gnodes_whc[i][j  ][1] - gnodes_whc[i+1][j+1][1])) \
               -((gnodes_whc[i][j  ][0] - gnodes_whc[i+1][j+1][0])  \
                *(gnodes_whc[i][j+1][1] - gnodes_whc[i+1][j  ][1]))
            );
        }
    }
    cout << "--2D array of grid cell areas computed." << endl;
}

/* Method to compute physical coordinate metrics, ie.
 * A = [[x_xi, x_eta], [y_xi, y_eta]];
 *
 * Parameters:
 * -----------
 *  int        gnxi     : size of xi dimension;
 *  int        gneta    : size of eta dimension;
 *  int        nhc      : number of halo cells;
 *  double***  gxc_whc  : 3d array of cell centers with halo cells;
 *  double***& gpcm_whc : 3d array of physical coordinate metrics considering 
 *                        halo cells;
 */
void computePhysicalCoordinateMetrics(
    int gnxi, int gneta, int nhc,
    double*** gxc_whc, double***& gpcm_whc
) 
{
    cout << "--Computing physical coordinate metrics:" << endl;
    // Compute physical coordinate metrics by second order central differences 
    // inside the domain and second order forward/backward at the boundaries.
    // 1. Inside - central in xi and central in eta:
    for (int i=1; i<(gnxi+2*nhc-2); i++) {
        for (int j=1; j<(gneta+2*nhc-2); j++) {
            gpcm_whc[i][j][0] = 0.5*(gxc_whc[i+1][j][0] - gxc_whc[i-1][j][0]);
            gpcm_whc[i][j][1] = 0.5*(gxc_whc[i][j+1][0] - gxc_whc[i][j-1][0]);
            gpcm_whc[i][j][2] = 0.5*(gxc_whc[i+1][j][1] - gxc_whc[i-1][j][1]);
            gpcm_whc[i][j][3] = 0.5*(gxc_whc[i][j+1][1] - gxc_whc[i][j-1][1]);
        }
    }
    cout << "----INSIDE ok." << endl;
    // 2. Faces: 
    // the indexed gnxi+2*nhc-2, gneta+2*nhc-2 exist.
    //cout << gxc_whc[gnxi+2*nhc-2][gneta+2*nhc-2][0] << ","
    //     << gxc_whc[gnxi+2*nhc-2][gneta+2*nhc-2][1] << endl;
    for (int j=1; j<(gneta+2*nhc-2); j++) {
        // a) West - forward in xi and central in eta:
        gpcm_whc[0][j][0] = 0.5*(
            -3*gxc_whc[0][j][0] + 4*gxc_whc[1][j][0] - gxc_whc[2][j][0]
        );
        gpcm_whc[0][j][1] = 0.5*(gxc_whc[0][j+1][0] - gxc_whc[0][j-1][0]);
        gpcm_whc[0][j][2] = 0.5*(
            -3*gxc_whc[0][j][1] + 4*gxc_whc[1][j][1] - gxc_whc[2][j][1]
        );
        gpcm_whc[0][j][3] = 0.5*(gxc_whc[0][j+1][1] - gxc_whc[0][j-1][1]);
        // b) East - backward in xi and central in eta:
        gpcm_whc[gnxi+2*nhc-2][j][0] = 0.5*(
             3*gxc_whc[gnxi+2*nhc-2][j][0] \
           - 4*gxc_whc[gnxi+2*nhc-3][j][0] \
           +   gxc_whc[gnxi+2*nhc-4][j][0]
        );
        gpcm_whc[gnxi+2*nhc-2][j][1] = 0.5*(
            gxc_whc[gnxi+2*nhc-2][j+1][0] - gxc_whc[gnxi+2*nhc-2][j-1][0]
        );
        gpcm_whc[gnxi+2*nhc-2][j][2] = 0.5*(
             3*gxc_whc[gnxi+2*nhc-2][j][1] \
           - 4*gxc_whc[gnxi+2*nhc-3][j][1] \
           +   gxc_whc[gnxi+2*nhc-4][j][1]
        );
        gpcm_whc[gnxi+2*nhc-2][j][3] = 0.5*(
            gxc_whc[gnxi+2*nhc-2][j+1][1] - gxc_whc[gnxi+2*nhc-2][j-1][1]
        );
    }
    cout << "----WE FACES ok." << endl;
    for (int i=1; i<(gnxi+2*nhc-2); i++) {
        // c. South - central in xi and forward in eta:
        gpcm_whc[i][0][0] = 0.5*(gxc_whc[i+1][0][0] - gxc_whc[i-1][0][0]);
        gpcm_whc[i][0][1] = 0.5*(
            - 3*gxc_whc[i][0][0] + 4*gxc_whc[i][1][0] - gxc_whc[i][2][0]
        ); 
        gpcm_whc[i][0][2] = 0.5*(gxc_whc[i+1][0][1] - gxc_whc[i-1][0][1]); 
        gpcm_whc[i][0][3] = 0.5*(
            - 3*gxc_whc[i][0][1] + 4*gxc_whc[i][1][1] - gxc_whc[i][2][1]
        ); 
        // d. North - central in xi and backward in eta:
        gpcm_whc[i][gneta+2*nhc-2][0] = 0.5*(
            gxc_whc[i+1][gneta+2*nhc-2][0] - gxc_whc[i-1][gneta+2*nhc-2][0]
        );
        gpcm_whc[i][gneta+2*nhc-2][1] = 0.5*(
              3*gxc_whc[i][gneta+2*nhc-2][0] \
            - 4*gxc_whc[i][gneta+2*nhc-3][0] \
            +   gxc_whc[i][gneta+2*nhc-4][0]
        ); 
        gpcm_whc[i][gneta+2*nhc-2][2] = 0.5*(
            gxc_whc[i+1][gneta+2*nhc-2][1] - gxc_whc[i-1][gneta+2*nhc-2][1]
        ); 
        gpcm_whc[i][gneta+2*nhc-2][3] = 0.5*(
              3*gxc_whc[i][gneta+2*nhc-2][1] \
            - 4*gxc_whc[i][gneta+2*nhc-3][1] \
            +   gxc_whc[i][gneta+2*nhc-4][1]
        ); 
    }
    cout << "----SN FACES ok." << endl;
    // 3. Corners:
    // a) SW - forward in xi and forward in eta:
    gpcm_whc[0][0][0] = 0.5*(
        - 3*gxc_whc[0][0][0] \
        + 4*gxc_whc[1][0][0] \
        -   gxc_whc[2][0][0]
    );
    gpcm_whc[0][0][1] = 0.5*(
        - 3*gxc_whc[0][0][0] \
        + 4*gxc_whc[0][1][0] \
        -   gxc_whc[0][2][0]
    );
    gpcm_whc[0][0][2] = 0.5*(
        - 3*gxc_whc[0][0][1] \
        + 4*gxc_whc[1][0][1] \
        -   gxc_whc[2][0][1]
    );
    gpcm_whc[0][0][3] = 0.5*(
        - 3*gxc_whc[0][0][1] \
        + 4*gxc_whc[0][1][1] \
        -   gxc_whc[0][2][1]
    );
    // b) SE - backward in xi and forward in eta:
    gpcm_whc[gnxi+2*nhc-2][0][0] = 0.5*(
        + 3*gxc_whc[gnxi+2*nhc-2][0][0] \
        - 4*gxc_whc[gnxi+2*nhc-3][0][0] \
        +   gxc_whc[gnxi+2*nhc-4][0][0]
    );
    gpcm_whc[gnxi+2*nhc-2][0][1] = 0.5*(
        - 3*gxc_whc[gnxi+2*nhc-2][0][0] \
        + 4*gxc_whc[gnxi+2*nhc-2][1][0] \
        -   gxc_whc[gnxi+2*nhc-2][2][0]
    );
    gpcm_whc[gnxi+2*nhc-2][0][2] = 0.5*(
        + 3*gxc_whc[gnxi+2*nhc-2][0][1] \
        - 4*gxc_whc[gnxi+2*nhc-3][0][1] \
        +   gxc_whc[gnxi+2*nhc-4][0][1]
    );
    gpcm_whc[gnxi+2*nhc-2][0][3] = 0.5*(
        - 3*gxc_whc[gnxi+2*nhc-2][0][1] \
        + 4*gxc_whc[gnxi+2*nhc-2][1][1] \
        -   gxc_whc[gnxi+2*nhc-2][2][1]
    );
    // c) NW - forward in xi and backward in eta:
    gpcm_whc[0][gneta+2*nhc-2][0] = 0.5*(
        - 3*gxc_whc[0][gneta+2*nhc-2][0] \
        + 4*gxc_whc[1][gneta+2*nhc-2][0] \
        -   gxc_whc[2][gneta+2*nhc-2][0]
    );
    gpcm_whc[0][gneta+2*nhc-2][1] = 0.5*(
        + 3*gxc_whc[0][gneta+2*nhc-2][0] \
        - 4*gxc_whc[0][gneta+2*nhc-3][0] \
        +   gxc_whc[0][gneta+2*nhc-4][0]
    );
    gpcm_whc[0][gneta+2*nhc-2][2] = 0.5*(
        - 3*gxc_whc[0][gneta+2*nhc-2][1] \
        + 4*gxc_whc[1][gneta+2*nhc-2][1] \
        -   gxc_whc[2][gneta+2*nhc-2][1]
    );
    gpcm_whc[0][gneta+2*nhc-2][3] = 0.5*(
        + 3*gxc_whc[0][gneta+2*nhc-2][1] \
        - 4*gxc_whc[0][gneta+2*nhc-3][1] \
        +   gxc_whc[0][gneta+2*nhc-4][1]
    );
    // d) NE - backward in xi and backward in eta:
    gpcm_whc[gnxi+2*nhc-2][gneta+2*nhc-2][0] = 0.5*(
        + 3*gxc_whc[gnxi+2*nhc-2][gneta+2*nhc-2][0] \
        - 4*gxc_whc[gnxi+2*nhc-3][gneta+2*nhc-2][0] \
        +   gxc_whc[gnxi+2*nhc-4][gneta+2*nhc-2][0]
    );
    gpcm_whc[gnxi+2*nhc-2][gneta+2*nhc-2][1] = 0.5*(
        + 3*gxc_whc[gnxi+2*nhc-2][gneta+2*nhc-2][0] \
        - 4*gxc_whc[gnxi+2*nhc-2][gneta+2*nhc-3][0] \
        +   gxc_whc[gnxi+2*nhc-2][gneta+2*nhc-4][0]
    );
    gpcm_whc[gnxi+2*nhc-2][gneta+2*nhc-2][2] = 0.5*(
        + 3*gxc_whc[gnxi+2*nhc-2][gneta+2*nhc-2][1] \
        - 4*gxc_whc[gnxi+2*nhc-3][gneta+2*nhc-2][1] \
        +   gxc_whc[gnxi+2*nhc-4][gneta+2*nhc-2][1]
    );
    gpcm_whc[gnxi+2*nhc-2][gneta+2*nhc-2][3] = 0.5*(
        + 3*gxc_whc[gnxi+2*nhc-2][gneta+2*nhc-2][1] \
        - 4*gxc_whc[gnxi+2*nhc-2][gneta+2*nhc-3][1] \
        +   gxc_whc[gnxi+2*nhc-2][gneta+2*nhc-4][1]
    );
    cout << "----CORNERS ok." << endl;
}

/* Method to compute generalized coordinate metrics.
 * 1. calls method to compute A = [[x_xi, x_eta], [y_xi, y_eta]];
 * 2. compute 1/J  = det([[x_xi, x_eta],[y_xi, y_eta]])
 *                 = x_xi*y_eta - x_eta*y_xi;
 * 3. compute A^-1 = 1/J [[y_eta, -x_eta],[-y_xi, x_xi]]
 *
 * Parameters:
 * -----------
 *  int        gnxi     : size of xi dimension;
 *  int        gneta    : size of eta dimension;
 *  int        nhc      : number of halo cells;
 *  double***  gxc_whc  : 3d array of cell centers with halo cells;
 *  double***& ggcm_whc : 3d array of generalized coordinate metrics considering 
 *                        halo cells;
 *  double**&  ginvj_whc: 2d array of the inverse of the jacobian. 
 */
void computeGeneralCoordinateMetrics(
    int gnxi, int gneta, int nhc,
    double*** gxc_whc, double***& ggcm_whc, double**& ginvj_whc
)
{
    double*** gpcm_whc;
    // Allocate memory for gpcm_whc (physical coordinate metrics), ggcm_whc (ge-
    // neralized coordinate metrics), ginvj_whc (jacobian inverse).
    allocate3D(gnxi+2*nhc-1, gneta+2*nhc-1, 4, gpcm_whc );
    allocate3D(gnxi+2*nhc-1, gneta+2*nhc-1, 4, ggcm_whc );
    allocate2D(gnxi+2*nhc-1, gneta+2*nhc-1   , ginvj_whc);

    // Calls method to compute physical coordinate metrics, separate because 
    // needs proper treatments near the boundaries.
    computePhysicalCoordinateMetrics(gnxi, gneta, nhc, gxc_whc, gpcm_whc);
    //print3Darray(gpcm_whc, gnxi+2*nhc-1, gneta+2*nhc-1, 4);
    cout << "--3D array of physical coordinate metrics computed." << endl;

    // Compute the inverse of the Jacobian:
    for (int i=0; i<(gnxi+2*nhc-1); i++) {
        for (int j=0; j<(gneta+2*nhc-1); j++) {
            ginvj_whc[i][j] = gpcm_whc[i][j][0]*gpcm_whc[i][j][3] \
                            - gpcm_whc[i][j][1]*gpcm_whc[i][j][2];
        }
    }

    // Compute the generalized coordinate metrics:
    for (int i=0; i<(gnxi+2*nhc-1); i++) {
        for (int j=0; j<(gneta+2*nhc-1); j++) {
            ggcm_whc[i][j][0] =   ginvj_whc[i][j]*gpcm_whc[i][j][3];
            ggcm_whc[i][j][1] = - ginvj_whc[i][j]*gpcm_whc[i][j][1];
            ggcm_whc[i][j][2] = - ginvj_whc[i][j]*gpcm_whc[i][j][2];
            ggcm_whc[i][j][3] =   ginvj_whc[i][j]*gpcm_whc[i][j][0];
        }
    }
    
    cout << "--2D array of generalized coordinate metrics computed." << endl;
    deallocate3D(gnxi+2*nhc-1, gneta+2*nhc-1, gpcm_whc);
}


/* Driver method to compute all grid metrics:
 * 1. cell centers;
 * 2. cell areas;
 * 3. projected cell face areas; - TODO
 * 4. grid metrics for generalized coordinates;
 */
void Grid::computeMetrics() 
{
    computeCellCenters             (gnxi, gneta, nhc, gnodes_whc, gxc_whc  );
    //print3Darray(gxc_whc, gnxi+2*nhc-1, gneta+2*nhc-1, 2);
    computeCellAreas               (gnxi, gneta, nhc, gnodes_whc, ga_whc   );
    //print2Darray(ga_whc, gnxi+2*nhc-1, gneta+2*nhc-1);
    computeGeneralCoordinateMetrics(
        gnxi, gneta, nhc, 
        gxc_whc, ggcm_whc, ginvj_whc
    );
    //print3Darray(gpcm_whc, gnxi+2*nhc-1, gneta+2*nhc-1, 4);
    //print3Darray(gpcm_whc, gnxi+2*nhc-1, gneta+2*nhc-1, 4);
}

