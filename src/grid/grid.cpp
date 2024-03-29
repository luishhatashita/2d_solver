#include "grid.h"

#include <iostream>
#include <vector>
#include <array>
#include <string>

// Ideally, I would like to include as io/readers.h and malloc/allocate.h, 
// without it looking for the directories io and malloc in the current folder.
#include "readers.h"
#include "writers.h"
#include "allocate.h"
#include "metrics_gc.h"
#include "metrics_fv.h"

//using namespace std;

/* Constructor method that reads the Oefelein formatted grid file based on spe-
 * cified path.                                                                
 *                                                                             
 * Parameters:                                                                 
 * -----------                                                                 
 *  string cfpath : relative path + filename of grid file.                     
 */
Grid::Grid(std::string file, int nhc, bool write) 
{
    // Implementation using the pointers properly;
    m_file  = file;
    // Reader function for Oefelein specific grid file;
    readGrid(m_file, m_coords, m_nx, m_ny);
    
    // Allocate memory;
    allocate3D(m_nx, m_ny, 2, m_x_woh);

    // Populate 2D array from list of coordinates;
    for (int j=0; j<m_ny; j++) {
        for (int i=0; i<m_nx; i++) {
            m_x_woh[i][j][0] = m_coords[i+m_nx*j][0];
            m_x_woh[i][j][1] = m_coords[i+m_nx*j][1];
        }
    }

    std::cout << "--2D array of grid nodes arranged." << std::endl;
    
    std::string wfpath = m_file + ".bin";
    writeBinGridWithHalos(wfpath, m_nx, m_ny, 0, m_x_woh);

    Grid::addHaloCells(nhc, write);
    Grid::computeMetrics();
}

// Destructor of Grid
Grid::~Grid() 
{
    deallocate3D(m_nx          , m_ny          , m_x_woh);
    deallocate3D(m_nx+2*m_nhc  , m_ny+2*m_nhc  , m_x    );
    deallocate3D(m_nx+2*m_nhc-1, m_ny+2*m_nhc-1, m_xc   );
    deallocate2D(m_nx+2*m_nhc-1,                 m_a    );
    deallocate3D(m_nx+2*m_nhc-1, m_ny+2*m_nhc-1, m_mgc  );
    deallocate2D(m_nx+2*m_nhc-1,                 m_invj );
}

// Method to get computational dimension and size.
//int* Grid::getCompDim()
//{
//    // Return array from private object variables;
//    int* acompdim = new int[2];
//    acompdim[0] = m_nx;
//    acompdim[1] = m_ny;
//    return acompdim;
//}

// Method to get nodes as 3d vector.
double*** Grid::getNodes()
{
    // Return array of nodes;
    return m_x_woh;
}

// Method to get generalized metrics matrix as 3d vector.
double*** Grid::getGeneralizedCoordinateMetrics()
{
    // Return array of nodes;
    return m_mgc;
}

// Method to get inverse jacobian as 2d vector.
double** Grid::getInverseJacobian()
{
    // Return array of nodes;
    return m_invj;
}

/* Method to add n halo cells around the grid.
 *
 * Parameters:
 * -----------
 *  int  nhc   : number of halos cells as input;
 *  bool write : flag to write grid nodes into list of points in csv and binary 
 *               file;
 * TODO:
 * -----
 *  - improve reflection and implement extrapolation of halo nodes;
 */
void Grid::addHaloCells(int nhc, bool write)
{
    // Set grid halo cells whilst initializing the padded domain;
    m_nhc = nhc;
    
    // Allocate memory;
    allocate3D(m_nx+2*m_nhc, m_ny+2*m_nhc, 2, m_x);

    // Copy previous node values first;
    for (int j=0; j<m_ny; j++) {
        for (int i=0; i<m_nx; i++) {
            m_x[m_nhc+i][m_nhc+j][0] = m_x_woh[i][j][0];
            m_x[m_nhc+i][m_nhc+j][1] = m_x_woh[i][j][1];
        }
    }

    // Update padded cells with information from boundary nodes;
    // 1. west and east;
    double dxw, dxe;
    for (int i=0; i<m_nhc; i++) {
        for (int j=0; j<m_ny; j++) {
            dxw                           = m_x_woh[1][j][0] - m_x_woh[0][j][0];
            m_x[m_nhc-(i+1)][m_nhc+j][0]  = m_x_woh[0][j][0] - (i+1)*dxw;
            m_x[m_nhc-(i+1)][m_nhc+j][1]  = m_x_woh[0][j][1];
            dxe                           = m_x_woh[m_nx-1][j][0] - m_x_woh[m_nx-2][j][0];
            m_x[m_nhc+m_nx+i][m_nhc+j][0] = m_x_woh[m_nx-1][j][0] + (i+1)*dxe;
            m_x[m_nhc+m_nx+i][m_nhc+j][1] = m_x_woh[m_nx-1][j][1];
        }
    }

    // 2. south and north;
    double dys, dyn;
    for (int j=0; j<m_nhc; j++) {
        for (int i=0; i<m_nx; i++) {
            dys                           = m_x_woh[i][1][1] - m_x_woh[i][0][1];
            m_x[m_nhc+i][m_nhc-(j+1)][0]  = m_x_woh[i][0][0];
            m_x[m_nhc+i][m_nhc-(j+1)][1]  = m_x_woh[i][0][1] - (j+1)*dys;
            dyn                           = m_x_woh[i][m_ny-1][1] - m_x_woh[i][m_ny-2][1];
            m_x[m_nhc+i][m_nhc+m_ny+j][0] = m_x_woh[i][m_ny-1][0];
            m_x[m_nhc+i][m_nhc+m_ny+j][1] = m_x_woh[i][m_ny-1][1] + (j+1)*dyn;
        }
    }

    // 3. corners - SW, SE, NW, NE;
    for (int i=0; i<m_nhc; i++) {
        for (int j=0; j<m_nhc; j++) {
            // SW - march left and down, copying x from above and y from the right
            m_x[m_nhc-(i+1)][m_nhc-(j+1)][0]   = m_x[m_nhc-(i+1)][m_nhc][0];
            m_x[m_nhc-(i+1)][m_nhc-(j+1)][1]   = m_x[m_nhc][m_nhc-(j+1)][1];
            // SE - march right and down, copying x from above and y from the left
            m_x[m_nhc+m_nx+i][m_nhc-(j+1)][0]  = m_x[m_nhc+m_nx+i][m_nhc][0];
            m_x[m_nhc+m_nx+i][m_nhc-(j+1)][1]  = m_x[m_nhc+m_nx-1][m_nhc-(j+1)][1];
            // NW - march left and up, copying x from below and y from the right
            m_x[m_nhc-(i+1)][m_nhc+m_ny+j][0]  = m_x[m_nhc-(i+1)][m_nhc+m_ny-1][0];
            m_x[m_nhc-(i+1)][m_nhc+m_ny+j][1]  = m_x[m_nhc][m_nhc+m_ny+j][1];
            // NE - march right and up, copying x from below and y from the left
            m_x[m_nhc+m_nx+i][m_nhc+m_ny+j][0] = m_x[m_nhc+m_nx+i][m_nhc+m_ny-1][0];
            m_x[m_nhc+m_nx+i][m_nhc+m_ny+j][1] = m_x[m_nhc+m_nx-1][m_nhc+m_ny+j][1];
        }
    }


    // Write file of grid nodes with halo cells;
    if (write) {
        std::string wfpath = m_file + "_whc.dat";
        writeCSVGridWithHalos(wfpath, m_nx, m_ny, m_nhc, m_x);
        wfpath = m_file + "_whc.bin";
        writeBinGridWithHalos(wfpath, m_nx, m_ny, m_nhc, m_x);
    }
    std::cout << "--2D array of grid nodes with halo cells computed." << std::endl;
}

/* Driver method to compute all grid metrics:
 * 1. cell centers;
 * 2. cell areas;
 * 3. projected cell face areas; - TODO
 * 4. grid metrics for generalized coordinates;
 */
void Grid::computeMetrics() 
{
    computeCellCenters             (m_nx, m_ny, m_nhc, m_x, m_xc);
    //print3Darray(gxc_whc, gm_nx+2*m_nhc-1, gm_ny+2*m_nhc-1, 2);
    computeCellAreas               (m_nx, m_ny, m_nhc, m_x, m_a );
    //print2Darray(ga_whc, gm_nx+2*m_nhc-1, gm_ny+2*m_nhc-1);
    computeGeneralCoordinateMetrics(
        m_nx, m_ny, m_nhc, 
        m_xc, m_mgc, m_invj
    );
    //print3Darray(gpcm_whc, gm_nx+2*m_nhc-1, gm_ny+2*m_nhc-1, 4);
    //print3Darray(gpcm_whc, gm_nx+2*m_nhc-1, gm_ny+2*m_nhc-1, 4);
}
