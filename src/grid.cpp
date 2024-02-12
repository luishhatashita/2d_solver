#include <iostream>
#include <vector>
#include <array>
#include <string>
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
    readGrid(cfpath, &gcoords, &gnxi, &gneta);
    // Initialized vector, however had to "allocate" it in some way;
    //gnodes = vector<vector<vector<float>>> (
    //        gnxi, vector<vector<float>>(gneta, vector<float>(2))
    //        );
    //double gnodes[gnxi][gneta][2];
    //double gnodes_whc[gnxi+2*nhc][gneta+2*nhc][2];
    //double gxc_whc[gnxi+2*nhc-1][gneta+2*nhc-1][2],
    //       ga_whc[gnxi+2*nhc-1][gneta+2*nhc-1][2];
    allocate3D(gnxi, gneta, 2, gnodes);
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
    deallocate3D(gnxi, gneta, gnodes);
    deallocate3D(gnxi+2*nhc, gneta+2*nhc, gnodes_whc);
    deallocate3D(gnxi+2*nhc-1, gneta+2*nhc-1, gxc_whc);
}

// Method to get computational dimension and size.
array<int, 2> Grid::getCompDim()
{
    array<int, 2> acompdim;
    acompdim[0] = gnxi;
    acompdim[1] = gneta;
    return acompdim;
}

// Method to get nodes as 3d vector.
//vector<vector<vector<double>>> Grid::getNodes()
double*** Grid::getNodes()
{
    return gnodes;
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

// Method to compute cell center based on padded nodes.
void Grid::computeCellCenters() 
{
    // Allocate memory;
    allocate3D(gnxi+2*nhc-1, gneta+2*nhc-1, 2, gxc_whc);
    for (int j=0; j<(gneta+2*nhc-1); j++) {
        for (int i=0; i<(gnxi+2*nhc-1); i++) {
            gxc_whc[i][j][0] = 0.5*(gnodes_whc[i][j][0] + gnodes_whc[i+1][j][0]);
            gxc_whc[i][j][1] = 0.5*(gnodes_whc[i][j][1] + gnodes_whc[i][j+1][1]);
        }
    }
    cout << "--2D array of grid cells centers computed." << endl;
}
