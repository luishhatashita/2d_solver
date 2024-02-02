#include <iostream>
#include <vector>
#include <string>
#include "readers.h"
#include "grid.h"

/* Constructor method that reads the Oefelein formatted grid file based on spe-
 * cified path.
 *
 * Parameters:
 * -----------
 *  string cfpath : relative path + filename of grid file.
 */
Grid::Grid(std::string cfpath) 
{
    // Implementation using the pointers properly;
    gfpath  = cfpath;
    readGrid(cfpath, &gcoords, &gnxi, &gneta);
    // Initialized vector, however had to "allocate" it in some way;
    gnodes = std::vector<std::vector<std::vector<float>>> (
            gnxi, std::vector<std::vector<float>>(gneta, std::vector<float>(2))
            );
    for (int j=0; j<gneta; j++) {
        for (int i=0; i<gnxi; i++) {
            gnodes[i][j][0] = gcoords[i+gnxi*j][0];
            gnodes[i][j][1] = gcoords[i+gnxi*j][1];
        }
    }
}

// Method to get computational dimension and size.
std::vector<int> Grid::getCompDim()
{
    std::vector<int> acompdim;
    acompdim.push_back(gnxi);
    acompdim.push_back(gneta);
    return acompdim;
}

// Method to get nodes as 3d vector.
std::vector<std::vector<std::vector<float>>> Grid::getNodes()
{
    return gnodes;
}
