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
    // Implementation with temporary variables;
    //std::vector<std::vector<float>> ccoords; 
    //int cnxi, cneta;
    //readGrid(cfpath, &ccoords, &cnxi, &cneta);
    //gcoords = ccoords;
    //gnxi    = cnxi;
    //gneta   = cneta;
}

// Method to get computational dimension and size.
std::vector<int> Grid::getCompDim()
{
    std::vector<int> acompdim;
    acompdim.push_back(gnxi);
    acompdim.push_back(gneta);
    return acompdim;
}
