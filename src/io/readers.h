#ifndef CFD_IO_READERS_H_
#define CFD_IO_READERS_H_

#include <vector>
#include <string>

// Method to read .dat file from Oefelein specific grid.
void readGrid(
    std::string fpath, 
    std::vector<std::vector<double>>& coords, 
    int& nxi, 
    int& neta
);

#endif // CFD_IO_READERS_H_
