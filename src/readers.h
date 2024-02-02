#ifndef READERS
#define READERS
#include <vector>
#include <string>
// Method to read .dat file from Oefelein specific grid.
void readGrid(
        std::string fpath, 
        std::vector<std::vector<float>>* coords, 
        int* nxi, 
        int* neta
);
void writeGridWithHalos(
        std::string* wfpath, 
        std::vector<std::vector<std::vector<float>>>* nodes_whc, 
        int* nhc, 
        int* nxi, 
        int* neta
);
#endif

