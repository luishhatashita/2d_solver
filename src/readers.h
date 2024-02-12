#ifndef READERS
#define READERS
#include <vector>
#include <string>
using namespace std;
// Method to read .dat file from Oefelein specific grid.
void readGrid(
        string fpath, 
        vector<vector<double>>* coords, 
        int* nxi, 
        int* neta
);
void writeGridWithHalos(
        string* wfpath, 
        double*** nodes_whc, 
        int* nhc, 
        int* nxi, 
        int* neta
);
#endif

