// 2d CFD solver

#include <iostream>
#include "grid/grid.h"

// using directive to std identifiers for brevity;
//using namespace std;

int main(int argc, char* argv[0]) 
{
    // Parameters:
    int nhc = 1;
    //int nhc = stoi(argv[1]);
    bool write = true;
    // Grid:
    std::string fname = "g33x25u";
    //string fname = argv[2];
    // Ideally add everything implicitly in the constructor.
    Grid grid(fname, nhc, write);
    //grid.addHaloCells(nhc, writecsv);
    //grid.computeMetrics();
    std::cout << "ok" << std::endl;
    return 0;
}
