// 2d CFD solver

#include <iostream>
#include "grid.h"

// using directive to std identifiers for brevity;
using namespace std;

int main(int argc, char* argv[0]) 
{
    // Parameters:
    int nhc = 1;
    //int nhc = stoi(argv[1]);
    bool writecsv = true;
    // Grid:
    string fname = "g33x22u";
    //string fname = argv[2];
    Grid grid(fname);
    grid.addHaloCells(nhc, writecsv);
    grid.computeCellCenters();
    cout << "ok" << endl;
    return 0;
}

