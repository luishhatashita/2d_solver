// 2d CFD solver

#include <iostream>
#include "grid.h"

int main(int argc, char* argv[0]) 
{
    // Parameters:
    int nhc = std::stoi(argv[1]);
    bool writecsv = true;
    // Grid:
    Grid grid("grid.dat");
    grid.addHaloCells(nhc, writecsv);
    grid.computeCellCenters();
    std::cout << "ok" << std::endl;
    return 0;
}

