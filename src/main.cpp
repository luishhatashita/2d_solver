// 2d CFD solver

#include <iostream>
#include "grid.h"

int main() 
{
    // Parameters:
    int nhc = 2;
    bool writecsv = true;
    // Grid:
    Grid gcoarse("grid.dat");
    gcoarse.addHaloCells(nhc, writecsv);
    std::cout << "ok" << std::endl;
    return 0;
}

