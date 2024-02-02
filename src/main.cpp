// 2d CFD solver

#include <iostream>
#include "grid.h"

int main() 
{
    //std::vector<std::vector<float>> coords;
    //int nxi, neta;
    //readGrid("grid.dat", &coords, &nxi, &neta);
    Grid gcoarse("grid.dat");
    std::cout << gcoarse.getNodes()[0][0][0] << std::endl;
    std::cout << "ok" << std::endl;
    return 0;
}

