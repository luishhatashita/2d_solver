// 2d CFD solver

#include <iostream>
#include "grid.h"

int main() 
{
    std::cout << "teste" << std::endl;
    //std::vector<std::vector<float>> coords;
    //int nxi, neta;
    //readGrid("grid.dat", &coords, &nxi, &neta);
    Grid gcoarse("grid.dat");
    std::cout << gcoarse.getCompDim()[0] << std::endl;
    return 0;
}

