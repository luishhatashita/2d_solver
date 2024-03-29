#ifndef CFD_IO_WRITERS_H_
#define CFD_IO_WRITERS_H_

#include <string>

void writeCSVGridWithHalos(
    std::string wfpath, 
    int nx, 
    int ny,
    int nhc, 
    double***& x_h 
);
void writeBinGridWithHalos(
    std::string wfpath, 
    int nx, 
    int ny,
    int nhc, 
    double***& x_h 
);
void print2Darray(
    double** arr,
    int dim1,
    int dim2
);
void print3Darray(
    double*** arr,
    int dim1,
    int dim2,
    int dim3
);

#endif
