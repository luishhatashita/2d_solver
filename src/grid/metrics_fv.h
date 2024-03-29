#ifndef CFD_GRID_METRICS_FV_H_
#define CFD_GRID_METRICS_FV_H_

#include <string>

void computeCellCenters(
    int nx, int ny, int nhc, 
    double***& x, double***& xc,
    std::string fpath, bool write
);
void computeFaceCenters(
    int nx, int ny, int nhc, 
    double***& x, double***& xu, double***& xv,
    std::string fpath, bool write
);
void computeProjectedFaceAreas(
    int nx, int ny, int nhc, 
    double***& x, double***& su, double***& sv,
    std::string fpath, bool write
);
void computeCellVolumes(
    int nx, int ny, int nhc, 
    double***& x, double**& v,
    std::string fpath, bool write
);

#endif
