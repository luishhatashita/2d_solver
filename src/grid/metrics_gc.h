#ifndef CFD_GRID_METRICS_GC_H_
#define CFD_GRID_METRICS_GC_H_

void computeCellAreas(
    int nx, int ny, int nhc,
    double***& x, double**& a
);
void computePhysicalCoordinateMetrics(
    int nx, int ny, int nhc,
    double***& xc, double***& mpc 
); 
void computeGeneralCoordinateMetrics(
    int nx, int ny, int nhc,
    double***& xc, double***& mgc, double**& invj 
);

#endif
