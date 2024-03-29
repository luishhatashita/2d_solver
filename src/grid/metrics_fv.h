#ifndef CFD_GRID_METRICS_FV_H_
#define CFD_GRID_METRICS_FV_H_

void computeCellCenters(
    int nx, int ny, int nhc, 
    double***& x, double***& xc
);

#endif
