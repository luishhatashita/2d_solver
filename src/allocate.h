#ifndef ALLOCATION
#define ALLOCATION
using namespace std;
void   allocate2D(int nrows, int ncols, double**& arr);
void   allocate3D(int nrows, int ncols, int ndeps, double***& arr);
void   allocate4D(int nrows, int ncols, int ndeps, int ncomp, double****& arr);
void deallocate2D(int nrows, double**& arr);
void deallocate3D(int nrows, int ncols, double***& arr);
void deallocate4D(int nrows, int ncols, int ndeps, double****& arr);
#endif
