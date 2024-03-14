#include <iostream>
#include "allocate.h"
#include "grid.h"

using namespace std;

/* Method to allocate memory for a 2D array.
 *
 * Parameters:
 * -----------
 *  int      nrows : number of rows;
 *  int      ncols : number of columns;
 *  double** arr   : 2d array of pointers;
 */
void allocate2D(int nrows, int ncols, double**& arr)
{

    // Allocate memory for the array of pointers (rows)
    double** newarr = new double*[nrows];
    // Contiguous memory allocation:
    //newarr[0] = new double[nrows*ncols];

    // Allocate memory for each row (columns)
    for (int i = 0; i < nrows; i++) {
        // Basic allocation
        newarr[i] = new double[ncols];
        // Contiguous memory allocation;
        //newarr[i] = &newarr[0][i*ncols];
    }
    
    // Initialize zeros
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            newarr[i][j] = 0;
        }
    }

    // Assining memory to array to be initialized;
    arr = newarr;
}

/* Method to allocate memory for a 3D array.
 *
 * Parameters:
 * -----------
 *  int       nrows : number of rows;
 *  int       ncols : number of columns;
 *  int       ndeps : number of "depths";
 *  double*** arr   : 3d array of pointers;
 */
void allocate3D(int nrows, int ncols, int ndeps, double***& arr)
{
    // Allocate memory for the array of arrays of pointers (rows)
    double*** newarr = new double**[nrows];
    // Contiguous memory allocation:
    //newarr[0][0] = new double[nrows*ncols*ndeps];

    // Allocate memory for the array of pointers (columns)
    for (int i = 0; i < nrows; i++) {
        newarr[i] = new double*[ncols];
        // Allocate memory for the each column (depths)
        for (int j = 0; j < ncols; j++) {
            // Basic allocation
            newarr[i][j] = new double[ndeps];
            // Contiguous memory allocation:
            //newarr[i][j] = &newarr[0][0][i*ncols*ndeps+j*ncols];
        }
    }

    // Initialize zeros
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            for (int k = 0; k < ndeps; k++) {
                newarr[i][j][k] = 0;
            }
        }
    }

    // Assining memory to array to be initialized;
    arr = newarr;
}

/* Method to allocate memory for a 4D array.
 *
 * Parameters:
 * -----------
 *  int        nrows : number of rows;
 *  int        ncols : number of columns;
 *  int        ndeps : number of "depths";
 *  int        ncomp : number of "complexity";
 *  double**** arr   : 4d array of pointers;
 */
void allocate4D(int nrows, int ncols, int ndeps, int ncomp, double****& arr)
{
    // Allocate memory for the array of arrays of arrays of pointers (rows)
    double**** newarr = new double***[nrows];

    // Allocate memory for the array of arrays of pointers (columns)
    for (int i = 0; i < nrows; i++) {
        newarr[i] = new double**[ncols];
        // Allocate memory for the array of pointers (depths)
        for (int j = 0; j < ncols; j++) {
            newarr[i][j] = new double*[ndeps];
            // Allocate memory for the each complexity (complexity)
            for (int k = 0; k < ndeps; k++) {
                newarr[i][j][k] = new double[ncomp];
            }
        }
    }
    
    // Initialize zeros
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            for (int k = 0; k < ndeps; k++) {
                for (int l = 0; l < ncomp; l++) {
                newarr[i][j][k][l] = 0;
                }
            }
        }
    }

    // Assining memory to array to be initialized;
    arr = newarr;
}

/* Method to deallocate memory for a 2D array.
 *
 * Parameters:
 * -----------
 *  int      nrows : number of rows;
 *  double** arr   : 2d array of pointers;
 */
void deallocate2D(int nrows, double**& arr)
{
    // Free memory for each row
    for (int i = 0; i < nrows; i++) {
        delete[] arr[i];
    }

    // Free memory for the array of pointers
    //delete[] arr[0];
    delete[] arr;
}

/* Method to deallocate memory for a 3D array.
 *
 * Parameters:
 * -----------
 *  int       nrows : number of rows;
 *  int       ncols : number of columns;
 *  double*** arr   : 3d array of pointers;
 */
void deallocate3D(int nrows, int ncols, double***& arr)
{
    // Free memory for each row
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            delete[] arr[i][j];
        }
        delete[] arr[i];
    }

    // Free memory for the array of pointers
    //delete[] arr[0][0];
    delete[] arr;
}

/* Method to deallocate memory for a 4D array.
 *
 * Parameters:
 * -----------
 *  int        nrows : number of rows;
 *  int        ncols : number of columns;
 *  int        ndeps : number of "depths";
 *  double**** arr   : 4d array of pointers;
 */
void deallocate4D(int nrows, int ncols, int ndeps, double****& arr)
{
    // Free memory for each row
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            for (int k = 0; k < ndeps; k++) {
                delete[] arr[i][j][k];
            }
            delete[] arr[i][j];
        }
        delete[] arr[i];
    }

    // Free memory for the array of pointers
    delete[] arr;
}
