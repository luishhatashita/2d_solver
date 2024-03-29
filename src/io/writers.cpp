#include "writers.h"

#include <iostream>
#include <fstream>
#include <string>

/*  Method to write grid file with hallo cells, with the same origina formating:
 *  - first line as dimension of the computational domain, nxi and neta, res-
 *  pectively;
 *  - afterwards, x and y locations, respectively, of all nodes.
 *
 * Parameters:
 * -----------
 *  string*   wfpath : filename for file to be written;
 *  int*      nxi    : size of xi dimension;
 *  int*      neta   : size of eta dimension;
 *  int*      nhc    : number of hallo cells;
 *  double*** x_h    : 3d array of x and y coordinates of all nodes organized
 *                     in xi by eta matrix;
 */
void writeCSVGridWithHalos(
    std::string wfpath, 
    int nx, int ny, int nhc, 
    double***& x_h
)
{
    std::ofstream out;
    out.open(wfpath);

    out << nx+2*nhc << "," << ny+2*nhc << "\n";
    // dimensions are now gnxi + 2*nhc by gneta + 2*nhc
    for (int j=0; j<(ny+2*nhc); j++) {
        for (int i=0; i<(nx+2*nhc); i++) {
            out << x_h[i][j][0] << "," << x_h[i][j][1] << ",0" << "\n";
        }
    }
    out.close();
}

/*  Method to write grid file with hallo cells as 3d array: nx, ny, 2 (for x 
 *  and y).
 *
 * Parameters:
 * -----------
 *  string     fname : reference of string literal of file name;
 *  int        nx    : number of rows;
 *  int        ny    : number of columns;
 *  int        ny    : number of halo cells;
 *  double***& x_h   : 3d array of x and y coordinates of all nodes organized
 *                     in xi by eta matrix;
 *
 * TODO:
 * -----
 *  - the way this is going to output is in form of a ravel row after row, I 
 *  would still need to reshape the array in python for example, however, that
 *  required prior knowledge of the dimensions of the array. Hence, ideally I
 *  should output in some other way that the metadata or the shape of the array
 *  is still written to the output.
 */
void writeBinGridWithHalos(
    std::string fname, 
    int nx, int ny, int nhc,
    double***& x_h
)
{
    //std::cout << "writeBinaryArray" << std::endl;
    
    // Open a binary file for writing
    std::ofstream out(fname, std::ios::binary);

    if (!out.is_open()) {
        std::cerr << "Error opening file for writing." << std::endl;
        return;
    }

    // Write the array to the binary file
    //out.write(reinterpret_cast<const char*>(arr), (nx*ny)*sizeof(double));

    // Write each row of the 2D array to the binary file
    for (int i=0; i<(nx+2*nhc); i++) {
        for (int j=0; j<(ny+2*nhc); j++) {
            out.write(reinterpret_cast<const char*>(x_h[i][j]), 2*sizeof(double));
        }
    }

    // Close the file
    out.close();
}

/*  Method to print 2D array on console.
 *
 * Parameters:
 * -----------
 *  double** arr  : 2d array of a quantity for all nodes organized in a xi by 
 *                  eta matrix;
 *  int      dim1 : size of first dimension;
 *  int      dim2 : size of second dimension;
 */
void print2Darray(double** arr, int dim1, int dim2)
{
    for (int j=0; j<dim2; j++) {
        for (int i=0; i<dim1; i++) {
            if (i % 5 == 0 && j % 5 == 0) {
                std::cout << arr[i][j] << ",";
            }
        }
        std::cout << std::endl;
    }
}

/*  Method to print 3D array on console.
 *
 * Parameters:
 * -----------
 *  double*** arr  : 3d array of a quantity for all nodes organized in a xi by 
 *                   eta matrix;
 *  int       dim1 : size of first dimension;
 *  int       dim2 : size of second dimension;
 *  int       dim3 : size of third dimension;
 */
void print3Darray(double*** arr, int dim1, int dim2, int dim3)
{
    for (int k=0; k<dim3; k++) {
        for (int j=0; j<dim2; j++) {
            for (int i=0; i<dim1; i++) {
                if (i % 5 == 0 && j % 5 == 0) {
                    std::cout << arr[i][j][k] << ",";
                }
            }
            std::cout << std::endl;
        }
        std::cout << std::endl << std::endl;
    }
}
