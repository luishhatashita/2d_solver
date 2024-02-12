#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <string>
#include "readers.h"
#include "grid.h"

using namespace std;

/*  Method to read Oefelein specific grid file which is formatted as follows:
 *  - first line as dimension of the computational domain, nxi and neta, res-
 *  pectively;
 *  - afterwards, x and y locations, respectively, of all nodes.
 *
 * Parameters:
 * -----------
 *  string fpath  : filename with relative path to the directory from where the 
 *                  executable is being run;
 *  vec*   coords : 2d vector of x and y coordinates of all nodes;
 *  int    nxi    : size of xi dimension;
 *  int    neta   : size of eta dimension;
 */
void readGrid(
    string fpath, 
    vector<vector<double>>* coords, 
    int* nxi, 
    int* neta
) 
{
    string row;
    string item;
    char separator = ',';

    string fname = fpath + ".dat";

    cout << "--Reading: " << fname << endl;
    ifstream in(fname);
    while(getline(in,row))
    {
        vector<double> R;
        stringstream ss(row);
        while (getline(ss, item, separator)) R.push_back(stof(item));
        // Two of them work;
        (*coords).push_back(R);
        //coords->push_back(R);
    }
    in.close();
    cout << "--Finished reading: " << fname << endl;
    *nxi  = (int)(*coords)[0][0];
    *neta = (int)(*coords)[0][1];
    // Can I improve this?
    (*coords).erase((*coords).begin());
    cout << "--Coords and dimensions obtained." << endl;
}

/*  Method to write grid file with hallo cells, with the same origina formating:
 *  - first line as dimension of the computational domain, nxi and neta, res-
 *  pectively;
 *  - afterwards, x and y locations, respectively, of all nodes.
 *
 * Parameters:
 * -----------
 *  string*   wfpath : filename for file to be written;
 *  double*** coords : 3d array of x and y coordinates of all nodes organized
 *                     in xi by eta matrix;
 *  int*      nhc    : number of hallo cells;
 *  int*      nxi    : size of xi dimension;
 *  int*      neta   : size of eta dimension;
 */
void writeGridWithHalos(
    string* wfpath, 
    double*** nodes_whc, 
    int* nhc, 
    int* nxi, 
    int* neta
)
{
    ofstream out;
    out.open((*wfpath));

    out << (*nxi)+2*(*nhc) << "," << (*neta)+2*(*nhc) << "\n";
    // dimensions are now gnxi + 2*nhc by gneta + 2*nhc
    for (int j=0; j<((*neta)+2*(*nhc)); j++) {
        for (int i=0; i<((*nxi)+2*(*nhc)); i++) {
            out << nodes_whc[i][j][0] << "," << nodes_whc[i][j][1] << ",0" << "\n";
        }
    }
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
                cout << arr[i][j] << ",";
            }
        }
        cout << endl;
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
                    cout << arr[i][j][k] << ",";
                }
            }
            cout << endl;
        }
        cout << endl << endl;
    }
}
