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
 *  string      fpath  : 
 *      filename with relative path to the directory from where the executable 
 *      is being run;
 *  vec pointer coords :
 *      2d vector of x and y coordinates of all nodes;
 *  int         nxi    : 
 *      size of xi dimension;
 *  int         neta   :
 *      size of eta dimension;
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

void writeGridWithHalos(
        string* wfpath, 
        //double (*nodes_whc)[][gneta][2],
        //array<array<array<double,2>, gneta>, gnxi>* nodes_whc, 
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
            //out << gnodes_whc[i][j][0] << "," << gnodes_whc[i][j][1] << ",0" << "\n";
            //out << (*nodes_whc)[i][j][0] << "," << (*nodes_whc)[i][j][1] << ",0" << "\n";
            out << nodes_whc[i][j][0] << "," << nodes_whc[i][j][1] << ",0" << "\n";
        }
    }
    out.close();
}
