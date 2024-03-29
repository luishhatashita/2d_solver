#include "readers.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

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
    std::string fpath, 
    std::vector<std::vector<double>>& coords, 
    int& nx, 
    int& ny
) 
{
    std::string row;
    std::string item;
    char separator = ',';

    std::string fname = fpath + ".dat";
    //double* R = new double[2];
    std::vector<double> R;

    std::cout << "--Reading: " << fname << std::endl;
    std::ifstream in(fname);
    // Read first line;
    std::getline(in,row);
    std::stringstream ss(row);
    while (std::getline(ss, item, separator)) R.push_back(stof(item));
    nx = (int) R[0];
    ny = (int) R[1];

    while(std::getline(in,row))
    {
        // Need to "redefine" R such that one can "overwrite" it;
        std::vector<double> R;
        std::stringstream ss(row);
        while (std::getline(ss, item, separator)) R.push_back(stof(item));
        // Two of them work;
        coords.push_back(R);
        //coords->push_back(R);
    }
    in.close();
    std::cout << "--Finished reading: " << fname << std::endl;
    std::cout << "--Coords and dimensions obtained." << std::endl;

    //delete[] R;
}
