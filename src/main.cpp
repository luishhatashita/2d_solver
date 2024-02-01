// 2d CFD solver

#include "readers.h"
#include <iostream>
//#include <fstream>
//#include <sstream>
//#include <vector>
//#include <string>
//#include <grid.h>

//void readGrid(std::string fpath, std::vector<std::vector<float>>* coords, int* nxi, int* neta) 
//{
//    std::string row;
//    std::string item;
//    char separator = ',';
//
//    std::ifstream in(fpath);
//    while(getline(in,row))
//    {
//        std::vector<float> R;
//        std::stringstream ss(row);
//        while (getline(ss, item, separator)) R.push_back(std::stof(item));
//        // Two of them work;
//        (*coords).push_back(R);
//        //coords->push_back(R);
//    }
//    in.close();
//    *nxi  = (int)(*coords)[0][0];
//    *neta = (int)(*coords)[0][1];
//    // Can I improve this?
//    (*coords).erase((*coords).begin());
//}

int main() 
{
    std::cout << "teste" << std::endl;
    std::vector<std::vector<float>> coords;
    int nxi, neta;
    readGrid("grid.dat", &coords, &nxi, &neta);
    std::cout << coords[0][0] << std::endl;
    return 0;
}

