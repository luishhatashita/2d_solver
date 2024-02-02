#ifndef READERS
#define READERS
#include <vector>
#include <string>
// Method to read .dat file from Oefelein specific grid.
void readGrid(std::string fpath, std::vector<std::vector<float>>* coords, int* nxi, int* neta);
#endif
