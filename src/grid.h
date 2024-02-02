#include <vector>
#include <string>

// Grid class - 2d with additional CFD methods and grid metrics calculations.
class Grid 
{
    private:
        std::string gfpath;
        std::vector<std::vector<float>> gcoords;
        std::vector<std::vector<std::vector<float>>> gnodes;
        int gnxi, gneta;
    public:
        // constructor
        Grid(std::string cfpath);//, std::vector<std::vector<float>> ccoords, int cnxi, int cneta);
        std::vector<int> getCompDim();
        std::vector<std::vector<std::vector<float>>> getNodes();
};
