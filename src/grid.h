#include <vector>
#include <string>

// Grid class - 2d with additional CFD methods and grid metrics calculations.
class Grid 
{
    private:
        std::string gfpath;
        std::vector<std::vector<float>> gcoords;
        std::vector<std::vector<std::vector<float>>> 
            gnodes, gnodes_whc, gxc_whc, ga_whc;
        int gnxi, gneta;
        int nhc;
    public:
        // constructor
        Grid(std::string cfpath);//, std::vector<std::vector<float>> ccoords, int cnxi, int cneta);
        std::vector<int> getCompDim();
        std::vector<std::vector<std::vector<float>>> getNodes();
        void addHaloCells(int inhc, bool writecsv = false, std::string = "grid_whc.dat");
        void computeCellCenters();
        //void computeCellAreas();
};
