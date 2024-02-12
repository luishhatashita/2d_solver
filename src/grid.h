#include <vector>
#include <array>
#include <string>

using namespace std;

// Grid class - 2d with additional CFD methods and grid metrics calculations.
class Grid 
{
    private:
        string gfpath;
        vector<vector<double>> gcoords;
        int gnxi, gneta;
        int nhc;
        //vector<vector<vector<float>>> 
        //    gnodes, gnodes_whc, gxc_whc, ga_whc;
        double ***gnodes, ***gnodes_whc, ***gxc_whc, ***ga_whc;
    public:
        // constructor
        Grid(string cfpath);//, std::vector<std::vector<float>> ccoords, int cnxi, int cneta);
        array<int, 2> getCompDim();
        double*** getNodes();
        //vector<vector<vector<double>>> getNodes();
        void addHaloCells(int inhc, bool writecsv = false);
        void computeCellCenters();
        //void computeCellAreas();
        ~Grid();
};
