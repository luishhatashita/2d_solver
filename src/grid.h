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
        double ***gnodes, ***gnodes_whc, ***gxc_whc, **ga_whc, ***ggcm_whc,
               **ginvj_whc;
    public:
        // constructor
        Grid(string cfpath);
        array<int, 2> getCompDim();
        double*** getNodes();
        void addHaloCells(int inhc, bool writecsv = false);
        void computeMetrics();
        double*** getGeneralizedCoordinateMetrics();
        double** getInverseJacobian();
        // destructor
        ~Grid();
};
