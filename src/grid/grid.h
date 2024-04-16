#ifndef CFD_GRID_GRID_H_
#define CFD_GRID_GRID_H_

#include <vector>
#include <array>
#include <string>

//using namespace std;

// Grid class - 2d with additional CFD methods and grid metrics calculations.
// TODO:
// -----
// - rename class variables;
// - compute FV metrics instead of generalized coordinate ones;
class Grid 
{
    private:
        std::string m_file;
        bool m_write;
        std::vector<std::vector<double>> m_coords;
        int m_nx, m_ny, m_nhc;
        double ***m_x_woh, ***m_x, ***m_xc, ***m_xu, ***m_xv,
               **m_a, ***m_mgc, **m_invj,
               **m_v, ***m_su, ***m_sv;
    public:
        // constructor
        Grid(std::string file, int nhc, bool write = false);
        //std::array<int, 2> getCompDim();
        void      getArraySizes(int*& dims) const;
        double*** getNodes();
        void      addHaloCells(int inhc);
        void      computeGeneralMetrics();
        double*** getGeneralizedCoordinateMetrics();
        double**  getInverseJacobian();
        void      computeFiniteVolumeMetrics();
        double*** getProjectedFaceAreasXi() const;
        double*** getProjectedFaceAreasEta() const;
        double**  getCellVolumes() const;
        // destructor
        ~Grid();
};

#endif // CFD_GRID_GRID_H_
