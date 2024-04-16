#ifndef CFD_SCHEMES_SOLUTION_H
#define CFD_SCHEMES_SOLUTION_H

#include "grid.h"
#include "parameters.h"

class Solution
{
    private:
        bool   m_write;
        int    m_nx, m_ny, m_nhc, s_nit;
        double ***m_Qn,   ***m_Qvn,    ***m_Qnp1, ***m_Qvnp1,
               ***m_QnuL, ***m_QnuR,   ***m_Qnu,
               ***m_QnvL, ***m_QnvR,   ***m_Qnv,
               ***m_Ehn,  ***m_EhnL,   ***m_EhnR, ***m_AQu,
               ***m_Fhn,  ***m_FhnL,   ***m_FhnR, ***m_BQv,
                **m_dt,
                **m_L2,    **m_Linfty;
        //Grid       *s_grid;
        //Parameters *s_par;
    public:
        // constructor
        Solution(const Grid& grid, const Parameters& par, int nrest, bool write=false);
        // restart reader;
        // restart writer;
        void writePrimitives(const Parameters& par, const Grid& grid, int it);
        void writeConservatives(const Grid& grid, int it);
        // Shock-capturing schemes
        void constructLeftRightStates(double kappa, double epsilon);
        void computeFluxes(const Parameters& par, const Grid& grid);
        // Time integration
        void computeTimeSteps(const Parameters& par, const Grid& grid, double& dt_max);
        void checkTimeStep(const Parameters& par, const Grid& grid, double dt);
        void integrateLocalTime(const Parameters& par, const Grid& grid);
        void integrateFixedTime(const Parameters& par, const Grid& grid);
        void updateVectors();
        // Error norms
        void computeErrorNorms(int it);
        void logErrorNorms(int it);
        void writeErrorNorms();
        // destructor
        ~Solution();
};

#endif // CFD_SCHEMES_SOLUTION_H
