// 2d CFD solver

#include <iostream>
#include "grid/grid.h"
#include "parameters/parameters.h"
#include "schemes/solution.h"
//#include "schemes/schemes.h"

int main(int argc, char* argv[0]) 
{
    // Parameters:
    int nhc = 1;
    bool write = true;
    double dt_max;
    Parameters sim_par = {
        //.rt    = {.nrest = 1,     .nit     = 11, .nlog=1, .CFL   = 0.01,},
        //.rt    = {.nrest = 10,     .nit     = 101, .nlog=10, .CFL   = 0.5, .dt = 1.0e-7},
        //.rt    = {.nrest = 100,     .nit     = 1001, .nlog=10, .CFL   = 0.01, .dt = 1.0e-7},
        .rt    = {.nrest = 1000,     .nit     = 10001, .nlog  = 100, .CFL   = 0.25, .dt = 4.0e-5},
        //.rt    = {.nrest = 5000,     .nit     = 100001, .nlog  = 100, .CFL   = 0.0001, .dt = 1.0e-6},
        .ref   = {.pref  = 101325.0, .uref    = 694.4, .Tref  = 300.0, 
                  .lref  = 1.0, .Mref = 2.0},
        .bcs   = {.s = 5, .e = 4, .n = 5, .w = 2,},
        .td    = {.cp    = 1005.0,   .R       = 287.0, .gamma = 1.4,
                  .c = 347.2},
        .muscl = {.kappa = 1.0,      .epsilon = 0.0}, // first order upwind;
        //.muscl = {.kappa = -1.0,     .epsilon = 1.0}, // second order accurate;
        //.muscl = {.kappa = 0.5,      .epsilon = 1.0}, // third order accurate;
        .ausm  = {.beta  = 0.125,    .kp      = 0.25,  .ku    = 0.75,
                  .sigma = 1.0},
    }; 

    // Grid:
    //std::string fname = "g33x25u";
    std::string fname = "g65x49u";
    Grid grid(fname, nhc, write);

    // Solution variable for schemes:
    Solution solution(grid, sim_par, -1, write);

    // Temporal loop:
    for (int i=1; i<sim_par.rt.nit; i++) {
        solution.setIteration(i);
        /*
        // Flux construction at cell faces;                
        solution.constructLeftRightStates(               // Will replace for
            sim_par.muscl.kappa, sim_par.muscl.epsilon   // construction of
        );                                               // xi and eta faces
        // Shock capturing scheme;                       // separately.    
        solution.computeFluxes(sim_par, grid);
        */
        //solution.computeXiFirstOrderFluxes(sim_par, grid);
        //solution.computeEtaFirstOrderFluxes(sim_par, grid);
        //solution.computeXiRoeFluxes(sim_par, grid);
        //solution.computeEtaRoeFluxes(sim_par, grid);
        solution.computeXiAUSMFluxes(sim_par, grid);
        solution.computeEtaAUSMFluxes(sim_par, grid);
        // Time step computation/check;
        solution.computeTimeSteps(sim_par, grid, dt_max);
        //solution.checkTimeStep(sim_par, grid, sim_par.rt.dt);
        // Time integration;
        solution.integrateLocalTime(sim_par, grid);
        //solution.integrateFixedTime(sim_par, grid);
        // Compute error metrics:
        solution.computeErrorNorms(sim_par);
        if (i % sim_par.rt.nlog == 0) {
            solution.log();
        }
        // Update values and write if specified.
        solution.updateVectors();
        if (i % sim_par.rt.nrest == 0) {
            solution.writePrimitives(sim_par, grid, i);
            solution.writeConservatives(grid, i);
        }
        //if (i == 1){
        //    break;
        //}
    }
    solution.writeErrorNorms();
    std::cout << "end" << std::endl;
    return 0;
}
