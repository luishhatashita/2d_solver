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
        //.rt    = {.nrest = 1,     .nit     = 11, .CFL   = 0.01,},
        //.rt    = {.nrest = 10,     .nit     = 101, .CFL   = 0.01, .dt = 1.0e-7},
        //.rt    = {.nrest = 100,     .nit     = 1001, .CFL   = 1.0, .dt = 1.0e-7},
        .rt    = {.nrest = 1000,     .nit     = 10001, .CFL   = 0.61, .dt = 4.0e-5},
        //.rt    = {.nrest = 10000,     .nit     = 100001, .CFL   = 0.000001, .dt = 1.0e-6},
        .ref   = {.pref  = 101325.0, .uref    = 694.4, .Tref  = 300.0, 
                  .lref  = 1.0,},
        .bcs   = {.s = 5, .e = 4, .n = 5, .w = 2,},
        .td    = {.cp    = 1005.0,   .R       = 287.0, .gamma = 1.4,
                  .c = 347.2},
        .muscl = {.kappa = 1.0,      .epsilon = 0.0},
    }; 

    // Grid:
    std::string fname = "g33x25u";
    //std::string fname = "g65x49u";
    Grid grid(fname, nhc, write);

    // Solution variable for schemes:
    Solution solution(grid, sim_par, -1, write);

    // Temporal loop:
    for (int i=1; i<sim_par.rt.nit; i++) {
        solution.setIteration(i);
        // Flux construction at cell faces;
        solution.constructLeftRightStates(
            sim_par.muscl.kappa, sim_par.muscl.epsilon
        );
        // Shock capturing scheme;
        solution.computeFluxes(sim_par, grid);
        // Time step computation/check;
        solution.computeTimeSteps(sim_par, grid, dt_max);
        //solution.checkTimeStep(sim_par, grid, sim_par.rt.dt);
        // Time integration;
        solution.integrateLocalTime(sim_par, grid);
        //solution.integrateFixedTime(sim_par, grid);
        // Compute error metrics:
        solution.computeErrorNorms(sim_par);
        if (i % 50 == 0) {
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
