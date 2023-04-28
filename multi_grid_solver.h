#pragma once
#include "grid.h"
#include "smoothing.h"

class MultiGridSolver
{
    public:
    MultiGridSolver(Lattice3D& grid, int levels, double tolerance);

    void Solve(int current_mesh_level);

    private:
    Lattice3D& grid_;

    int levels_;
    double tolerance_;
    int pre_smooth_iter_ = 100;
    int post_smooth_iter_ = 100;
    int max_steps_ = 1000;

    bool init_guess_ = 1;
    
};