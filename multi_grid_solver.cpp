#include "multi_grid_solver.h"
#include <iostream>
MultiGridSolver::MultiGridSolver(Lattice3D& grid, int levels, double tolerance) : grid_(grid), levels_(levels), tolerance_(tolerance){};

void MultiGridSolver::Solve(int current_mesh_level)
{

    Gauss_Siedel pre_gs(grid_, pre_smooth_iter_);  //pre-smoothing
    if (init_guess_)
    {
        pre_gs.InitGuess();
        init_guess_ = 0;
    }
    pre_gs.Smoothify();
    std::cout << "Pre-smoothing on level " << current_mesh_level << " done" << std::endl;
    if (current_mesh_level < levels_)
    {
        grid_.Coarsify(); //restrict
        std::cout << "Coarsification on level " << current_mesh_level << " done : " << grid_.GetDimensions().XDim << " points in current grid" << std::endl;

        Solve(current_mesh_level + 1);
        grid_.Finify(); //prolongate
        std::cout << "Finification on level " << current_mesh_level << " done : " << grid_.GetDimensions().XDim << " points in current grid" << std::endl;
    }
    Gauss_Siedel post_gs(grid_, post_smooth_iter_);
    post_gs.Smoothify();
    std::cout << "Post-smoothing on level " << current_mesh_level << " done" << std::endl;

}