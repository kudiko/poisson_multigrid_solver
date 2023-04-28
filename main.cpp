#include "tests.h"

#include "multi_grid_solver.h"
#include "finite_diff.h"
#include <fstream>

void PrintGridToFile(Lattice3D& lat)
{
    std::ofstream out("output.txt");
    Dimensions3D dims = lat.GetDimensions();

    for (int k = 0; k < dims.ZDim; ++k)
    {
        for (int j = 0; j < dims.YDim; ++j)
        {
            for (int i = 0; i < dims.XDim; ++i)
            {
                Point3D point = lat.GetCoordinates(i, j, k);
                out << point.x << ' ' << point.y << ' ' << point.z << ' ' << lat(i, j, k) << std::endl;
            }
        }
    }
}

void CalculateError(Lattice3D& lat)
{
    std::ofstream out("error.txt");
    Dimensions3D dims = lat.GetDimensions();

    for (int k = 1; k < dims.ZDim - 1; ++k)
    {
        for (int j = 1; j < dims.YDim - 1; ++j)
        {
            for (int i = 1; i < dims.XDim - 1; ++i)
            {
                Point3D point = lat.GetCoordinates(i, j, k);
                double err = Calc_Error(lat, i, j, k);
                out << point.x << ' ' << point.y << ' ' << point.z << ' ' << std::abs(err) << ' ' << std::abs(err/lat(i,j,k)) << std::endl;
            }
        }
    }
}

int main()
{
    Tests tests;
    tests.RunTests();


    Domain3D dom({0, 0, 0}, 5, 5, 5);
    Lattice3D lat(dom, {100, 100, 100});
    lat.SetBoundaryCondition();
    // 200 -> 100 -> 50 -> 25 -> 13 -> 7 -> 4
    MultiGridSolver solve(lat, 5, 0.000001);

    solve.Solve(0);
    PrintGridToFile(lat);

    CalculateError(lat);

    return 0;
}