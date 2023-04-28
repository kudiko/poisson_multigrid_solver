#include "finite_diff.h"
#include <cmath>
#include "grid.h"

#define x_0 2
#define y_0 0.5
#define z_0 0.5

double Eq_rhs (double a, double b, double c)
{
    return 4 * M_PI * std::exp(- a * a / x_0 / x_0 - b * b / y_0 / y_0 - c * c / z_0 / z_0);
}

double Boundary_function (double a, double b, double c)
{
    return -1.0 * x_0 * y_0 * z_0 * std::pow(M_PI, 3.0 / 2)/ std::sqrt(a * a + b * b + c * c);
   //return 0;
}

double Laplace_operator(Lattice3D& lat, int i, int j, int k, double h)
{
    return 1 / h / h * ( (lat(i + 1, j, k ) - 2 * lat(i, j, k) + lat(i - 1, j, k)) + 
    (lat(i, j + 1, k) - 2 * lat(i, j, k) + lat(i, j - 1, k)) + 
    (lat(i, j, k + 1) - 2 * lat(i, j, k) + lat(i, j, k - 1)) );
}

double Calc_Error(Lattice3D& lat, int i, int j, int k)
{
    Point3D point = lat.GetCoordinates(i, j, k);
    Point3D another_point = lat.GetCoordinates(i + 1, j, k);
    return Laplace_operator(lat, i, j, k, another_point.x - point.x) - Eq_rhs(point.x, point.y, point.z);
}