#pragma once
#include "grid.h"
double Eq_rhs (double x, double y, double z);
double Boundary_function (double a, double b, double c);
double Laplace_operator(int i, int j, int k);
double Calc_Error(Lattice3D& lat, int i, int j, int k);