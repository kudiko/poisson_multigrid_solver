#include "smoothing.h"

#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>

#include "grid.h"
#include "finite_diff.h"

Gauss_Siedel::Gauss_Siedel(Lattice3D &lat, int steps) : lat_(lat), steps_(steps)
{
	Point3D deltas = lat.GetDeltas();
	assert(std::abs(deltas.x - deltas.y) < std::pow(10, -6) && std::abs(deltas.x - deltas.z) < std::pow(10, -6));
	h_ = deltas.x;
}

void Gauss_Siedel::InitGuess()
{
	Dimensions3D dims = lat_.GetDimensions();

	for (int i = 1; i < dims.XDim - 1; ++i)
	{
		for (int j = 1; j < dims.YDim - 1; ++j)
		{
			for (int k = 1; k < dims.ZDim - 1; ++k)
			{
				lat_(i, j, k) = 0;
			}
		}
	}
}

void Gauss_Siedel::Smoothify()
{
	for (int i = 0; i < steps_; ++i)
	{
		Step();
	}
}

void Gauss_Siedel::Step()
{
	Dimensions3D dims = lat_.GetDimensions();
	for (int i = 1; i < dims.XDim - 1; ++i)
	{
		for (int j = 1; j < dims.YDim - 1; ++j)
		{
			for (int k = 1; k < dims.ZDim - 1; ++k)
			{
				lat_(i, j, k) = (-h_ * h_ / 6.0) * Eq_rhs(lat_.GetCoordinates(i, j, k).x,
				lat_.GetCoordinates(i, j, k).y, lat_.GetCoordinates(i, j, k).z) +
				1 / 6.0 * (lat_(i + 1, j, k) + lat_(i - 1, j, k) + lat_(i, j + 1, k) + 
				lat_(i, j - 1, k) + lat_(i, j, k + 1) + lat_(i, j, k - 1));
			}
		}
	}
}