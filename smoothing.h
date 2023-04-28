#pragma once
#include "grid.h"

class Gauss_Siedel
{
	public:
	Gauss_Siedel(Lattice3D &lat, int steps);
	void InitGuess();
	void Smoothify();

    private:
	Lattice3D& lat_;
	int steps_;
	double h_;

	void Step();
};