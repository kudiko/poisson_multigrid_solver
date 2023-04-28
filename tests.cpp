#include "tests.h"

#include <cassert>
#include <cmath>
#include <iostream>

#include "grid.h"

void Tests::TestGrid()
{

    {
        Domain3D dom({0,0,0}, 3,3,3);
        Lattice3D lat(dom, {50,50,50});

        lat(0,0,0) = 2;
        assert(lat(0,0,0) == 2);

        lat(25, 25, 25) = 14;
        assert(lat(25, 25, 25) == 14);

        lat(49, 49, 49) = 49;
        assert(lat(49, 49, 49) == 49);

        lat(20, 30, 40) = 42;
        assert(lat(20, 30, 40) == 42);


    }
    /*
    {
        Domain3D dom({0,0,0}, 1,1,1);
        Lattice3D lat(dom, {3,3,3});

        lat.Fill(1);
        //lat.SetBoundaryCondition();
       // lat.print();
       // lat.Finify();
        lat.Finify();
       assert(lat.points_.GetDimensions().XDim == 5);
       assert(lat.points_.GetDimensions().YDim == 5);
       assert(lat.points_.GetDimensions().ZDim == 5);

       //lat.print();
        assert(lat(0, 0, 0) == 1);
        assert(lat(2, 2, 2) == 1);
        assert(lat(0, 2, 4) == 1);
        
    }
    */
    /*
    // boundary set test
    {
        Domain3D dom({0,0,0}, 1,1,1);
        Lattice3D lat(dom, {3,3,3});

        lat.Fill(0.);
        lat.SetBoundaryCondition();
        assert(lat.h_ == 1);
        assert(std::abs(lat(2, 1, 1) - 1) < std::pow(10, -6));
        assert(std::abs(lat(2, 2, 2) - std::sqrt(3)) < std::pow(10, -6));

    }
    */
    {
        Domain3D dom({0,0,0}, 1,1,1);
        Lattice3D lat(dom, {4,4,4});
        lat.Fill_Interior(1.);
        lat.SetBoundaryCondition();
        //lat.Fill_Interior(0);
        lat.Finify();
        //lat.print();
       // lat.InterpolationProcedure();
        //lat.print();
    }

}