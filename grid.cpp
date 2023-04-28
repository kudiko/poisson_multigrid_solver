#include "grid.h"

#include <utility>
#include <cmath>
#include <iostream>
#include <cassert>

#include "finite_diff.h"

Array3D::Array3D(Dimensions3D dimensions) : dimensions_(dimensions)
{
    points_ = new double[dimensions_.XDim * dimensions_.YDim * dimensions_.ZDim];
}

double& Array3D::operator()(int i, int j, int k)
{
    return *(points_ + (i * (dimensions_.YDim * dimensions_.ZDim) + j * dimensions_.ZDim + k));
}

Dimensions3D Array3D::GetDimensions() const
{
    return dimensions_;
}

void Array3D::swap(Array3D& other)
{
    std::swap(points_, other.points_);
    std::swap(dimensions_, other.dimensions_);
}

Array3D::~Array3D()
{
    delete[] points_;
}


Lattice3D::Lattice3D(Domain3D domain, Dimensions3D dimensions) : 
    domain_(domain), points_(dimensions)
    {
        assert(points_.GetDimensions().XDim == points_.GetDimensions().YDim &&
        points_.GetDimensions().XDim == points_.GetDimensions().ZDim);

        h_ = (domain_.x_dom_.re - domain_.x_dom_.le) / (dimensions.XDim - 1);
        double check_y = (domain_.y_dom_.re - domain_.y_dom_.le) / (dimensions.YDim - 1);
        double check_z = (domain_.z_dom_.re - domain_.z_dom_.le) / (dimensions.ZDim - 1);

        assert(std::abs(h_ - check_y) < std::pow(10, -6) && 
        std::abs(h_ - check_z) < std::pow(10, -6));
    };

void Lattice3D::Finify() // coarse to fine step
{
    Dimensions3D dims = points_.GetDimensions();
    Lattice3D new_lattice(domain_, {dims.XDim * 2 - 1, dims.YDim * 2 - 1, dims.ZDim * 2 - 1});
    new_lattice.SetBoundaryCondition();
    new_lattice.InterpolationProcedure(*this);
    swap(new_lattice);
}

void Lattice3D::Coarsify() // fine to coarse step
{
    Dimensions3D dims = points_.GetDimensions();
    Lattice3D new_lattice(domain_, {(dims.XDim + 1) / 2, (dims.YDim + 1) / 2, (dims.ZDim + 1) / 2});
    new_lattice.SetBoundaryCondition();
    for (int i = 1; i < (dims.XDim + 1) / 2 - 1; ++i)
    {
        for (int j = 1; j < (dims.YDim + 1) / 2 - 1; ++j)
        {
            for (int k = 1; k < (dims.ZDim + 1) / 2 - 1; ++k)
            {
                new_lattice(i, j, k) = 1 / 2.0 * (this->operator()(2 * i, 2 * j, 2 * k)) +
                1 / 12.0 * ((this->operator()(2 * i + 1, 2 * j, 2 * k)) + (this->operator()(2 * i - 1, 2 * j, 2 * k)) 
                + (this->operator()(2 * i, 2 * j + 1, 2 * k)) + (this->operator()(2 * i, 2 * j - 1, 2 * k))
                + (this->operator()(2 * i, 2 * j, 2 * k + 1)) + (this->operator()(2 * i, 2 * j, 2 * k - 1)));
            }
        }
    }
    swap(new_lattice);
}

void Lattice3D::InterpolationProcedure(Lattice3D& old)
{
    /*
    Dimensions3D dims = points_.GetDimensions();
    for (int k = 2; k < dims.ZDim - 1 ; k += 2)
    {
        for (int j = 2; j < dims.YDim - 1; j += 2)
        {
            for (int i = 1; i < dims.XDim - 1; ++i)
            {
                this->operator()(i, j, k) = 1 / 2.0 * (this->operator()(i - 1, j, k) + this->operator()(i - 1, j, k));
            }
        }

        for (int i = 2; i < dims.XDim - 1; i += 2)
        {
            for (int j = 1; j < dims.YDim - 1; ++j)
            {
                this->operator()(i, j, k) = 1 / 2.0 * (this->operator()(i, j - 1, k) + this->operator()(i, j + 1, k));
            }
        }

        for (int i = 1; i < dims.XDim - 1; i += 2)
        {
            for (int j = 1; j < dims.YDim - 1; j += 2)
            {
                this->operator()(i, j, k) = 1 / 4.0 * (this->operator()(i - 1, j - 1, k) + this->operator()(i - 1, j + 1, k) +
                this->operator()(i + 1, j - 1, k) + this->operator()(i + 1, j + 1, k)) +
                1 / 2.0 * (this->operator()(i, j - 1, k) + this->operator()(i, j + 1, k) +
                this->operator()(i - 1, j, k) + this->operator()(i + 1, j, k));
            }
        }
        
    }

    for (int k = 1; k < dims.ZDim - 1; k += 2)
    {
        for (int i = 1; i < dims.XDim - 1; ++i)
        {
            for (int j = 1; j < dims.YDim - 1; ++j)
            {
                this->operator()(i, j, k) = 1 / 2.0 * (this->operator()(i, j, k + 1) + this->operator()(i, j, k - 1));
            }
        }
    }
    */
    Dimensions3D dims = points_.GetDimensions();
    for (int k = 1; k < dims.ZDim - 1 ; ++k)
    {
        for (int j = 1; j < dims.YDim - 1; ++j)
        {
            for (int i = 1; i < dims.XDim - 1; ++i)
            {
                this->operator()(i, j, k) = TrilinearInterpolation(old, i / 2., j / 2., k / 2.,
                                                                    i / 2, j / 2, k / 2,
                                                                   i / 2 + 1, j / 2 + 1, k / 2 + 1);
            }
        }
    }
}

double Lattice3D::TrilinearInterpolation(Lattice3D& old, double i, double j, double k, int i1, int j1, int k1, int i2, int j2, int k2)
{
    return (
    old(i1, j1, k1) * (i2 - i) * (j2 - j) * (k2 - k) + 
    old(i1, j1, k2) * (i2 - i) * (j2 - j) * (k - k1) + 
    old(i1, j2, k1) * (i2 - i) * (j - j1) * (k2 - k) +
    old(i1, j2, k2) * (i2 - i) * (j - j1) * (k - k1) +
    old(i2, j1, k1) * (i - i1) * (j2 - j) * (k2 - k) +
    old(i2, j1, k2) * (i - i1) * (j2 - j) * (k - k1) +
    old(i2, j2, k1) * (i - i1) * (j - j1) * (k2 - k) +
    old(i2, j2, k2) * (i - i1) * (j - j1) * (k - k1)
    );
}

double& Lattice3D::operator()(int i, int j, int k)
{
    return points_(i, j, k);
}

Point3D Lattice3D::GetCoordinates(int i, int j, int k) const
{
    double x = domain_.x_dom_.le + i * h_;
    double y = domain_.y_dom_.le + j * h_;
    double z = domain_.z_dom_.le + k * h_;
    return {x, y, z};
}

Point3D Lattice3D::GetDeltas() const
{
    return {h_, h_, h_};
}

double Lattice3D::Get_h() const
{
    return h_;
}

Dimensions3D Lattice3D::GetDimensions() const
{
    return points_.GetDimensions();
}

void Lattice3D::SetBoundaryCondition()
{
    Dimensions3D dims = points_.GetDimensions();
    for (int i = 0; i < dims.XDim; ++i)
    {
        for (int j = 0; j < dims.YDim; ++j)
        {
            for (int k = 0; k < dims.ZDim ; ++k)
            {
                if (i == 0 || i == dims.XDim - 1 || j == 0 || j == dims.YDim - 1 || k == 0 || k == dims.ZDim - 1)
                {
                    this->operator()(i, j, k) = Boundary_function(GetCoordinates(i, j, k).x, 
                    GetCoordinates(i, j, k).y, GetCoordinates(i, j, k).z);
                }
            }
        }
    }

    
}

Lattice3D::~Lattice3D()
{
    
}

void Lattice3D::Fill(double value)
{
    Dimensions3D dims = points_.GetDimensions();
    for (int i = 0; i < dims.XDim; ++i)
    {
        for (int j = 0; j < dims.YDim; ++j)
        {
            for (int k = 0; k < dims.ZDim ; ++k)
            {
                this->operator()(i, j, k) = value;
            }
        }
    }
}

void Lattice3D::Fill_Interior(double value)
{
    Dimensions3D dims = points_.GetDimensions();
    for (int i = 1; i < dims.XDim - 1; ++i)
    {
        for (int j = 1; j < dims.YDim - 1; ++j)
        {
            for (int k = 1; k < dims.ZDim - 1; ++k)
            {
                this->operator()(i, j, k) = value;
            }
        }
    }
}

void Lattice3D::swap(Lattice3D& other)
{
    points_.swap(other.points_);
    std::swap(domain_, other.domain_);
    std::swap(h_, other.h_);
}

void Lattice3D::print()
{
    Dimensions3D dims = points_.GetDimensions();
    for (int k = 0; k < dims.ZDim; ++k)
    {
        for (int j = 0; j < dims.YDim; ++j)
        {
            for (int i = 0; i < dims.XDim ; ++i)
            {
                std::cout << this->operator()(i, j, k) << ' ';
            }
            std::cout << std::endl;
        }
        std::cout << std::endl << std::endl;
    }
}