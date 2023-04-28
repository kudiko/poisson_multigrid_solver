#pragma once
#include <utility>
#include <vector>

struct Domain1D
{
    double le; //left edge
    double re; //right edge
};

struct Point3D
{
    double x;
    double y;
    double z;
};

struct Domain3D
{
    Domain3D(Point3D center, double x_len, double y_len, double z_len) : center_(center)
    {
        x_dom_.le = center_.x - x_len;
        x_dom_.re = center_.x + x_len;
        y_dom_.le = center_.y - y_len;
        y_dom_.re = center_.y + y_len;
        z_dom_.le = center_.z - z_len;
        z_dom_.re = center_.z + z_len;
    }
    Point3D center_;
    Domain1D x_dom_;
    Domain1D y_dom_;
    Domain1D z_dom_;
};

struct Dimensions3D
{
    std::size_t XDim;
    std::size_t YDim;
    std::size_t ZDim;
};

class Array3D
{
    public:
    Array3D(Dimensions3D dimensions);
    double& operator()(int i, int j, int k);
    ~Array3D();
    Dimensions3D GetDimensions() const;
    void swap(Array3D& other);

    private:
    double* points_ = nullptr;
    Dimensions3D dimensions_;

};

class Lattice3D
{
    public:
    Lattice3D(Domain3D domain, Dimensions3D dimensions);

    void Finify();
    void Coarsify();
    double& operator()(int i, int j, int k);
    Point3D GetCoordinates(int i, int j, int k) const;
    Point3D GetDeltas() const;
    double Get_h() const;
    void SetBoundaryCondition();
    Dimensions3D GetDimensions() const;
    void Fill(double value);
    void Fill_Interior(double value);
    ~Lattice3D();
    
    private:
    friend class Tests;

    Domain3D domain_;

    Array3D points_;

    double h_;

    

    void swap(Lattice3D& other);
    void InterpolationProcedure(Lattice3D& old_lat);
    double TrilinearInterpolation(Lattice3D& old, double i, double j, double k, int i1, int j1, int k1, int i2, int j2, int k2);
    void print();

};