#pragma once
#include <vector>
#include <array>
#include <string>
#include <iostream>

class Mesh
{
public:
    Mesh() = default;

    static Mesh getOrthogonalMesh(int nx, int ny, double lengthX, double lengthY,
                                  double x0 = 0, double y0 = 0); // just an easy way to generate a simple mesh

    int nx() const;
    int ny() const;

    double x(int i, int j) const;
    double y(int i, int j) const;

    std::vector<std::array<double, 2>> &coords() { return _coords; }
    std::vector<std::vector<double>> dx() const
    {
        std::vector<std::vector<double>> dx(_nx, std::vector<double>(_ny, 0.0));
        for (int i = 0; i < _nx; ++i)
        {
            for (int j = 0; j < _ny; ++j)
            {
                dx[i][j] = x(i + 1, j) - x(i, j);
            }
        }
        return dx;
    }

    std::vector<std::vector<double>> dy() const
    {
        std::vector<std::vector<double>> dy(_nx, std::vector<double>(_ny, 0.0));
        for (int i = 0; i < _nx; ++i)
        {
            for (int j = 0; j < _ny; ++j)
            {
                dy[i][j] = y(i, j + 1) - y(i, j);
            }
        }
        return dy;
    }

private:
    int _nx;
    int _ny;
    std::vector<std::array<double, 2>> _coords; // flatten vector (_nx x _ny ), each element is an array<double, 2> (x, y)
};

class MeshRegion
{
public:
    MeshRegion() = default;
    MeshRegion(std::string name, Mesh &mesh, std::vector<std::array<int, 2>> nodes)
        : _name(std::move(name)), _mesh(mesh), _nodes(std::move(nodes)) {}
    MeshRegion(std::string name, const Mesh &mesh, std::vector<std::array<int, 2>> nodes)
        : _name(std::move(name)), _mesh(mesh), _nodes(std::move(nodes)) {}
    MeshRegion(std::string name, Mesh &mesh, std::array<int, 2> node0, std::array<int, 2> node1);

    const std::string &name() const { return _name; }
    std::vector<std::array<int, 2>> nodes() const { return _nodes; }

    const Mesh &mesh() const { return _mesh; }


private:
    void sortNodes();
    const std::string _name;
    std::vector<std::array<int, 2>> _nodes;
    const Mesh &_mesh;
};
