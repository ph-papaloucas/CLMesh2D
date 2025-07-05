#pragma once

#include "BoundaryCondition.hpp"
#include "Mesh.hpp"

class DirichletBC : public BoundaryCondition
{
public:
    DirichletBC(Mesh &mesh, MeshRegion &region, double value)
        : _mesh(mesh),
          _region(region),
          _value(value)
    {
        // check if MeshRegion is on a boundary
        int i0 = region.corner0()[0];
        int j0 = region.corner0()[1];
        int i1 = region.corner1()[0];
        int j1 = region.corner1()[1];

        bool same_i = (i0 == i1);
        bool same_j = (j0 == j1);

        if (same_i && !(i0 == 0 || i0 == _mesh.nx() - 1))
        {
            throw std::runtime_error("BoundaryRegion '" + region.name() + "' is aligned vertically but not on left/right boundary.");
        }
        if (same_j && !(j0 == 0 || j0 == _mesh.ny() - 1))
        {
            throw std::runtime_error("BoundaryRegion '" + region.name() + "' is aligned horizontally but not on top/bottom boundary.");
        }
    }

    void apply(ScalarField &field) const override
    {
        // filed = rhs
        // Warning: corners are getting overwritten every time !
        int nx = _mesh.nx();
        int ny = _mesh.ny();

        for (const auto &node : _region.nodes())
        {
            int i = node[0];
            int j = node[1];

            // not really elegant..
            // if we are sure region is just a line we can do if only once outside for loop
            // We do not elseifs in order to handle corner cells

            if (i == 0)
                field(i - 1, j) = _value; // left ghost
            if (i == nx - 1)
                field(i + 1, j) = _value; // right ghost
            if (j == 0)
                field(i, j - 1) = _value; // bottom ghost
            if (j == ny - 1)
                field(i, j + 1) = _value; // top ghosts
        }
    }

    double value() const
    {
        return _value;
    }

private:
    const MeshRegion &_region;
    const Mesh &_mesh;
    double _value;
};