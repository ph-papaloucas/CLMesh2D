#pragma once

#include "BoundaryCondition.hpp"
#include "Mesh.hpp"

class DirichletBC : public BoundaryCondition
{
public:
    DirichletBC(Mesh &mesh, MeshRegion &region, ScalarFunction func)
        : BoundaryCondition(mesh, region, func)
    {
        checkIfIsBoundaryRegion(_region);
    }

    DirichletBC(Mesh &mesh, MeshRegion &region, double value)
        : BoundaryCondition(mesh, region,
                            [value](double, double)
                            { return value; })
    {
        checkIfIsBoundaryRegion(_region);
    }

    void applyBCsToStencilField(StencilField &stencilField) const
    {
        for (auto &node : _region.nodes())
        {
            int i = node[0];
            int j = node[1];
            stencilField.center(i, j) = 1.0;
            // stencilField.east(i, j) = 0.0;
            // stencilField.west(i, j) = 0.0;
            // stencilField.north(i, j) = 0.0;
            // stencilField.south(i, j) = 0.0;
        }
    }
};