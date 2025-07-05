#pragma once

#include "BoundaryCondition.hpp"
#include "Mesh.hpp"

class DirichletBC : public BoundaryCondition
{
public:
    DirichletBC(Mesh &mesh, MeshRegion &region, ScalarFunction func)
        : BoundaryCondition(mesh, region, func)
    {
        chehckIfIsBoundaryRegion(_region);
    }

    DirichletBC(Mesh &mesh, MeshRegion &region, double value)
        : BoundaryCondition(mesh, region,
                            [value](double, double)
                            { return value; })
    {
        chehckIfIsBoundaryRegion(_region);
    }

    void applyBCsToStencilField(ScalarField &u) const override
    {
        std::cout << "test" << std::endl;
    }
};