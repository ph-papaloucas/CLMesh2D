#pragma once

#include "BoundaryCondition.hpp"
#include "Mesh.hpp"

class NeumannBC : public BoundaryCondition
{
public:
    NeumannBC(Mesh &mesh, MeshRegion &region, ScalarFunction func)
        : BoundaryCondition(mesh, region, func)
    {
        checkIfIsBoundaryRegion(_region);
    }

    NeumannBC(Mesh &mesh, MeshRegion &region, double value)
        : BoundaryCondition(mesh, region, [value](double, double)
                            { return value; })
    {
        checkIfIsBoundaryRegion(_region);
    }

    void applyBCsToStencilField(StencilField &f) const override
    {
        std::cout << "test" << std::endl;
    }
};