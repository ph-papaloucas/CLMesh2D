#pragma once

#include "BoundaryCondition.hpp"
#include "Mesh.hpp"

class NeumannBC : public BoundaryCondition
{
public:
    NeumannBC(Mesh &mesh, MeshRegion &region, double value)
        : BoundaryCondition(mesh),
          _region(region),
          _value(value)
    {
        chehckIfIsBoundaryRegion(_region);
    }

    void applyBCsToStencilField(ScalarField &u) const override
    {
        std::cout << "test" << std::endl;
    }

    double value() const
    {
        return _value;
    }

private:
    const MeshRegion &_region;
    double _value;
};