#pragma once
#include "ScalarField.hpp"
#include "Mesh.hpp"

class BoundaryCondition
{
public:
    BoundaryCondition(Mesh &mesh, MeshRegion &region, ScalarFunction func) : _mesh(mesh), _region(region), _func(func) {}

    virtual void applyBCsToStencilField(ScalarField &u) const = 0;

    void chehckIfIsBoundaryRegion(const MeshRegion &region) const
    {
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

    const MeshRegion &region() const
    {
        return _region;
    }

    const std::vector<double> values() const
    {
        std::vector<std::array<int, 2>> nodes = _region.nodes();
        std::vector<double> values(nodes.size());
        for (size_t i = 0; i < nodes.size(); ++i)
        {
            std::array<int, 2> node = nodes[i];
            values[i] = _func(_mesh.x(node[0], node[1]), _mesh.y(node[0], node[1]));
        }
        return values;
    }

protected:
    const Mesh &_mesh;
    const MeshRegion &_region;
    ScalarFunction _func;
};