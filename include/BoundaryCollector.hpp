#pragma once
#include <vector>
#include <memory>
#include <BoundaryCondition.hpp>
#include <ScalarField.hpp>
#include "StencilField.hpp"
#include "Mesh.hpp"
#include <optional>

#include <unordered_map>

struct NodeHash
{
    size_t operator()(const std::array<int, 2> &node) const
    {
        return std::hash<int>()(node[0]) ^ (std::hash<int>()(node[1]) << 1);
    }
};

class BoundaryCollector
{
public:
    BoundaryCollector(std::vector<std::unique_ptr<BoundaryCondition>> &&bcs)
        : _bcs(std::move(bcs))
    {

        for (const auto &bc : _bcs)
        {
            const auto &region = bc->region();
            auto bcNodes = region.nodes();
            auto bcValues = bc->values();

            for (size_t k = 0; k < bcNodes.size(); ++k)
            {
                _nodeValueMap[bcNodes[k]] = bcValues[k];
            }
        }
    }

    void checkNodes(const std::vector<std::array<int, 2>> &nodes) const
    {
        // check if we have dulpicated nodes (if they are corners we do not care, just remove them)
        std::cout << " Not Implemented yet: BoundaryCollector::checkNodes\n";
    }

    void applyBCsToStencilField(StencilField &stencilField) const
    {
        for (const auto &bc : _bcs)
            bc->applyBCsToStencilField(stencilField);
    }

    void applyBCsToSourceField(ScalarField &sourceField) const
    {
        for (const auto &bc : _bcs)
            bc->applyBCsToSourceField(sourceField);
    }

    // The static factory method
    template <typename... BCTypes>
    static BoundaryCollector makeCollector(BCTypes &&...bcs)
    {
        std::vector<std::unique_ptr<BoundaryCondition>> vec;
        // Emplace new unique_ptrs by forwarding the args to the constructors of BCTypes
        (vec.emplace_back(std::make_unique<std::decay_t<BCTypes>>(std::forward<BCTypes>(bcs))), ...);
        return BoundaryCollector(std::move(vec));
    }

    double value(int i, int j) const
    {
        auto it = _nodeValueMap.find({i, j});
        if (it != _nodeValueMap.end())
            return it->second;
        return 0; // or throw, or some default
    }

private:
    std::vector<std::unique_ptr<BoundaryCondition>> _bcs;
    std::vector<std::array<int, 2>> _nodes; // all nodes of the boundary conditions
    std::vector<double> _values;            // this could be ScalarFunction(t) instead of double in the future
    std::unordered_map<std::array<int, 2>, double, NodeHash> _nodeValueMap;
    // std::optional<MeshRegion> _region; // change this.. _region is not optional, but MeshRegion does not have a default constructor
};