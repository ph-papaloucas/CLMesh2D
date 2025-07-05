#pragma once
#include <vector>
#include <memory>
#include <BoundaryCondition.hpp>
#include <ScalarField.hpp>

class BoundaryCollector
{
public:
    BoundaryCollector(std::vector<std::unique_ptr<BoundaryCondition>> &&bcs)
        : _bcs(std::move(bcs))
    {
    }

    void applyBCsToStencilField(ScalarField &u) const
    {
        for (const auto &bc : _bcs)
            bc->applyBCsToStencilField(u);
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
        return 0;
    }

private:
    std::vector<std::unique_ptr<BoundaryCondition>> _bcs;
};