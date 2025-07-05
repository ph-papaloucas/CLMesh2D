#pragma once
#include "ScalarField.hpp"

class BoundaryCondition
{
public:
    virtual void apply(ScalarField &field) const = 0;

};