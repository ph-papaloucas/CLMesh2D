#pragma once
#include <iostream>
#include <vector>
#include <ScalarField.hpp>

class VectorField {
public:
    ScalarField& u();  // x-component
    ScalarField& v();  // y-component
};