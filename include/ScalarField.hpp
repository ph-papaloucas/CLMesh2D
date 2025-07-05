#pragma once
#include <iostream>
#include <vector>
#include <assert.h>
#include "types.hpp"
#include "Mesh.hpp"

// always has ghost cells
class ScalarField
{
public:
    // Constructors
    ScalarField(Mesh &mesh, int nghostLayers = 1);

    // Getters
    int nx() const;
    int ny() const;
    int nghostLayers() const;
    std::vector<double> &data() { return _data; }
    const std::vector<double> &data() const { return _data; }

    Mesh &mesh() { return _mesh; }
    const Mesh &mesh() const { return _mesh; }

    // Operator overloading
    double &operator()(int i, int j);             // allows writing
    const double &operator()(int i, int j) const; // read-only
    ScalarField &operator=(const ScalarField &other);
    ScalarField operator-(const ScalarField &other);
    friend std::ostream &operator<<(std::ostream &os, const ScalarField &u);
    // other
    void applySourceFunction(ScalarFunction func);


private:
    Mesh &_mesh;
    int _nx;
    int _ny;
    int _nghostLayers;         // num of ghost layers
    std::vector<double> _data; // size = (_nx + 2*nghostLayers) x  (_ny + 2*nghostLayers)
};
