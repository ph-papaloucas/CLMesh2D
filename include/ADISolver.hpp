#pragma once
#include "ScalarField.hpp"
#include "Mesh.hpp"
#include "BoundaryCollector.hpp"
#include "StencilField.hpp"

std::vector<double> tridiagonalSolver(const std::vector<double> &a,
                                      const std::vector<double> &b,
                                      const std::vector<double> &c,
                                      const std::vector<double> &rhs);

class ADISolver
// for now this class is solving only the poison Laplace(u) = f
// It mixes discretization logic with the solver logic
// It should be refactored later to separate the discretization and the solver logic
{
public:
    ADISolver(ScalarField &u,
              ScalarField &f,
              BoundaryCollector &bcs)
        : _u(u), _f(f), _bcs(bcs), _nx(_u.nx()), _ny(_u.ny()), _nghostLayers(_u.nghostLayers()), _stencilField(_u.mesh()) {}

    void solve(int numSweeps = 1000);

    const StencilField &stencilField() const;

    void sweepRows();

    void applyBoundaryConditionsRow(std::vector<double> &a, std::vector<double> &b, std::vector<double> &c, std::vector<double> &rhs, int j);

    void applyBoundaryConditionsColumn(std::vector<double> &a, std::vector<double> &b, std::vector<double> &c, std::vector<double> &rhs, int i);

    void sweepColumns();

    ScalarField &u();

    void updateRow(std::vector<double> urow, int j);

    void updateColumn(std::vector<double> ucol, int i);

private:
    ScalarField &_u;                  // not const, because we will apply boundary conditions on ghost cells
    const ScalarField &_f;            // not const, because we will apply boundary conditions on ghost cells
    const StencilField _stencilField; // holds the stencil for each cell
    const BoundaryCollector &_bcs;
    const int _nx;           // number of grid points in one direction (including ghost cells)
    const int _ny;           // number of grid points in the other direction (including ghost cells)
    const int _nghostLayers; // number of ghost layers
};