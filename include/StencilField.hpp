#pragma once
#include "ScalarField.hpp"

class StencilField
{
public:
    ScalarField center; // i, j
    ScalarField east;   // i + 1, j
    ScalarField west;   // i - 1, j
    ScalarField north;  // i, j + 1
    ScalarField south;  // i, j - 1

    StencilField(Mesh &mesh)
        : center(mesh, 0), east(mesh, 0), west(mesh, 0), north(mesh, 0), south(mesh, 0)
    {
        initialize();
    }

private:
    void initialize()
    {
        // mixes inner point and boundary conditions logic
        // it assumes dirichlet bc
        const Mesh &mesh = center.mesh();
        const int nx = mesh.nx();
        const int ny = mesh.ny();

        double dx0, dx1, dy0, dy1;
        for (int i = 1; i < nx - 1; ++i)
        {
            for (int j = 1; j < ny - 1; ++j)
            {
                // Yikes, everywhere we use ghost cells to avoid ifs except for here...
                dx0 = mesh.x(i, j) - mesh.x(i - 1, j);
                dx1 = mesh.x(i + 1, j) - mesh.x(i, j);
                dy0 = mesh.y(i, j) - mesh.y(i, j - 1);
                dy1 = mesh.y(i, j + 1) - mesh.y(i, j);

                east(i, j) = 2 / (dx1 * (dx0 + dx1));
                west(i, j) = 2 / (dx0 * (dx0 + dx1));
                north(i, j) = 2 / (dy1 * (dy0 + dy1));
                south(i, j) = 2 / (dy0 * (dy0 + dy1));
                center(i, j) = -(east(i, j) + west(i, j) +
                                 north(i, j) + south(i, j));
                // _stencilField.center(i, j) = -4 / ((dx0 * dx1) * (dx0 + dx1)) - 4 / ((dy0 * dy1) * (dy0 + dy1));
            }   
        }
    }
};