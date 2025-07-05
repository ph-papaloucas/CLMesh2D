#pragma once
#include "ScalarField.hpp"
#include "StencilField.hpp"
#include "DirichletBC.hpp"
#include "StencilField.hpp"
#include <vector>
#include <cassert>

class DumbSolver
{
public:
    DumbSolver(ScalarField &u, ScalarField &f, DirichletBC &bcs)
        : _u(u), _f(f), _bcs(bcs),
          _nx(u.nx()), _ny(u.ny()), _nghostLayers(u.nghostLayers()),
          _stencilField(u.mesh()) {}

    void solve()
    {
        int N = _nx * _ny;
        std::vector<std::vector<double>> A(N, std::vector<double>(N, 0.0));
        std::vector<double> b(N);
        auto index = [this](int i, int j)
        { return j * _nx + i; };

        // Fill A and b
        for (int i = 0; i < _nx; ++i)
        {
            for (int j = 0; j < _ny; ++j)
            {
                int k = index(i, j);

                if (i == 0 || i == _nx - 1 || j == 0 || j == _ny - 1)
                {
                    A[k][k] = 1.0;
                    b[k] = _bcs.value();
                    continue;
                }

                A[k][index(i, j)] = _stencilField.center(i, j);
                A[k][index(i + 1, j)] = _stencilField.east(i, j);
                A[k][index(i - 1, j)] = _stencilField.west(i, j);
                A[k][index(i, j + 1)] = _stencilField.north(i, j);
                A[k][index(i, j - 1)] = _stencilField.south(i, j);

                b[k] = _f(i, j);
            }
        }

        // Solve Ax = b using naive LU
        std::vector<double> x = solveLU(A, b);

        // Unpack result
        for (int i = 0; i < _nx; ++i)
            for (int j = 0; j < _ny; ++j)
                _u(i, j) = x[index(i, j)];
    }

    ScalarField &u() { return _u; }
    StencilField &stencilField() { return _stencilField; }

private:
    ScalarField &_u;
    ScalarField &_f;
    DirichletBC &_bcs;
    StencilField _stencilField;
    int _nx, _ny, _nghostLayers;

    std::vector<double> solveLU(const std::vector<std::vector<double>> &A, const std::vector<double> &b)
    {
        int n = A.size();
        std::vector<std::vector<double>> L(n, std::vector<double>(n, 0.0));
        std::vector<std::vector<double>> U = A;

        // LU decomposition
        for (int i = 0; i < n; ++i)
        {
            L[i][i] = 1.0;
            for (int j = i + 1; j < n; ++j)
            {
                double factor = U[j][i] / U[i][i];
                L[j][i] = factor;
                for (int k = i; k < n; ++k)
                {
                    U[j][k] -= factor * U[i][k];
                }
            }
        }

        // Forward substitution: Ly = b
        std::vector<double> y(n);
        for (int i = 0; i < n; ++i)
        {
            y[i] = b[i];
            for (int j = 0; j < i; ++j)
                y[i] -= L[i][j] * y[j];
        }

        // Back substitution: Ux = y
        std::vector<double> x(n);
        for (int i = n - 1; i >= 0; --i)
        {
            x[i] = y[i];
            for (int j = i + 1; j < n; ++j)
                x[i] -= U[i][j] * x[j];
            x[i] /= U[i][i];
        }

        return x;
    }
};
