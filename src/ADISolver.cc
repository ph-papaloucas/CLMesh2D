#include "../include/ADISolver.hpp"

void ADISolver::solve(int numSweeps)
{
    for (int i = 0; i < numSweeps; ++i)
    {
        sweepRows();
        sweepColumns();
    }

    // ADI logic.. settings must be applied for residual etc. for now we will do fixed runs
}

void ADISolver::sweepRows()
{
    // Tridiagonal System to setup and solve:
    // (Cpx + Cpy)*u_{i,j} + Cpw*u_{i-1,j} + Cpe*u_{i+1,j} = f_i,j - Cpn*u_{i,j+1} - Cps*u_{i,j-1}
    std::vector<double> b(_nx);
    std::vector<double> a(_nx - 1), c(_nx - 1);
    std::vector<double> rhs(_nx);

    for (int j = 1; j < _ny - 1; ++j)
    {
        for (int i = 1; i < _nx - 1; ++i)
        {
            a[i - 1] = _stencilField.west(i, j);
            b[i] = _stencilField.center(i, j);
            c[i] = _stencilField.east(i, j);
            rhs[i] = _f(i, j) - _stencilField.north(i, j) * _u(i, j + 1) - _stencilField.south(i, j) * _u(i, j - 1);
        }

        applyBoundaryConditionsRow(a, b, c, rhs, j);
        std::vector<double> urow = tridiagonalSolver(a, b, c, rhs);
        updateRow(urow, j);
    }
}

void ADISolver::applyBoundaryConditionsRow(std::vector<double> &a, std::vector<double> &b, std::vector<double> &c, std::vector<double> &rhs, int j)
{
    // hardcoded dirhclet
    //  boundary 1
    b[0] = 1;
    c[0] = 0;
    rhs[0] = _bcs.value();

    // boundary 2
    a[_nx - 2] = 0;
    b[_nx - 1] = 1;
    rhs[_nx - 1] = _bcs.value();
}

void ADISolver::applyBoundaryConditionsColumn(std::vector<double> &a, std::vector<double> &b, std::vector<double> &c, std::vector<double> &rhs, int i)
{
    // hardcoded dirhclet
    // boundary 1
    b[0] = 1;
    c[0] = 0;
    rhs[0] = _bcs.value();

    // boundary 2
    a[_ny - 2] = 0;
    b[_ny - 1] = 1;
    rhs[_ny - 1] = _bcs.value();
}

void ADISolver::sweepColumns()
{
    // Tridiagonal System to setup and solve:
    // (Cpx + Cpy)*u_{i,j} + Cps*u_{i,j-1} + Cpn*u_{i,j+1} = f_i,j - Cpw*u_{i-1,j} - Cpe*u_{i+1,j}
    std::vector<double> b(_ny);
    std::vector<double> a(_ny - 1), c(_ny - 1);
    std::vector<double> rhs(_ny);

    for (int i = 1; i < _nx - 1; ++i)
    {
        for (int j = 1; j < _ny - 1; ++j)
        {
            a[j - 1] = _stencilField.south(i, j);
            b[j] = _stencilField.center(i, j);
            c[j] = _stencilField.north(i, j);
            rhs[j] = _f(i, j) - _stencilField.south(i, j) * _u(i + 1, j) - _stencilField.west(i, j) * _u(i - 1, j);
        }

        applyBoundaryConditionsColumn(a, b, c, rhs, i);
        std::vector<double> urow = tridiagonalSolver(a, b, c, rhs);
        updateColumn(urow, i);
    }
}

ScalarField &ADISolver::u()
{
    return _u;
}
void ADISolver::updateRow(std::vector<double> urow, int j)
{
    for (int i = 0; i < _nx; ++i)
        _u(i, j) = urow[i];
}

void ADISolver::updateColumn(std::vector<double> ucol, int i)
{
    for (int j = 0; j < _ny; ++j)
    {
        _u(i, j) = ucol[j];
    }
}

std::vector<double> tridiagonalSolver(const std::vector<double> &a,
                                      const std::vector<double> &b,
                                      const std::vector<double> &c,
                                      const std::vector<double> &rhs)
{
    int n = b.size();
    std::vector<double> x(n, 0);
    std::vector<double> c_prime(n - 1, 0);
    std::vector<double> d_prime(n, 0);

    // Forward elimination
    c_prime[0] = c[0] / b[0];
    d_prime[0] = rhs[0] / b[0];

    for (int i = 1; i < n - 1; ++i)
    {
        double m = b[i] - a[i - 1] * c_prime[i - 1];
        c_prime[i] = c[i] / m;
        d_prime[i] = (rhs[i] - a[i - 1] * d_prime[i - 1]) / m;
    }

    d_prime[n - 1] = (rhs[n - 1] - a[n - 2] * d_prime[n - 2]) / (b[n - 1] - a[n - 2] * c_prime[n - 2]);

    // Back substitution
    x[n - 1] = d_prime[n - 1];
    for (int i = n - 2; i >= 0; --i)
    {
        x[i] = d_prime[i] - c_prime[i] * x[i + 1];
    }

    return x;

} // end of tridiagonal_solver

const StencilField &ADISolver::stencilField() const
{
    return _stencilField;
}