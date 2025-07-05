#include "../include/ScalarField.hpp"
#include "../include/Mesh.hpp"
#include "../include/DirichletBC.hpp"
#include "../include/ADISolver.hpp"
#include "../include/DumbSolver.hpp"
#include <iostream>
#include <math.h>

double sourceEquation(double x, double y)
{
    return -2 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y); // Poisson equation source term for Laplace equation
}

double solutionEquation(double x, double y)
{
    return sin(M_PI * x) * sin(M_PI * y); // Solution to the Laplace equation with the given source term
}

int main()
{
    // Generate a simple orthogonal mesh
    int nx = 8;
    int ny = 8;
    double lengthX = 1;                                            // M_PI;                                         // nx - 1; // 1D because ny = 1
    double lengthY = 1;                                            // M_PI;                                         // ny - 1;
    Mesh mesh = Mesh::getOrthogonalMesh(nx, ny, lengthX, lengthY); // 1D because ny = 1

    // Function space (assumes central differences)
    ScalarField phi(mesh);
    ScalarField f(mesh);
    f.applySourceFunction(sourceEquation);

    // BCs
    MeshRegion bottom("bottom", {0, 0}, {nx - 1, 0});
    DirichletBC bc1(mesh, bottom, 0);

    ADISolver solver(phi, f, bc1);
    solver.solve(100);

    std::cout << "Calculated: \n"
              << solver.u() << std::endl;

    ScalarField solution(mesh);
    solution.applySourceFunction(solutionEquation);

    std::cout << "Expected Solution:\n"
              << solution << std::endl;
    std::cout << "Difference:\n"
              << solver.u() - solution << std::endl;

    // Dumb Solver
    DumbSolver dumbSolver(phi, f, bc1);
    dumbSolver.solve();
    std::cout << "Dumb Solver Result:\n"
              << dumbSolver.u() << std::endl;
}