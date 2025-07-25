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
    // prepare boundaries (this logic could be a static method in MeshRegion)
    std::vector<std::array<int, 2>> boundaryNodes(2 * (nx + ny));
    for (int i = 0; i < nx; ++i)
    {
        boundaryNodes[i] = {i, 0};           // bottom
        boundaryNodes[i + nx] = {i, ny - 1}; // top
    }
    for (int j = 0; j < ny; ++j)
    {
        boundaryNodes[j + 2 * nx] = {0, j};           // left
        boundaryNodes[j + 2 * nx + ny] = {nx - 1, j}; // right
    }
    MeshRegion boundaryRegion("boundary", mesh, boundaryNodes);
    double dirichletValue = 0.0;                           // Dirichlet boundary condition value
    DirichletBC bc1(mesh, boundaryRegion, dirichletValue); // Dirichlet BC with the solution function
    BoundaryCollector bcs = BoundaryCollector::makeCollector(bc1);

    ADISolver solver(phi, f, bcs);
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
    DumbSolver dumbSolver(phi, f, bcs);
    dumbSolver.solve();
    std::cout << "Dumb Solver Result:\n"
              << dumbSolver.u() << std::endl;
}