#include "../include/ScalarField.hpp"
#include "../include/Mesh.hpp"
#include "../include/DirichletBC.hpp"
#include <iostream>

double sourceEquation(double x, double y)
{
    return 2;
}

int main()
{
    // Generate a simple orthogonal mesh
    int nx = 4;
    int ny = 4;
    double lengthX = 1;
    double lengthY = 1;
    Mesh mesh = Mesh::getOrthogonalMesh(nx, ny, lengthX, lengthY); // 1D because ny = 1

    // Function space (assumes central differences)
    ScalarField phi(mesh);
    ScalarField rhs(mesh);

    // BCs
    MeshRegion bottom("bottom", {0, 0}, {nx - 1, 0});
    MeshRegion top("top", {0, ny - 1}, {nx - 1, ny - 1});
    MeshRegion left("left", {0, 0}, {0, ny - 1});
    MeshRegion right("right", {nx - 1, 0}, {nx - 1, ny - 1});

    std::cout << " rhs before applying source\n";
    std::cout << rhs << std::endl;

    rhs.applySourceFunction(sourceEquation);
    std::cout << "rhs after applying source\n";
    std::cout << rhs << std::endl;

    DirichletBC bc1(mesh, bottom, 11);
    DirichletBC bc2(mesh, top, 22);
    DirichletBC bc3(mesh, left, 33);
    DirichletBC bc4(mesh, right, 44);
    bc1.apply(rhs);
    bc2.apply(rhs);
    bc3.apply(rhs);
    bc4.apply(rhs);
    std::cout << "rhs after applying boundary conditions\n";
    std::cout << rhs << std::endl;

    // // Create scalar field quantities for the equations
    // ScalarField T;
    // ScalarField source(sourceEquation);

    // // Write the Partial Differential Equation
    // PDE eq(sourceEquation);
    // eq.addLaplacian(&T);

    // // // Create a solver object (discretization and initialization of fields will happen here)
    // // //solver is a class because you can have other settings such us timesteps etc
    // FDSolver solver(eq, mesh, bcs);
    // // Actual solve and get a results object?
    // solver.solve();
}