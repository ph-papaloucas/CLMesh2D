#include "../include/ScalarField.hpp"
#include "../include/Mesh.hpp"
#include "../include/DirichletBC.hpp"
#include "../include/ADISolver.hpp"
#include "../include/BoundaryCollector.hpp"
#include "../include/NeumannBC.hpp"
#include <iostream>
#include <math.h>

int main()
{

    // CHECK STENCILS
    // CHECK IF RHS IS CORRECT
    // CHECK EQUATIONS AGAIN

    // Generate a simple orthogonal mesh
    int nx = 6;
    int ny = 6;
    double lengthX = nx - 1;                                       // M_PI;                                         // nx - 1; // 1D because ny = 1
    double lengthY = ny - 1;                                       // M_PI;                                         // ny - 1;
    Mesh mesh = Mesh::getOrthogonalMesh(nx, ny, lengthX, lengthY); // 1D because ny = 1

    // Prepare MeshRegions (should find an easier way to do this)
    MeshRegion topBoundaryRegion("top", mesh, {0, ny - 1}, {nx - 1, ny - 1});
    std::vector<std::array<int, 2>> otherBoundaryNodes(2 * ny + nx);
    for (int j = 0; j < ny; ++j)
    {
        otherBoundaryNodes[j] = {0, j};           // left
        otherBoundaryNodes[j + ny] = {nx - 1, j}; // right
    }
    for (int i = 0; i < nx; ++i)
    {
        otherBoundaryNodes[i + 2 * ny] = {i, 0};
    }
    MeshRegion otherBoundaryRegions("otherBoundary", mesh, otherBoundaryNodes);

    // Prepare Boundary Conditions
    // Collector could be hidden inside solver
    DirichletBC dirichletBC1(mesh, topBoundaryRegion, 0.0);
    NeumannBC neumannBC1(mesh, otherBoundaryRegions, 1.0);
    BoundaryCollector bcs = BoundaryCollector::makeCollector(neumannBC1, dirichletBC1);

    // Prepare StencilField and rhs
    StencilField stencilField(mesh);
    ScalarField f(mesh, 0); // just a dummy source term
    f.applySourceFunction([](double x, double y) { return 3; }); // just a dummy source term

    bcs.applyBCsToStencilField(stencilField); // this logic is handled internaly by the solver
    bcs.applyBCsToSourceField(f);             // because e.g. ADI is implemented in a way that it needs to reconstruct the system a lot of times

    std::cout << "StencilField.center:\n"
              << stencilField.center << std::endl;
    std::cout << "StencilField.east:\n"
              << stencilField.east << std::endl;
    std::cout << "StencilField.west:\n"
              << stencilField.west << std::endl;
    std::cout << "StencilField.north:\n"
              << stencilField.north << std::endl;
    std::cout << "StencilField.south:\n"
              << stencilField.south << std::endl;
    std::cout << "f:\n"
              << f << std::endl;
}