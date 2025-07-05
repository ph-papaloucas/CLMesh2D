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
   
    // prepare boundaries
    MeshRegion bottom("bottom", {0, 0}, {nx - 1, 0});
    MeshRegion top("top", {0, ny - 1}, {nx - 1, ny - 1});
    MeshRegion left("left", {0, 0}, {0, ny - 1});
    MeshRegion right("right", {nx - 1, 0}, {nx - 1, ny - 1});
    StencilField stencilField(mesh);

    double dirichletValue = 0.0;
    DirichletBC bc1(mesh, bottom, dirichletValue);
    DirichletBC bc2(mesh, top, dirichletValue);
    DirichletBC bc3(mesh, left, dirichletValue);
    DirichletBC bc4(mesh, right, dirichletValue);


    NeumannBC neumannBC1(mesh, left, 1.0);
    DirichletBC dirichletBC1(mesh, right, 0.0);
    
    int i =0, j = 0;
    BoundaryCollector bcs = BoundaryCollector::makeCollector(neumannBC1, dirichletBC1);
    bcs.applyBCsToStencilField(stencilField.center);
    double test = bcs.value(i, j);

    //BoundaryCollector bcs = BoundaryCollector::makeCollector(std::move(bc1), std::move(bc2), std::move(bc3), std::move(bc4));

    BoundaryCollector bcs = BoundaryCollector::makeCollector(bc1, bc2, bc3, bc4);
    bcs.applyBCsToStencilField(stencilField.center);

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
}