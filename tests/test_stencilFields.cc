#include "../include/ScalarField.hpp"
#include "../include/Mesh.hpp"
#include "../include/DirichletBC.hpp"
#include "../include/ADISolver.hpp"
#include <iostream>
#include <math.h>

int main()
{

    // CHECK STENCILS
    // CHECK IF RHS IS CORRECT
    // CHECK EQUATIONS AGAIN

    // Generate a simple orthogonal mesh
    int nx = 3;
    int ny = 3;
    double lengthX = nx - 1;                                       // M_PI;                                         // nx - 1; // 1D because ny = 1
    double lengthY = ny - 1;                                       // M_PI;                                         // ny - 1;
    Mesh mesh = Mesh::getOrthogonalMesh(nx, ny, lengthX, lengthY); // 1D because ny = 1

    StencilField stencilField(mesh);

    // for (int i = 0; i < mesh.coords().size(); ++i)
    // {
    //     std::cout << i << std::endl;
    // }

    ScalarField phi(mesh, 0);
    for (int i = 0; i < phi.data().size(); ++i)
    {
        std::cout << i << std::endl;
    }
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