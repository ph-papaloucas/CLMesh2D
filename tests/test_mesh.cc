#include "../include/Mesh.hpp"

int main()
{
    // Generate a simple orthogonal mesh
    int nx = 6;
    int ny = 6;
    double lengthX = nx - 1;                                       // M_PI;                                         // nx - 1; // 1D because ny = 1
    double lengthY = ny - 1;                                       // M_PI;                                         // ny - 1;
    Mesh mesh = Mesh::getOrthogonalMesh(nx, ny, lengthX, lengthY); // 1D because ny = 1

    // just a region
    MeshRegion orthogonalRegion("box", mesh, {1, 1}, {nx - 2, ny - 2});

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
    std::cout << "Mesh size: " << mesh.nx() << " x " << mesh.ny() << "\n";
    std::cout << "Boundary Nodes:\n";
    for (const auto &node : boundaryRegion.nodes())
    {
        std::cout << "{" << node[0] << ", " << node[1] << "}\n";
    }

    return 0;
}