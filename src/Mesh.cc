#include "../include/Mesh.hpp"

Mesh Mesh::getOrthogonalMesh(int nx, int ny, double lengthX, double lengthY, double x0, double y0)
{
    Mesh mesh;
    mesh._nx = nx;
    mesh._ny = ny;

    double dx = lengthX / (nx - 1);
    double dy = lengthY / (ny - 1);

    mesh._coords.reserve(nx * ny);

    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            double x = x0 + i * dx;
            double y = y0 + j * dy;
            mesh._coords.push_back({x, y});
        }
    }

    return mesh;
}

int Mesh::nx() const
{
    return _nx;
}

int Mesh::ny() const
{
    return _ny;
}

double Mesh::x(int i, int j) const
{
    return _coords[i + j * _nx][0];
}

double Mesh::y(int i, int j) const
{
    return _coords[i + j * _nx][1];
}

// void Mesh::addBoundaryRegion(MeshRegion boundary)
// {
//     // check if MeshRegion is on a boundary
//     int i0 = boundary.corner0()[0];
//     int j0 = boundary.corner0()[1];
//     int i1 = boundary.corner1()[0];
//     int j1 = boundary.corner1()[1];

//     bool same_i = (i0 == i1);
//     bool same_j = (j0 == j1);

//     if (same_i && !(i0 == 0 || i0 == _nx - 1))
//     {
//         throw std::runtime_error("BoundaryRegion '" + boundary.name() + "' is aligned vertically but not on left/right boundary.");
//     }
//     if (same_j && !(j0 == 0 || j0 == _ny - 1))
//     {
//         throw std::runtime_error("BoundaryRegion '" + boundary.name() + "' is aligned horizontally but not on top/bottom boundary.");
//     }

//     _boundaryRegions.push_back(boundary);
// }

#include <algorithm>
MeshRegion::MeshRegion(std::string name, Mesh &mesh, std::array<int, 2> node0, std::array<int, 2> node1)
    : _name(std::move(name)), _mesh(mesh)
{
    // Extract i, j indices
    int i0 = node0[0];
    int j0 = node0[1];
    int i1 = node1[0];
    int j1 = node1[1];

    // Ensure bottom-left (min i, min j) and top-right (max i, max j)
    int imin = std::min(i0, i1);
    int imax = std::max(i0, i1);
    int jmin = std::min(j0, j1);
    int jmax = std::max(j0, j1);

    _nodes.reserve((imax - imin) * (jmax - jmin));
    // Fill nodes in ascending order (j outer, i inner)
    for (int j = jmin; j <= jmax; ++j)
    {
        for (int i = imin; i <= imax; ++i)
        {
            _nodes.emplace_back(std::array<int, 2>{i, j});
        }
    }
}


void MeshRegion::sortNodes()
{
    // Sort nodes in ascending order (j outer, i inner)
    std::sort(_nodes.begin(), _nodes.end(), [](const std::array<int, 2> &a, const std::array<int, 2> &b)
              {
        if (a[1] == b[1])
            return a[0] < b[0];
        return a[1] < b[1]; });
}
