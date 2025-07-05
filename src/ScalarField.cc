#include "../include/ScalarField.hpp"

ScalarField::ScalarField(Mesh &mesh, int nghostLayers)
    : _mesh(mesh),
      _nghostLayers(nghostLayers),
      _nx(_mesh.nx()),
      _ny(_mesh.ny()),
      _data((_nx + 2 * _nghostLayers) * (_ny + 2 * _nghostLayers), 0.0) {}

// Getters
int ScalarField::nx() const
{
    return this->_nx;
}
int ScalarField::ny() const
{
    return this->_ny;
}
int ScalarField::nghostLayers() const
{
    return this->_nghostLayers;
}

// Operator overloading
const double &ScalarField::operator()(int i, int j) const
{ // read - only
    return _data[(j + _nghostLayers) * (_nx + 2 * _nghostLayers) + (i + _nghostLayers)];
    // without ghost cells: _data[j*_nx + i]
    // with ghost cells: j -> j + _nghostLayers
    //                   i -> i + _nghostLayers
    //                   _nx -> _nx + 2*_nghostLayers
}

double &ScalarField::operator()(int i, int j)
{
    assert(i >= -1 && i <= _nx); // assertions can be disabled in release mode for performance !
    assert(j >= -1 && j <= _ny);
    return _data[(j + _nghostLayers) * (_nx + 2 * _nghostLayers) + (i + _nghostLayers)];
}

ScalarField &ScalarField::operator=(const ScalarField &other)
{
    if (this != &other) // Avoid self-assignment
    {
        _nx = other._nx;
        _ny = other._ny;
        _nghostLayers = other._nghostLayers;
        _data = other._data; // std::vector has its own deep copy
    }
    return *this;
}

ScalarField ScalarField::operator-(const ScalarField &other)
{
    assert(_nx == other._nx && _ny == other._ny && _nghostLayers == other._nghostLayers);
    ScalarField result(_mesh, _nghostLayers);
    for (int r = 0; r < _data.size(); ++r)
    {
        result._data[r] = _data[r] - other._data[r];
    }
    return result;
}

void ScalarField::applySourceFunction(ScalarFunction func)
{
    for (int i = 0; i < _nx; ++i)
    {
        for (int j = 0; j < _ny; ++j)
        {
            double x = _mesh.x(i, j);
            double y = _mesh.y(i, j);
            _data[(j + _nghostLayers) * (_nx + 2 * _nghostLayers) + (i + _nghostLayers)] = func(x, y);
        }
    }
}

#include <iomanip> // for std::fixed and std::setprecision
#include <cmath>   // for std::log10, std::floor
std::ostream &operator<<(std::ostream &os, const ScalarField &u)
{
    int nx = u.nx();
    int ny = u.ny();
    double max = *std::max_element(u.data().begin(), u.data().end());

    int precision = 6;
    os << std::fixed << std::setprecision(precision); // Set fixed format with 1 decimal
    int intPartDigits = (max > 0) ? static_cast<int>(std::floor(std::log10(max))) + 1 : precision;
    int maxChars = intPartDigits + 1 /*decimal point*/ + 1 /*digit after decimal*/;

    for (int j = ny - 1 + u.nghostLayers(); j >= -u.nghostLayers(); --j)
    {
        for (int i = -u.nghostLayers(); i < nx + u.nghostLayers(); ++i)
        {
            if (i < 0 || j < 0 || i > nx - 1 || j > ny - 1)
                os << "[" << std::setw(maxChars) << std::setfill('0') << u(i, j) << "]";
            else
                os << " " << std::setw(maxChars) << std::setfill('0') << u(i, j) << " ";
        }
        os << "\n";
    }
    return os;
}

