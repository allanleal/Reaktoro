// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#include "BilinearInterpolator.hpp"

// C++ includes
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Math/BilinearInterpolator.hpp>

namespace Reaktoro {
namespace {

auto binarySearchHelper(real p, const std::vector<real>& coordinates, unsigned begin, unsigned end) -> unsigned
{
    if(end - begin == 1)
        return begin;

    unsigned mid = (begin + end)/2;

    if(p < coordinates[mid])
        return binarySearchHelper(p, coordinates, begin, mid);
    else
        return binarySearchHelper(p, coordinates, mid, end);
}

auto binarySearch(real p, const std::vector<real>& coordinates) -> unsigned
{
    return binarySearchHelper(p, coordinates, 0, coordinates.size());
}

auto interpolationOutOfBoundsError(real x, real xA, real xB, real y, real yA, real yB) -> void
{
    Exception exception;
    exception.error << "Unable to perform an interpolation at the coordinate pair (" << x << ", " << y << ").";
    exception.reason << "Either the x- or y-coordinate is out of bound, where " << xA << " < x < " << xB << " and " << yA << " < y < " << yB << ".";
    RaiseError(exception);
}

} // namespace

BilinearInterpolator::BilinearInterpolator()
{}

BilinearInterpolator::BilinearInterpolator(
    const std::vector<real>& xcoordinates,
    const std::vector<real>& ycoordinates,
    const std::vector<real>& data)
: m_xcoordinates(xcoordinates),
  m_ycoordinates(ycoordinates),
  m_data(data)
{}

BilinearInterpolator::BilinearInterpolator(
    const std::vector<real>& xcoordinates,
    const std::vector<real>& ycoordinates,
    const std::function<real(real, real)>& function)
: m_xcoordinates(xcoordinates),
  m_ycoordinates(ycoordinates),
  m_data(xcoordinates.size() * ycoordinates.size())
{
    unsigned k = 0;
    for(unsigned j = 0; j < ycoordinates.size(); ++j)
        for(unsigned i = 0; i < xcoordinates.size(); ++i, ++k)
            m_data[k] = function(xcoordinates[i], ycoordinates[j]);
}

auto BilinearInterpolator::setCoordinatesX(const std::vector<real>& xcoordinates) -> void
{
    m_xcoordinates = xcoordinates;
}

auto BilinearInterpolator::setCoordinatesY(const std::vector<real>& ycoordinates) -> void
{
    m_ycoordinates = ycoordinates;
}

auto BilinearInterpolator::setData(const std::vector<real>& data) -> void
{
    m_data = data;
}

auto BilinearInterpolator::xCoodinates() const -> const std::vector<real>&
{
    return m_xcoordinates;
}

auto BilinearInterpolator::yCoodinates() const -> const std::vector<real>&
{
    return m_ycoordinates;
}

auto BilinearInterpolator::data() const -> const std::vector<real>&
{
    return m_data;
}

auto BilinearInterpolator::empty() const -> bool
{
    return m_data.empty();
}

auto BilinearInterpolator::operator()(real x, real y) const -> real
{
    // Check if the interpolation data contains only one point
    if(m_data.size() == 1) return m_data[0];

    const real xA = m_xcoordinates.front();
    const real xB = m_xcoordinates.back();
    const real yA = m_ycoordinates.front();
    const real yB = m_ycoordinates.back();

    x = std::max(xA, std::min(x, xB));
    y = std::max(yA, std::min(y, yB));

    const unsigned sizex = m_xcoordinates.size();
    const unsigned sizey = m_ycoordinates.size();

    const real i = binarySearch(x, m_xcoordinates);
    const real j = binarySearch(y, m_ycoordinates);

    const auto k = [=](unsigned i, unsigned j) { return i + j*sizex; };

    if(i == sizex || j == sizey)
        interpolationOutOfBoundsError(x, xA, xB, y, yA, yB);

    const real x1 = m_xcoordinates[i];
    const real x2 = m_xcoordinates[i + 1];

    const real y1 = m_ycoordinates[j];
    const real y2 = m_ycoordinates[j + 1];

    const real z11 = m_data[k(i  , j  )]; // z at (x1, y1)
    const real z21 = m_data[k(i+1, j  )]; // z at (x2, y1)
    const real z12 = m_data[k(i  , j+1)]; // z at (x1, y2)
    const real z22 = m_data[k(i+1, j+1)]; // z at (x2, y2)

    const real f11 =  z11*(x2 - x)*(y2 - y);
    const real f12 = -z12*(x2 - x)*(y1 - y);
    const real f21 = -z21*(x1 - x)*(y2 - y);
    const real f22 =  z22*(x1 - x)*(y1 - y);

    return (f11 + f12 + f21 + f22)/((x2 - x1)*(y2 - y1));
}

auto operator<<(std::ostream& out, const BilinearInterpolator& interpolator) -> std::ostream&
{
    const auto& xcoordinates = interpolator.xCoodinates();
    const auto& ycoordinates = interpolator.yCoodinates();
    const auto& data         = interpolator.data();

    const unsigned sizex = xcoordinates.size();
    const unsigned sizey = ycoordinates.size();

    out << std::setw(15) << std::right << "y/x";
    for(auto x : xcoordinates)
        out << std::setw(15) << std::right << x;
    out << std::endl;

    for(unsigned j = 0; j < sizey; ++j)
    {
        out << std::setw(15) << std::right << ycoordinates[j];
        for(unsigned i = 0; i < sizex; ++i)
            out << std::setw(15) << std::right << data[i + j*sizex];
        out << std::endl;
    }

    return out;
}

} // namespace Reaktoro
