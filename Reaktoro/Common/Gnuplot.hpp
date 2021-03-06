// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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

#pragma once

// C++ includes
#include <sstream>
#include <string>

namespace Reaktoro {

class Gnuplot
{
public:
    Gnuplot();

    virtual ~Gnuplot();

    auto operator<<(std::string str) -> Gnuplot&;

    auto operator<<(const std::stringstream& str) -> Gnuplot&;

private:
    FILE* pipe;
};

} // namespace Reaktoro
