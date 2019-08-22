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

#include "GaseousMixture.hpp"

namespace Reaktoro {

GaseousMixture::GaseousMixture()
: GeneralMixture<GaseousSpecies>()
{}

GaseousMixture::GaseousMixture(const std::vector<GaseousSpecies>& species)
: GeneralMixture<GaseousSpecies>(species)
{}

GaseousMixture::~GaseousMixture()
{}

auto GaseousMixture::state(const real& T, const real& P, VectorXrConstRef n) const -> GaseousMixtureState
{
    GaseousMixtureState res;
    res.T = T;
    res.P = P;
    res.x = moleFractions(n);
    return res;
}

} // namespace Reaktoro
