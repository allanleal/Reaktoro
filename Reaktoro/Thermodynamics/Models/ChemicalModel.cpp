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

#include "ChemicalModel.hpp"

namespace Reaktoro {

ChemicalModelResult::ChemicalModelResult()
{}

ChemicalModelResult::ChemicalModelResult(Index nphases, Index nspecies)
: ln_activity_coefficients(nspecies),
  ln_activities(nspecies),
  phase_molar_volumes(nphases),
  phase_residual_molar_gibbs_energies(nphases),
  phase_residual_molar_enthalpies(nphases),
  phase_residual_molar_heat_capacities_cp(nphases),
  phase_residual_molar_heat_capacities_cv(nphases)
{}

auto ChemicalModelResult::resize(Index nphases, Index nspecies) -> void
{
    ln_activity_coefficients.resize(nspecies);
    ln_activities.resize(nspecies);
    phase_molar_volumes.resize(nphases);
    phase_residual_molar_gibbs_energies.resize(nphases);
    phase_residual_molar_enthalpies.resize(nphases);
    phase_residual_molar_heat_capacities_cp.resize(nphases);
    phase_residual_molar_heat_capacities_cv.resize(nphases);
}

auto ChemicalModelResult::phaseProperties(Index iphase, Index ispecies, Index nspecies) -> PhaseChemicalModelResult
{
    return {
        ln_activity_coefficients.segment(ispecies, nspecies),
        ln_activities.segment(ispecies, nspecies),
        phase_molar_volumes[iphase],
        phase_residual_molar_gibbs_energies[iphase],
        phase_residual_molar_enthalpies[iphase],
        phase_residual_molar_heat_capacities_cp[iphase],
        phase_residual_molar_heat_capacities_cv[iphase]
    };
}

auto ChemicalModelResult::phaseProperties(Index iphase, Index ispecies, Index nspecies) const -> PhaseChemicalModelResultConst
{
    return {
        ln_activity_coefficients.segment(ispecies, nspecies),
        ln_activities.segment(ispecies, nspecies),
        phase_molar_volumes[iphase],
        phase_residual_molar_gibbs_energies[iphase],
        phase_residual_molar_enthalpies[iphase],
        phase_residual_molar_heat_capacities_cp[iphase],
        phase_residual_molar_heat_capacities_cv[iphase]
    };
}

} // namespace Reaktoro
