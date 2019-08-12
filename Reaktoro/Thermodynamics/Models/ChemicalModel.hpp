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

#pragma once

// C++ includes
#include <functional>

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Models/PhaseChemicalModel.hpp>

namespace Reaktoro {

/// The result of a chemical model function that calculates the chemical properties of species.
class ChemicalModelResult
{
public:
    /// Construct a default ChemicalModelResultBase instance.
    ChemicalModelResult();

    /// Construct a ChemicalModelResultBase instance with allocated memory
    /// @param nphases The number of phases in the chemical system.
    /// @param nspecies The number of species in the chemical system.
    ChemicalModelResult(Index nphases, Index nspecies);

    /// Resize this ChemicalModelResultBase with a given number of species.
    /// @param nphases The number of phases in the chemical system.
    /// @param nspecies The number of species in the chemical system.
    auto resize(Index nphases, Index nspecies) -> void;

    /// Return a view of the chemical properties of a phase.
    /// @param iphase The index of the phase.
    /// @param ispecies The index of the first species in the phase.
    /// @param nspecies The number of species in the phase.
    auto phaseProperties(Index iphase, Index ispecies, Index nspecies) -> PhaseChemicalModelResult;

    /// Return a view of the chemical properties of a phase.
    /// @param iphase The index of the phase.
    /// @param ispecies The index of the first species in the phase.
    /// @param nspecies The number of species in the phase.
    auto phaseProperties(Index iphase, Index ispecies, Index nspecies) const -> PhaseChemicalModelResultConst;

    /// Return the natural log of the activity coefficients of the species.
    inline auto lnActivityCoefficients() -> VectorXrRef { return ln_activity_coefficients; }

    /// Return the natural log of the activity coefficients of the species.
    inline auto lnActivityCoefficients() const -> VectorXrConstRef { return ln_activity_coefficients; }

    /// Return the natural log of the activities of the species.
    inline auto lnActivities() -> VectorXrRef { return ln_activities; }

    /// Return the natural log of the activities of the species.
    inline auto lnActivities() const -> VectorXrConstRef { return ln_activities; }

    /// Return the molar volumes of the phases (in units of m3/mol).
    inline auto phaseMolarVolumes() -> VectorXrRef { return phase_molar_volumes; }

    /// Return the molar volumes of the phases (in units of m3/mol).
    inline auto phaseMolarVolumes() const -> VectorXrConstRef { return phase_molar_volumes; }

    /// Return the residual molar Gibbs energies of the phases w.r.t. to its ideal state (in units of J/mol).
    inline auto phaseResidualMolarGibbsEnergies() -> VectorXrRef { return phase_residual_molar_gibbs_energies; }

    /// Return the residual molar Gibbs energies of the phases w.r.t. to its ideal state (in units of J/mol).
    inline auto phaseResidualMolarGibbsEnergies() const -> VectorXrConstRef { return phase_residual_molar_gibbs_energies; }

    /// Return the residual molar enthalpies of the phases w.r.t. to its ideal state (in units of J/mol).
    inline auto phaseResidualMolarEnthalpies() -> VectorXrRef { return phase_residual_molar_enthalpies; }

    /// Return the residual molar enthalpies of the phases w.r.t. to its ideal state (in units of J/mol).
    inline auto phaseResidualMolarEnthalpies() const -> VectorXrConstRef { return phase_residual_molar_enthalpies; }

    /// Return the residual molar isobaric heat capacities of the phases w.r.t. to its ideal state (in units of J/(mol*K)).
    inline auto phaseResidualMolarHeatCapacitiesCp() -> VectorXrRef { return phase_residual_molar_heat_capacities_cp; }

    /// Return the residual molar isobaric heat capacities of the phases w.r.t. to its ideal state (in units of J/(mol*K)).
    inline auto phaseResidualMolarHeatCapacitiesCp() const -> VectorXrConstRef { return phase_residual_molar_heat_capacities_cp; }

    /// Return the residual molar isochoric heat capacities of the phases w.r.t. to its ideal state (in units of J/(mol*K)).
    inline auto phaseResidualMolarHeatCapacitiesCv() -> VectorXrRef { return phase_residual_molar_heat_capacities_cv; }

    /// Return the residual molar isochoric heat capacities of the phases w.r.t. to its ideal state (in units of J/(mol*K)).
    inline auto phaseResidualMolarHeatCapacitiesCv() const -> VectorXrConstRef { return phase_residual_molar_heat_capacities_cv; }

private:
    /// The natural log of the activity coefficients of the species.
    VectorXr ln_activity_coefficients;

    /// The natural log of the activities of the species.
    VectorXr ln_activities;

    /// The molar volumes of the phases (in units of m3/mol).
    VectorXr phase_molar_volumes;

    /// The residual molar Gibbs energies of the phases w.r.t. to its ideal state (in units of J/mol).
    VectorXr phase_residual_molar_gibbs_energies;

    /// The residual molar enthalpies of the phases w.r.t. to its ideal state (in units of J/mol).
    VectorXr phase_residual_molar_enthalpies;

    /// The residual molar isobaric heat capacities of the phases w.r.t. to its ideal state (in units of J/(mol*K)).
    VectorXr phase_residual_molar_heat_capacities_cp;

    /// The residual molar isochoric heat capacities of the phases w.r.t. to its ideal state (in units of J/(mol*K)).
    VectorXr phase_residual_molar_heat_capacities_cv;
};

/// The signature of the chemical model function that calculates the chemical properties of the species in a chemical system.
using ChemicalModel = std::function<void(ChemicalModelResult&, const real&, const real&, VectorXrConstRef)>;

} // namespace Reaktoro
