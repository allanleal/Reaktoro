// Reaktoro is a unified framework for modeling chemically reactive phases.
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

// Reaktoro includes
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ThermoPropsPhase.hpp>

namespace Reaktoro {

/// The standard thermodynamic properties of the species in a chemical system.
class ThermoProps
{
public:
    /// Construct a ThermoProps object.
    explicit ThermoProps(const ChemicalSystem& system);

    /// Construct a copy of a ThermoProps object.
    ThermoProps(const ThermoProps& other);

    /// Destroy this ThermoProps object.
    ~ThermoProps();

    /// Assign a ThermoProps object to this.
    auto operator=(ThermoProps other) -> ThermoProps&;

    /// Update the standard thermodynamic properties of the species in the chemical system.
    /// @param T The temperature condition (in K)
    /// @param P The pressure condition (in Pa)
    auto update(const real& T, const real& P) -> void;

    /// Update the standard thermodynamic properties of the species in the chemical system.
    /// @param T The temperature condition (in K)
    /// @param P The pressure condition (in Pa)
    /// @param wrtvar The variable with respect to automatic differentiation should be carried out.
    auto update(const real& T, const real& P, Wrt<real&> wrtvar) -> void;

    /// Return the chemical system associated with these standard thermodynamic properties.
    auto system() const -> const ChemicalSystem&;

    /// Return the standard thermodynamic properties of a phase with given index.
    auto phaseProps(Index idx) const -> ThermoPropsPhaseConstRef;

    /// Return the temperature of the system (in K).
    auto temperature() const -> const real&;

    /// Return the pressure of the system (in Pa).
    auto pressure() const -> const real&;

    /// Return the standard partial molar Gibbs energies of the species in the system (in J/mol).
    auto standardGibbsEnergies() const -> ArrayXrConstRef;

    /// Return the standard partial molar enthalpies of the species in the system (in J/mol).
    auto standardEnthalpies() const -> ArrayXrConstRef;

    /// Return the standard partial molar volumes of the species in the system (in m3/mol).
    auto standardVolumes() const -> ArrayXrConstRef;

    /// Return the standard partial molar entropies of the species in the system (in J/(mol*K)).
    auto standardEntropies() const -> ArrayXr;

    /// Return the standard partial molar internal energies of the species in the system (in J/mol).
    auto standardInternalEnergies() const -> ArrayXr;

    /// Return the standard partial molar Helmholtz energies of the species in the system (in J/mol).
    auto standardHelmholtzEnergies() const -> ArrayXr;

    /// Return the standard partial molar isobaric heat capacities of the species in the system (in J/(mol*K)).
    auto standardHeatCapacitiesConstP() const -> ArrayXrConstRef;

    /// Return the standard partial molar isochoric heat capacities of the species in the system (in J/(mol*K)).
    auto standardHeatCapacitiesConstV() const -> ArrayXrConstRef;

private:
    struct Impl;

    Ptr<Impl> pimpl;
};

} // namespace Reaktoro
