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

// Reaktoro includes
#include <Reaktoro/Common/Real.hpp>

namespace Reaktoro {

struct WaterThermoState
{
	/// The temperature of water (in units of K)
	real temperature = 0.0;

	/// The specific volume of water (in units of m3/kg)
	real volume = 0.0;

	/// The specific entropy of water (in units of J/(kg*K))
	real entropy = 0.0;

	/// The specific Helmholtz free energy of water (in units of J/kg)
	real helmholtz = 0.0;

	/// The specific internal energy of water (in units of J/kg)
	real internal_energy = 0.0;

	/// The specific enthalpy of water (in units of J/kg)
	real enthalpy = 0.0;

	/// The specific Gibbs free energy of water (in units of J/kg)
	real gibbs = 0.0;

	/// The specific isochoric heat capacity of water (in units of J/(kg*K))
	real cv = 0.0;

	/// The specific isobaric heat capacity of water (in units of J/(kg*K))
	real cp = 0.0;

	/// The specific density of water (in units of kg/m3)
	real density = 0.0;

	/// The first-order partial derivative of density with respect to temperature (in units of (kg/m3)/K)
	real densityT = 0.0;

	/// The first-order partial derivative of density with respect to pressure (in units of (kg/m3)/Pa)
	real densityP = 0.0;

	/// The second-order partial derivative of density with respect to temperature (in units of (kg/m3)/(K*K))
	real densityTT = 0.0;

	/// The second-order partial derivative of density with respect to temperature and pressure (in units of (kg/m3)/(K*Pa))
	real densityTP = 0.0;

	/// The second-order partial derivative of density with respect to pressure (in units of (kg/m3)/(Pa*Pa))
	real densityPP = 0.0;

	/// The pressure of water (in units of Pa)
	real pressure = 0.0;

	/// The first-order partial derivative of pressure with respect to temperature (in units of Pa/K)
	real pressureT = 0.0;

	/// The first-order partial derivative of pressure with respect to density (in units of Pa/(kg/m3))
	real pressureD = 0.0;

	/// The second-order partial derivative of pressure with respect to temperature (in units of Pa/(K*K))
	real pressureTT = 0.0;

	/// The second-order partial derivative of pressure with respect to temperature and density (in units of Pa/(K*kg/m3))
	real pressureTD = 0.0;

	/// The second-order partial derivative of pressure with respect to density (in units of Pa/((kg/m3)*(kg/m3)))
	real pressureDD = 0.0;
};

} // namespace Reaktoro
