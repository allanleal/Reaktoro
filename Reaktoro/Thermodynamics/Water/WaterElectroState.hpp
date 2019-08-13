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

struct WaterElectroState
{
	/// The dielectric constant of water
	real epsilon = 0.0;

	/// The first-order partial derivative of the dielectric constant with respect to temperature
	real epsilonT = 0.0;

	/// The first-order partial derivative of the dielectric constant with respect to pressure
	real epsilonP = 0.0;

	/// The second-order partial derivative of the dielectric constant with respect to temperature
	real epsilonTT = 0.0;

	/// The second-order partial derivative of the dielectric constant with respect to temperature and pressure
	real epsilonTP = 0.0;

	/// The second-order partial derivative of the dielectric constant with respect to pressure
	real epsilonPP = 0.0;

	/// The Born function \f$ Z\equiv-\frac{1}{\epsilon} \f$ (see Helgeson and Kirkham, 1974)
	real bornZ = 0.0;

	/// The Born function \f$ Y\equiv\left[\frac{\partial Z}{\partial T}\right]_{P} \f$ (see Helgeson and Kirkham, 1974)
	real bornY = 0.0;

	/// The Born function \f$ Q\equiv\left[\frac{\partial Z}{\partial P}\right]_{T} \f$ (see Helgeson and Kirkham, 1974)
	real bornQ = 0.0;

	/// The Born function \f$ N\equiv\left[\frac{\partial Q}{\partial P}\right]_{T} \f$ (see Helgeson and Kirkham, 1974)
	real bornN = 0.0;

	/// The Born function \f$ U\equiv\left[\frac{\partial Q}{\partial T}\right]_{P} \f$ (see Helgeson and Kirkham, 1974)
	real bornU = 0.0;

	/// The Born function \f$ X\equiv\left[\frac{\partial Y}{\partial T}\right]_{P} \f$ (see Helgeson and Kirkham, 1974)
	real bornX = 0.0;
};

} // namespace Reaktoro
