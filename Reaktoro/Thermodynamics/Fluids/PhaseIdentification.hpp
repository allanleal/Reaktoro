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
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Real.hpp>

namespace Reaktoro {

/// The options for the fluid types supported by the cubic equation of state.
enum class CubicEOSFluidType
{
    Liquid,
    Vapor
};

/// Defines the enumeration of available phase identification methods.
enum class PhaseIdentificationMethod
{
    None,
    VolumeMethod,
    IsothermalCompressibilityMethods,
    GibbsEnergyAndEquationOfStateMethod,
};

/// Return a CubicEOSFluidType that says if the phase is a Liquid or Gas based on Volume Method
/// @param Temperature Phase temperature
/// @param real Phase pressure
/// @param Z Phase compressibility factor
/// @return The type of the phase
///
/// Reference: Bennett, J. and Schmidt, K.A., 2016. Comparison of Phase Identification Methods Used in Oil Industry Flow Simulations. Energy & Fuels, 31(4), pp.3370-3379.
auto identifyPhaseUsingVolume(
    const real& temperature,
    const real& pressure,
    const real& Z,
    const real& b) -> CubicEOSFluidType;

/// Return a CubicEOSFluidType that says if the phase is a Liquid or Gas based on Isothermal Compressibility
/// @param Temperature Phase temperature
/// @param real Phase pressure
/// @param Z Phase compressibility factor
/// @return The type of the phase
///
/// Reference: Bennett, J. and Schmidt, K.A., 2016. Comparison of Phase Identification Methods Used in Oil Industry Flow Simulations. Energy & Fuels, 31(4), pp.3370-3379.
auto identifyPhaseUsingIsothermalCompressibility(
    const real& temperature,
    const real& pressure,
    const real& Z) -> CubicEOSFluidType;

/// Return a CubicEOSFluidType that says if the phase is a Liquid or Gas based on gibbs residual energy and
/// equation of state
/// @param pressure Phase pressure
/// @param temperature Phase temperature
/// @param amix attractive parameter
/// @param Zs Z-roots, as calculated by the cubic EOS. Must have size 1 or 2 here.
///     If size(Z) == 2, the values od Gibbs residual energy are compared. It is a liquid phase if
///     Gibbs residual energy of Z_min is the smallest and gaseous if Gibbs residual energy of Z_max
///     is the smallest.
///     If size(Z) == 1, the pressue is compared with the local P_min and local P_max of the EoS.
///     It is a liquid phase if P \textgreater P_min and gaseous if P \textless Pmax.
/// @return The type of the phase
///
/// Reference: Bennett, J. and Schmidt, K.A., 2016. Comparison of Phase Identification Methods Used in Oil Industry Flow Simulations. Energy & Fuels, 31(4), pp.3370-3379.
auto identifyPhaseUsingGibbsEnergyAndEos(
    const real& pressure,
    const real& temperature,
    const real& amix,
    const real& bmix,
    const real& A,
    const real& B,
    const real& C,
    std::vector<real> Zs,
    const real epsilon,
    const real sigma) -> CubicEOSFluidType;

}
