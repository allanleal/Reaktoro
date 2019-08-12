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

#include "CubicEOS.hpp"

// C++ includes
#include <algorithm>

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Common/TableUtils.hpp>
#include <Reaktoro/Math/Roots.hpp>

namespace Reaktoro {
namespace internal {

using AlphaResult = std::tuple<real, real, real>;

/// A high-order function that return an `alpha` function that calculates alpha, alphaT and alphaTT (temperature derivatives) for a given EOS.
auto alpha(CubicEOS::Model type) -> std::function<AlphaResult(const real&, const real&)>
{
    // The alpha function for van der Waals EOS
    auto alphaVDW = [](const real& Tr, const real& omega) -> AlphaResult
    {
        const auto val = 1.0;
        const auto ddT = 0.0;
        const auto d2dT2 = 0.0;
        return std::make_tuple(val, ddT, d2dT2);
    };

    // The alpha function for Redlich-Kwong EOS
    auto alphaRK = [](const real& T, const real& omega) -> AlphaResult
    {
        const auto sqrtT = std::sqrt(T);
        const auto val = 1.0/sqrtT;
        const auto ddT = -0.5 / T * val;
        const auto d2dT2 = 0.5/(T*T) * val - 0.5/T * ddT;
        return std::make_tuple(val, ddT, d2dT2);
    };

    // The alpha function for Soave-Redlich-Kwong EOS
    auto alphaSRK = [](const real& T, const real& omega) -> AlphaResult
    {
        const auto m = 0.480 + 1.574*omega - 0.176*omega*omega;
        const auto sqrtT = std::sqrt(T);
        const auto aux_val = 1.0 + m*(1.0 - sqrtT);
        const auto aux_ddT = -0.5*m/sqrtT;
        const auto aux_d2dT2 = 0.25*m/(T*sqrtT);
        const auto val = aux_val*aux_val;
        const auto ddT = 2.0*aux_val*aux_ddT;
        const auto d2dT2 = 2.0*(aux_ddT*aux_ddT + aux_val*aux_d2dT2);
        return std::make_tuple(val, ddT, d2dT2);
    };

    // The alpha function for Peng-Robinson (1978) EOS
    auto alphaPR = [](const real& T, const real& omega) -> AlphaResult
    {
        // Jaubert, J.-N., Vitu, S., Mutelet, F. and Corriou, J.-P., 2005.
        // Extension of the PPR78 model (predictive 1978, Peng–Robinson EOS
        // with temperature dependent kij calculated through a group
        // contribution method) to systems containing aromatic compounds.
        // Fluid Phase Equilibria, 237(1-2), pp.193–211.
        const auto m = omega <= 0.491 ?
            0.374640 + 1.54226*omega - 0.269920*omega*omega :
            0.379642 + 1.48503*omega - 0.164423*omega*omega + 0.016666*omega*omega*omega;
        const auto sqrtT = std::sqrt(T);
        const auto aux_val = 1.0 + m*(1.0 - sqrtT);
        const auto aux_ddT = -0.5*m/sqrtT;
        const auto aux_d2dT2 = 0.25*m/(T*sqrtT);
        const auto val = aux_val*aux_val;
        const auto ddT = 2.0*aux_val*aux_ddT;
        const auto d2dT2 = 2.0*(aux_ddT*aux_ddT + aux_val*aux_d2dT2);
        return std::make_tuple(val, ddT, d2dT2);
    };

    switch(type)
    {
        case CubicEOS::VanDerWaals: return alphaVDW;
        case CubicEOS::RedlichKwong: return alphaRK;
        case CubicEOS::SoaveRedlichKwong: return alphaSRK;
        case CubicEOS::PengRobinson: return alphaPR;
        default: return alphaPR;
    }
}

auto sigma(CubicEOS::Model type) -> double
{
    switch(type)
    {
        case CubicEOS::VanDerWaals: return 0.0;
        case CubicEOS::RedlichKwong: return 1.0;
        case CubicEOS::SoaveRedlichKwong: return 1.0;
        case CubicEOS::PengRobinson: return 1.0 + 1.4142135623730951;
        default: return 1.0 + 1.4142135623730951;
    }
}

auto epsilon(CubicEOS::Model type) -> double
{
    switch(type)
    {
        case CubicEOS::VanDerWaals: return 0.0;
        case CubicEOS::RedlichKwong: return 0.0;
        case CubicEOS::SoaveRedlichKwong: return 0.0;
        case CubicEOS::PengRobinson: return 1.0 - 1.4142135623730951;
        default: return 1.0 - 1.4142135623730951;
    }
}

auto Omega(CubicEOS::Model type) -> double
{
    switch(type)
    {
        case CubicEOS::VanDerWaals: return 1.0/8.0;
        case CubicEOS::RedlichKwong: return 0.08664;
        case CubicEOS::SoaveRedlichKwong: return 0.08664;
        case CubicEOS::PengRobinson: return 0.0777960739;
        default: return 0.0777960739;
    }
}

auto Psi(CubicEOS::Model type) -> double
{
    switch(type)
    {
        case CubicEOS::VanDerWaals: return 27.0/64.0;
        case CubicEOS::RedlichKwong: return 0.42748;
        case CubicEOS::SoaveRedlichKwong: return 0.42748;
        case CubicEOS::PengRobinson: return 0.457235529;
        default: return 0.457235529;
    }
}

} // namespace internal

struct CubicEOS::Impl
{
    /// The number of species in the phase.
    unsigned nspecies;

    /// The flag that indicates if the phase is vapor (false means liquid instead).
    bool isvapor = true;

    /// The type of the cubic equation of state.
    CubicEOS::Model model = CubicEOS::PengRobinson;

    /// The critical temperatures of the species (in units of K).
    std::vector<real> critical_temperatures;

    /// The critical pressures of the species (in units of Pa).
    std::vector<real> critical_pressures;

    /// The acentric factors of the species.
    std::vector<real> acentric_factors;

    /// The function that calculates the interaction parameters kij and its temperature derivatives.
    InteractionParamsFunction calculate_interaction_params;

    /// The result with thermodynamic properties calculated from the cubic equation of state
    Result result;

    /// Construct a CubicEOS::Impl instance.
    Impl(unsigned nspecies)
    : nspecies(nspecies)
    {
        // Initialize the dimension of the chemical vector quantities
        VectorXr vec(nspecies);
        result.partial_molar_volumes = vec;
        result.residual_partial_molar_enthalpies = vec;
        result.residual_partial_molar_gibbs_energies = vec;
        result.ln_fugacity_coefficients = vec;
    }

    auto operator()(const real& T, const real& P, const VectorXr& x) -> Result
    {
        // Check if the mole fractions are zero or non-initialized
        if(x.size() == 0 || min(x) <= 0.0)
            return Result(nspecies); // result with zero values

        // Auxiliary variables
        const auto R = universalGasConstant;
        const auto Psi = internal::Psi(model);
        const auto Omega = internal::Omega(model);
        const auto epsilon = internal::epsilon(model);
        const auto sigma = internal::sigma(model);
        const auto alpha = internal::alpha(model);

        // Calculate the parameters `a` of the cubic equation of state for each species
        VectorXr a(nspecies);
        VectorXr aT(nspecies);
        VectorXr aTT(nspecies);
        for(unsigned i = 0; i < nspecies; ++i)
        {
            const auto Tc = critical_temperatures[i];
            const auto Pc = critical_pressures[i];
            const auto omega = acentric_factors[i];
            const auto factor = Psi*R*R*(Tc*Tc)/Pc;
            real alpha_val, alpha_ddT, alpha_d2dT2;
            std::tie(alpha_val, alpha_ddT, alpha_d2dT2) = alpha(T, omega);
            a[i] = factor * alpha_val;
            aT[i] = factor * alpha_ddT;
            aTT[i] = factor * alpha_d2dT2;
        };

        // Calculate the parameters `b` of the cubic equation of state for each species
        Vector b(nspecies);
        for(unsigned i = 0; i < nspecies; ++i)
        {
            const double Tci = critical_temperatures[i];
            const double Pci = critical_pressures[i];
            b[i] = Omega*R*Tci/Pci;
        }

        // Calculate the table of binary interaction parameters and its temperature derivatives
        InteractionParamsResult kres;
        InteractionParamsArgs kargs{T, a, aT, aTT, b};

        if(calculate_interaction_params)
            kres = calculate_interaction_params(kargs);

        // Calculate the parameter `amix` of the phase and the partial molar parameters `abar` of each species
        real amix = 0.0;
        real amixT = 0.0;
        real amixTT = 0.0;
        VectorXr abar = zeros(nspecies);
        VectorXr abarT = zeros(nspecies);
        for(unsigned i = 0; i < nspecies; ++i)
        {
            for(unsigned j = 0; j < nspecies; ++j)
            {
                const auto r = kres.k.empty() ? 1.0 : 1.0 - kres.k[i][j];
                const auto rT = kres.kT.empty() ? 0.0 : -kres.kT[i][j];
                const auto rTT = kres.kTT.empty() ? 0.0 : -kres.kTT[i][j];

                const auto s = std::sqrt(a[i]*a[j]);
                const auto sT = 0.5*s/(a[i]*a[j]) * (aT[i]*a[j] + a[i]*aT[j]);
                const auto sTT = 0.5*s/(a[i]*a[j]) * (aTT[i]*a[j] + 2*aT[i]*aT[j] + a[i]*aTT[j]) - sT*sT/s;

                const auto aij = r*s;
                const auto aijT = rT*s + r*sT;
                const auto aijTT = rTT*s + 2.0*rT*sT + r*sTT;

                amix += x[i] * x[j] * aij;
                amixT += x[i] * x[j] * aijT;
                amixTT += x[i] * x[j] * aijTT;

                abar[i] += 2 * x[j] * aij;
                abarT[i] += 2 * x[j] * aijT;
            }
        }

        // Finalize the calculation of `abar` and `abarT`
        for(unsigned i = 0; i < nspecies; ++i)
        {
            abar[i] -= amix;
            abarT[i] -= amixT;
        }

        // Calculate the parameter `bmix` of the cubic equation of state
        real bmix = 0.0;
        Vector bbar = zeros(nspecies);
        for(unsigned i = 0; i < nspecies; ++i)
        {
            const auto Tci = critical_temperatures[i];
            const auto Pci = critical_pressures[i];
            bbar[i] = Omega*R*Tci/Pci;
            bmix += x[i] * bbar[i];
        }

        // Calculate the temperature derivative of `bmix`
        const auto bmixT = 0.0; // no temperature dependence

        // Calculate auxiliary quantities `beta` and `q`
        const auto beta = P*bmix/(R*T);
        const auto betaT = beta * (bmixT/bmix - 1.0/T);

        const auto q = amix/(bmix*R*T);
        const auto qT = q*(amixT/amix - 1.0/T);
        const auto qTT = qT*qT/q + q*(1.0/(T*T) + amixTT/amix - amixT*amixT/(amix*amix));

        // Calculate the coefficients A, B, C of the cubic equation of state
        const auto A = (epsilon + sigma - 1)*beta - 1;
        const auto B = (epsilon*sigma - epsilon - sigma)*beta*beta - (epsilon + sigma - q)*beta;
        const auto C = -epsilon*sigma*beta*beta*beta - (epsilon*sigma + q)*beta*beta;

        // Calculate the partial temperature derivative of the coefficients A, B, C
        const auto AT = (epsilon + sigma - 1)*betaT;
        const auto BT = 2*(epsilon*sigma - epsilon - sigma)*beta*betaT + qT*beta - (epsilon + sigma - q)*betaT;
        const auto CT = -3*epsilon*sigma*beta*beta*betaT - qT*beta*beta - 2*(epsilon*sigma + q)*beta*betaT;

        // Define the non-linear function and its derivative for calculation of its root
        const auto f = [&](const real& Z) -> std::tuple<real, real>
        {
            const auto val = Z*Z*Z + A*Z*Z + B*Z + C;
            const auto grad = 3*Z*Z + 2*A*Z + B;
            return std::make_tuple(val, grad);
        };

        // Define the parameters for Newton's method
        const auto tolerance = 1e-6;
        const auto maxiter = 100;

        // Determine the appropriate initial guess for the cubic equation of state
        const real Z0 = isvapor ? 1.0 : beta;

        // Calculate the compressibility factor Z using Newton's method
        real Z = newton(f, Z0, tolerance, maxiter);

        // Calculate the partial derivatives of Z (dZdT, dZdP, dZdn)
        const auto factor = -1.0/(3*Z*Z + 2*A*Z + B);

        // Calculate the partial temperature derivative of Z
        const real ZT = -(AT*Z*Z + BT*Z + CT)/(3*Z*Z + 2*A*Z + B);

        // Calculate the integration factor I and its temperature derivative IT
        real I;
        if(epsilon != sigma) I = std::log((Z + sigma*beta)/(Z + epsilon*beta))/(sigma - epsilon);
                        else I = beta/(Z + epsilon*beta);

        // Calculate the temperature derivative IT of the integration factor I
        real IT;
        if(epsilon != sigma) IT = ((ZT + sigma*betaT)/(Z + sigma*beta) - (ZT + epsilon*betaT)/(Z + epsilon*beta))/(sigma - epsilon);
                        else IT = I*(betaT/beta - (ZT + epsilon*betaT)/(Z + epsilon*beta));

        real& V = result.molar_volume;
        real& G_res = result.residual_molar_gibbs_energy;
        real& H_res = result.residual_molar_enthalpy;
        real& Cp_res = result.residual_molar_heat_capacity_cp;
        real& Cv_res = result.residual_molar_heat_capacity_cv;
        VectorXr& Vi = result.partial_molar_volumes;
        VectorXr& Gi_res = result.residual_partial_molar_gibbs_energies;
        VectorXr& Hi_res = result.residual_partial_molar_enthalpies;
        VectorXr& ln_phi = result.ln_fugacity_coefficients;

        // Calculate the partial molar Zi for each species
        V = Z*R*T/P;
        G_res = R*T*(Z - 1 - std::log(Z - beta) - q*I);
        H_res = R*T*(Z - 1 + T*qT*I);
        Cp_res = R*T*(ZT + qT*I + T*qTT + T*qT*IT) + H_res/T;

        const real dPdT = P*(1.0/T + ZT/Z);
        const real dVdT = V*(1.0/T + ZT/Z);

        Cv_res = Cp_res - T*dPdT*dVdT + R;

        for(unsigned i = 0; i < nspecies; ++i)
        {
            const real bi = bbar[i];
            const real betai = P*bi/(R*T);
            const real ai = abar[i];
            const real aiT = abarT[i];
            const real qi = q*(1 + ai/amix - bi/bmix);
            const real qiT = qi*qT/q + q*(aiT - ai*amixT/amix)/amix;
            const real Ai = (epsilon + sigma - 1.0)*betai - 1.0;
            const real Bi = (epsilon*sigma - epsilon - sigma)*(2*beta*betai - beta*beta) - (epsilon + sigma - q)*(betai - beta) - (epsilon + sigma - qi)*beta;
            const real Ci = -3*sigma*epsilon*beta*beta*betai + 2*epsilon*sigma*beta*beta*beta - (epsilon*sigma + qi)*beta*beta - 2*(epsilon*sigma + q)*(beta*betai - beta*beta);
            const real Zi = -(Ai*Z*Z + (Bi + B)*Z + Ci + 2*C)/(3*Z*Z + 2*A*Z + B);
            real Ii;
            if(epsilon != sigma) Ii = I + ((Zi + sigma*betai)/(Z + sigma*beta) - (Zi + epsilon*betai)/(Z + epsilon*beta))/(sigma - epsilon);
                            else Ii = I * (1 + betai/beta - (Zi + epsilon*betai)/(Z + epsilon*beta));

            Vi[i] = R*T*Zi/P;
            Gi_res[i] = R*T*(Zi - (Zi - betai)/(Z - beta) - std::log(Z - beta) - qi*I - q*Ii + q*I);
            Hi_res[i] = R*T*(Zi - 1 + T*(qiT*I + qT*Ii - qT*I));
            ln_phi[i] = Gi_res[i]/(R*T);
        }

        return result;
    }
};

CubicEOS::Result::Result()
{}

CubicEOS::Result::Result(unsigned nspecies)
: molar_volume(nspecies),
  residual_molar_gibbs_energy(nspecies),
  residual_molar_enthalpy(nspecies),
  residual_molar_heat_capacity_cp(nspecies),
  residual_molar_heat_capacity_cv(nspecies),
  partial_molar_volumes(nspecies),
  residual_partial_molar_gibbs_energies(nspecies),
  residual_partial_molar_enthalpies(nspecies),
  ln_fugacity_coefficients(nspecies)
{}

CubicEOS::CubicEOS(unsigned nspecies)
: pimpl(new Impl(nspecies))
{}

CubicEOS::CubicEOS(const CubicEOS& other)
: pimpl(new Impl(*other.pimpl))
{}

CubicEOS::~CubicEOS()
{}

auto CubicEOS::operator=(CubicEOS other) -> CubicEOS&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto CubicEOS::numSpecies() const -> unsigned
{
    return pimpl->nspecies;
}

auto CubicEOS::setModel(Model model) -> void
{
    pimpl->model = model;
}

auto CubicEOS::setPhaseAsLiquid() -> void
{
    pimpl->isvapor = false;
}

auto CubicEOS::setPhaseAsVapor() -> void
{
    pimpl->isvapor = true;
}

auto CubicEOS::setCriticalTemperatures(const std::vector<double>& values) -> void
{
    Assert(values.size() == numSpecies(), "Cannot set the critical "
        "temperatures of the species in the CubicEOS object.", "Expecting " +
        std::to_string(numSpecies()) + " values, but only " +
        std::to_string(values.size()) + " were given.");

    Assert(minValue(values) > 0, "Cannot set the critical temperatures of the species "
        "in the CubicEOS object.", "Expecting non-zero critical "
        "temperatures of the gases.");

    pimpl->critical_temperatures = values;
}

auto CubicEOS::setCriticalPressures(const std::vector<double>& values) -> void
{
    Assert(values.size() == numSpecies(), "Cannot set the critical "
        "pressures of the species in the CubicEOS object.", "Expecting " +
        std::to_string(numSpecies()) + " values, but only " +
        std::to_string(values.size()) + " were given.");

    Assert(minValue(values) > 0, "Cannot set the critical pressures of the species "
        "in the CubicEOS object.", "Expecting non-zero critical "
        "pressures of the gases.");

    pimpl->critical_pressures = values;
}

auto CubicEOS::setAcentricFactors(const std::vector<double>& values) -> void
{
    Assert(values.size() == numSpecies(), "Cannot set the acentric "
        "factors of the species in CubicEOS.", "Expecting " +
        std::to_string(numSpecies()) + " values, but only " +
        std::to_string(values.size()) + " values were given.");

    pimpl->acentric_factors = values;
}

auto CubicEOS::setInteractionParamsFunction(const InteractionParamsFunction& func) -> void
{
    pimpl->calculate_interaction_params = func;
}

auto CubicEOS::operator()(const real& T, const real& P, const VectorXr& x) -> Result
{
    return pimpl->operator()(T, P, x);
}

} // namespace Reaktoro
