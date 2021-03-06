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

// pybind11 includes
#include <Reaktoro/pybind11.hxx>

// Optima includes
#include <Optima/State.hpp>

// Reaktoro includes
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
using namespace Reaktoro;

void exportChemicalState(py::module& m)
{
    const auto return_internal_ref = py::return_value_policy::reference_internal;

    py::class_<ChemicalState>(m, "ChemicalState")
        .def(py::init<const ChemicalSystem&>())
        .def("setTemperature", py::overload_cast<real>(&ChemicalState::setTemperature))
        .def("setTemperature", py::overload_cast<real, String>(&ChemicalState::setTemperature))
        .def("setPressure", py::overload_cast<real>(&ChemicalState::setPressure))
        .def("setPressure", py::overload_cast<real, String>(&ChemicalState::setPressure))
        .def("setSpeciesAmounts", [](ChemicalState& s, double val) { s.setSpeciesAmounts(val); })
        .def("setSpeciesAmounts", [](ChemicalState& s, const real& val) { s.setSpeciesAmounts(val); })
        .def("setSpeciesAmounts", [](ChemicalState& s, ArrayXrConstRef vals) { s.setSpeciesAmounts(vals); })
        .def("setSpeciesAmounts", [](ChemicalState& s, const py::array_t<double>& vals) { s.setSpeciesAmounts(ArrayXd::Map(vals.data(), vals.size())); })
        // .def("setSpeciesAmounts", [](ChemicalState& s, VectorXdConstRef vals) { s.setSpeciesAmounts(vals.array()); })
        // .def("setSpeciesAmounts", py::overload_cast<real>(&ChemicalState::setSpeciesAmounts))
        // .def("setSpeciesAmounts", py::overload_cast<ArrayXrConstRef>(&ChemicalState::setSpeciesAmounts))
        // .def("setSpeciesAmounts", py::overload_cast<ArrayXdConstRef>(&ChemicalState::setSpeciesAmounts))
        .def("setSpeciesAmount", py::overload_cast<Index, real>(&ChemicalState::setSpeciesAmount))
        .def("setSpeciesAmount", py::overload_cast<Index, real, String>(&ChemicalState::setSpeciesAmount))
        .def("setSpeciesAmount", py::overload_cast<String, real>(&ChemicalState::setSpeciesAmount))
        .def("setSpeciesAmount", py::overload_cast<String, real, String>(&ChemicalState::setSpeciesAmount))
        .def("setSpeciesMass", py::overload_cast<Index, real>(&ChemicalState::setSpeciesMass))
        .def("setSpeciesMass", py::overload_cast<Index, real, String>(&ChemicalState::setSpeciesMass))
        .def("setSpeciesMass", py::overload_cast<String, real>(&ChemicalState::setSpeciesMass))
        .def("setSpeciesMass", py::overload_cast<String, real, String>(&ChemicalState::setSpeciesMass))
        .def("system", &ChemicalState::system, return_internal_ref)
        .def("temperature", &ChemicalState::temperature)
        .def("pressure", &ChemicalState::pressure)
        .def("speciesAmounts", &ChemicalState::speciesAmounts, return_internal_ref)
        .def("elementAmounts", &ChemicalState::elementAmounts)
        .def("speciesAmount", py::overload_cast<Index>(&ChemicalState::speciesAmount, py::const_))
        .def("speciesAmount", py::overload_cast<Index, String>(&ChemicalState::speciesAmount, py::const_))
        .def("speciesAmount", py::overload_cast<String>(&ChemicalState::speciesAmount, py::const_))
        .def("speciesAmount", py::overload_cast<String, String>(&ChemicalState::speciesAmount, py::const_))
        .def("speciesMass", py::overload_cast<Index>(&ChemicalState::speciesMass, py::const_))
        .def("speciesMass", py::overload_cast<Index, String>(&ChemicalState::speciesMass, py::const_))
        .def("speciesMass", py::overload_cast<String>(&ChemicalState::speciesMass, py::const_))
        .def("speciesMass", py::overload_cast<String, String>(&ChemicalState::speciesMass, py::const_))
        .def("equilibrium", py::overload_cast<>(&ChemicalState::equilibrium, py::const_), return_internal_ref)
        .def("equilibrium", py::overload_cast<>(&ChemicalState::equilibrium), return_internal_ref)
        .def("props", py::overload_cast<>(&ChemicalState::props, py::const_), return_internal_ref)
        .def("props", py::overload_cast<>(&ChemicalState::props), return_internal_ref)
        ;

    py::class_<ChemicalState::Equilibrium>(m, "_ChemicalStateEquilibrium")
        .def("setInputNames", &ChemicalState::Equilibrium::setInputNames)
        .def("setInputValues", &ChemicalState::Equilibrium::setInputValues)
        .def("setInitialComponentAmounts", &ChemicalState::Equilibrium::setInitialComponentAmounts)
        .def("setControlVariablesP", &ChemicalState::Equilibrium::setControlVariablesP)
        .def("setControlVariablesQ", &ChemicalState::Equilibrium::setControlVariablesQ)
        .def("setOptimaState", &ChemicalState::Equilibrium::setOptimaState)
        .def("setIndicesPrimarySecondarySpecies", &ChemicalState::Equilibrium::setIndicesPrimarySecondarySpecies)
        .def("setIndicesStrictlyUnstableElements", &ChemicalState::Equilibrium::setIndicesStrictlyUnstableElements)
        .def("setIndicesStrictlyUnstableSpecies", &ChemicalState::Equilibrium::setIndicesStrictlyUnstableSpecies)
        .def("numPrimarySpecies", &ChemicalState::Equilibrium::numPrimarySpecies)
        .def("numSecondarySpecies", &ChemicalState::Equilibrium::numSecondarySpecies)
        .def("indicesPrimarySpecies", &ChemicalState::Equilibrium::indicesPrimarySpecies, return_internal_ref)
        .def("indicesSecondarySpecies", &ChemicalState::Equilibrium::indicesSecondarySpecies, return_internal_ref)
        .def("indicesStrictlyUnstableElements", &ChemicalState::Equilibrium::indicesStrictlyUnstableElements, return_internal_ref)
        .def("indicesStrictlyUnstableSpecies", &ChemicalState::Equilibrium::indicesStrictlyUnstableSpecies, return_internal_ref)
        .def("elementChemicalPotentials", &ChemicalState::Equilibrium::elementChemicalPotentials, return_internal_ref)
        .def("speciesStabilities", &ChemicalState::Equilibrium::speciesStabilities, return_internal_ref)
        .def("explicitTitrantAmounts", &ChemicalState::Equilibrium::explicitTitrantAmounts, return_internal_ref)
        .def("implicitTitrantAmounts", &ChemicalState::Equilibrium::implicitTitrantAmounts, return_internal_ref)
        .def("inputNames", &ChemicalState::Equilibrium::inputNames, return_internal_ref)
        .def("inputValues", &ChemicalState::Equilibrium::inputValues, return_internal_ref)
        .def("initialComponentAmounts", &ChemicalState::Equilibrium::initialComponentAmounts, return_internal_ref)
        .def("controlVariablesP", &ChemicalState::Equilibrium::controlVariablesP, return_internal_ref)
        .def("controlVariablesQ", &ChemicalState::Equilibrium::controlVariablesQ, return_internal_ref)
        .def("p", &ChemicalState::Equilibrium::p, return_internal_ref)
        .def("q", &ChemicalState::Equilibrium::q, return_internal_ref)
        .def("w", &ChemicalState::Equilibrium::w, return_internal_ref)
        .def("b", &ChemicalState::Equilibrium::b, return_internal_ref)
        .def("optimaState", &ChemicalState::Equilibrium::optimaState, return_internal_ref)
        ;
}
