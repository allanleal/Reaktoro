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

// Catch includes
#include <catch2/catch.hpp>

// Reaktoro includes
#include <Reaktoro/Extensions/ThermoFun/ThermoFunDatabase.hpp>
using namespace Reaktoro;

TEST_CASE("Testing ThermoFunDatabase", "[ThermoFunDatabase]")
{
    ThermoFunDatabase aq17("aq17");
    ThermoFunDatabase cemdata18("cemdata18");
    ThermoFunDatabase heracles("heracles");
    ThermoFunDatabase mines16("mines16");
    ThermoFunDatabase psinagra("psinagra-12-07");
    ThermoFunDatabase slop98organic("slop98-organic");
    ThermoFunDatabase slop98("slop98");

    const auto T = 298.15;
    const auto P = 1.0e5;

    Species species;
    StandardThermoProps props;

    //-------------------------------------------------------------------
    // Testing attributes and thermodynamic properties of H2O@
    //-------------------------------------------------------------------
    species = aq17.species().get("H2O@");
    CHECK( species.formula().equivalent("H2O") );
    CHECK( species.substance() == "Water HGK" );
    CHECK( species.aggregateState() == AggregateState::Aqueous );
    CHECK( species.charge() == 0 );
    CHECK( species.molarMass() == Approx(0.0180153) );

    props = species.props(T, P);
    CHECK( props.G0  == Approx(-2.371817e+05) );
    CHECK( props.H0  == Approx(-2.858310e+05) );
    CHECK( props.V0  == Approx( 1.806862e-05) );
    CHECK( props.Cp0 == Approx( 7.532758e+01) );
    CHECK( props.Cv0 == Approx( 7.453942e+01) );

    //-------------------------------------------------------------------
    // Testing attributes and thermodynamic properties of CO3-2
    //-------------------------------------------------------------------
    species = aq17.species().get("CO3-2");
    CHECK( species.formula().equivalent("CO3-2") );
    CHECK( species.substance() == "CO3-2 carbonate ion" );
    CHECK( species.aggregateState() == AggregateState::Aqueous );
    CHECK( species.charge() == -2 );
    CHECK( species.molarMass() == Approx(0.0600100979) );

    props = species.props(T, P);
    CHECK( props.G0  == Approx(-5.279830e+05) );
    CHECK( props.H0  == Approx(-6.752359e+05) );
    CHECK( props.V0  == Approx(-6.063738e-06) );
    CHECK( props.Cp0 == Approx(-3.228612e+02) );
    CHECK( props.Cv0 == Approx(-3.228612e+02) );

    //-------------------------------------------------------------------
    // Testing attributes and thermodynamic properties of Ca+2
    //-------------------------------------------------------------------
    species = aq17.species().get("Ca+2");
    CHECK( species.formula().equivalent("Ca+2") );
    CHECK( species.substance() == "Ca+2 ion" );
    CHECK( species.aggregateState() == AggregateState::Aqueous );
    CHECK( species.charge() == +2 );
    CHECK( species.molarMass() == Approx(0.040076902) );

    props = species.props(T, P);
    CHECK( props.G0  == Approx(-5.528210e+05) );
    CHECK( props.H0  == Approx(-5.431003e+05) );
    CHECK( props.V0  == Approx(-1.844093e-05) );
    CHECK( props.Cp0 == Approx(-3.099935e+01) );
    CHECK( props.Cv0 == Approx(-3.099935e+01) );

    //-------------------------------------------------------------------
    // Testing attributes and thermodynamic properties of CO2
    //-------------------------------------------------------------------
    species = aq17.species().get("CO2");
    CHECK( species.formula().equivalent("CO2") );
    CHECK( species.substance() == "Carbon dioxide (CO2)" );
    CHECK( species.aggregateState() == AggregateState::Gas );
    CHECK( species.charge() == 0 );
    CHECK( species.molarMass() == Approx(0.0440096006) );

    props = species.props(T, P);
    CHECK( props.G0  == Approx(-3.943510e+05) );
    CHECK( props.H0  == Approx(-3.935472e+05) );
    CHECK( props.V0  == Approx( 2.466434e-02) );
    CHECK( props.Cp0 == Approx( 3.710812e+01) );
    CHECK( props.Cv0 == Approx( 0.000000e+00) );

    //-------------------------------------------------------------------
    // Testing attributes and thermodynamic properties of Calcite
    //-------------------------------------------------------------------
    species = aq17.species().get("Calcite");
    CHECK( species.formula().equivalent("CaCO3") );
    CHECK( species.substance() == "Calcite (cc)" );
    CHECK( species.aggregateState() == AggregateState::CrystallineSolid );
    CHECK( species.charge() == 0 );
    CHECK( species.molarMass() == Approx(0.1000869999) );

    props = species.props(T, P);
    CHECK( props.G0  == Approx(-1.129195e+06) );
    CHECK( props.H0  == Approx(-1.207470e+06) );
    CHECK( props.V0  == Approx( 3.689000e-05) );
    CHECK( props.Cp0 == Approx( 8.337073e+01) );
    CHECK( props.Cv0 == Approx( 0.000000e+00) );

    CHECK_THROWS( ThermoFunDatabase("not-a-valid-file-name") );
    CHECK_THROWS( ThermoFunDatabase::withName("not-a-valid-file-name") );
    CHECK_THROWS( ThermoFunDatabase::fromFile("not-a-valid-file-name") );
}
