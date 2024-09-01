#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/reaktoro-pse/"
#################################################################################
import pytest
import reaktoro as rkt
import reaktoro_pse.core.reaktoro_state as rktState

from pyomo.environ import ConcreteModel, Var, units as pyunits

__author__ = "Alexander V. Dudchenko (SLAC)"


@pytest.fixture
def build_rkt_state_with_species():
    m = ConcreteModel()
    m.temp = Var(initialize=293.15, units=pyunits.K)
    m.temp.fix()
    # temp.construct()
    m.pressure = Var(initialize=1e5, units=pyunits.Pa)
    m.pressure.fix()
    # pressure.construct()
    m.pH = Var(initialize=7, units=pyunits.dimensionless)
    m.pH.fix()
    # pH.construct()
    m.composition = Var(
        ["H2O", "Mg", "Na", "Cl", "Ca", "HCO3", "CO2"],
        initialize=1,
        units=pyunits.mol / pyunits.s,
    )
    m.composition.construct()
    m.composition["H2O"].fix(50)
    m.composition["Mg"].fix(0.1)
    m.composition["Na"].fix(0.5)
    m.composition["Cl"].fix(0.5)
    m.composition["Ca"].fix(0.01)
    m.composition["HCO3"].fix(0.01)
    m.composition["CO2"].fix(0.001)

    rkt_state = rktState.ReaktoroState()
    rkt_state.set_database()
    rkt_state.set_input_options(
        "aqueous_phase",
    )
    rkt_state.register_system_inputs(temperature=m.temp, pressure=m.pressure)
    rkt_state.register_aqueous_inputs(composition=m.composition, pH=m.pH)
    rkt_state.set_aqueous_phase_activity_model()
    return m, rkt_state


@pytest.fixture
def build_rkt_state_with_species_mass_basis():
    m = ConcreteModel()
    m.temp = Var(initialize=293.15, units=pyunits.K)
    m.temp.fix()
    # temp.construct()
    m.pressure = Var(initialize=1e5, units=pyunits.Pa)
    m.pressure.fix()
    # pressure.construct()
    m.pH = Var(initialize=7, units=pyunits.dimensionless)
    m.pH.fix()
    # pH.construct()
    m.composition = Var(
        ["H2O", "Mg", "Na", "Cl", "Ca", "HCO3", "CO2"],
        initialize=1,
        units=pyunits.kg / pyunits.s,
    )
    m.composition.construct()
    m.composition["H2O"].fix(50 * 0.018016)
    m.composition["Mg"].fix(0.1 * 0.02430390284018224)
    m.composition["Na"].fix(0.5 * 0.022989251420091117)
    m.composition["Cl"].fix(0.5 * 0.035453548579908886)
    m.composition["Ca"].fix(0.01 * 0.040078902840182236)
    m.composition["HCO3"].fix(0.01 * 0.06001219715981776)
    m.composition["CO2"].fix(0.001 * 0.0440111)
    rkt_state = rktState.ReaktoroState()

    rkt_state.set_database()
    rkt_state.set_input_options(
        "aqueous_phase",
    )
    rkt_state.register_system_inputs(temperature=m.temp, pressure=m.pressure)
    rkt_state.register_aqueous_inputs(composition=m.composition, pH=m.pH)
    rkt_state.set_aqueous_phase_activity_model()
    return m, rkt_state


@pytest.fixture
def build_rkt_state_with_species_no_ph():
    m = ConcreteModel()
    m.temp = Var(initialize=293.15, units=pyunits.K)
    m.temp.fix()
    # temp.construct()
    m.pressure = Var(initialize=1e5, units=pyunits.Pa)
    m.pressure.fix()
    # pH.construct()
    m.composition = Var(
        ["H2O", "Mg", "Na", "Cl", "Ca", "HCO3"],
        initialize=1,
        units=pyunits.mol / pyunits.s,
    )
    m.composition.construct()
    m.composition["H2O"].fix(50)
    m.composition["Mg"].fix(0.1)
    m.composition["Na"].fix(0.5)
    m.composition["Cl"].fix(0.5)
    m.composition["Ca"].fix(0.01)
    m.composition["HCO3"].fix(0.01)
    rkt_state = rktState.ReaktoroState()

    rkt_state.set_database()
    rkt_state.set_input_options(
        "aqueous_phase",
    )
    rkt_state.register_system_inputs(temperature=m.temp, pressure=m.pressure)
    rkt_state.register_aqueous_inputs(composition=m.composition)
    rkt_state.set_aqueous_phase_activity_model()
    return m, rkt_state


@pytest.fixture
def build_rkt_state_with_elements():
    temp = Var(initialize=293.15, units=pyunits.K)
    temp.construct()
    pressure = Var(initialize=1e5, units=pyunits.Pa)
    pressure.construct()
    pH = Var(initialize=7, units=pyunits.dimensionless)
    pH.construct()
    composition = Var(
        ["H", "O", "Na", "Cl", "Ca", "C", "O"],
        initialize=1,
        units=pyunits.mol / pyunits.s,
    )
    composition.construct()
    composition["H"].fix(100.01)
    composition["O"].fix(50)
    composition["Na"].fix(0.5)
    composition["Cl"].fix(0.5)
    composition["Ca"].fix(0.01)
    composition["C"].fix(0.01)
    rkt_state = rktState.ReaktoroState()

    rkt_state.set_database()
    rkt_state.set_input_options(
        "aqueous_phase",
        convert_to_rkt_species=False,
        composition_is_elements=True,
    )
    rkt_state.register_system_inputs(temperature=temp, pressure=pressure)
    rkt_state.register_aqueous_inputs(composition=composition, pH=pH)
    rkt_state.set_aqueous_phase_activity_model()
    return rkt_state


def test_state_with_species(build_rkt_state_with_species):
    m, rkt_state = build_rkt_state_with_species
    rkt_state.build_state()
    rkt_state.equilibrate_state()
    print(rkt_state.state)
    assert (
        pytest.approx(float(rkt_state.state.props().elementAmount("O")), 1e-3) == 50.03
    )
    assert (
        pytest.approx(float(rkt_state.state.props().elementAmount("H")), 0.1) == 100.0
    )
    assert (
        pytest.approx(float(rkt_state.state.props().elementAmount("Na")), 1e-3) == 0.5
    )
    assert (
        pytest.approx(float(rkt_state.state.props().elementAmount("Mg")), 1e-3) == 0.1
    )
    assert (
        pytest.approx(float(rkt_state.state.props().elementAmount("Cl")), 1e-3) == 0.5
    )
    assert (
        pytest.approx(float(rkt_state.state.props().elementAmount("C")), 1e-3) == 0.011
    )


def test_state_with_species_mass_basis(build_rkt_state_with_species_mass_basis):
    m, rkt_state = build_rkt_state_with_species_mass_basis
    rkt_state.build_state()
    rkt_state.equilibrate_state()
    print(rkt_state.state)
    assert (
        pytest.approx(float(rkt_state.state.props().elementAmount("O")), 1e-3) == 50.03
    )
    assert (
        pytest.approx(float(rkt_state.state.props().elementAmount("H")), 0.1) == 100.0
    )
    assert (
        pytest.approx(float(rkt_state.state.props().elementAmount("Na")), 1e-3) == 0.5
    )
    assert (
        pytest.approx(float(rkt_state.state.props().elementAmount("Mg")), 1e-3) == 0.1
    )
    assert (
        pytest.approx(float(rkt_state.state.props().elementAmount("Cl")), 1e-3) == 0.5
    )
    assert (
        pytest.approx(float(rkt_state.state.props().elementAmount("C")), 1e-3) == 0.011
    )
    test_dict = {
        "H2O": 0.018016,
        "Mg+2": 0.02430390284018224,
        "Mg": 0.02430390284018224,
        "Na+": 0.022989251420091117,
        "Na": 0.022989251420091117,
        "Cl-": 0.035453548579908886,
        "Cl": 0.035453548579908886,
        "Ca+2": 0.040078902840182236,
        "Ca": 0.040078902840182236,
        "CO3-2": 0.06001219715981776,
        "HCO3": 0.06001219715981776,
        "CO2": 0.0440111,
        "temperature": None,
        "pressure": None,
        "pH": None,
    }
    for key, inp in rkt_state.inputs.items():
        # print(key, inp.get_unit_conversion_value())
        assert test_dict[key] == inp.get_unit_conversion_value()
    # print(test_dict)


def test_chian_activity(build_rkt_state_with_species):
    m, rkt_state = build_rkt_state_with_species

    rkt_state.set_aqueous_phase_activity_model(
        ["ActivityModelHKF", rkt.ActivityModelSetschenow("H2O", 0.123)]
    )

    rkt_state.build_state()
    rkt_state.equilibrate_state()
    props = rkt_state.state.props()
    for species in rkt_state.system.species():
        print(f"{species.name():<20}{props.speciesActivityCoefficient(species.name())}")
    assert (
        pytest.approx(float(props.speciesActivityCoefficient("H2O")), 1e-2) == 1.25736
    )


def test_state_with_solids(build_rkt_state_with_species):
    m, rkt_state = build_rkt_state_with_species
    rkt_state.register_mineral_phases("Calcite")
    rkt_state.set_mineral_phase_activity_model()
    rkt_state.build_state()
    rkt_state.equilibrate_state()
    print(rkt_state.state)
    assert (
        pytest.approx(float(rkt_state.state.props().speciesAmount("Calcite")), 1e-5)
        == 0.008835542753970915
    )


def test_state_with_gas(build_rkt_state_with_species):
    m, rkt_state = build_rkt_state_with_species
    rkt_state.register_gas_phase("CO2(g)")
    rkt_state.set_gas_phase_activity_model()
    rkt_state.build_state()
    rkt_state.equilibrate_state()
    print(rkt_state.state.props())
    assert (
        pytest.approx(float(rkt_state.state.props().speciesAmount("CO2(g)")), 1e-5)
        == 1e-16
    )
    assert (
        pytest.approx(float(rkt_state.state.props().speciesActivity("CO2(g)")), 1e-5)
        == 1.0
    )


def test_state_with_elements(build_rkt_state_with_elements):
    rkt_state = build_rkt_state_with_elements
    rkt_state.build_state()
    """ when user provides element amounts, they do not get set in itial condiitons"""
    assert pytest.approx(float(rkt_state.state.props().elementAmount("O")), 1e-3) == 0
    assert pytest.approx(float(rkt_state.state.props().elementAmount("H")), 0.1) == 0
    assert pytest.approx(float(rkt_state.state.props().elementAmount("Na")), 1e-3) == 0
    assert pytest.approx(float(rkt_state.state.props().elementAmount("Cl")), 1e-3) == 0
    assert pytest.approx(float(rkt_state.state.props().elementAmount("C")), 1e-3) == 0
