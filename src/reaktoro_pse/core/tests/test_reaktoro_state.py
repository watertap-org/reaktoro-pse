import pytest
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
    rkt_state = rktState.reaktoroState()
    rkt_state.register_inputs(
        composition=m.composition, temperature=m.temp, pressure=m.pressure, pH=m.pH
    )
    rkt_state.set_database()
    rkt_state.set_activity_model()
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
    rkt_state = rktState.reaktoroState()
    rkt_state.register_inputs(
        composition=m.composition, temperature=m.temp, pressure=m.pressure
    )
    rkt_state.set_database()
    rkt_state.set_activity_model()
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
    rkt_state = rktState.reaktoroState()
    rkt_state.register_inputs(
        composition=composition,
        temperature=temp,
        pressure=pressure,
        pH=pH,
        convert_to_rkt_species=False,
        composition_is_elements=True,
    )
    rkt_state.set_database()
    rkt_state.set_activity_model()
    return rkt_state


def test_state_with_species(build_rkt_state_with_species):
    m, rkt_state = build_rkt_state_with_species
    rkt_state.build_state()
    rkt_state.equilibrate_state()
    print(rkt_state.rktState)
    assert (
        pytest.approx(float(rkt_state.rktState.props().elementAmount("O")), 1e-3)
        == 50.03
    )
    assert (
        pytest.approx(float(rkt_state.rktState.props().elementAmount("H")), 0.1)
        == 100.0
    )
    assert (
        pytest.approx(float(rkt_state.rktState.props().elementAmount("Na")), 1e-3)
        == 0.5
    )
    assert (
        pytest.approx(float(rkt_state.rktState.props().elementAmount("Mg")), 1e-3)
        == 0.1
    )
    assert (
        pytest.approx(float(rkt_state.rktState.props().elementAmount("Cl")), 1e-3)
        == 0.5
    )
    assert (
        pytest.approx(float(rkt_state.rktState.props().elementAmount("C")), 1e-3)
        == 0.011
    )


def test_state_with_solids(build_rkt_state_with_species):
    m, rkt_state = build_rkt_state_with_species
    rkt_state.register_mineral_phases("Calcite")
    rkt_state.build_state()
    rkt_state.equilibrate_state()
    print(rkt_state.rktState)
    assert (
        pytest.approx(float(rkt_state.rktState.props().speciesAmount("Calcite")), 1e-5)
        == 0.00861693780618118
    )


def test_state_with_gas(build_rkt_state_with_species):
    m, rkt_state = build_rkt_state_with_species
    rkt_state.register_gas_phases("CO2")
    rkt_state.build_state()
    rkt_state.equilibrate_state()
    print(rkt_state.rktState)
    assert (
        pytest.approx(float(rkt_state.rktState.props().speciesAmount("CO2(g)")), 1e-5)
        == 1e-16
    )


def test_state_with_elements(build_rkt_state_with_elements):
    rkt_state = build_rkt_state_with_elements
    rkt_state.build_state()
    """ when user provides element amounts, they do not get set in itial condiitons"""
    assert (
        pytest.approx(float(rkt_state.rktState.props().elementAmount("O")), 1e-3) == 0
    )
    assert pytest.approx(float(rkt_state.rktState.props().elementAmount("H")), 0.1) == 0
    assert (
        pytest.approx(float(rkt_state.rktState.props().elementAmount("Na")), 1e-3) == 0
    )
    assert (
        pytest.approx(float(rkt_state.rktState.props().elementAmount("Cl")), 1e-3) == 0
    )
    assert (
        pytest.approx(float(rkt_state.rktState.props().elementAmount("C")), 1e-3) == 0
    )
