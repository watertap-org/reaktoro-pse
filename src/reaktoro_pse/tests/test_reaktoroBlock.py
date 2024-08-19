import pytest
from reaktoro_pse.reaktoroBlock import reaktorBlock

from pyomo.environ import (
    ConcreteModel,
    Var,
    assert_optimal_termination,
    units as pyunits,
)
from watertap.core.solvers import get_solver


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
    m.outputs = Var([("scalingTendency", "Calcite"), ("pH", None)], initialize=1)
    return m


@pytest.fixture
def build_rkt_state_with_indexed_species():
    m = ConcreteModel()
    m.temp = Var([0, 1], initialize=293.15, units=pyunits.K)
    m.temp.fix()
    # temp.construct()
    m.pressure = Var([0, 1], initialize=1e5, units=pyunits.Pa)
    m.pressure.fix()
    # pressure.construct()
    m.pH = Var([0, 1], initialize=7, units=pyunits.dimensionless)
    m.pH.fix()
    # pH.construct()
    m.composition = Var(
        [0, 1],
        ["H2O", "Mg", "Na", "Cl", "Ca", "HCO3"],
        initialize=1,
        units=pyunits.mol / pyunits.s,
    )
    for idx in [0, 1]:
        m.composition[(idx, "H2O")].fix(50)
        m.composition[(idx, "Mg")].fix(0.1 * (1 + idx))
        m.composition[(idx, "Na")].fix(0.5 * (1 + idx))
        m.composition[(idx, "Cl")].fix(0.5 * (1 + idx))
        m.composition[(idx, "Ca")].fix(0.01 * (1 + idx))
        m.composition[(idx, "HCO3")].fix(0.01 * (1 + idx))

    m.outputs = Var(
        [0, 1], [("scalingTendency", "Calcite"), ("pH", None)], initialize=1
    )
    return m


def test_blockBuild(build_rkt_state_with_species):
    m = build_rkt_state_with_species
    m.outputs.display()
    m.reaktoroBlock = reaktorBlock(
        composition=m.composition,
        temperature=m.temp,
        pressure=m.pressure,
        activity_model="ActivityModelPitzer",
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        numerical_jac_type="average",
        numerical_jac_order=2,
        numerical_jac_step=1e-8,
        outputs=m.outputs,
        convert_to_rkt_species=True,
    )
    print("rkt block")
    m.reaktoroBlock.reaktoro_model.display()
    print("rkt block")
    m.reaktoroBlock.initialize()
    cy_solver = get_solver(solver="cyipopt-watertap")
    cy_solver.options["max_iter"] = 20
    m.pH.fix()
    m.composition["H2O"].unfix()
    m.composition["H2O"].setlb(30)
    m.outputs[("scalingTendency", "Calcite")].fix(5)
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.display()
    assert pytest.approx(m.composition["H2O"].value) == 1601.999


def test_blockBuild_with_speciation_block(build_rkt_state_with_species):
    m = build_rkt_state_with_species
    m.CaO = Var(["CaO"], initialize=0.001, units=pyunits.mol / pyunits.s)
    m.CaO.fix()
    m.outputs.display()
    m.reaktoroBlock = reaktorBlock(
        composition=m.composition,
        temperature=m.temp,
        pressure=m.pressure,
        chemical_addition=m.CaO,
        pH=m.pH,
        activity_model="ActivityModelPitzer",
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        numerical_jac_type="average",
        numerical_jac_order=2,
        numerical_jac_step=1e-8,
        outputs=m.outputs,
        convert_to_rkt_species=True,
        build_speciation_block=True,
    )
    m.reaktoroBlock.initialize()
    cy_solver = get_solver(solver="cyipopt-watertap")
    cy_solver.options["max_iter"] = 20
    m.pH.unfix()
    m.outputs[("scalingTendency", "Calcite")].fix(5)
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.display()
    assert pytest.approx(m.pH.value, 1e-2) == 6.401

    m.reaktoroBlock.display()
    m.reaktoroBlock.speciation_block.outputs.display()
    m.reaktoroBlock.speciation_block.reaktoro_model.display()
    m.reaktoroBlock.reaktoro_model.display()


def test_indexed_blockBuild(build_rkt_state_with_indexed_species):
    m = build_rkt_state_with_indexed_species
    m.outputs.display()
    m.reaktoroBlock = reaktorBlock(
        [0, 1],
        composition=m.composition,
        temperature=m.temp,
        pressure=m.pressure,
        pH=m.pH,
        activity_model="ActivityModelPitzer",
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        numerical_jac_type="average",
        numerical_jac_order=2,
        numerical_jac_step=1e-8,
        outputs=m.outputs,
        convert_to_rkt_species=True,
    )
    for blk in m.reaktoroBlock:
        m.reaktoroBlock[blk].initialize()
    m.reaktoroBlock[0].reaktoro_model.display()
    cy_solver = get_solver(solver="cyipopt-watertap")
    cy_solver.options["max_iter"] = 20
    m.pH.unfix()
    m.outputs[0, ("scalingTendency", "Calcite")].fix(5)
    m.outputs[1, ("scalingTendency", "Calcite")].fix(2.5)
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.display()
    assert pytest.approx(m.pH[0].value) == 6.782103
    assert pytest.approx(m.pH[1].value) == 6.161411138058621


def test_indexed_blockBuild_with_speciation_block(
    build_rkt_state_with_indexed_species,
):
    m = build_rkt_state_with_indexed_species
    m.CaO = Var([0, 1], ["CaO"], initialize=0.01, units=pyunits.mol / pyunits.s)
    m.CaO.fix()
    m.outputs.display()
    m.reaktoroBlock = reaktorBlock(
        [0, 1],
        composition=m.composition,
        temperature=m.temp,
        pressure=m.pressure,
        pH=m.pH,
        chemical_addition=m.CaO,
        activity_model="ActivityModelPitzer",
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        numerical_jac_type="average",
        numerical_jac_order=2,
        numerical_jac_step=1e-8,
        outputs=m.outputs,
        convert_to_rkt_species=True,
        build_speciation_block=True,
    )
    for blk in m.reaktoroBlock:
        m.reaktoroBlock[blk].initialize()
    m.reaktoroBlock.display()
    cy_solver = get_solver(solver="cyipopt-watertap")
    cy_solver.options["max_iter"] = 20
    m.CaO.unfix()
    m.outputs[(0, "pH", None)].fix(11.5)
    m.outputs[(1, "pH", None)].fix(10)
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.display()
    assert pytest.approx(m.CaO[(0, "CaO")].value) == 0.01732553618254949
    assert pytest.approx(m.CaO[(1, "CaO")].value) == 0.01205362158984656
