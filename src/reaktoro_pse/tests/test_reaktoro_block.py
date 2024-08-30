###############################################################################
# #################################################################################
# # WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# # through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# # National Renewable Energy Laboratory, and National Energy Technology
# # Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# # of Energy). All rights reserved.
# #
# # Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# # information, respectively. These files are also available online at the URL
# # "https://github.com/watertap-org/reaktoro-pse/"
# #################################################################################
###############################################################################
import pytest

from reaktoro_pse.reaktoro_block import ReaktoroBlock

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
    m.pressure = Var(initialize=1e5, units=pyunits.Pa)
    m.pressure.fix()
    m.pH = Var(initialize=7, units=pyunits.dimensionless)
    m.pH.fix()
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
    m.pressure = Var([0, 1], initialize=1e5, units=pyunits.Pa)
    m.pressure.fix()
    m.pH = Var([0, 1], initialize=7, units=pyunits.dimensionless)
    m.pH.fix()
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
    m.property_block = ReaktoroBlock(
        composition=m.composition,
        temperature=m.temp,
        pressure=m.pressure,
        pH=m.pH,
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        numerical_jac_type="average",
        numerical_jac_order=2,
        numerical_jac_step=1e-8,
        outputs=m.outputs,
        convert_to_rkt_species=True,
    )
    print("rkt block")
    m.property_block.reaktoro_model.display()

    print("rkt block")
    m.property_block.initialize()
    cy_solver = get_solver(solver="cyipopt-watertap")
    cy_solver.options["max_iter"] = 20
    m.pH.fix()
    m.composition["H2O"].unfix()
    m.composition["H2O"].setlb(30)
    m.outputs[("scalingTendency", "Calcite")].fix(5)
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.display()
    assert pytest.approx(m.composition["H2O"].value, 1e-3) == 68.0601837


def test_blockBuild_solids_gas(build_rkt_state_with_species):
    m = build_rkt_state_with_species
    m.outputs.display()
    m.solid_gas_outputs = Var(
        [
            ("speciesAmount", "Calcite"),
            ("vaporPressure", "H2O(g)"),
            # ("speciesActivityLn", "H2O(g)"),
        ],
        initialize=0.5,
    )
    m.property_block = ReaktoroBlock(
        composition=m.composition,
        temperature=m.temp,
        pressure=m.pressure,
        pH=m.pH,
        mineral_phases=["Calcite"],
        gas_phases=["H2O(g)"],
        gas_phase_activity_model="ActivityModelRedlichKwong",
        aqueous_phase_activity_model="ActivityModelPitzer",
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        numerical_jac_type="average",
        numerical_jac_order=2,
        numerical_jac_step=1e-8,
        outputs=m.solid_gas_outputs,
        convert_to_rkt_species=True,
    )
    m.display()
    m.property_block.initialize()
    cy_solver = get_solver(solver="cyipopt-watertap")
    cy_solver.options["max_iter"] = 20
    m.temp.fix(273.15 + 50)
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.display()
    assert (
        pytest.approx(m.solid_gas_outputs[("vaporPressure", "H2O(g)")].value, 1e-1)
        == 49382.90
    )
    # assert (
    #     pytest.approx(m.solid_gas_outputs[("speciesAmount", "Calcite")].value) == 0.0001
    # )


def test_blockBuild_with_speciation_block(build_rkt_state_with_species):
    m = build_rkt_state_with_species
    m.CaO = Var(["CaO"], initialize=0.001, units=pyunits.mol / pyunits.s)
    m.CaO.fix()
    m.outputs.display()
    m.property_block = ReaktoroBlock(
        composition=m.composition,
        temperature=m.temp,
        pressure=m.pressure,
        chemistry_modifier=m.CaO,
        pH=m.pH,
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        numerical_jac_type="average",
        numerical_jac_order=2,
        numerical_jac_step=1e-8,
        outputs=m.outputs,
        convert_to_rkt_species=True,
        build_speciation_block=True,
        # dissolve_species_in_reaktoro=False,
        # presolve_during_initialization=True,
        # presolve_tolerance=1e-16,
    )
    m.property_block.initialize()
    cy_solver = get_solver(solver="cyipopt-watertap")
    cy_solver.options["max_iter"] = 20
    m.pH.unfix()
    m.outputs[("scalingTendency", "Calcite")].fix(5)
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.display()
    assert pytest.approx(m.outputs[("pH", None)].value, 1e-2) == 6.7496301
    assert pytest.approx(m.pH.value, 1e-2) == 6.401

    m.property_block.display_jacobian_outputs()

    scaling_result = m.property_block.display_jacobian_scaling()
    expected_scaling = {
        "speciation_block": {
            ("speciesAmount", "H+"): 9.007999999999993e-08,
            ("speciesAmount", "H2O"): 50.0,
            ("speciesAmount", "CO3-2"): 3.2175702176273733e-06,
            ("speciesAmount", "CO2"): 0.00189035577659813,
            ("speciesAmount", "Ca+2"): 0.01,
            ("speciesAmount", "Cl-"): 0.7116050981506346,
            ("speciesAmount", "HCO3-"): 0.007825323588838813,
            ("speciesAmount", "Mg+2"): 0.09971792990850152,
            ("speciesAmount", "MgCO3"): 0.0002811030643454316,
            ("speciesAmount", "MgOH+"): 9.670271530541402e-07,
            ("speciesAmount", "Na+"): 0.5,
            ("speciesAmount", "OH-"): 6.004424745615723e-08,
        },
        "property_block": {
            ("saturationIndex", "Calcite"): 1.554873983061197,
            ("pH", None): 7.520409745594153,
        },
    }
    assert "speciation_block" in scaling_result
    assert "property_block" in scaling_result
    new_scaling = {}
    for key in scaling_result["speciation_block"]:
        new_scaling[key] = 1
        assert (
            pytest.approx(scaling_result["speciation_block"][key], 1e-3)
            == expected_scaling["speciation_block"][key]
        )
    m.property_block.set_jacobian_scaling(new_scaling, speciation_block=True)
    new_scaling = {}
    for key in scaling_result["property_block"]:
        new_scaling[key] = 1
        assert (
            pytest.approx(scaling_result["property_block"][key], 1e-3)
            == expected_scaling["property_block"][key]
        )
    m.property_block.set_jacobian_scaling(new_scaling, speciation_block=False)
    scaling_result = m.property_block.display_jacobian_scaling()
    assert "speciation_block" in scaling_result
    assert "property_block" in scaling_result
    for key in scaling_result["speciation_block"]:
        assert scaling_result["speciation_block"][key] == 1
    for key in scaling_result["property_block"]:
        assert scaling_result["property_block"][key] == 1


def test_blockBuild_with_speciation_block_no_chem_addition(
    build_rkt_state_with_species,
):
    m = build_rkt_state_with_species
    m.outputs.display()
    m.property_block = ReaktoroBlock(
        composition=m.composition,
        temperature=m.temp,
        pressure=m.pressure,
        pH=m.pH,
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        numerical_jac_type="average",
        numerical_jac_order=2,
        numerical_jac_step=1e-8,
        outputs=m.outputs,
        convert_to_rkt_species=True,
        build_speciation_block=True,
        presolve_during_initialization=True,
    )
    m.property_block.initialize()
    cy_solver = get_solver(solver="cyipopt-watertap")
    cy_solver.options["max_iter"] = 20
    m.pH.unfix()
    m.outputs[("scalingTendency", "Calcite")].fix(5)
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.display()
    assert pytest.approx(m.outputs[("pH", None)].value, 1e-2) == m.pH.value

    m.property_block.display()
    m.property_block.speciation_block.outputs.display()
    m.property_block.speciation_block.reaktoro_model.display()
    m.property_block.reaktoro_model.display()


def test_blockBuild_with_speciation_block_no_chem_super_critical_db(
    build_rkt_state_with_species,
):
    translation_dict = {
        "H2O": "H2O(aq)",
        "Mg": "Mg+2",
        "Na": "Na+",
        "Cl": "Cl-",
        "SO4": "SO4-2",
        "Ca": "Ca+2",
        "HCO3": "HCO3-",
    }
    m = build_rkt_state_with_species
    m.outputs.display()
    m.CaO = Var(["CaO"], initialize=0.002, units=pyunits.mol / pyunits.s)
    m.CaO.fix()
    m.property_block = ReaktoroBlock(
        composition=m.composition,
        temperature=m.temp,
        pressure=m.pressure,
        pH=m.pH,
        chemistry_modifier=m.CaO,
        database="SupcrtDatabase",
        database_file="supcrtbl",
        numerical_jac_type="average",
        numerical_jac_order=2,
        numerical_jac_step=1e-8,
        outputs=m.outputs,
        convert_to_rkt_species=True,
        species_to_rkt_species_dict=translation_dict,
        build_speciation_block=True,
        presolve_during_initialization=True,
    )
    for e, con in m.property_block.rkt_inputs.constraint_dict.items():
        print(e, con)
    m.property_block.initialize()

    m.display()
    cy_solver = get_solver(solver="cyipopt-watertap")
    cy_solver.options["max_iter"] = 20
    m.pH.unfix()
    m.outputs[("scalingTendency", "Calcite")].fix(5)
    result = cy_solver.solve(m, tee=True)

    m.display()
    assert_optimal_termination(result)
    assert pytest.approx(m.outputs[("pH", None)].value, 1e-2) == 6.899783669305352
    assert pytest.approx(m.pH.value, 1e-2) == 6.0677628977

    m.property_block.display()
    m.property_block.speciation_block.outputs.display()
    m.property_block.speciation_block.reaktoro_model.display()
    m.property_block.reaktoro_model.display()


def test_indexed_blockBuild(build_rkt_state_with_indexed_species):
    m = build_rkt_state_with_indexed_species
    m.outputs.display()
    m.property_block = ReaktoroBlock(
        [0, 1],
        composition=m.composition,
        temperature=m.temp,
        pressure=m.pressure,
        pH=m.pH,
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        numerical_jac_type="average",
        numerical_jac_order=2,
        numerical_jac_step=1e-8,
        outputs=m.outputs,
        convert_to_rkt_species=True,
    )
    for blk in m.property_block:
        m.property_block[blk].initialize()
    m.property_block[0].reaktoro_model.display()
    cy_solver = get_solver(solver="cyipopt-watertap")
    cy_solver.options["max_iter"] = 20
    m.pH.unfix()
    m.outputs[0, ("scalingTendency", "Calcite")].fix(5)
    m.outputs[1, ("scalingTendency", "Calcite")].fix(2.5)
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.display()
    assert pytest.approx(m.pH[0].value, 1e-3) == 6.78206
    assert pytest.approx(m.pH[1].value, 1e-3) == 6.161411138058621


def test_indexed_blockBuild_with_speciation_block(
    build_rkt_state_with_indexed_species,
):
    m = build_rkt_state_with_indexed_species
    m.CaO = Var([0, 1], ["CaO"], initialize=0.01, units=pyunits.mol / pyunits.s)
    m.CaO.fix()
    m.outputs.display()
    m.property_block = ReaktoroBlock(
        [0, 1],
        composition=m.composition,
        temperature=m.temp,
        pressure=m.pressure,
        pH=m.pH,
        chemistry_modifier=m.CaO,
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        numerical_jac_type="average",
        numerical_jac_order=2,
        numerical_jac_step=1e-8,
        outputs=m.outputs,
        convert_to_rkt_species=True,
        build_speciation_block=True,
    )
    for blk in m.property_block:
        m.property_block[blk].initialize()
    m.property_block.display()
    cy_solver = get_solver(solver="cyipopt-watertap")
    cy_solver.options["max_iter"] = 20
    m.CaO.unfix()
    m.outputs[(0, "pH", None)].fix(11.5)
    m.outputs[(1, "pH", None)].fix(10)
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.display()
    assert pytest.approx(m.CaO[(0, "CaO")].value, 1e-3) == 0.01732553618254949
    assert pytest.approx(m.CaO[(1, "CaO")].value, 1e-3) == 0.01205362158984656
