#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/reaktoro-pse/"
#################################################################################
import logging

import pytest

from reaktoro_pse.reaktoro_block import ReaktoroBlock

from pyomo.environ import (
    ConcreteModel,
    Block,
    Var,
    Constraint,
    assert_optimal_termination,
    value,
    units as pyunits,
)
import idaes.core.util.scaling as iscale

from pyomo.contrib.pynumero.interfaces.external_grey_box import (
    ExternalGreyBoxModel,
)
from reaktoro_pse.core.util_classes.cyipopt_solver import (
    get_cyipopt_watertap_solver,
)
from watertap_solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


def build_comp(blk):
    blk.temp = Var(initialize=293.15, units=pyunits.K)
    blk.temp.fix()
    blk.pressure = Var(initialize=1e5, units=pyunits.Pa)
    blk.pressure.fix()
    blk.pH = Var(initialize=7, units=pyunits.dimensionless)
    blk.pH.fix()
    blk.composition = Var(
        ["H2O", "Mg", "Na", "Cl", "Ca", "HCO3", "SO4"],
        initialize=1,
        units=pyunits.mol / pyunits.s,
    )
    blk.composition["H2O"].fix(50)
    blk.composition["Mg"].fix(0.1)
    blk.composition["Na"].fix(0.5)
    blk.composition["Cl"].fix(0.5)
    blk.composition["Ca"].fix(0.01)
    blk.composition["HCO3"].fix(0.01)
    blk.composition["SO4"].fix(0.01)
    blk.outputs = Var(
        [("scalingTendency", "Calcite"), ("pH", None), ("pE", None)],
        initialize=1,
    )


def add_simple_test_constraints(blk):
    blk.recovery = Var(initialize=0.9, bounds=(0, 1))
    blk.outlet_composition = Var(
        blk.composition.index_set(), initialize=1, units=pyunits.mol / pyunits.s
    )

    @blk.Constraint(blk.composition.index_set())
    def mass_balance_rule(b, i):
        return b.outlet_composition[i] == b.composition[i] * b.recovery


@pytest.fixture
def build_rkt_state_with_species():
    m = ConcreteModel()
    build_comp(m)
    return m


@pytest.fixture
def build_rkt_state_with_species_and_pE():
    m = ConcreteModel()
    build_comp(m)
    m.pE = Var(initialize=0, units=pyunits.dimensionless)
    m.pE.fix()
    return m


@pytest.fixture
def build_rkt_state_with_species_and_mixing():
    m = ConcreteModel()
    m.comp_a = Block()
    build_comp(m.comp_a)
    m.comp_b = Block()
    build_comp(m.comp_b)
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
        if idx == 0:
            m.composition[(idx, "H2O")].fix(50)
        else:
            m.composition[(idx, "H2O")].fix(20)
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
        aqueous_phase={
            "composition": m.composition,
            "convert_to_rkt_species": True,
        },
        system_state={
            "temperature": m.temp,
            "pressure": m.pressure,
            "pH": m.pH,
        },
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        outputs=m.outputs,
    )
    print("rkt block")
    m.property_block.reaktoro_model.display()
    print("rkt block")
    m.property_block.initialize()

    m.property_block.display_jacobian_scaling()
    cy_solver = get_cyipopt_watertap_solver()
    cy_solver.options["max_iter"] = 20
    m.pH.fix()
    m.composition["H2O"].unfix()
    m.composition["H2O"].setlb(30)
    m.outputs[("scalingTendency", "Calcite")].fix(5)
    m.property_block.output_constraints.pprint()
    print(degrees_of_freedom(m))
    assert degrees_of_freedom(m) == 0
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.display()
    assert pytest.approx(m.composition["H2O"].value, 1e-3) == 68.0601837


def test_blockBuild_with_specie_balance_warning(build_rkt_state_with_species, caplog):
    m = build_rkt_state_with_species
    m.outputs.display()
    with caplog.at_level(logging.WARNING):
        m.property_block = ReaktoroBlock(
            aqueous_phase={
                "composition": m.composition,
                "convert_to_rkt_species": True,
            },
            system_state={
                "temperature": m.temp,
                "pressure": m.pressure,
                "pH": m.pH,
            },
            database="PhreeqcDatabase",
            database_file="pitzer.dat",
            outputs=m.outputs,
            charge_neutrality_ion="SO4-2",
        )
    print(caplog.text)
    assert f"""The charge neutrality ion SO4-2 is not an element,
                    and inexact speciation is provided. Ignore this warning, if you want to only adjust specie ratios to
                    achieve charge neutrality, otherwise supply an element (such as "S" instead of "SO4-2")
                    to find the amount of element that should be added or removed to achieve charge neutrality. Adjustment
                    of specie ratios should in general be done when exact speciation is provided ('exact_speciation=True')
                    as otherwise Reaktoro solver might not converge, as shifting specie ratios might be insufficient to
                    achieve a charge neutral solution.""".replace(
        " ", ""
    ).replace(
        "\n", ""
    ) in caplog.text.replace(
        " ", ""
    ).replace(
        "\n", ""
    )


def test_blockBuild_with_pE(build_rkt_state_with_species_and_pE):
    m = build_rkt_state_with_species_and_pE
    m.outputs.display()
    m.property_block = ReaktoroBlock(
        aqueous_phase={
            "composition": m.composition,
            "convert_to_rkt_species": True,
        },
        system_state={
            "temperature": m.temp,
            "pressure": m.pressure,
            "pH": m.pH,
            "pE": m.pE,
        },
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        outputs=m.outputs,
    )
    print("rkt block")
    m.property_block.reaktoro_model.display()
    print("rkt block")
    m.property_block.initialize()

    m.property_block.display_jacobian_scaling()
    cy_solver = get_cyipopt_watertap_solver()
    cy_solver.options["max_iter"] = 20
    m.pH.fix()
    m.composition["H2O"].unfix()
    m.composition["H2O"].setlb(30)
    m.outputs[("scalingTendency", "Calcite")].fix(5)
    m.property_block.output_constraints.pprint()
    print(degrees_of_freedom(m))
    assert degrees_of_freedom(m) == 0
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.display()
    assert pytest.approx(m.composition["H2O"].value, 1e-3) == 68.0601837
    assert pytest.approx(m.outputs[("pE", None)].value, 1e-3) == -9.675807527465942


def build_O2_outputs(m):
    m.O2_outputs = Var(
        [("pH", None), ("pE", None), ("speciesAmount", "O2")],
        initialize=1,
    )
    m.dissolved_O2 = Var(initialize=1e-6, units=pyunits.mol / pyunits.s)
    m.dissolved_O2.fix()
    m.O2_fugacity = Var(initialize=0.21, units=pyunits.dimensionless)
    m.O2_fugacity.fix()


def test_blockBuild_with_dissolved_O2(build_rkt_state_with_species):
    """dissolved O2 defines redox state of block with open H and O"""
    m = build_rkt_state_with_species
    build_O2_outputs(m)
    m.property_block = ReaktoroBlock(
        aqueous_phase={
            "composition": m.composition,
            "convert_to_rkt_species": True,
        },
        system_state={
            "temperature": m.temp,
            "pressure": m.pressure,
            "pH": m.pH,
            "dissolved_O2": m.dissolved_O2,
        },
        database="PhreeqcDatabase",
        database_file="phreeqc.dat",
        outputs=m.O2_outputs,
    )
    assert (
        "dissolved_O2_constraint"
        in m.property_block.rkt_solver.input_specs.equilibrium_specs.namesConstraints()
    )
    m.property_block.initialize()
    cy_solver = get_cyipopt_watertap_solver()
    cy_solver.options["max_iter"] = 20
    assert degrees_of_freedom(m) == 0
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    assert pytest.approx(m.O2_outputs[("speciesAmount", "O2")].value, 1e-6) == 1e-6
    pE_low_O2 = m.O2_outputs[("pE", None)].value
    # more O2 is more oxidizing
    m.dissolved_O2.fix(1e-4)
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    assert pytest.approx(m.O2_outputs[("speciesAmount", "O2")].value, 1e-6) == 1e-4
    # O2 + 4H+ + 4e- = 2H2O, 100x O2 raises pE by log10(100)/4
    assert pytest.approx(m.O2_outputs[("pE", None)].value - pE_low_O2, 1e-3) == 0.5


@pytest.mark.parametrize(
    "concentration_input",
    ["dissolved_O2", "composition"],
)
def test_blockBuild_concentration_input_raises(
    build_rkt_state_with_species, concentration_input
):
    """concentrations (e.g. mg/L) would otherwise be silently treated as amounts"""
    m = build_rkt_state_with_species
    build_O2_outputs(m)
    m.DO_concentration = Var(initialize=8, units=pyunits.mg / pyunits.L)
    m.DO_concentration.fix()
    m.composition_concentration = Var(
        ["H2O", "Na", "Cl"], initialize=1, units=pyunits.mol / pyunits.L
    )
    m.composition_concentration.fix()
    if concentration_input == "dissolved_O2":
        composition = m.composition
        dissolved_O2 = m.DO_concentration
    else:
        composition = m.composition_concentration
        dissolved_O2 = m.dissolved_O2
    with pytest.raises(ValueError, match="is a concentration"):
        m.property_block = ReaktoroBlock(
            aqueous_phase={
                "composition": composition,
                "convert_to_rkt_species": True,
            },
            system_state={
                "temperature": m.temp,
                "pressure": m.pressure,
                "pH": m.pH,
                "dissolved_O2": dissolved_O2,
            },
            database="PhreeqcDatabase",
            database_file="phreeqc.dat",
            outputs=m.O2_outputs,
        )


def test_blockBuild_multiple_redox_specs_raise(build_rkt_state_with_species_and_pE):
    m = build_rkt_state_with_species_and_pE
    build_O2_outputs(m)
    with pytest.raises(ValueError, match="Only one redox specification"):
        m.property_block = ReaktoroBlock(
            aqueous_phase={
                "composition": m.composition,
                "convert_to_rkt_species": True,
            },
            system_state={
                "temperature": m.temp,
                "pressure": m.pressure,
                "pH": m.pH,
                "pE": m.pE,
                "dissolved_O2": m.dissolved_O2,
            },
            database="PhreeqcDatabase",
            database_file="phreeqc.dat",
            outputs=m.O2_outputs,
        )


@pytest.mark.parametrize(
    "coupling_type",
    [True, False],
)
def test_blockBuild_with_O2_dose_in_property_block(
    build_rkt_state_with_species, coupling_type
):
    """O2 added as chemistry modifier dissolves into water, dissolved O2 leaving
    block is found from equilibrium, and O2 dose can be solved for target dissolved O2
    """
    m = build_rkt_state_with_species
    build_O2_outputs(m)
    m.O2_dose = Var(initialize=1e-4, units=pyunits.mol / pyunits.s)
    m.O2_dose.fix()
    m.property_block = ReaktoroBlock(
        aqueous_phase={
            "composition": m.composition,
            "convert_to_rkt_species": True,
        },
        system_state={
            "temperature": m.temp,
            "pressure": m.pressure,
            "pH": m.pH,
            "dissolved_O2": m.dissolved_O2,
        },
        database="PhreeqcDatabase",
        database_file="phreeqc.dat",
        chemistry_modifier={"O2": m.O2_dose},
        outputs=m.O2_outputs,
        direct_speciation_to_property_block_coupling=coupling_type,
        build_speciation_block=True,
    )
    m.property_block.initialize()
    cy_solver = get_cyipopt_watertap_solver()
    cy_solver.options["max_iter"] = 20
    assert degrees_of_freedom(m) == 0
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    # nothing in water consumes O2, so dose adds to inlet dissolved O2
    assert (
        pytest.approx(m.O2_outputs[("speciesAmount", "O2")].value, 1e-4) == 1e-6 + 1e-4
    )
    # solve for O2 dose that gives target dissolved O2
    m.O2_dose.unfix()
    m.O2_outputs[("speciesAmount", "O2")].fix(5e-5)
    assert degrees_of_freedom(m) == 0
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    assert pytest.approx(m.O2_dose.value, 1e-4) == 5e-5 - 1e-6


@pytest.mark.parametrize(
    "coupling_type",
    [True, False],
)
def test_blockBuild_with_aeration_in_property_block(
    build_rkt_state_with_species, coupling_type
):
    """water with low dissolved O2 is aerated in property block (open to O2 gas)"""
    m = build_rkt_state_with_species
    build_O2_outputs(m)
    m.dissolved_O2.fix(1e-8)
    m.CaO = Var(["CaO"], initialize=0.001, units=pyunits.mol / pyunits.s)
    m.CaO.fix()
    m.property_block = ReaktoroBlock(
        aqueous_phase={
            "composition": m.composition,
            "convert_to_rkt_species": True,
        },
        system_state={
            "temperature": m.temp,
            "pressure": m.pressure,
            "pH": m.pH,
            "dissolved_O2": m.dissolved_O2,
        },
        system_state_modifier={"O2_fugacity": m.O2_fugacity},
        database="PhreeqcDatabase",
        database_file="phreeqc.dat",
        chemistry_modifier=m.CaO,
        outputs=m.O2_outputs,
        direct_speciation_to_property_block_coupling=coupling_type,
        build_speciation_block=True,
    )
    assert (
        "f[O2]"
        in m.property_block.rkt_solver.input_specs.equilibrium_specs.namesConstraints()
    )
    m.property_block.initialize()
    cy_solver = get_cyipopt_watertap_solver()
    cy_solver.options["max_iter"] = 20
    assert degrees_of_freedom(m) == 0
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    o2_saturated = m.O2_outputs[("speciesAmount", "O2")].value
    # ~50 mol of water (0.9 kg) at 20 C in equilibrium with air holds ~2.5e-4 mol O2
    assert 2e-4 < o2_saturated < 3e-4
    # inlet dissolved O2 does not change the aerated state
    m.dissolved_O2.fix(1e-6)
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    assert (
        pytest.approx(m.O2_outputs[("speciesAmount", "O2")].value, 1e-4) == o2_saturated
    )
    # lower O2 fugacity (less O2 in gas) lowers dissolved O2
    m.O2_fugacity.fix(0.021)
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    assert (
        pytest.approx(m.O2_outputs[("speciesAmount", "O2")].value, 1e-2)
        == o2_saturated / 10
    )


@pytest.mark.parametrize(
    "scaling_type",
    [
        "no_scaling",
        "variable_output_scaling",
        "jacobian_matrix_square_sum",
        "jacobian_matrix_inverse_sum",
        "variable_oi_scaling_square_sum",
        "variable_oi_scaling_inverse_sum",
    ],
)
def test_block_jacobian_scaling(build_rkt_state_with_species, scaling_type):
    m = build_rkt_state_with_species
    m.property_block = ReaktoroBlock(
        aqueous_phase={
            "composition": m.composition,
            "convert_to_rkt_species": True,
        },
        system_state={
            "temperature": m.temp,
            "pressure": m.pressure,
            "pH": m.pH,
        },
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        outputs=m.outputs,
        jacobian_options={"scaling_type": scaling_type},
    )
    m.property_block.initialize()
    scaling_factors = m.property_block.display_jacobian_scaling()

    def dict_test(dict_a, dict_b, tol=1e-4):
        for blk in dict_a.keys():
            for key in dict_a[blk].keys():
                assert key in dict_b[blk]
                assert pytest.approx(dict_a[blk][key], tol) == dict_b[blk][key]

    if scaling_type == "no_scaling":
        dict_test(
            scaling_factors,
            {
                "property_block": {
                    ("scalingTendency", "Calcite"): 1.0,
                    ("pH", None): 1.0,
                    ("pE", None): 1,
                }
            },
        )
    elif scaling_type == "variable_output_scaling":
        dict_test(
            scaling_factors,
            {
                "property_block": {
                    ("scalingTendency", "Calcite"): 0.1088858058987964,
                    ("pH", None): 0.14285714285714285,
                    ("pE", None): 0.10361133064195724,
                }
            },
        )
    elif scaling_type == "jacobian_matrix_square_sum":
        dict_test(
            scaling_factors,
            {
                "property_block": {
                    ("scalingTendency", "Calcite"): 0.0007698149159929651,
                    ("pH", None): 0.9999999999999998,
                    ("pE", None): 0.1359733544810525,
                }
            },
        )
    elif scaling_type == "jacobian_matrix_inverse_sum":
        dict_test(
            scaling_factors,
            {
                "property_block": {
                    ("scalingTendency", "Calcite"): 9.627321324829857e-08,
                    ("pH", None): 1e-08,
                    ("pE", None): 1e-08,
                }
            },
        )
    elif scaling_type == "variable_oi_scaling_square_sum":
        dict_test(
            scaling_factors,
            {
                "property_block": {
                    ("scalingTendency", "Calcite"): 0.0006275654371006124,
                    ("pH", None): 0.0008233598912186447,
                    ("pE", None): 0.0005971658974846666,
                }
            },
        )

    elif scaling_type == "variable_oi_scaling_inverse_sum":
        dict_test(
            scaling_factors,
            {
                "property_block": {
                    ("scalingTendency", "Calcite"): 100.0,
                    ("pH", None): 100.0,
                    ("pE", None): 100,
                }
            },
        )


def test_activate_deactivate(build_rkt_state_with_species):
    m = build_rkt_state_with_species
    m.property_block = ReaktoroBlock(
        aqueous_phase={
            "composition": m.composition,
            "convert_to_rkt_species": True,
        },
        system_state={
            "temperature": m.temp,
            "pressure": m.pressure,
            "pH": m.pH,
        },
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        outputs=m.outputs,
    )
    m.property_block.initialize()

    m.pH.fix()
    m.composition["H2O"].unfix()
    m.composition["H2O"].setlb(30)
    m.outputs[("scalingTendency", "Calcite")].fix(5)
    m.property_block.deactivate()

    assert m.property_block.active == False
    for v in m.property_block.component_data_objects(Constraint):
        print(v.name)
        assert v.active == False
    for v in m.property_block.component_data_objects(ExternalGreyBoxModel):
        assert v.active == False

    cy_solver = get_cyipopt_watertap_solver()
    cy_solver.options["max_iter"] = 20

    # this should fail run solve, raising value error
    with pytest.raises(ValueError):
        result = cy_solver.solve(m, tee=True)
    add_simple_test_constraints(m)
    m.composition["H2O"].fix()
    water_solver = get_solver()
    results = water_solver.solve(m, tee=True)
    assert_optimal_termination(results)
    assert pytest.approx(m.composition["H2O"].value, 1e-3) == 50
    m.property_block.activate()
    assert m.property_block.active == True
    for v in m.property_block.component_data_objects(Constraint):
        print(v.name)
        assert v.active == True
    for v in m.property_block.component_data_objects(ExternalGreyBoxModel):
        assert v.active == True

    # this solve should solve
    m.composition["H2O"].unfix()
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
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
        aqueous_phase={
            "composition": m.composition,
            "convert_to_rkt_species": True,
            "activity_model": "ActivityModelPitzer",
        },
        system_state={
            "temperature": m.temp,
            "pressure": m.pressure,
            "pH": m.pH,
        },
        mineral_phase={"phase_components": "Calcite"},
        gas_phase={
            "phase_components": ["H2O(g)"],
            "activity_model": "ActivityModelRedlichKwong",
        },
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        outputs=m.solid_gas_outputs,
    )
    m.display()
    m.property_block.initialize()
    cy_solver = get_cyipopt_watertap_solver()
    cy_solver.options["max_iter"] = 20
    m.temp.fix(273.15 + 50)
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.display()
    assert (
        pytest.approx(m.solid_gas_outputs[("vaporPressure", "H2O(g)")].value, 1e-1)
        == 49382.90
    )


def _offset_value(rkt_value, a, offset_type):
    if offset_type == "smooth_floor":
        return (rkt_value**2 + a**2) ** 0.5
    return rkt_value + a


@pytest.mark.parametrize(
    "offset, offset_type",
    [(1e-8, "additive"), (1e-4, "smooth_floor"), (0.0, "additive"), (1e-4, None)],
)
def test_species_amount_zero_offset(build_rkt_state_with_species, offset, offset_type):
    """species that are not present (Brucite at pH 7) are returned by Reaktoro as epsilon,
    the offset keeps the lower bounded output var away from 0, other outputs are untouched
    """
    m = build_rkt_state_with_species
    m.mineral_outputs = Var(
        [("speciesAmount", "Calcite"), ("speciesAmount", "Brucite"), ("pH", None)],
        initialize=1e-3,
    )
    sf = 1e3
    for key in [("speciesAmount", "Calcite"), ("speciesAmount", "Brucite")]:
        m.mineral_outputs[key].setlb(0)
        iscale.set_scaling_factor(m.mineral_outputs[key], sf)
    m.property_block = ReaktoroBlock(
        aqueous_phase={
            "composition": m.composition,
            "convert_to_rkt_species": True,
            "activity_model": "ActivityModelPitzer",
        },
        system_state={
            "temperature": m.temp,
            "pressure": m.pressure,
            "pH": m.pH,
        },
        mineral_phase={"phase_components": ["Calcite", "Brucite"]},
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        outputs=m.mineral_outputs,
        reaktoro_solve_options={
            "species_amount_zero_offset": offset,
            "species_amount_offset_type": offset_type,
        },
    )
    m.property_block.initialize()
    builder = m.property_block.rkt_block_builder
    offset_index = builder.species_amount_offset_index
    gray_box = m.property_block.reaktoro_model.get_external_model()
    # gray box outputs (floored) and raw reaktoro values
    rkt_out = m.property_block.reaktoro_model.outputs
    raw = m.property_block.rkt_solver.output_specs.user_outputs
    disabled = offset == 0 or offset_type is None
    if disabled:
        assert offset_index == {}
        assert gray_box._output_floors == {}
    else:
        # only speciesAmount outputs get an offset
        assert set(offset_index.keys()) == {
            ("speciesAmount", "Calcite"),
            ("speciesAmount", "Brucite"),
        }
    # offset is applied in the gray box, output constraints stay linear
    for con in m.property_block.output_constraints.values():
        assert con.body.polynomial_degree() == 1
    a = 0.0 if disabled else offset / sf
    for key in offset_index:
        assert pytest.approx(builder.get_species_amount_offset(key), rel=1e-12) == a
        # initialization is consistent with the gray box floor
        expected = _offset_value(raw[key].value, a, offset_type)
        assert (
            pytest.approx(m.mineral_outputs[key].value, rel=1e-12, abs=1e-30)
            == expected
        )
        assert pytest.approx(rkt_out[key].value, rel=1e-12, abs=1e-30) == expected
    cy_solver = get_cyipopt_watertap_solver()
    cy_solver.options["max_iter"] = 20
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    brucite = ("speciesAmount", "Brucite")
    # reaktoro returns epsilon for brucite, the output var sits offset away from its bound
    assert raw[brucite].value < 1e-30
    # (abs tolerance covers ipopt constraint violation)
    assert pytest.approx(m.mineral_outputs[brucite].value, rel=1e-3, abs=1e-13) == a
    for key in [("speciesAmount", "Calcite"), ("pH", None)]:
        expected = raw[key].value
        if key in offset_index:
            expected = _offset_value(expected, a, offset_type)
        assert (
            pytest.approx(m.mineral_outputs[key].value, rel=1e-9, abs=1e-30) == expected
        )
    if offset_type == "smooth_floor":
        # present species are shifted by ~a**2/(2n), far less than a
        calcite = ("speciesAmount", "Calcite")
        shift = m.mineral_outputs[calcite].value - raw[calcite].value
        assert 0 < shift < 1e-3 * a
    # report compares reaktoro amounts (from the gray box) with floored outputs
    report = m.property_block.get_species_amount_offset_report()
    status = m.property_block.check_species_amount_offsets()
    assert set(report) == set(offset_index)
    if disabled:
        assert status == {"species": {}, "elements": {}}
        return
    calcite = ("speciesAmount", "Calcite")
    for key, r in report.items():
        assert r["reaktoro_value"] == pytest.approx(raw[key].value, rel=1e-12)
        assert r["offset"] == pytest.approx(a, rel=1e-12)
        assert r["floored_value"] == pytest.approx(
            _offset_value(raw[key].value, a, offset_type), rel=1e-12
        )
        assert r["output_var_value"] == m.mineral_outputs[key].value
    assert status["species"] == {brucite: "absent", calcite: "ok"}
    assert status["elements"] == {e: "ok" for e in ["Mg", "O", "H", "Ca", "C"]}
    # calcite is shifted by more than a tiny tolerance
    calcite_rel_error = report[calcite]["relative_error"]
    assert calcite_rel_error > 0
    with pytest.raises(RuntimeError, match="Calcite"):
        m.property_block.check_species_amount_offsets(
            relative_tolerance=calcite_rel_error / 2, raise_error=True
        )
    # offset amount is compared with the system amount of each element of the species
    props = m.property_block.rkt_solver.state.state.props()
    mg_total = float(props.elementAmount("Mg"))
    brucite_elements = report[brucite]["elements"]
    assert set(brucite_elements) == {"Mg", "O", "H"}
    assert brucite_elements["Mg"]["coefficient"] == 1
    assert brucite_elements["Mg"]["system_amount"] == pytest.approx(mg_total)
    assert set(report[calcite]["elements"]) == {"Ca", "C", "O"}
    # per element, offsets of all species with the element are summed (O: 3 in
    # calcite, 2 in brucite)
    elements = m.property_block.get_species_amount_offset_element_report()
    assert elements["Mg"]["offset_fraction"] == pytest.approx(
        report[brucite]["error"] / mg_total, rel=1e-9
    )
    o_amount = 3 * report[calcite]["error"] + 2 * report[brucite]["error"]
    assert elements["O"]["offset_amount"] == pytest.approx(o_amount, rel=1e-12)
    assert elements["O"]["offset_fraction"] == pytest.approx(
        o_amount / float(props.elementAmount("O")), rel=1e-9
    )
    assert set(elements["O"]["species"]) == {brucite, calcite}
    # an element that the offsets take a noticeable part of is flagged
    mg_fraction = elements["Mg"]["offset_fraction"]
    assert 0 < mg_fraction < 1e-2
    status = m.property_block.check_species_amount_offsets(
        element_tolerance=mg_fraction / 2
    )
    assert status["elements"]["Mg"] == "depleted"
    with pytest.raises(RuntimeError, match="Mg.*from Brucite"):
        m.property_block.check_species_amount_offsets(
            element_tolerance=mg_fraction / 2, raise_error=True
        )
    # value based rescaling keeps the user scaling factor of offset outputs
    m.property_block.update_block_scaling()
    assert iscale.get_scaling_factor(m.mineral_outputs[brucite]) == sf
    # user rescaling after initialization changes the offset on the next solve
    iscale.set_scaling_factor(m.mineral_outputs[brucite], 10 * sf)
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    assert (
        pytest.approx(m.mineral_outputs[brucite].value, rel=1e-3, abs=1e-14) == a / 10
    )


def test_species_amount_offset_value_based_scaling(build_rkt_state_with_species):
    """offset outputs without a user scaling factor are scaled from their values, and
    rescaling (update_block_scaling) keeps following values, while a user scaling
    factor is kept"""
    m = build_rkt_state_with_species
    calcite = ("speciesAmount", "Calcite")
    brucite = ("speciesAmount", "Brucite")
    m.mineral_outputs = Var([calcite, brucite], initialize=1e-3, bounds=(0, None))
    sf = 1e3
    # only brucite has a user scaling factor
    iscale.set_scaling_factor(m.mineral_outputs[brucite], sf)
    m.property_block = ReaktoroBlock(
        aqueous_phase={
            "composition": m.composition,
            "convert_to_rkt_species": True,
            "activity_model": "ActivityModelPitzer",
        },
        system_state={
            "temperature": m.temp,
            "pressure": m.pressure,
            "pH": m.pH,
        },
        mineral_phase={"phase_components": ["Calcite", "Brucite"]},
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        outputs=m.mineral_outputs,
    )
    m.property_block.initialize()
    raw = m.property_block.rkt_solver.output_specs.user_outputs
    calcite_sf = iscale.get_scaling_factor(m.mineral_outputs[calcite])
    assert pytest.approx(calcite_sf, rel=1e-12) == 1 / raw[calcite].value
    # chemistry changes, rescaling follows the new value for calcite only
    raw[calcite].value *= 10
    m.property_block.update_block_scaling()
    assert (
        pytest.approx(iscale.get_scaling_factor(m.mineral_outputs[calcite]), rel=1e-12)
        == calcite_sf / 10
    )
    assert iscale.get_scaling_factor(m.mineral_outputs[brucite]) == sf
    # a value based scaling factor of an absent species (1/epsilon, capped at 1e32) puts
    # its scaled amount at ~1e-8, it is still classified as absent, not distorted
    iscale.set_scaling_factor(m.mineral_outputs[brucite], 1e32)
    report = m.property_block.get_species_amount_offset_report()
    assert report[brucite]["reaktoro_value"] < 1e-30
    status = m.property_block.check_species_amount_offsets(
        absent_scaled_amount=report[brucite]["scaled_amount"] / 2, raise_error=True
    )
    assert status["species"][brucite] == "absent"


@pytest.mark.parametrize("coupling_type", [True, False])
def test_species_amount_offset_not_applied_to_speciation(
    build_rkt_state_with_species, coupling_type
):
    """species amounts passed from speciation block to property block are never offset,
    only user speciesAmount outputs of the property block are"""
    m = build_rkt_state_with_species
    m.CaO = Var(["CaO"], initialize=0.001, units=pyunits.mol / pyunits.s)
    m.CaO.fix()
    m.mineral_outputs = Var([("speciesAmount", "Ca+2"), ("pH", None)], initialize=1e-3)
    iscale.set_scaling_factor(m.mineral_outputs[("speciesAmount", "Ca+2")], 1e3)
    m.property_block = ReaktoroBlock(
        aqueous_phase={
            "composition": m.composition,
            "convert_to_rkt_species": True,
        },
        system_state={
            "temperature": m.temp,
            "pressure": m.pressure,
            "pH": m.pH,
        },
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        chemistry_modifier=m.CaO,
        outputs=m.mineral_outputs,
        direct_speciation_to_property_block_coupling=coupling_type,
        build_speciation_block=True,
        reaktoro_solve_options={
            "species_amount_zero_offset": 1e-4,
            "species_amount_offset_type": "smooth_floor",
        },
    )
    m.property_block.initialize()
    builder = m.property_block.rkt_block_builder
    assert set(builder.species_amount_offset_index) == {("speciesAmount", "Ca+2")}
    # terminal property block output is floored in the gray box
    floors = m.property_block.reaktoro_model.get_external_model()._output_floors
    assert len(floors) == 1
    if coupling_type:
        # speciation outputs stay inside the coupled solver, never in pyomo
        assert m.property_block.coupled_solver.outputs
        assert not hasattr(m.property_block.speciation_block, "rkt_block_builder")
    else:
        spec_block = m.property_block.speciation_block
        assert spec_block.rkt_block_builder.species_amount_offset_index == {}
        assert spec_block.reaktoro_model.get_external_model()._output_floors == {}
        # property block inputs equal speciation outputs (an offset would be 1e-7 here)
        for key in spec_block.reaktoro_model.outputs:
            assert value(spec_block.outputs[key]) == pytest.approx(
                value(spec_block.reaktoro_model.outputs[key]), rel=1e-6, abs=1e-30
            )


@pytest.mark.parametrize(
    "coupling_type",
    [True, False],
)
def test_blockBuild_with_speciation_block(build_rkt_state_with_species, coupling_type):
    m = build_rkt_state_with_species
    m.CaO = Var(["CaO"], initialize=0.001, units=pyunits.mol / pyunits.s)
    m.CaO.fix()

    m.property_block = ReaktoroBlock(
        aqueous_phase={
            "composition": m.composition,
            "convert_to_rkt_species": True,
        },
        system_state={
            "temperature": m.temp,
            "pressure": m.pressure,
            "pH": m.pH,
        },
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        chemistry_modifier=m.CaO,
        outputs=m.outputs,
        direct_speciation_to_property_block_coupling=coupling_type,
        build_speciation_block=True,
    )
    m.property_block.initialize()
    cy_solver = get_cyipopt_watertap_solver(limited_memory=False)
    cy_solver.options["max_iter"] = 50
    m.pH.unfix()
    new_scaling = m.property_block.display_jacobian_scaling()
    m.property_block.display_jacobian_scaling()
    m.outputs[("scalingTendency", "Calcite")].fix(5)
    print(degrees_of_freedom(m))
    assert degrees_of_freedom(m) == 0
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)

    assert pytest.approx(m.outputs[("pH", None)].value, 1e-3) == 6.7496301
    assert pytest.approx(m.pH.value, 1e-3) == 6.401


def test_coupled_solver_solves_without_initialize(build_rkt_state_with_species):
    """coupled solver must work when initialize was skipped (e.g. state loaded from file)"""
    m = build_rkt_state_with_species
    m.CaO = Var(["CaO"], initialize=0.001, units=pyunits.mol / pyunits.s)
    m.CaO.fix()
    m.property_block = ReaktoroBlock(
        aqueous_phase={
            "composition": m.composition,
            "convert_to_rkt_species": True,
        },
        system_state={
            "temperature": m.temp,
            "pressure": m.pressure,
            "pH": m.pH,
        },
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        chemistry_modifier=m.CaO,
        outputs=m.outputs,
        direct_speciation_to_property_block_coupling=True,
        build_speciation_block=True,
    )
    coupled_solver = m.property_block.coupled_solver
    assert coupled_solver.prop_block_not_equilibrated
    jac, out = coupled_solver.solve_reaktoro_block()
    assert not coupled_solver.prop_block_not_equilibrated
    assert len(out) > 0


def test_blockBuild_with_speciation_and_mixing(
    build_rkt_state_with_species_and_mixing,
):
    m = build_rkt_state_with_species_and_mixing
    m.comp_a.CaO = Var(["CaO"], initialize=0.001, units=pyunits.mol / pyunits.s)
    m.comp_a.CaO.fix()
    m.comp_a.outputs.display()
    m.mixing_prop = ReaktoroBlock(
        aqueous_phase={
            "composition": m.comp_b.composition,
            "convert_to_rkt_species": True,
        },
        system_state={
            "temperature": m.comp_b.temp,
            "pressure": m.comp_b.pressure,
            "pH": m.comp_b.pH,
        },
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        outputs={"speciesAmount": True},
        build_graybox_model=False,
    )
    m.property_block = ReaktoroBlock(
        aqueous_phase={
            "composition": m.comp_a.composition,
            "convert_to_rkt_species": True,
        },
        system_state={
            "temperature": m.comp_a.temp,
            "pressure": m.comp_a.pressure,
            "pH": m.comp_a.pH,
        },
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        chemistry_modifier=m.comp_a.CaO,
        outputs=m.comp_a.outputs,
        build_speciation_block=True,
        external_speciation_reaktoro_blocks=[m.mixing_prop],
    )
    m.property_block.display()
    m.property_block.initialize()
    cy_solver = get_cyipopt_watertap_solver()
    cy_solver.options["max_iter"] = 100
    m.comp_a.pH.fix()
    m.comp_b.pH.unfix()
    m.comp_a.outputs[("scalingTendency", "Calcite")].fix(5)

    result = cy_solver.solve(m, tee=True)
    m.property_block.display_jacobian_scaling()
    m.display()
    assert_optimal_termination(result)
    assert (
        pytest.approx(m.comp_a.outputs[("pH", None)].value, 1e-3) == 6.765353236871279
    )
    assert pytest.approx(m.comp_b.pH.value, 1e-3) == 6.222585017249812


def test_blockBuild_with_speciation_block_no_chem_addition(
    build_rkt_state_with_species,
):
    m = build_rkt_state_with_species
    m.outputs.display()
    m.property_block = ReaktoroBlock(
        aqueous_phase={
            "composition": m.composition,
            "convert_to_rkt_species": True,
        },
        system_state={
            "temperature": m.temp,
            "pressure": m.pressure,
            "pH": m.pH,
        },
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        outputs=m.outputs,
        build_speciation_block=True,
    )
    m.property_block.initialize()
    cy_solver = get_cyipopt_watertap_solver()
    cy_solver.options["max_iter"] = 20
    m.pH.unfix()
    m.outputs[("scalingTendency", "Calcite")].fix(5)
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.display()
    assert pytest.approx(m.outputs[("pH", None)].value, 1e-2) == m.pH.value


def test_blockBuild_with_temp_and_pressure_modification_in_speciation_block(
    build_rkt_state_with_species,
):
    m = build_rkt_state_with_species
    m.CaO = Var(["CaO"], initialize=0.001, units=pyunits.mol / pyunits.s)
    m.CaO.fix()
    m.outputs.display()
    m.temp_mod = Var(initialize=333.15, units=pyunits.K)
    m.temp_mod.fix()
    m.pressure_mod = Var(initialize=5e5, units=pyunits.Pa)
    m.pressure_mod.fix()

    m.outputs_mod = Var(
        [
            ("scalingTendency", "Calcite"),
            ("pH", None),
            ("temperature", None),
            ("pressure", None),
        ],
        initialize=1,
    )
    m.property_block = ReaktoroBlock(
        aqueous_phase={
            "composition": m.composition,
            "convert_to_rkt_species": True,
        },
        system_state={
            "temperature": m.temp,
            "pressure": m.pressure,
            "pH": m.pH,
        },
        system_state_modifier={
            "temperature": m.temp_mod,
            "pressure": m.pressure_mod,
        },
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        chemistry_modifier=m.CaO,
        outputs=m.outputs_mod,
        build_speciation_block=True,
    )
    m.property_block.initialize()
    cy_solver = get_cyipopt_watertap_solver()
    cy_solver.options["max_iter"] = 20
    m.pH.unfix()
    m.outputs_mod[("scalingTendency", "Calcite")].fix(5)
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.display()
    assert pytest.approx(m.outputs_mod[("pH", None)].value, 1e-3) == 6.250981308052
    assert pytest.approx(m.pH.value, 1e-3) == 5.995934005877454
    assert pytest.approx(m.outputs_mod[("temperature", None)].value, 1e-3) == 333.15
    assert pytest.approx(m.outputs_mod[("pressure", None)].value, 1e-3) == 5e5


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
        aqueous_phase={
            "composition": m.composition,
            "convert_to_rkt_species": True,
            "species_to_rkt_species_dict": translation_dict,
        },
        system_state={
            "temperature": m.temp,
            "pressure": m.pressure,
            "pH": m.pH,
        },
        chemistry_modifier=m.CaO,
        database="SupcrtDatabase",
        database_file="supcrtbl",
        reaktoro_solve_options={"solver_tolerance": 1e-8, "epsilon": 1e-64},
        outputs=m.outputs,
        build_speciation_block=True,
    )

    m.property_block.initialize()
    m.property_block.display_jacobian_outputs()
    cy_solver = get_cyipopt_watertap_solver(limited_memory=False)

    cy_solver.options["max_iter"] = 50
    result = cy_solver.solve(m, tee=True)
    m.pH.unfix()
    m.outputs[("scalingTendency", "Calcite")].fix(5)
    result = cy_solver.solve(m, tee=True)

    m.display()
    assert_optimal_termination(result)
    m.property_block.display()
    m.property_block.display_jacobian_scaling()
    m.outputs.display()
    m.pH.display()
    assert pytest.approx(m.outputs[("pH", None)].value, 1e-3) == 6.903711162478472


def test_indexed_blockBuild(build_rkt_state_with_indexed_species):
    m = build_rkt_state_with_indexed_species
    m.outputs.display()
    m.property_block = ReaktoroBlock(
        [0, 1],
        aqueous_phase={
            "composition": m.composition,
            "convert_to_rkt_species": True,
        },
        system_state={
            "temperature": m.temp,
            "pressure": m.pressure,
            "pH": m.pH,
        },
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        outputs=m.outputs,
    )
    for blk in m.property_block:
        m.property_block[blk].initialize()
    m.property_block[0].reaktoro_model.display()
    cy_solver = get_cyipopt_watertap_solver()
    cy_solver.options["max_iter"] = 20
    m.pH.unfix()
    m.outputs[0, ("scalingTendency", "Calcite")].fix(5)
    m.outputs[1, ("scalingTendency", "Calcite")].fix(2.5)
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.display()
    assert pytest.approx(m.pH[0].value, 1e-3) == 6.78206
    assert pytest.approx(m.pH[1].value, 1e-3) == 5.719012533419923


def test_indexed_blockBuild_with_speciation_block(
    build_rkt_state_with_indexed_species,
):
    m = build_rkt_state_with_indexed_species
    m.CaO = Var([0, 1], ["CaO"], initialize=0.01, units=pyunits.mol / pyunits.s)
    m.CaO.fix()
    m.outputs.display()
    m.property_block = ReaktoroBlock(
        [0, 1],
        aqueous_phase={
            "composition": m.composition,
            "convert_to_rkt_species": True,
        },
        system_state={
            "temperature": m.temp,
            "pressure": m.pressure,
            "pH": m.pH,
        },
        chemistry_modifier=m.CaO,
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        outputs=m.outputs,
        build_speciation_block=True,
    )
    for blk in m.property_block:
        m.property_block[blk].initialize()
    m.property_block.display()
    cy_solver = get_cyipopt_watertap_solver()
    cy_solver.options["max_iter"] = 20
    m.CaO.unfix()
    m.outputs[(0, "pH", None)].fix(11.5)
    m.outputs[(1, "pH", None)].fix(10)
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.display()
    assert pytest.approx(m.CaO[(0, "CaO")].value, 1e-3) == 0.01732553618254949
    assert pytest.approx(m.CaO[(1, "CaO")].value, 1e-3) == 0.011351679127420139


def test_blockBuild_element_amount_in_phase(build_rkt_state_with_species):
    m = build_rkt_state_with_species
    m.outputs_element = Var(
        [
            ("speciesAmount", "Na+", None),
            ("elementAmountInPhase", "Na", "AqueousPhase"),
            ("elementAmountInPhase", "Ca", "AqueousPhase"),
            ("pH", None, None),
        ],
        initialize=0,
        units=pyunits.mol / pyunits.s,
    )
    m.CaO = Var(["CaO"], initialize=0.002, units=pyunits.mol / pyunits.s)
    m.CaO.fix()
    m.property_block = ReaktoroBlock(
        aqueous_phase={
            "composition": m.composition,
            "convert_to_rkt_species": True,
        },
        system_state={
            "temperature": m.temp,
            "pressure": m.pressure,
            "pH": m.pH,
        },
        chemistry_modifier=m.CaO,
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        outputs=m.outputs_element,
        build_speciation_block=True,
    )
    m.property_block.initialize()
    m.outputs_element.display()
    m.outputs_element["pH", None, None].fix(11.5)
    m.CaO.unfix()
    cy_solver = get_cyipopt_watertap_solver()
    cy_solver.options["max_iter"] = 20
    result = cy_solver.solve(m, tee=True)
    m.display()
    assert_optimal_termination(result)
    assert (
        pytest.approx(
            m.outputs_element["elementAmountInPhase", "Na", "AqueousPhase"].value, 1e-3
        )
        == 0.5
    )
    assert (
        pytest.approx(
            m.outputs_element["elementAmountInPhase", "Ca", "AqueousPhase"].value, 1e-3
        )
        == m.CaO["CaO"].value + m.composition["Ca"].value
    )


def _check_scaling_tendencies(database_file, composition, minerals, pH, pE):
    """build property block with scalingTendency outputs, check values against
    Reaktoro saturation index and jacobian against central finite difference.
    returns jacobian types of each output"""
    import numpy as np
    import reaktoro as rkt

    m = ConcreteModel()
    m.temp = Var(initialize=293.15, units=pyunits.K)
    m.temp.fix()
    m.pressure = Var(initialize=1e5, units=pyunits.Pa)
    m.pressure.fix()
    m.pH = Var(initialize=pH)
    m.pH.fix()
    m.pE = Var(initialize=pE)
    m.pE.fix()
    m.composition = Var(list(composition), initialize=1, units=pyunits.mol / pyunits.s)
    for k, v in composition.items():
        m.composition[k].fix(v)
    m.st = Var(minerals, initialize=1)
    m.blk = ReaktoroBlock(
        aqueous_phase={
            "composition": m.composition,
            "convert_to_rkt_species": False,
            "activity_model": "ActivityModelPhreeqc",
        },
        system_state={
            "temperature": m.temp,
            "pressure": m.pressure,
            "pH": m.pH,
            "pE": m.pE,
        },
        database="PhreeqcDatabase",
        database_file=database_file,
        outputs={("scalingTendency", s): m.st[s] for s in minerals},
        reaktoro_solve_options={"solver_tolerance": 1e-12},
    )
    m.blk.initialize()
    jac_types = {
        s: m.blk.rkt_outputs.rkt_outputs[("scalingTendency", s)].jacobian_type
        for s in minerals
    }
    gb = m.blk.reaktoro_model.get_external_model()
    names_in, names_out = gb.input_names(), gb.output_names()
    x0 = np.array([m.blk.reaktoro_model.inputs[n].value for n in names_in])
    gb.set_input_values(x0)
    y0 = np.array(gb.evaluate_outputs())
    jac = gb.evaluate_jacobian_outputs().toarray()
    aq_props = rkt.AqueousProps(m.blk.rkt_state.state)
    for i, (_, s) in enumerate(names_out):
        assert pytest.approx(y0[i], rel=1e-6) == 10 ** float(
            aq_props.saturationIndex(s)
        ), f"{database_file} {s}"
    # compare with central finite difference, scaled by input magnitude
    jac_fd = np.zeros_like(jac)
    for j in range(len(names_in)):
        if x0[j] == 0:
            continue
        step = abs(x0[j]) * 1e-5
        x_up, x_down = x0.copy(), x0.copy()
        x_up[j] += step
        x_down[j] -= step
        gb.set_input_values(x_up)
        y_up = np.array(gb.evaluate_outputs())
        gb.set_input_values(x_down)
        y_down = np.array(gb.evaluate_outputs())
        jac_fd[:, j] = (y_up - y_down) / (2 * step)
    for i, (_, s) in enumerate(names_out):
        exact = jac[i] * np.abs(x0) / y0[i]
        numerical = jac_fd[i] * np.abs(x0) / y0[i]
        assert (
            np.linalg.norm(exact - numerical) / np.linalg.norm(numerical) < 1e-2
        ), f"{database_file} {s}"
    return jac_types


@pytest.mark.parametrize(
    "pH, pE, Mn",
    [(9.5, 11.54, 4e-6), (8.5, 12.54, 1e-14), (7.57, 13.0, 4e-6)],
)
def test_redox_mineral_scaling_tendency_exact_derivatives(pH, pE, Mn):
    """Birnessite and Manganite reactions contain e-, which is not a system species,
    e- is replaced using the O2/H2O couple so scaling tendency has exact derivatives
    and matches Reaktoro saturation index (including trace Mn)"""
    comp = {
        "H2O": 55.5,
        "Na+": 7.6e-3,
        "Cl-": 6.4e-3,
        "Ca+2": 1.1e-3,
        "Mg+2": 1.9e-3,
        "HCO3-": 7.8e-3,
        "SO4-2": 8.8e-4,
        "Mn+2": Mn,
    }
    jac_types = _check_scaling_tendencies(
        "wateq4f.dat", comp, ["Birnessite", "Manganite", "Calcite"], pH, pE
    )
    # exact (calculated from species activities), not numerical saturation index
    assert all(t == "calculated" for t in jac_types.values())


@pytest.mark.parametrize(
    "database_file",
    [
        "phreeqc.dat",
        "wateq4f.dat",
        "minteq.dat",
        "minteq.v4.dat",
        "Amm.dat",
        "Tipping_Hurley.dat",
        "llnl.dat",
    ],
)
@pytest.mark.parametrize("pH, pE", [(8.5, 12.5), (7.0, -2.0)])
def test_scaling_tendency_matches_saturation_index_all_databases(database_file, pH, pE):
    """for each database, scaling tendency of every redox mineral (e- in reaction) that
    can form from system elements, plus Calcite, matches Reaktoro saturation index and
    has correct derivatives, under oxidizing and reducing conditions"""
    import json
    import reaktoro as rkt

    db = rkt.PhreeqcDatabase(database_file)
    species = {s.name() for s in db.species()}
    candidates = {
        "H2O": 55.5,
        "Na+": 7.6e-3,
        "Cl-": 6.4e-3,
        "Ca+2": 1.1e-3,
        "Mg+2": 1.9e-3,
        "HCO3-": 7.8e-3,
        "SO4-2": 8.8e-4,
        "Mn+2": 4e-6,
        "Fe+2": 4e-6,
    }
    comp = {k: v for k, v in candidates.items() if k in species}
    elements = set()
    for k in comp:
        elements.update(e.symbol() for e, _ in db.species().get(k).elements())

    def supported(s):
        js = json.loads(s.standardThermoModel().params().dumpJson())
        return (
            isinstance(js, list)
            and js[0] is not None
            and ("PhreeqcLgK" in js[0] or "VantHoff" in js[0])
        )

    redox_minerals = []
    for s in db.species():
        if s.aggregateState() in (rkt.AggregateState.Aqueous, rkt.AggregateState.Gas):
            continue
        reactants = [x.name() for x, _ in s.reaction().reactants()]
        if "e-" not in reactants or not supported(s):
            continue
        if {e.symbol() for e, _ in s.elements()} <= elements:
            redox_minerals.append(s.name())
    minerals = redox_minerals + ["Calcite"]
    jac_types = _check_scaling_tendencies(database_file, comp, minerals, pH, pE)
    # all minerals (redox ones through O2 or H2 couple) get exact derivatives
    assert all(t == "calculated" for t in jac_types.values()), jac_types


@pytest.mark.parametrize(
    "database_file",
    [
        "phreeqc.dat",
        "wateq4f.dat",
        "minteq.dat",
        "minteq.v4.dat",
        "Amm.dat",
        "Tipping_Hurley.dat",
        "llnl.dat",
        "core10.dat",
    ],
)
@pytest.mark.parametrize(
    "pH, pE, temperature, pressure",
    [
        (8.5, 12.5, 293.15, 1e5),
        (7.0, 2.0, 293.15, 1e5),
        (7.0, 6.0, 293.15, 1e5),
        (7.5, 12.0, 313.15, 10e5),
    ],
)
def test_pE_matches_aqueous_props_all_databases(
    database_file, pH, pE, temperature, pressure
):
    """pE output is calculated from O2/H2 couple activities (exact derivatives), check
    against Reaktoro aqueous props pE in a property block (after speciation at given
    pH and pE and HCl addition, so pE is not an input) for oxidizing, reducing, mixed
    (both couples above epsilon floor), and elevated temperature/pressure conditions
    (sit.dat is excluded as it needs SIT activity model)"""
    if database_file == "core10.dat" and pE < 4:
        # Reaktoro fails to converge on finite difference perturbations here
        pytest.skip("Reaktoro solve fails for core10.dat reducing state perturbations")
    import numpy as np
    import reaktoro as rkt

    db = rkt.PhreeqcDatabase(database_file)
    species = {s.name() for s in db.species()}
    candidates = {
        "H2O": 55.5,
        "Na+": 7.6e-3,
        "Cl-": 6.4e-3,
        "Ca+2": 1.1e-3,
        "Mg+2": 1.9e-3,
        "HCO3-": 7.8e-3,
        "SO4-2": 8.8e-4,
        "Mn+2": 4e-6,
        "Fe+2": 4e-6,
    }
    comp = {k: v for k, v in candidates.items() if k in species}
    m = ConcreteModel()
    m.temp = Var(initialize=temperature, units=pyunits.K)
    m.temp.fix()
    m.pressure = Var(initialize=pressure, units=pyunits.Pa)
    m.pressure.fix()
    m.pH = Var(initialize=pH)
    m.pH.fix()
    m.pE = Var(initialize=pE)
    m.pE.fix()
    m.composition = Var(list(comp), initialize=1, units=pyunits.mol / pyunits.s)
    for k, v in comp.items():
        m.composition[k].fix(v)
    m.HCl = Var(["HCl"], initialize=1e-3, units=pyunits.mol / pyunits.s)
    m.HCl.fix()
    m.outputs = Var([("pE", None), ("pH", None)], initialize=1)
    m.blk = ReaktoroBlock(
        aqueous_phase={
            "composition": m.composition,
            "convert_to_rkt_species": False,
            "activity_model": "ActivityModelPhreeqc",
        },
        system_state={
            "temperature": m.temp,
            "pressure": m.pressure,
            "pH": m.pH,
            "pE": m.pE,
        },
        database="PhreeqcDatabase",
        database_file=database_file,
        chemistry_modifier=m.HCl,
        outputs=m.outputs,
        build_speciation_block=True,
    )
    m.blk.initialize()
    # exact (calculated from couple activities), not numerical aqueous props pE
    assert m.blk.rkt_outputs.rkt_outputs[("pE", None)].jacobian_type == "calculated"
    gb = m.blk.reaktoro_model.get_external_model()
    names_in, names_out = gb.input_names(), gb.output_names()
    x0 = np.array([m.blk.reaktoro_model.inputs[n].value for n in names_in])
    gb.set_input_values(x0)
    y0 = np.array(gb.evaluate_outputs())
    jac = gb.evaluate_jacobian_outputs().toarray()
    i_pE = names_out.index(("pE", None))
    rkt_pE = float(rkt.AqueousProps(m.blk.rkt_state.state).pE())
    # pressure correction uses species volumes at 25 C and 1 atm (as in scaling
    # tendency), Reaktoro uses volumes at T and P, small difference away from 1 atm
    tolerance = 1e-6 if pressure == 1e5 else 1e-4
    assert pytest.approx(y0[i_pE], abs=tolerance) == rkt_pE
    # HCl addition changes pE (pE is not just passed through)
    assert abs(rkt_pE - pE) > 1e-4
    jac_fd = np.zeros(len(names_in))
    for j in range(len(names_in)):
        if x0[j] == 0:
            continue
        # pE is sensitive, larger step keeps Reaktoro solver tolerance noise small
        step = abs(x0[j]) * 1e-4
        x_up, x_down = x0.copy(), x0.copy()
        x_up[j] += step
        x_down[j] -= step
        gb.set_input_values(x_up)
        y_up = gb.evaluate_outputs()[i_pE]
        gb.set_input_values(x_down)
        y_down = gb.evaluate_outputs()[i_pE]
        jac_fd[j] = (y_up - y_down) / (2 * step)
    exact = jac[i_pE] * np.abs(x0)
    numerical = jac_fd * np.abs(x0)
    assert np.linalg.norm(exact - numerical) / np.linalg.norm(numerical) < 1e-2
