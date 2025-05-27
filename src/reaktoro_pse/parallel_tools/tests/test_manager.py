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

from reaktoro_pse.parallel_tools.reaktoro_block_manager import ReaktoroBlockManager

from reaktoro_pse.reaktoro_block import ReaktoroBlock
from reaktoro_pse.tests.test_reaktoro_block import (
    build_rkt_state_with_species,
)
from pyomo.environ import (
    ConcreteModel,
    Var,
    assert_optimal_termination,
    units as pyunits,
)
from watertap_solvers import get_solver


def test_blockBuild_with_speciation_block(build_rkt_state_with_species):
    m = build_rkt_state_with_species
    m.CaO = Var(["CaO"], initialize=0.001, units=pyunits.mol / pyunits.s)
    m.CaO.fix()
    m.reaktoro_manager = ReaktoroBlockManager()
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
        jacobian_options={
            "numerical_type": "average",
            "numerical_order": 2,
            "numerical_step": 1e-8,
        },
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        chemistry_modifier=m.CaO,
        outputs=m.outputs,
        build_speciation_block=True,
        reaktoro_block_manager=m.reaktoro_manager,
    )
    m.reaktoro_manager.build_reaktoro_blocks()
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
    m.property_block.update_block_scaling()
    m.property_block.update_jacobian_scaling()
    scaling_result = m.property_block.display_jacobian_scaling()
    print(scaling_result)
    expected_scaling = {
        "speciation_block": {
            ("speciesAmount", "H+"): 9.008000000000023e-08,
            ("speciesAmount", "H2O"): 49.999999999999964,
            ("speciesAmount", "CO3-2"): 3.2176918698800403e-06,
            ("speciesAmount", "CO2"): 0.001890049191646547,
            ("speciesAmount", "Ca+2"): 0.01000000000000001,
            ("speciesAmount", "Cl-"): 0.6916048636482265,
            ("speciesAmount", "HCO3-"): 0.00782561945440841,
            ("speciesAmount", "SO4-2"): 0.009999916992739206,
            ("speciesAmount", "HSO4-"): 8.300726079741266e-08,
            ("speciesAmount", "Mg+2"): 0.09971791911744696,
            ("speciesAmount", "MgCO3"): 0.0002811136620751731,
            ("speciesAmount", "MgOH+"): 9.672204778991323e-07,
            ("speciesAmount", "Na+"): 0.5000000000000001,
            ("speciesAmount", "OH-"): 6.005625780099479e-08,
        },
        "property_block": {
            ("saturationIndex", "Calcite"): 1.0044028388381383,
            ("pH", None): 7.000000001868148,
            ("elementAmount", "H"): 100.00508467563753,
            ("elementAmount", "O"): 50.066130674180215,
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
    m.property_block.update_jacobian_scaling(new_scaling)
    scaling_result = m.property_block.display_jacobian_scaling()

    assert "speciation_block" in scaling_result
    for key in scaling_result["speciation_block"]:
        assert scaling_result["speciation_block"][key] == 1
    new_scaling = {}
    for key in scaling_result["property_block"]:
        new_scaling[key] = 1
        assert (
            pytest.approx(scaling_result["property_block"][key], 1e-3)
            == expected_scaling["property_block"][key]
        )
    m.property_block.update_jacobian_scaling(new_scaling)
    scaling_result = m.property_block.display_jacobian_scaling()

    assert "property_block" in scaling_result
    for key in scaling_result["property_block"]:
        assert scaling_result["property_block"][key] == 1
    m.property_block.display_reaktoro_state()
    m.reaktoro_manager.terminate_workers()


def test_blockBuild_with_wateqf_data_base(build_rkt_state_with_species):
    m = build_rkt_state_with_species
    m.CaO = Var(["CaO"], initialize=0.001, units=pyunits.mol / pyunits.s)
    m.CaO.fix()
    m.reaktoro_manager = ReaktoroBlockManager()
    m.property_block = ReaktoroBlock(
        aqueous_phase={
            "composition": m.composition,
            "convert_to_rkt_species": True,
            "activity_model": "ActivityModelPhreeqc",
        },
        system_state={
            "temperature": m.temp,
            "pressure": m.pressure,
            "pH": m.pH,
        },
        database="PhreeqcDatabase",
        database_file="wateq4f.dat",
        chemistry_modifier=m.CaO,
        outputs=m.outputs,
        build_speciation_block=True,
        reaktoro_block_manager=m.reaktoro_manager,
        exclude_species_list=[
            "CH4",
            "O2",
            "S2-2",
            "S4-2",
            "S3-2",
            "S5-2",
            "S6-2",
        ],
    )
    m.reaktoro_manager.build_reaktoro_blocks()
    m.property_block.initialize()
    cy_solver = get_solver(solver="cyipopt-watertap")
    cy_solver.options["max_iter"] = 20
    m.pH.unfix()
    m.outputs[("scalingTendency", "Calcite")].fix(5)
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.display()
    assert pytest.approx(m.outputs[("pH", None)].value, 1e-2) == 7.49301431889365
    assert pytest.approx(m.pH.value, 1e-2) == 6.669409618808208
    m.reaktoro_manager.terminate_workers()
