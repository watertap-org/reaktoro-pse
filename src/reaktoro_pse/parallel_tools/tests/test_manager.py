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


@pytest.mark.parametrize(
    "hess_type",
    [
        "ZeroHessian",
        "GaussNewton",
        "LBFGS",
        "BFGS",
        "BFGS_mod",  # Does not work on this example
        "BFGS_damp",
        "BFGS_ipopt",
    ],
)
def test_blockBuild_with_speciation_block(build_rkt_state_with_species, hess_type):
    m = build_rkt_state_with_species
    m.CaO = Var(["CaO"], initialize=0.001, units=pyunits.mol / pyunits.s)
    m.CaO.fix()
    m.reaktoro_manager = ReaktoroBlockManager(
        hessian_options={"hessian_type": hess_type}
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
        jacobian_options={
            "numerical_type": "average",
            "numerical_order": 2,
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
    assert pytest.approx(m.outputs[("pH", None)].value, 1e-3) == 6.7496301
    assert pytest.approx(m.pH.value, 1e-3) == 6.401
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
    assert pytest.approx(m.outputs[("pH", None)].value, 1e-3) == 8.027346120955238
    assert pytest.approx(m.pH.value, 1e-3) == 7.2526416924401556
    m.reaktoro_manager.terminate_workers()
