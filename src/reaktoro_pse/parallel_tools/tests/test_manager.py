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

import pytest

from reaktoro_pse.parallel_tools.reaktoro_block_manager import ReaktoroBlockManager

from reaktoro_pse.reaktoro_block import ReaktoroBlock
from reaktoro_pse.tests.test_reaktoro_block import (
    build_rkt_state_with_species,
    add_simple_test_constraints,
)
from pyomo.environ import (
    ConcreteModel,
    Constraint,
    Var,
    assert_optimal_termination,
    units as pyunits,
)
from watertap_solvers import get_solver

from reaktoro_pse.core.util_classes.cyipopt_solver import (
    get_cyipopt_watertap_solver,
)


from pyomo.contrib.pynumero.interfaces.external_grey_box import (
    ExternalGreyBoxModel,
)


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
        "SR1",
        "PSB",
        "SR1_lagrangian",
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
    m.reaktoro_manager.initialize()

    cy_solver = get_cyipopt_watertap_solver()
    cy_solver.options["max_iter"] = 20
    m.pH.unfix()
    m.outputs[("scalingTendency", "Calcite")].fix(5)
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.display()
    assert pytest.approx(m.outputs[("pH", None)].value, 1e-3) == 6.7496301
    assert pytest.approx(m.pH.value, 1e-3) == 6.401

    m.reaktoro_manager.terminate_workers()


def test_deactivate(build_rkt_state_with_species):
    m = build_rkt_state_with_species
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
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        outputs=m.outputs,
        reaktoro_block_manager=m.reaktoro_manager,
    )
    m.reaktoro_manager.build_reaktoro_blocks()
    m.property_block.initialize()
    m.pH.fix()
    m.composition["H2O"].unfix()
    m.composition["H2O"].setlb(30)
    m.outputs[("scalingTendency", "Calcite")].fix(5)
    m.reaktoro_manager.deactivate()

    assert m.property_block.active == False
    assert m.reaktoro_manager.active == False
    for v in m.reaktoro_manager.component_data_objects(Constraint):
        print(v.name)
        assert v.active == False
    for v in m.reaktoro_manager.component_data_objects(ExternalGreyBoxModel):
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

    m.reaktoro_manager.activate()
    assert m.property_block.active == True
    assert m.reaktoro_manager.active == True
    for v in m.reaktoro_manager.component_data_objects(Constraint):
        print(v.name)
        assert v.active == True
    for v in m.reaktoro_manager.component_data_objects(ExternalGreyBoxModel):
        assert v.active == True

    # this solve should solve
    m.composition["H2O"].unfix()
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    assert pytest.approx(m.composition["H2O"].value, 1e-3) == 68.0601837

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
    assert pytest.approx(m.outputs[("pH", None)].value, 1e-3) == 8.02734612095523
    assert pytest.approx(m.pH.value, 1e-3) == 7.2526416924401556
    m.reaktoro_manager.terminate_workers()


@pytest.mark.parametrize("parallel_mode", [False, True])
def test_species_amount_offset_with_manager(
    build_rkt_state_with_species, parallel_mode
):
    """managed blocks floor species amounts in the shared gray box like unmanaged blocks"""
    import idaes.core.util.scaling as iscale

    m = build_rkt_state_with_species
    keys = [("speciesAmount", "Calcite"), ("speciesAmount", "Brucite"), ("pH", None)]
    m.mineral_outputs = Var(keys, initialize=1e-3)
    sf = 1e3
    for key in keys[:2]:
        m.mineral_outputs[key].setlb(0)
        iscale.set_scaling_factor(m.mineral_outputs[key], sf)
    m.reaktoro_manager = ReaktoroBlockManager(use_parallel_mode=parallel_mode)
    m.property_block = ReaktoroBlock(
        aqueous_phase={
            "composition": m.composition,
            "convert_to_rkt_species": True,
            "activity_model": "ActivityModelPitzer",
        },
        system_state={"temperature": m.temp, "pressure": m.pressure, "pH": m.pH},
        mineral_phase={"phase_components": ["Calcite", "Brucite"]},
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        outputs=m.mineral_outputs,
        reaktoro_solve_options={
            "species_amount_zero_offset": 1e-4,
            "species_amount_offset_type": "smooth_floor",
        },
        reaktoro_block_manager=m.reaktoro_manager,
    )
    m.reaktoro_manager.build_reaktoro_blocks()
    try:
        m.reaktoro_manager.initialize()
        gray_box = m.reaktoro_manager.reaktoro_model.get_external_model()
        floored = {gray_box.outputs[i] for i in gray_box._output_floors}
        assert floored == {(0, keys[0]), (0, keys[1])}
        cy_solver = get_cyipopt_watertap_solver()
        cy_solver.options["max_iter"] = 20
        result = cy_solver.solve(m, tee=True)
        assert_optimal_termination(result)
        a = 1e-4 / sf
        brucite = m.mineral_outputs[keys[1]].value
        assert pytest.approx(brucite, rel=1e-3, abs=1e-13) == a
        calcite_raw = m.property_block.rkt_solver.output_specs.user_outputs[keys[0]]
        if not parallel_mode:
            # raw value is only tracked in the main process solver in serial mode
            assert (
                pytest.approx(m.mineral_outputs[keys[0]].value, rel=1e-9)
                == (calcite_raw.value**2 + a**2) ** 0.5
            )
    finally:
        m.reaktoro_manager.terminate_workers()
