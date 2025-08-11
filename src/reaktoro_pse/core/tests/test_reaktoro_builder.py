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
from reaktoro_pse.core.reaktoro_jacobian import (
    ReaktoroJacobianSpec,
)
from reaktoro_pse.core.reaktoro_outputs import (
    ReaktoroOutputSpec,
)

from reaktoro_pse.core.reaktoro_inputs import (
    ReaktoroInputSpec,
)
from reaktoro_pse.core.reaktoro_solver import (
    ReaktoroSolver,
)
from reaktoro_pse.core.reaktoro_block_builder import (
    ReaktoroBlockBuilder,
)

from pyomo.environ import Var, units as pyunits

from reaktoro_pse.core.tests.test_reaktoro_state import (
    build_rkt_state_with_species,
    build_rkt_state_with_species_mass_basis,
)
from pyomo.environ import Block, assert_optimal_termination
from idaes.core.util.model_statistics import degrees_of_freedom
from reaktoro_pse.core.util_classes.cyipopt_solver import (
    get_cyipopt_watertap_solver,
)


@pytest.fixture
def build_with_dissolve_in_rkt(build_rkt_state_with_species):
    m, rkt_state = build_rkt_state_with_species
    # rkt_state.register_solid_phases("Calcite")
    rkt_state.build_state()
    rkt_state.equilibrate_state()
    rkt_inputs = ReaktoroInputSpec(rkt_state)
    m.lime = Var(initialize=0.01, units=pyunits.mol / pyunits.s)
    m.lime.fix()
    rkt_inputs.register_chemistry_modifier("CaO", m.lime)

    rkt_inputs.configure_specs(dissolve_species_in_rkt=True)
    rkt_inputs.build_input_specs()
    rkt_outputs = ReaktoroOutputSpec(rkt_state)
    rkt_outputs.register_output("scalingTendencySaturationIndex", "Calcite")

    rkt_outputs.register_output("scalingTendencySaturationIndex", "Brucite")
    rkt_outputs.register_output("saturationIndex", "Calcite")
    rkt_outputs.register_output("scalingTendency", "Calcite")

    rkt_outputs.register_output("scalingTendency", "Brucite")
    rkt_outputs.register_output("scalingTendencyPyomo", "Calcite")
    rkt_outputs.register_output("scalingTendencyPyomo", "Brucite")
    rkt_outputs.register_output("osmoticPressure", "H2O")
    rkt_outputs.register_output("osmoticPressurePyomo", "H2O")
    rkt_outputs.register_output("pH")
    rkt_jacobian = ReaktoroJacobianSpec(rkt_state, rkt_outputs)
    rkt_solver = ReaktoroSolver(rkt_state, rkt_inputs, rkt_outputs, rkt_jacobian)
    return m, rkt_solver


@pytest.fixture
def build_with_dissolve_in_pyomo(build_rkt_state_with_species):
    m, rkt_state = build_rkt_state_with_species
    rkt_state.build_state()
    rkt_state.equilibrate_state()
    rkt_inputs = ReaktoroInputSpec(rkt_state)
    m.lime = Var(initialize=0.01, units=pyunits.mol / pyunits.s)
    m.lime.fix()
    rkt_inputs.register_chemistry_modifier("CaO", m.lime)
    rkt_inputs.configure_specs(dissolve_species_in_rkt=False)
    rkt_inputs.build_input_specs()
    rkt_outputs = ReaktoroOutputSpec(rkt_state)
    rkt_outputs.register_output("scalingTendencySaturationIndex", "Calcite")
    rkt_outputs.register_output("speciesAmount", get_all_indexes=True)
    rkt_outputs.register_output("scalingTendency", "Calcite")
    rkt_outputs.register_output("scalingTendencyPyomo", "Calcite")
    rkt_outputs.register_output("pH")
    rkt_jacobian = ReaktoroJacobianSpec(rkt_state, rkt_outputs)
    rkt_solver = ReaktoroSolver(rkt_state, rkt_inputs, rkt_outputs, rkt_jacobian)
    return m, rkt_solver


@pytest.fixture
def build_with_dissolve_in_rkt_mass_basis(build_rkt_state_with_species_mass_basis):
    m, rkt_state = build_rkt_state_with_species_mass_basis
    # rkt_state.register_solid_phases("Calcite")
    rkt_state.build_state()
    rkt_state.equilibrate_state()
    rkt_inputs = ReaktoroInputSpec(rkt_state)
    m.lime = Var(initialize=0.01 * 0.05608, units=pyunits.kg / pyunits.s)
    m.lime.fix()
    rkt_inputs.register_chemistry_modifier("CaO", m.lime)

    rkt_inputs.configure_specs(dissolve_species_in_rkt=True)
    rkt_inputs.build_input_specs()
    rkt_outputs = ReaktoroOutputSpec(rkt_state)
    rkt_outputs.register_output("saturationIndex", "Calcite")
    rkt_outputs.register_output("scalingTendency", "Calcite")
    rkt_outputs.register_output("scalingTendency", "Brucite")
    rkt_outputs.register_output("scalingTendencyPyomo", "Calcite")
    rkt_outputs.register_output("scalingTendencyPyomo", "Brucite")
    rkt_outputs.register_output("scalingTendencySaturationIndex", "Calcite")
    rkt_outputs.register_output("scalingTendencySaturationIndex", "Brucite")
    rkt_outputs.register_output("pH")
    # rkt_outputs.register_output("pHDirect")
    rkt_jacobian = ReaktoroJacobianSpec(rkt_state, rkt_outputs)
    rkt_solver = ReaktoroSolver(rkt_state, rkt_inputs, rkt_outputs, rkt_jacobian)
    return m, rkt_solver


@pytest.fixture
def build_with_dissolve_in_pyomo_mass_basis(build_rkt_state_with_species_mass_basis):
    m, rkt_state = build_rkt_state_with_species_mass_basis
    rkt_state.build_state()
    rkt_state.equilibrate_state()
    rkt_inputs = ReaktoroInputSpec(rkt_state)
    m.lime = Var(initialize=0.01 * 0.05608, units=pyunits.kg / pyunits.s)
    m.lime.fix()
    rkt_inputs.register_chemistry_modifier("CaO", m.lime)
    rkt_inputs.configure_specs(dissolve_species_in_rkt=False)
    rkt_inputs.build_input_specs()
    rkt_outputs = ReaktoroOutputSpec(rkt_state)
    rkt_outputs.register_output("scalingTendencySaturationIndex", "Calcite")

    rkt_outputs.register_output("speciesAmount", get_all_indexes=True)
    rkt_outputs.register_output("scalingTendency", "Calcite")
    rkt_outputs.register_output("scalingTendencyPyomo", "Calcite")
    rkt_outputs.register_output("pH")
    rkt_jacobian = ReaktoroJacobianSpec(rkt_state, rkt_outputs)
    rkt_solver = ReaktoroSolver(rkt_state, rkt_inputs, rkt_outputs, rkt_jacobian)
    return m, rkt_solver


def test_build_with_rkt_dissolution(build_with_dissolve_in_rkt):
    m, rkt_solver = build_with_dissolve_in_rkt
    m.rkt_block = Block()
    builder = ReaktoroBlockBuilder(m.rkt_block, rkt_solver)
    builder.initialize()
    m.display()
    m.rkt_block.reaktoro_model.display()
    # will have as many DOFs as outputs due to pyomo not
    # knowing tha graybox exists.
    assert len(m.rkt_block.reaktoro_model.inputs) == len(
        rkt_solver.input_specs.rkt_inputs
    )
    assert len(m.rkt_block.outputs) == len(rkt_solver.output_specs.user_outputs)
    assert len(m.rkt_block.reaktoro_model.outputs) == len(
        rkt_solver.output_specs.rkt_outputs
    )
    assert degrees_of_freedom(m) == 0
    cy_solver = get_cyipopt_watertap_solver()
    cy_solver.options["max_iter"] = 40
    m.pH.unfix()
    m.lime.value = 0.1
    m.rkt_block.outputs[("scalingTendency", "Calcite")].fix(1000)
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.rkt_block.outputs.display()
    assert pytest.approx(m.pH.value, 1e-3) == 7.9832838874247205
    for scalant in ["Calcite", "Brucite"]:
        assert (
            pytest.approx(m.rkt_block.outputs[("scalingTendency", scalant)].value, 1e-5)
            == m.rkt_block.outputs[("scalingTendencyPyomo", scalant)].value
        )
        assert (
            pytest.approx(m.rkt_block.outputs[("scalingTendency", scalant)].value, 1e-5)
            == m.rkt_block.outputs[("scalingTendencySaturationIndex", scalant)].value
        )
        assert (
            pytest.approx(
                m.rkt_block.outputs[("scalingTendencySaturationIndex", scalant)].value,
                1e-5,
            )
            == m.rkt_block.outputs[("scalingTendencyPyomo", scalant)].value
        )
    assert (
        pytest.approx(
            m.rkt_block.outputs[("osmoticPressure", "H2O")].value,
            1e-5,
        )
        == m.rkt_block.outputs[("osmoticPressurePyomo", "H2O")].value
    )


def test_build_with_pyomo_dissolution(build_with_dissolve_in_pyomo):
    m, rkt_solver = build_with_dissolve_in_pyomo
    m.rkt_block = Block()
    builder = ReaktoroBlockBuilder(m.rkt_block, rkt_solver)
    builder.initialize()
    # will have as many DOFs as outputs due to pyomo not
    # knowing tha graybox exists.
    print(rkt_solver.output_specs.rkt_outputs)
    assert degrees_of_freedom(m) == 0
    cy_solver = get_cyipopt_watertap_solver()
    cy_solver.options["max_iter"] = 20
    m.pH.unfix()
    m.display()
    m.rkt_block.outputs[("scalingTendency", "Calcite")].fix(5)
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    assert pytest.approx(m.pH.value, 1e-3) == 6.5257440
    assert (
        pytest.approx(m.rkt_block.outputs[("scalingTendency", "Calcite")].value, 1e-3)
        == m.rkt_block.outputs[("scalingTendencyPyomo", "Calcite")].value
    )
    assert (
        pytest.approx(
            m.rkt_block.outputs[("scalingTendencySaturationIndex", "Calcite")].value,
            1e-3,
        )
        == m.rkt_block.outputs[("scalingTendencyPyomo", "Calcite")].value
    )


def test_build_with_rkt_dissolution_mass_basis(build_with_dissolve_in_rkt_mass_basis):
    m, rkt_solver = build_with_dissolve_in_rkt_mass_basis
    m.rkt_block = Block()
    builder = ReaktoroBlockBuilder(m.rkt_block, rkt_solver)
    builder.initialize()
    # will have as many DOFs as outputs due to pyomo not
    # knowing tha graybox exists.
    assert len(m.rkt_block.reaktoro_model.inputs) == len(
        rkt_solver.input_specs.rkt_inputs
    )
    assert len(m.rkt_block.outputs) == len(rkt_solver.output_specs.user_outputs)
    assert len(m.rkt_block.reaktoro_model.outputs) == len(
        rkt_solver.output_specs.rkt_outputs
    )
    assert degrees_of_freedom(m) == 0
    cy_solver = get_cyipopt_watertap_solver()
    cy_solver.options["max_iter"] = 20
    m.pH.unfix()
    m.rkt_block.outputs[("scalingTendency", "Calcite")].fix(5)

    # m.display()
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    assert pytest.approx(m.pH.value, 1e-3) == 6.5257440
    assert (
        pytest.approx(m.rkt_block.outputs[("scalingTendency", "Calcite")].value, 1e-3)
        == m.rkt_block.outputs[("scalingTendencyPyomo", "Calcite")].value
    )
    assert (
        pytest.approx(
            m.rkt_block.outputs[("scalingTendencySaturationIndex", "Calcite")].value,
            1e-3,
        )
        == m.rkt_block.outputs[("scalingTendencyPyomo", "Calcite")].value
    )


def test_build_with_pyomo_dissolution_mass_basis(
    build_with_dissolve_in_pyomo_mass_basis,
):
    m, rkt_solver = build_with_dissolve_in_pyomo_mass_basis
    m.rkt_block = Block()
    builder = ReaktoroBlockBuilder(m.rkt_block, rkt_solver)
    builder.initialize()
    assert degrees_of_freedom(m) == 0
    cy_solver = get_cyipopt_watertap_solver()
    cy_solver.options["max_iter"] = 20
    m.pH.unfix()
    m.display()
    m.rkt_block.outputs[("scalingTendency", "Calcite")].fix(5)

    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    assert pytest.approx(m.pH.value, 1e-3) == 6.5257440
    assert (
        pytest.approx(m.rkt_block.outputs[("scalingTendency", "Calcite")].value, 1e-3)
        == m.rkt_block.outputs[("scalingTendencyPyomo", "Calcite")].value
    )
    assert (
        pytest.approx(
            m.rkt_block.outputs[("scalingTendencySaturationIndex", "Calcite")].value,
            1e-3,
        )
        == m.rkt_block.outputs[("scalingTendencyPyomo", "Calcite")].value
    )
