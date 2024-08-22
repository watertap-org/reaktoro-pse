import pytest
from reaktoro_pse.core.reaktoro_jacobian import (
    reaktoroJacobianSpec,
)
from reaktoro_pse.core.reaktoro_outputs import (
    reaktoroOutputSpec,
)

from reaktoro_pse.core.reaktoro_inputs import (
    reaktoroInputSpec,
)
from reaktoro_pse.core.reaktoro_solver import (
    reaktoroSolver,
)
from reaktoro_pse.core.reaktoro_block_builder import (
    reaktoroBlockBuilder,
)

from pyomo.environ import Var, units as pyunits

from reaktoro_pse.core.tests.test_reaktoro_state import (
    build_rkt_state_with_species,
)
from pyomo.environ import Block, assert_optimal_termination
from idaes.core.util.model_statistics import degrees_of_freedom
from watertap.core.solvers import get_solver


@pytest.fixture
def build_with_dissolve_in_rkt(build_rkt_state_with_species):
    m, rkt_state = build_rkt_state_with_species
    # rkt_state.register_solid_phases("Calcite")
    rkt_state.build_state()
    rkt_state.equilibrate_state()
    rkt_inputs = reaktoroInputSpec(rkt_state)
    m.lime = Var(initialize=0.01, units=pyunits.mol / pyunits.s)
    m.lime.fix()
    rkt_inputs.register_chemical_addition("CaO", m.lime)

    rkt_inputs.configure_specs(dissolve_species_in_rkt=True)
    rkt_outputs = reaktoroOutputSpec(rkt_state)
    rkt_outputs.register_output("saturationIndex", "Calcite")
    rkt_outputs.register_output("scalingTendency", "Calcite")
    rkt_outputs.register_output("scalingTendencyDirect", "Calcite")

    # rkt_outputs.register_output("scalingTendencyDirect", "Brucite")
    # rkt_outputs.register_output("osmoticPressure", "H2O")
    rkt_outputs.register_output("pH")
    # rkt_outputs.register_output("pHDirect")
    rkt_jacobian = reaktoroJacobianSpec(rkt_state, rkt_outputs)
    rkt_solver = reaktoroSolver(rkt_state, rkt_inputs, rkt_outputs, rkt_jacobian)
    return m, rkt_solver


@pytest.fixture
def build_with_dissolve_in_pyomo(build_rkt_state_with_species):
    m, rkt_state = build_rkt_state_with_species
    rkt_state.build_state()
    rkt_state.equilibrate_state()
    rkt_inputs = reaktoroInputSpec(rkt_state)
    m.lime = Var(initialize=0.01, units=pyunits.mol / pyunits.s)
    m.lime.fix()
    rkt_inputs.register_chemical_addition("CaO", m.lime)
    rkt_inputs.configure_specs(dissolve_species_in_rkt=False)
    rkt_outputs = reaktoroOutputSpec(rkt_state)

    rkt_outputs.register_output("speciesAmount", get_all_indexes=True)
    rkt_outputs.register_output("scalingTendency", "Calcite")
    rkt_outputs.register_output("scalingTendencyDirect", "Calcite")
    rkt_outputs.register_output("pH")
    rkt_jacobian = reaktoroJacobianSpec(rkt_state, rkt_outputs)
    rkt_solver = reaktoroSolver(rkt_state, rkt_inputs, rkt_outputs, rkt_jacobian)
    return m, rkt_solver


def test_build_with_rkt_dissolution(build_with_dissolve_in_rkt):
    m, rkt_solver = build_with_dissolve_in_rkt
    m.rktBlock = Block()
    builder = reaktoroBlockBuilder(m.rktBlock, rkt_solver)
    builder.initialize()
    # will have as many DOFs as outputs due to pyomo not
    # knowing tha graybox exists.
    assert len(m.rktBlock.reaktoro_model.inputs) == len(
        rkt_solver.rktInputSpec.rktInputs
    )
    assert len(m.rktBlock.outputs) == len(rkt_solver.rktOutputSpec.userOutputs)
    assert len(m.rktBlock.reaktoro_model.outputs) == len(
        rkt_solver.rktOutputSpec.rktOutputs
    )
    assert degrees_of_freedom(m) == len(rkt_solver.rktOutputSpec.rktOutputs)
    cy_solver = get_solver(solver="cyipopt-watertap")
    cy_solver.options["max_iter"] = 20
    m.pH.unfix()
    m.rktBlock.outputs[("scalingTendency", "Calcite")].fix(5)
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    assert pytest.approx(m.pH.value, 1e-3) == 6.5257440
    assert (
        pytest.approx(m.rktBlock.outputs[("scalingTendency", "Calcite")].value, 1e-3)
        == m.rktBlock.outputs[("scalingTendencyDirect", "Calcite")].value
    )


def test_build_with_pyomo_dissolution(build_with_dissolve_in_pyomo):
    m, rkt_solver = build_with_dissolve_in_pyomo
    m.rktBlock = Block()
    builder = reaktoroBlockBuilder(m.rktBlock, rkt_solver)
    builder.initialize()
    # will have as many DOFs as outputs due to pyomo not
    # knowing tha graybox exists.
    assert degrees_of_freedom(m) == len(rkt_solver.rktOutputSpec.rktOutputs)
    cy_solver = get_solver(solver="cyipopt-watertap")
    cy_solver.options["max_iter"] = 20
    m.pH.unfix()
    m.rktBlock.outputs[("scalingTendency", "Calcite")].fix(5)
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    assert pytest.approx(m.pH.value, 1e-3) == 6.5257440
    assert (
        pytest.approx(m.rktBlock.outputs[("scalingTendency", "Calcite")].value, 1e-3)
        == m.rktBlock.outputs[("scalingTendencyDirect", "Calcite")].value
    )
    m.display()
