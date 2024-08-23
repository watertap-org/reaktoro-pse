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
from reaktoro_pse.core.tests.test_reaktoro_state import (
    build_rkt_state_with_species,
)

__author__ = "Alexander V. Dudchenko (SLAC)"


@pytest.fixture
def build_standard_state(build_rkt_state_with_species):
    m, rkt_state = build_rkt_state_with_species
    rkt_state.register_mineral_phases("Calcite")
    rkt_state.build_state()
    rkt_state.equilibrate_state()
    rkt_inputs = ReaktoroInputSpec(rkt_state)
    rkt_inputs.configure_specs(dissolve_species_in_rkt=True)
    rkt_outputs = ReaktoroOutputSpec(rkt_state)

    rkt_outputs.register_output("speciesAmount", get_all_indexes=True)
    rkt_outputs.register_output("scalingTendency", "Calcite")
    rkt_outputs.register_output("pH")
    rkt_jacobian = ReaktoroJacobianSpec(rkt_state, rkt_outputs)
    rkt_solver = ReaktoroSolver(rkt_state, rkt_inputs, rkt_outputs, rkt_jacobian)
    return rkt_solver


def test_solver(build_standard_state):
    rkt_solver = build_standard_state
    rkt_inputs = rkt_solver.input_specs.rkt_inputs.rkt_input_list
    rkt_outputs = list(rkt_solver.output_specs.rkt_outputs.keys())
    print(rkt_inputs, rkt_outputs)
    jacobian, outputs = rkt_solver.solve_reaktoro_block()
    print(outputs)
    expected_input_names = [
        "temperature",
        "pressure",
        "pH",
        "CO3-2",
        "CO2",
        "Na",
        "Mg",
        "Ca",
        "H2O",
    ]
    for i, ein in enumerate(expected_input_names):
        assert ein in rkt_inputs[i]
    assert len(expected_input_names) == len(rkt_inputs)
    expected_output_keys = [
        ("speciesAmount", "H+"),
        ("speciesAmount", "H2O"),
        ("speciesAmount", "CO3-2"),
        ("speciesAmount", "CO2"),
        ("speciesAmount", "Ca+2"),
        ("speciesAmount", "Cl-"),
        ("speciesAmount", "HCO3-"),
        ("speciesAmount", "Mg+2"),
        ("speciesAmount", "MgCO3"),
        ("speciesAmount", "MgOH+"),
        ("speciesAmount", "Na+"),
        ("speciesAmount", "OH-"),
        ("speciesAmount", "Calcite"),
        ("saturationIndex", "Calcite"),
        ("pH", None),
    ]
    for i, eon in enumerate(expected_output_keys):
        assert eon == rkt_outputs[i]
    assert len(expected_output_keys) == len(rkt_outputs)
    expected_output_values = [
        9.00799999999998e-08,
        50.0,
        1.2347717588732544e-06,
        0.0007251176324639623,
        0.0028374543608618453,
        0.7024523350480891,
        0.0030030389116422994,
        0.09989096781756118,
        0.00010806304499670852,
        9.69137442114817e-07,
        0.5,
        6.007103894733061e-08,
        0.007162545639138155,
        -5.1857360357099634e-15,
        7.0,
    ]
    for i, eov in enumerate(expected_output_values):
        assert pytest.approx(eov, 1e-3) == outputs[i]
    assert len(expected_output_values) == len(outputs)

    expected_jac_shape = (len(expected_output_keys), len(expected_input_names))
    assert jacobian.shape[0] == expected_jac_shape[0]
    assert jacobian.shape[1] == expected_jac_shape[1]
