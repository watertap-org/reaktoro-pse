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
    rkt_inputs = reaktoroInputSpec(rkt_state)
    rkt_inputs.configure_specs(dissolve_species_in_rkt=True)
    rkt_outputs = reaktoroOutputSpec(rkt_state)

    rkt_outputs.register_output("speciesAmount", get_all_indexes=True)
    rkt_outputs.register_output("scalingTendency", "Calcite")
    rkt_outputs.register_output("pH")
    rkt_jacobian = reaktoroJacobianSpec(rkt_state, rkt_outputs)
    rkt_solver = reaktoroSolver(rkt_state, rkt_inputs, rkt_outputs, rkt_jacobian)
    return rkt_solver


def test_solver(build_standard_state):
    rkt_solver = build_standard_state
    rkt_inputs = rkt_solver.rktInputSpec.rktInputs.rktInputList
    rkt_outputs = list(rkt_solver.rktOutputSpec.rktOutputs.keys())
    print(rkt_inputs, rkt_outputs)
    jacobian, outputs = rkt_solver.solve_reaktoro_block()
    print(outputs)
    expected_input_names = [
        "temperature",
        "pressure",
        "pH",
        "CO3-2",
        "CO2",
        "H2O",
        "Na",
        "Mg",
        "Ca",
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
        1.1685140780757134e-07,
        50.000118490980476,
        2.9234940748061945e-05,
        0.001118822642133669,
        0.008160604938364287,
        0.7082105028633766,
        0.00797307238937095,
        0.09996030471842361,
        3.9474966111604523e-05,
        2.203154647963093e-07,
        0.5,
        1.1134620475464961e-07,
        0.0018393950616357136,
        -2.5928680178549815e-13,
        6.999999999999981,
    ]
    for i, eov in enumerate(expected_output_values):
        assert pytest.approx(eov, 1e-3) == outputs[i]
    assert len(expected_output_values) == len(outputs)

    expected_jac_shape = (len(expected_output_keys), len(expected_input_names))
    assert jacobian.shape[0] == expected_jac_shape[0]
    assert jacobian.shape[1] == expected_jac_shape[1]
