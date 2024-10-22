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
import enum
from matplotlib.font_manager import json_load
import pytest

from reaktoro_pse.core import reaktoro_jacobian
from reaktoro_pse.core.reaktoro_state import (
    ReaktoroState,
)
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
import pickle

__author__ = "Alexander V. Dudchenko (SLAC)"


@pytest.fixture
def build_standard_state(build_rkt_state_with_species):
    m, rkt_state = build_rkt_state_with_species
    rkt_state.register_mineral_phases("Calcite")
    rkt_state.set_mineral_phase_activity_model()
    rkt_state.build_state()
    rkt_state.equilibrate_state()
    rkt_inputs = ReaktoroInputSpec(rkt_state)
    rkt_inputs.configure_specs(dissolve_species_in_rkt=True)
    rkt_inputs.build_input_specs()
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
    print(rkt_solver.input_specs.constraint_dict)
    print(outputs)

    print([list(l) for l in jacobian])
    expected_jacboian = [
        [
            7.65801494995876e-29,
            -3.2593494210913563e-36,
            -2.0741686517690312e-07,
            -1.2444007911747241e-25,
            -1.2444007911747241e-25,
            -6.59081971194882e-25,
            -4.207269823736039e-25,
            -1.49021692379296e-25,
            1.8015999999999957e-09,
        ],
        [
            1.0164395367051604e-20,
            6.866245319043687e-27,
            -4.336808689942018e-19,
            0.0,
            0.0,
            7.572474548453445e-18,
            2.021901526413905e-16,
            3.608224830031759e-16,
            1.0000000000000002,
        ],
        [
            8.113696485592391e-09,
            2.259849376471579e-14,
            1.8969408354645065e-06,
            0.00018499351646543246,
            0.00018499351646543246,
            -5.365768096869442e-09,
            -2.079524546710911e-07,
            -0.00018499351646543243,
            2.1465127438467364e-08,
        ],
        [
            -2.2978260821996717e-05,
            -3.5929224848216754e-12,
            -0.0022255415992898475,
            0.10861973878704173,
            0.10861973878704173,
            2.5854592687684853e-05,
            -7.862397271260802e-05,
            -0.1086197387870417,
            1.2228659892086806e-05,
        ],
        [
            -6.97588777427779e-05,
            1.752980296255229e-11,
            -0.004359091473550756,
            -0.4251074388881366,
            -0.4251074388881366,
            1.23303128504726e-05,
            0.00047786612798534213,
            0.4251074388881366,
            6.417220061052413e-05,
        ],
        [
            -9.506954286039246e-05,
            1.2343205845216623e-11,
            -0.0067549379919742926,
            -1.332845072372407,
            -1.332845072372407,
            1.0000387004805447,
            1.9993273694083085,
            1.332845072372407,
            7.666185878716834e-05,
        ],
        [
            -4.902071817094089e-05,
            1.9528117516004683e-11,
            -0.002301278845038802,
            0.44991531783507444,
            0.44991531783507444,
            -1.3049869556825355e-05,
            -0.0005057528313723928,
            -0.44991531783507444,
            5.220447623445768e-05,
        ],
        [
            -2.319471918733151e-06,
            -1.5705239528364755e-12,
            -0.000168062223965199,
            -0.01617237731406427,
            -0.01617237731406427,
            5.078063204763274e-07,
            0.9989279156972547,
            0.01617237731406427,
            2.8189443979070547e-07,
        ],
        [
            2.2319875536741277e-06,
            1.5720094376045971e-12,
            0.0001658320299424294,
            0.016172510973281984,
            0.016172510973281984,
            -4.6904451234755533e-07,
            0.0010624508845250356,
            -0.016172510973281984,
            -2.824006434580601e-07,
        ],
        [
            8.748436505902249e-08,
            -1.4854847681209652e-15,
            2.2301940227696188e-06,
            -1.3365921772339677e-07,
            -1.3365921772339677e-07,
            -3.876180812891046e-08,
            9.633418220193148e-06,
            1.3365921772339675e-07,
            5.062036675536108e-10,
        ],
        [
            -4.2305648105828074e-23,
            2.1210189671414445e-28,
            7.711300019601966e-21,
            -3.607486416211351e-18,
            -3.607486416211351e-18,
            1.0,
            4.975628662739501e-19,
            3.469447455833694e-18,
            2.8879777636005514e-18,
        ],
        [
            4.818680399137594e-09,
            5.521859127979695e-16,
            1.3833746784683956e-07,
            1.4407827958482247e-09,
            1.4407827958482247e-09,
            -2.402918241875753e-09,
            -3.603326006163468e-09,
            -1.4407827958482247e-09,
            1.2326277977217798e-09,
        ],
        [
            6.975887774277788e-05,
            -1.7529802962552286e-11,
            0.004359091473550756,
            0.4251074388881365,
            0.4251074388881365,
            -1.2330312850461427e-05,
            -0.0004778661279853019,
            0.5748925611118636,
            -6.417220061052484e-05,
        ],
        [
            -8.844850965738781e-18,
            -6.482170044637456e-13,
            0.0,
            -2.592868017854982e-05,
            -0.0001296434008927491,
            -5.185736035709962e-07,
            0.0,
            2.5928680178549815e-05,
            -8.271806125530277e-25,
        ],
        [0.0, 0.0, 1.0000000002666605, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ]
    for i, jrow in enumerate(expected_jacboian):
        for k, jval in enumerate(jrow):
            assert pytest.approx(jval, 1e-3) == jacobian[i][k]

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


def test_pickled_solver(build_standard_state):
    old_rkt_solver = build_standard_state
    export_state = old_rkt_solver.state.export_config()
    export_inputs = old_rkt_solver.input_specs.export_config()
    export_outputs = old_rkt_solver.output_specs.export_config()
    export_jac = old_rkt_solver.jacobian_specs.export_config()
    export_solver = old_rkt_solver.export_config()

    export_data = [
        export_state,
        export_inputs,
        export_outputs,
        export_jac,
        export_solver,
    ]
    pickled_epxort = pickle.dumps(export_data)
    unpickeld_export = pickle.loads(pickled_epxort)

    rkt_state = ReaktoroState()
    rkt_state.load_from_export_object(unpickeld_export[0])
    rkt_state.build_state()
    rkt_state.equilibrate_state()
    rkt_inputs = ReaktoroInputSpec(rkt_state)
    rkt_inputs.load_from_export_object(unpickeld_export[1])
    rkt_inputs.build_input_specs()
    print(rkt_inputs.constraint_dict)
    rkt_outputs = ReaktoroOutputSpec(rkt_state)
    rkt_outputs.load_from_export_object(unpickeld_export[2])
    rkt_jacobian = ReaktoroJacobianSpec(rkt_state, rkt_outputs)
    rkt_jacobian.load_from_export_object(unpickeld_export[3])
    rkt_solver = ReaktoroSolver(rkt_state, rkt_inputs, rkt_outputs, rkt_jacobian)
    rkt_solver.load_from_export_object(unpickeld_export[4])

    rkt_inputs = rkt_solver.input_specs.rkt_inputs.rkt_input_list
    rkt_outputs = list(rkt_solver.output_specs.rkt_outputs.keys())

    jacobian, outputs = rkt_solver.solve_reaktoro_block()
    print(rkt_solver.state.state)
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
    expected_jacboian = [
        [
            7.65801494995876e-29,
            -3.2593494210913563e-36,
            -2.0741686517690312e-07,
            -1.2444007911747241e-25,
            -1.2444007911747241e-25,
            -6.59081971194882e-25,
            -4.207269823736039e-25,
            -1.49021692379296e-25,
            1.8015999999999957e-09,
        ],
        [
            1.0164395367051604e-20,
            6.866245319043687e-27,
            -4.336808689942018e-19,
            0.0,
            0.0,
            7.572474548453445e-18,
            2.021901526413905e-16,
            3.608224830031759e-16,
            1.0000000000000002,
        ],
        [
            8.113696485592391e-09,
            2.259849376471579e-14,
            1.8969408354645065e-06,
            0.00018499351646543246,
            0.00018499351646543246,
            -5.365768096869442e-09,
            -2.079524546710911e-07,
            -0.00018499351646543243,
            2.1465127438467364e-08,
        ],
        [
            -2.2978260821996717e-05,
            -3.5929224848216754e-12,
            -0.0022255415992898475,
            0.10861973878704173,
            0.10861973878704173,
            2.5854592687684853e-05,
            -7.862397271260802e-05,
            -0.1086197387870417,
            1.2228659892086806e-05,
        ],
        [
            -6.97588777427779e-05,
            1.752980296255229e-11,
            -0.004359091473550756,
            -0.4251074388881366,
            -0.4251074388881366,
            1.23303128504726e-05,
            0.00047786612798534213,
            0.4251074388881366,
            6.417220061052413e-05,
        ],
        [
            -9.506954286039246e-05,
            1.2343205845216623e-11,
            -0.0067549379919742926,
            -1.332845072372407,
            -1.332845072372407,
            1.0000387004805447,
            1.9993273694083085,
            1.332845072372407,
            7.666185878716834e-05,
        ],
        [
            -4.902071817094089e-05,
            1.9528117516004683e-11,
            -0.002301278845038802,
            0.44991531783507444,
            0.44991531783507444,
            -1.3049869556825355e-05,
            -0.0005057528313723928,
            -0.44991531783507444,
            5.220447623445768e-05,
        ],
        [
            -2.319471918733151e-06,
            -1.5705239528364755e-12,
            -0.000168062223965199,
            -0.01617237731406427,
            -0.01617237731406427,
            5.078063204763274e-07,
            0.9989279156972547,
            0.01617237731406427,
            2.8189443979070547e-07,
        ],
        [
            2.2319875536741277e-06,
            1.5720094376045971e-12,
            0.0001658320299424294,
            0.016172510973281984,
            0.016172510973281984,
            -4.6904451234755533e-07,
            0.0010624508845250356,
            -0.016172510973281984,
            -2.824006434580601e-07,
        ],
        [
            8.748436505902249e-08,
            -1.4854847681209652e-15,
            2.2301940227696188e-06,
            -1.3365921772339677e-07,
            -1.3365921772339677e-07,
            -3.876180812891046e-08,
            9.633418220193148e-06,
            1.3365921772339675e-07,
            5.062036675536108e-10,
        ],
        [
            -4.2305648105828074e-23,
            2.1210189671414445e-28,
            7.711300019601966e-21,
            -3.607486416211351e-18,
            -3.607486416211351e-18,
            1.0,
            4.975628662739501e-19,
            3.469447455833694e-18,
            2.8879777636005514e-18,
        ],
        [
            4.818680399137594e-09,
            5.521859127979695e-16,
            1.3833746784683956e-07,
            1.4407827958482247e-09,
            1.4407827958482247e-09,
            -2.402918241875753e-09,
            -3.603326006163468e-09,
            -1.4407827958482247e-09,
            1.2326277977217798e-09,
        ],
        [
            6.975887774277788e-05,
            -1.7529802962552286e-11,
            0.004359091473550756,
            0.4251074388881365,
            0.4251074388881365,
            -1.2330312850461427e-05,
            -0.0004778661279853019,
            0.5748925611118636,
            -6.417220061052484e-05,
        ],
        [
            -8.844850965738781e-18,
            -6.482170044637456e-13,
            0.0,
            -2.592868017854982e-05,
            -0.0001296434008927491,
            -5.185736035709962e-07,
            0.0,
            2.5928680178549815e-05,
            -8.271806125530277e-25,
        ],
        [0.0, 0.0, 1.0000000002666605, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ]
    for i, jrow in enumerate(expected_jacboian):
        for k, jval in enumerate(jrow):
            assert pytest.approx(jval, 1e-3) == jacobian[i][k]
