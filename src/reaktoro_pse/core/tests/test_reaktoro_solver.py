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

__author__ = "Alexander V. Dudchenko"


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
            1.801600000000002e-09,
            -1.4315700532822493e-28,
            5.794102885036005e-35,
            -2.074168651769038e-07,
            -8.965030771081448e-25,
            -8.965030771081448e-25,
            1.0644138100964667e-24,
            9.379863075151076e-25,
            -3.1454196417214056e-25,
        ],
        [
            1.0000000000000002,
            -2.0328790734103208e-20,
            2.3829910224916325e-26,
            -8.673617379884035e-19,
            -1.3877787807814457e-16,
            -1.3877787807814457e-16,
            1.796150305311689e-16,
            6.898236322439022e-17,
            -4.163336342344337e-16,
        ],
        [
            2.146512744557476e-08,
            8.113696482996199e-09,
            2.2598493767992473e-14,
            1.8969408352455002e-06,
            0.00018499351662456116,
            0.00018499351662456116,
            -5.365768080910646e-09,
            -2.079524548999669e-07,
            -0.0001849935166245611,
        ],
        [
            1.2228659844280275e-05,
            -2.2978260737355303e-05,
            -3.592922474492613e-12,
            -0.002225541591336931,
            0.10861973842012142,
            0.10861973842012142,
            2.5854592599210893e-05,
            -7.86239724961944e-05,
            -0.10861973842012139,
        ],
        [
            6.41722005797586e-05,
            -6.975887769188371e-05,
            1.7529802966639437e-11,
            -0.004359091468947705,
            -0.4251074388539886,
            -0.4251074388539886,
            1.2330312802203184e-05,
            0.00047786612806184994,
            0.42510743885398855,
        ],
        [
            7.666185870883717e-05,
            -9.506954272524157e-05,
            1.2343205859507627e-11,
            -0.0067549379793988215,
            -1.3328450727192371,
            -1.3328450727192371,
            1.0000387004804068,
            1.9993273694081843,
            1.3328450727192367,
        ],
        [
            5.220447625174331e-05,
            -4.902071820507287e-05,
            1.9528117509637544e-11,
            -0.0023012788483691717,
            0.4499153182220849,
            0.4499153182220849,
            -1.3049869518012584e-05,
            -0.0005057528319290325,
            -0.4499153182220848,
        ],
        [
            2.818944400422693e-07,
            -2.319471919120467e-06,
            -1.5705239529583894e-12,
            -0.00016806222394592018,
            -0.016172377327962704,
            -0.016172377327962704,
            5.078063190815328e-07,
            0.9989279156968379,
            0.0161723773279627,
        ],
        [
            -2.8240064371000217e-07,
            2.2319875540614784e-06,
            1.5720094377265117e-12,
            0.00016583202992315134,
            0.016172510987180554,
            0.016172510987180554,
            -4.690445109520167e-07,
            0.0010624508849418864,
            -0.01617251098718055,
        ],
        [
            5.062036675557344e-10,
            8.748436505898914e-08,
            -1.4854847681221633e-15,
            2.230194022768744e-06,
            -1.3365921785260104e-07,
            -1.3365921785260104e-07,
            -3.876180812890835e-08,
            9.633418220190593e-06,
            1.3365921785260098e-07,
        ],
        [
            -2.3566414962092144e-18,
            1.3643949762567108e-23,
            9.997123117458554e-30,
            1.213262760138025e-20,
            4.1389711976061096e-18,
            4.1389711976061096e-18,
            1.0000000000000002,
            9.143744537895742e-19,
            0.0,
        ],
        [
            1.232627797721995e-09,
            4.818680399138192e-09,
            5.52185912798044e-16,
            1.3833746784683872e-07,
            1.4407827961998218e-09,
            1.4407827961998218e-09,
            -2.4029182418759077e-09,
            -3.6033260061636596e-09,
            -1.440782796199822e-09,
        ],
        [
            -6.417220057975922e-05,
            6.975887769188371e-05,
            -1.7529802966639434e-11,
            0.004359091468947705,
            0.4251074388539888,
            0.4251074388539888,
            -1.2330312802241344e-05,
            -0.00047786612806182663,
            0.5748925611460114,
        ],
        [
            6.99010730437894e-18,
            8.211203626796692e-08,
            2.617088630456062e-10,
            7.669146299661466e-16,
            0.0,
            0.0,
            9.9858675776842e-19,
            3.195477624858944e-17,
            6.544338175711117e-14,
        ],
        [
            3.2817658674128153e-22,
            1.3504437466636273e-22,
            9.04216384611702e-31,
            0.9999999999999998,
            1.1997251576853272e-18,
            1.1997251576853272e-18,
            -3.7434118656825577e-19,
            8.715019409346677e-19,
            -1.82661571664871e-19,
        ],
    ]
    for i, jrow in enumerate(expected_jacboian):
        for k, jval in enumerate(jrow):
            assert pytest.approx(jval, 1e-3) == jacobian[i][k]

    expected_input_names = [
        "H2O",
        "temperature",
        "pressure",
        "pH",
        "HCO3-",
        "CO2",
        "Na+",
        "Mg+2",
        "Ca+2",
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
        ("scalingTendency", "Calcite"),
        ("pH", None),
    ]
    for i, eon in enumerate(expected_output_keys):
        assert eon == rkt_outputs[i]
    assert len(expected_output_keys) == len(rkt_outputs)
    expected_output_values = [
        9.008000000000008e-08,
        50.0,
        1.2347717593728471e-06,
        0.0007251176296841211,
        0.002837454359341226,
        0.7024523350437445,
        0.0030030389128573387,
        0.09989096781751751,
        0.00010806304504039303,
        9.691374421145397e-07,
        0.5,
        6.007103894734163e-08,
        0.0071625456406587745,
        1.0000000041068018,
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
        "H2O",
        "temperature",
        "pressure",
        "pH",
        "HCO3-",
        "CO2",
        "Na+",
        "Mg+2",
        "Ca+2",
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
        ("scalingTendency", "Calcite"),
        ("pH", None),
    ]
    for i, eon in enumerate(expected_output_keys):
        assert eon == rkt_outputs[i]
    assert len(expected_output_keys) == len(rkt_outputs)
    expected_output_values = [
        9.008000000000008e-08,
        50.0,
        1.2347717593728471e-06,
        0.0007251176296841211,
        0.002837454359341226,
        0.7024523350437445,
        0.0030030389128573387,
        0.09989096781751751,
        0.00010806304504039303,
        9.691374421145397e-07,
        0.5,
        6.007103894734163e-08,
        0.0071625456406587745,
        1.0000000041068018,
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
            1.801600000000002e-09,
            -1.4315700532822493e-28,
            5.794102885036005e-35,
            -2.074168651769038e-07,
            -8.965030771081448e-25,
            -8.965030771081448e-25,
            1.0644138100964667e-24,
            9.379863075151076e-25,
            -3.1454196417214056e-25,
        ],
        [
            1.0000000000000002,
            -2.0328790734103208e-20,
            2.3829910224916325e-26,
            -8.673617379884035e-19,
            -1.3877787807814457e-16,
            -1.3877787807814457e-16,
            1.796150305311689e-16,
            6.898236322439022e-17,
            -4.163336342344337e-16,
        ],
        [
            2.146512744557476e-08,
            8.113696482996199e-09,
            2.2598493767992473e-14,
            1.8969408352455002e-06,
            0.00018499351662456116,
            0.00018499351662456116,
            -5.365768080910646e-09,
            -2.079524548999669e-07,
            -0.0001849935166245611,
        ],
        [
            1.2228659844280275e-05,
            -2.2978260737355303e-05,
            -3.592922474492613e-12,
            -0.002225541591336931,
            0.10861973842012142,
            0.10861973842012142,
            2.5854592599210893e-05,
            -7.86239724961944e-05,
            -0.10861973842012139,
        ],
        [
            6.41722005797586e-05,
            -6.975887769188371e-05,
            1.7529802966639437e-11,
            -0.004359091468947705,
            -0.4251074388539886,
            -0.4251074388539886,
            1.2330312802203184e-05,
            0.00047786612806184994,
            0.42510743885398855,
        ],
        [
            7.666185870883717e-05,
            -9.506954272524157e-05,
            1.2343205859507627e-11,
            -0.0067549379793988215,
            -1.3328450727192371,
            -1.3328450727192371,
            1.0000387004804068,
            1.9993273694081843,
            1.3328450727192367,
        ],
        [
            5.220447625174331e-05,
            -4.902071820507287e-05,
            1.9528117509637544e-11,
            -0.0023012788483691717,
            0.4499153182220849,
            0.4499153182220849,
            -1.3049869518012584e-05,
            -0.0005057528319290325,
            -0.4499153182220848,
        ],
        [
            2.818944400422693e-07,
            -2.319471919120467e-06,
            -1.5705239529583894e-12,
            -0.00016806222394592018,
            -0.016172377327962704,
            -0.016172377327962704,
            5.078063190815328e-07,
            0.9989279156968379,
            0.0161723773279627,
        ],
        [
            -2.8240064371000217e-07,
            2.2319875540614784e-06,
            1.5720094377265117e-12,
            0.00016583202992315134,
            0.016172510987180554,
            0.016172510987180554,
            -4.690445109520167e-07,
            0.0010624508849418864,
            -0.01617251098718055,
        ],
        [
            5.062036675557344e-10,
            8.748436505898914e-08,
            -1.4854847681221633e-15,
            2.230194022768744e-06,
            -1.3365921785260104e-07,
            -1.3365921785260104e-07,
            -3.876180812890835e-08,
            9.633418220190593e-06,
            1.3365921785260098e-07,
        ],
        [
            -2.3566414962092144e-18,
            1.3643949762567108e-23,
            9.997123117458554e-30,
            1.213262760138025e-20,
            4.1389711976061096e-18,
            4.1389711976061096e-18,
            1.0000000000000002,
            9.143744537895742e-19,
            0.0,
        ],
        [
            1.232627797721995e-09,
            4.818680399138192e-09,
            5.52185912798044e-16,
            1.3833746784683872e-07,
            1.4407827961998218e-09,
            1.4407827961998218e-09,
            -2.4029182418759077e-09,
            -3.6033260061636596e-09,
            -1.440782796199822e-09,
        ],
        [
            -6.417220057975922e-05,
            6.975887769188371e-05,
            -1.7529802966639434e-11,
            0.004359091468947705,
            0.4251074388539888,
            0.4251074388539888,
            -1.2330312802241344e-05,
            -0.00047786612806182663,
            0.5748925611460114,
        ],
        [
            6.99010730437894e-18,
            8.211203626796692e-08,
            2.617088630456062e-10,
            7.669146299661466e-16,
            0.0,
            0.0,
            9.9858675776842e-19,
            3.195477624858944e-17,
            6.544338175711117e-14,
        ],
        [
            3.2817658674128153e-22,
            1.3504437466636273e-22,
            9.04216384611702e-31,
            0.9999999999999998,
            1.1997251576853272e-18,
            1.1997251576853272e-18,
            -3.7434118656825577e-19,
            8.715019409346677e-19,
            -1.82661571664871e-19,
        ],
    ]
    for i, jrow in enumerate(expected_jacboian):
        for k, jval in enumerate(jrow):
            assert pytest.approx(jval, 1e-3) == jacobian[i][k]
