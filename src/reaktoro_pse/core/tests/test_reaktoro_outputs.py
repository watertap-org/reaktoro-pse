###############################################################################
# #################################################################################
# # WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# # through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# # National Renewable Energy Laboratory, and National Energy Technology
# # Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# # of Energy). All rights reserved.
# #
# # Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# # information, respectively. These files are also available online at the URL
# # "https://github.com/watertap-org/reaktoro-pse/"
# #################################################################################
###############################################################################
import pytest
from reaktoro_pse.core.reaktoro_outputs import (
    ReaktoroOutputSpec,
    PropTypes,
    PyomoProperties,
)
from reaktoro_pse.core.tests.test_reaktoro_state import (
    build_rkt_state_with_species,
)


__author__ = "Alexander V. Dudchenko (SLAC)"


@pytest.fixture
def build_standard_state(build_rkt_state_with_species):
    m, rkt_state = build_rkt_state_with_species
    rkt_state.build_state()
    rkt_state.equilibrate_state()
    rkt_outputs = ReaktoroOutputSpec(rkt_state)
    return rkt_outputs


def test_pyomo_properties(build_standard_state):
    """just test that there are no errors during build"""
    rkt_outputs = build_standard_state
    props = PyomoProperties(
        rkt_outputs.state,
        rkt_outputs.supported_properties[PropTypes.chem_prop],
        rkt_outputs.supported_properties[PropTypes.aqueous_prop],
    )
    props.scalingTendency("Calcite")
    props.scalingTendencyDirect("Calcite")
    props.osmoticPressure("H2O")
    props.pHDirect()
    props.vaporPressure("H2O")


def test_output_setup(build_standard_state):
    """testing setting out outputs"""
    rkt_outputs = build_standard_state
    """ make sure we have all expected props"""
    expected_properties = ["chemProp", "aqueousProp", "pyomoBuiltProperties"]
    for ep in expected_properties:
        assert ep in rkt_outputs.supported_properties
    assert len(expected_properties) == len(rkt_outputs.supported_properties)

    """test registering chem props"""
    rkt_outputs.register_output("speciesAmount", "Na+")
    assert ("speciesAmount", "Na+") in rkt_outputs.rkt_outputs
    assert ("speciesAmount", "Na+") in rkt_outputs.user_outputs
    rkt_outputs.register_output("speciesAmount", get_all_indexes=True)
    for specie in rkt_outputs.species:
        assert (
            rkt_outputs.rkt_outputs[("speciesAmount", specie)].property_type
            == PropTypes.chem_prop
        )
    assert len(rkt_outputs.rkt_outputs.keys()) == len(rkt_outputs.species)
    rkt_outputs.register_output("saturationIndex", "Calcite")
    assert ("saturationIndex", "Calcite") in rkt_outputs.rkt_outputs
    assert (
        rkt_outputs.rkt_outputs["saturationIndex", "Calcite"].property_type
        == PropTypes.aqueous_prop
    )
    with pytest.raises(NotImplementedError) as a:
        rkt_outputs.register_output("NotRealProp")

    value = rkt_outputs.evaluate_property(
        rkt_outputs.rkt_outputs[("speciesAmount", "Na+")]
    )
    assert pytest.approx(value, 1e-3) == 0.5

    rkt_outputs.register_output("pH")
    assert ("pH", None) in rkt_outputs.rkt_outputs


def test_pyomo_constraints(build_standard_state):
    rkt_outputs = build_standard_state
    rkt_outputs.register_output("osmoticPressure", "H2O")
    assert ("osmoticPressure", "H2O") in rkt_outputs.user_outputs
    assert ("osmoticPressure", "H2O") not in rkt_outputs.rkt_outputs
    assert ("speciesStandardVolume", "H2O") in rkt_outputs.rkt_outputs
    assert ("speciesActivityLn", "H2O") in rkt_outputs.rkt_outputs
    assert ("temperature", None) in rkt_outputs.rkt_outputs
    assert (
        rkt_outputs.user_outputs[("osmoticPressure", "H2O")].property_type
        == PropTypes.pyomo_built_prop
    )
    assert (
        rkt_outputs.user_outputs[("osmoticPressure", "H2O")]
        .pyomo_build_options.properties[("speciesStandardVolume", "H2O")]
        .property_type
        == PropTypes.chem_prop
    )

    rkt_outputs.register_output("pHDirect")

    assert ("pHDirect", None) in rkt_outputs.user_outputs

    rkt_outputs.register_output("scalingTendency", "Calcite")
    assert ("scalingTendency", "Calcite") in rkt_outputs.user_outputs
    assert ("scalingTendency", "Calcite") not in rkt_outputs.rkt_outputs
    assert (
        rkt_outputs.user_outputs[("scalingTendency", "Calcite")].property_type
        == PropTypes.pyomo_built_prop
    )

    rkt_outputs.register_output("scalingTendencyDirect", "Calcite")
    assert ("scalingTendency", "Calcite") in rkt_outputs.user_outputs
    assert ("scalingTendency", "Calcite") not in rkt_outputs.rkt_outputs
    assert (
        rkt_outputs.user_outputs[("scalingTendencyDirect", "Calcite")].property_type
        == PropTypes.pyomo_built_prop
    )
    assert (
        rkt_outputs.user_outputs[
            ("scalingTendencyDirect", "Calcite")
        ].pyomo_build_options.options["logk_type"]
        == "Analytical"
    )
    assert (
        rkt_outputs.user_outputs[("scalingTendencyDirect", "Calcite")]
        .pyomo_build_options.properties[("speciesActivityLn", "Ca+2")]
        .stoichiometric_coeff
        == 1
    )
    assert (
        rkt_outputs.user_outputs[("scalingTendencyDirect", "Calcite")]
        .pyomo_build_options.properties[("speciesActivityLn", "CO3-2")]
        .stoichiometric_coeff
        == 1
    )
