import pytest
from reaktoro_pse.core.reaktoro_outputs import (
    reaktoroOutputSpec,
    propTypes,
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
    rkt_outputs = reaktoroOutputSpec(rkt_state)
    return rkt_outputs


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
    assert ("speciesAmount", "Na+") in rkt_outputs.rktOutputs
    assert ("speciesAmount", "Na+") in rkt_outputs.userOutputs
    rkt_outputs.register_output("speciesAmount", get_all_indexes=True)
    for specie in rkt_outputs.species:
        assert (
            rkt_outputs.rktOutputs[("speciesAmount", specie)].propertyType
            == propTypes.chemProp
        )
    assert len(rkt_outputs.rktOutputs.keys()) == len(rkt_outputs.species)
    rkt_outputs.register_output("saturationIndex", "Calcite")
    assert ("saturationIndex", "Calcite") in rkt_outputs.rktOutputs
    assert (
        rkt_outputs.rktOutputs["saturationIndex", "Calcite"].propertyType
        == propTypes.aqueousProp
    )
    with pytest.raises(NotImplementedError) as a:
        rkt_outputs.register_output("NotRealProp")

    value = rkt_outputs.evaluate_property(
        rkt_outputs.rktOutputs[("speciesAmount", "Na+")]
    )
    assert pytest.approx(value, 1e-3) == 0.5

    rkt_outputs.register_output("pH")
    assert ("pH", None) in rkt_outputs.rktOutputs


def test_pyomo_constraints(build_standard_state):
    rkt_outputs = build_standard_state
    rkt_outputs.register_output("osmoticPressure", "H2O")
    assert ("osmoticPressure", "H2O") in rkt_outputs.userOutputs
    assert ("osmoticPressure", "H2O") not in rkt_outputs.rktOutputs
    assert ("speciesStandardVolume", "H2O") in rkt_outputs.rktOutputs
    assert ("speciesActivityLn", "H2O") in rkt_outputs.rktOutputs
    assert ("temperature", None) in rkt_outputs.rktOutputs
    assert ("pressure", None) in rkt_outputs.rktOutputs
    assert (
        rkt_outputs.userOutputs[("osmoticPressure", "H2O")].propertyType
        == propTypes.pyomoBuiltProperties
    )
    assert (
        rkt_outputs.userOutputs[("osmoticPressure", "H2O")]
        .pyomoBuildOptions.properties[("speciesStandardVolume", "H2O")]
        .propertyType
        == propTypes.chemProp
    )

    rkt_outputs.register_output("pHDirect")

    assert ("pHDirect", None) in rkt_outputs.userOutputs

    rkt_outputs.register_output("scalingTendency", "Calcite")
    assert ("scalingTendency", "Calcite") in rkt_outputs.userOutputs
    assert ("scalingTendency", "Calcite") not in rkt_outputs.rktOutputs
    assert (
        rkt_outputs.userOutputs[("scalingTendency", "Calcite")].propertyType
        == propTypes.pyomoBuiltProperties
    )

    rkt_outputs.register_output("scalingTendencyDirect", "Calcite")
    assert ("scalingTendency", "Calcite") in rkt_outputs.userOutputs
    assert ("scalingTendency", "Calcite") not in rkt_outputs.rktOutputs
    assert (
        rkt_outputs.userOutputs[("scalingTendencyDirect", "Calcite")].propertyType
        == propTypes.pyomoBuiltProperties
    )
    assert (
        rkt_outputs.userOutputs[
            ("scalingTendencyDirect", "Calcite")
        ].pyomoBuildOptions.options["logk_type"]
        == "Analytical"
    )
    assert (
        rkt_outputs.userOutputs[("scalingTendencyDirect", "Calcite")]
        .pyomoBuildOptions.properties[("speciesActivityLn", "Ca+2")]
        .stoichiometricCoeff
        == 1
    )
    assert (
        rkt_outputs.userOutputs[("scalingTendencyDirect", "Calcite")]
        .pyomoBuildOptions.properties[("speciesActivityLn", "CO3-2")]
        .stoichiometricCoeff
        == 1
    )
