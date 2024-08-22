import pytest
from reaktoro_pse.examples import (
    simple_desalination,
    themal_precipitation,
)


def test_desal():
    m = simple_desalination.main()
    assert (
        pytest.approx(m.desal_properties[("scalingTendency", "Gypsum")].value, 1e-3)
        == 0.604051223942643
    )
    assert (
        pytest.approx(m.desal_properties[("osmoticPressure", "H2O")].value, 1e-1)
        == 1548396.415543
    )

    assert pytest.approx(m.desal_properties[("pH", None)].value, 1e-2) == 6.284055
    assert pytest.approx(m.water_recovery.value, 1e-3) == 0.899999
    assert pytest.approx(m.acid_addition.value, 1e-3) == 0.003043


def test_thermal_precipt():
    m = themal_precipitation.main()
    assert (
        pytest.approx(
            m.precipitation_properties[("speciesAmount", "Calcite")].value, 1e-3
        )
        == 0.0009549645
    )
    assert (
        pytest.approx(m.precipitation_properties[("molarEnthalpy", None)].value, 1e-3)
        == -283275.717
    )
    assert (
        pytest.approx(
            m.precipitation_properties[("vaporPressure", "H2O(g)")].value, 1e-3
        )
        == 69560.00913
    )

    assert pytest.approx(m.Q_heating.value, 1e-3) == 165000
    assert pytest.approx(m.Q_recoverable.value, 1e-3) == 82500
    assert pytest.approx(m.precipitator_temperature.value, 1e-3) == 273.15 + 58.273
    assert pytest.approx(m.cooled_treated_temperature.value, 1e-1) == 273.15 + 38.199
