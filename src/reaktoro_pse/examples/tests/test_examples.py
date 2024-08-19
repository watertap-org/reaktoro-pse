import pytest
from reaktoro_pse.examples import (
    simple_desalination,
    themal_precipitation,
)


def test_desal():
    m = simple_desalination.main()
    assert (
        pytest.approx(m.desal_properties[("scalingTendency", "Gypsum")].value, 1e-3)
        == 1
    )
    assert (
        pytest.approx(m.desal_properties[("osmoticPressure", "H2O")].value, 1e-3)
        == 323585.01770930213
    )

    assert pytest.approx(m.desal_properties[("pH", None)].value, 1e-3) == 4.000000
    assert pytest.approx(m.water_recovery.value, 1e-3) == 0.3691153516087952
    assert pytest.approx(m.acid_addition.value, 1e-3) == 0.00824853470415005


def test_thermal_precipt():
    m = themal_precipitation.main()
    assert (
        pytest.approx(
            m.precipitation_properties[("speciesAmount", "Calcite")].value, 1e-3
        )
        == 0.0016553088974723934
    )
    assert (
        pytest.approx(m.precipitation_properties[("molarEnthalpy", None)].value, 1e-3)
        == -280691.8099
    )

    assert pytest.approx(m.Q_heating.value, 1e-3) == 310000
    assert pytest.approx(m.Q_recoverable.value, 1e-3) == 155000
    assert pytest.approx(m.precipitator_temperature.value, 1e-3) == 371.0945
    assert pytest.approx(m.cooled_treated_temperature.value, 1e-1) == 337.95
