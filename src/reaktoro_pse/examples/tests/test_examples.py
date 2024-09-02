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
from reaktoro_pse.examples import (
    simple_desalination,
    simple_ion_exchange,
    thermal_precipitation,
    bigasificaiton_example,
)


def test_desal():
    m, m_open = simple_desalination.main()
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
    for key, obj in m.desal_properties.items():
        assert pytest.approx(obj.value, 1e-3) == m_open.desal_properties[key].value


def test_thermal_precipt():
    m = thermal_precipitation.main()
    assert (
        pytest.approx(
            m.precipitation_properties[("speciesAmount", "Calcite")].value, 1e-3
        )
        == 0.0009467511444760701
    )
    assert (
        pytest.approx(
            m.precipitation_properties[("vaporPressure", "H2O(g)")].value, 1e-3
        )
        == 17801.227149565908
    )
    assert (
        pytest.approx(m.precipitation_properties[("pH", None)].value, 1e-3)
        == 6.913772075650711
    )
    assert pytest.approx(m.Q_heating.value, 1e-3) == 165000
    assert pytest.approx(m.Q_recoverable.value, 1e-3) == 82500
    assert (
        pytest.approx(m.precipitator_temperature.value, 1e-3)
        == 273.15 + 57.90940479948432
    )
    assert (
        pytest.approx(m.cooled_treated_temperature.value, 1e-1)
        == 273.15 + 37.8773655288286
    )


def test_ion_exchange():
    m = simple_ion_exchange.main()

    assert pytest.approx(m.removal_percent["Mg"].value, 1e-1) == -5.000000005534077
    assert pytest.approx(m.removal_percent["Ca"].value, 1e-1) == -10.800409844332075
    assert pytest.approx(m.removal_percent["Na"].value, 1e-1) == 126.94292855983791
    assert pytest.approx(m.treated_pH.value, 1e-2) == 12.57850184575838
    assert pytest.approx(m.base_addition.value, 1e-2) == 0.09050676244974129
    assert pytest.approx(m.acid_addition.value, 1e-2) == 8.32016553459784e-11


def test_biogas():
    m = bigasificaiton_example.main()

    assert pytest.approx(m.air_to_fuel_ratio.value, 1e-1) == 3.8751662012681587
    assert pytest.approx(m.exhaust_temperature.value, 1e-1) == 2000
