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
from reaktoro_pse.examples import (
    simple_desalination,
    themal_precipitation,
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
    m = themal_precipitation.main()
    assert (
        pytest.approx(
            m.precipitation_properties[("speciesAmount", "Calcite")].value, 1e-3
        )
        == 0.0009467002694194435
    )
    assert (
        pytest.approx(
            m.precipitation_properties[("vaporPressure", "H2O(g)")].value, 1e-3
        )
        == 17832.261033720
    )
    assert (
        pytest.approx(m.precipitation_properties[("pH", None)].value, 1e-3)
        == 6.9136846633442275
    )
    assert pytest.approx(m.Q_heating.value, 1e-3) == 165000
    assert pytest.approx(m.Q_recoverable.value, 1e-3) == 82500
    assert (
        pytest.approx(m.precipitator_temperature.value, 1e-3)
        == 273.15 + 57.9465892050192
    )
    assert (
        pytest.approx(m.cooled_treated_temperature.value, 1e-1)
        == 273.15 + 37.870756668942306
    )
