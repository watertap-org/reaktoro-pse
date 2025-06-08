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
    biogas_combustion,
    simple_desalination,
    simple_ion_exchange,
    thermal_precipitation,
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
    m = thermal_precipitation.main()
    assert (
        pytest.approx(
            m.precipitation_properties[("speciesAmount", "Calcite")].value, 1e-2
        )
        == 0.0005126288369679213
    )
    assert (
        pytest.approx(
            m.precipitation_properties[("vaporPressure", "H2O(g)")].value, 1e-3
        )
        == 12162.679
    )
    assert (
        pytest.approx(m.precipitation_properties[("pH", None)].value, 1e-3)
        == 6.937058009543962
    )
    assert pytest.approx(m.Q_heating.value, abs=4e4) == 127043.64015901716
    assert pytest.approx(m.precipitator_temperature.value, 1e-3) == 273.15 + 50


def test_ion_exchange():
    m = simple_ion_exchange.main()

    assert pytest.approx(m.removal_percent["Mg"].value, 1e-1) == -35.54130924283
    assert pytest.approx(m.removal_percent["Ca"].value, 1e-1) == -79.15299911033
    assert pytest.approx(m.treated_pH.value, 1e-2) == 13.374349619456911
    assert pytest.approx(m.base_addition.value, abs=1e-1) == 0.31567192053040094


def test_biogas():
    m = biogas_combustion.main()

    assert pytest.approx(m.air_to_fuel_ratio.value, 1e-1) == 3.8751662012681587
    assert pytest.approx(m.exhaust_temperature.value, 1e-1) == 2000
