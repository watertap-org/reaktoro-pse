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
from reaktoro_pse.examples.reaktoro_pse_to_phreeqc_comparison import (
    water_removal_comparison,
    vapor_pressure_comparison,
    precipitation_comparison,
    acid_base_addition_comparison,
    solution_mixing_comparison,
)


def test_water_vapor_comp():
    result = vapor_pressure_comparison.main(False, False)
    print(result)

    expected_result = {"Vapor pressure": 0.0008019869790473859}
    for key in result:
        assert pytest.approx(result[key], 1e-3) == expected_result[key]


def test_water_removal_comp():
    result = water_removal_comparison.main(False, False)
    expected_result = {
        "Calcite": 0.25142049795489263,
        "pH": 0.08756483809897087,
        "Osmotic pressure": 0.045095730902661786,
    }
    for key in result:
        assert pytest.approx(result[key], 1e-3) == expected_result[key]
    print(result)


def test_water_removal_comp_wateq4f():
    result = water_removal_comparison.main_wate4qf(False, False)
    expected_result = {
        "Calcite": 0.03128858529075547,
        "pH": 0.0012100923124952804,
        "pE": 0.008019965754384479,
        "Osmotic pressure": 0.02607370346850101,
    }
    for key in result:
        assert pytest.approx(result[key], 1e-3) == expected_result[key]
    print(result)


def test_precip_comp():
    result = precipitation_comparison.main(False, False)
    print(result)
    expected_result = {
        "Calcite": 4.9384302203136643e-08,
        "pH": 0.03764096550592823,
        "formed phase Calcite": 0.12859845852994584,
    }
    for key in result:
        assert pytest.approx(result[key], abs=1e-3) == expected_result[key]


def test_precip_comp_wate4qf():
    result = precipitation_comparison.main_wate4qf(False, False)
    print(result)
    expected_result = {
        "Calcite": 9.398249956049654e-09,
        "pH": 0.00032505223102499676,
        "pE": 0.002168859127643667,
        "formed phase Calcite": 0.0010451446143509827,
    }
    for key in result:
        assert pytest.approx(result[key], abs=1e-3) == expected_result[key]


def test_acid_base_comp_wateq4f():
    result = acid_base_addition_comparison.main_wateq4f(False, False)
    print(result)
    hcl_result = {
        "Calcite": 0.4133777233123765,
        "pH": 0.03130413015969841,
        "pE": 0.08936653585804519,
        "Osmotic pressure": 0.033316159769095396,
    }
    naoh_result = {
        "Calcite": 0.09095827551791304,
        "pH": 0.004364515100135042,
        "pE": 0.4380635586124001,
        "Osmotic pressure": 0.032932650803246095,
    }
    for key in result[0]:
        assert pytest.approx(result[0][key], 1e-3) == hcl_result[key]
    for key in result[1]:
        assert pytest.approx(result[1][key], 1e-3) == naoh_result[key]


def test_acid_base_comp():
    result = acid_base_addition_comparison.main(False, False)
    print(result)
    hcl_result = {
        "Calcite": 0.5976297918466753,
        "pH": 0.007230990243316364,
        "Osmotic pressure": 0.05902701717931409,
    }
    naoh_result = {
        "Calcite": 0.7774901776376242,
        "pH": 0.016994006109543942,
        "Osmotic pressure": 0.05628415918558991,
    }
    for key in result[0]:
        assert pytest.approx(result[0][key], 1e-3) == hcl_result[key]
    for key in result[1]:
        assert pytest.approx(result[1][key], 1e-3) == naoh_result[key]


def test_mix_comp():
    result = solution_mixing_comparison.main(False, False)
    print(result)
    expected_result = {
        "Calcite": 0.6899966986598822,
        "pH": 0.00587538473081846,
        "Osmotic pressure": 0.05964941573467817,
    }
    for key in result:
        assert pytest.approx(result[key], abs=1e-3) == expected_result[key]


def test_mix_comp_wateq4f():
    result = solution_mixing_comparison.main_wateq4f(False, False)
    expected_result = {
        "Calcite": 0.45709052273393186,
        "pH": 0.037751350399219494,
        "pE": 0.6433578729052789,
        "Osmotic pressure": 0.03128788989509176,
    }
    for key in result:
        assert pytest.approx(result[key], abs=1e-3) == expected_result[key]
