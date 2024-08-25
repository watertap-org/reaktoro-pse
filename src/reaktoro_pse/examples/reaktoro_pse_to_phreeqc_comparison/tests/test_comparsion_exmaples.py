import pytest
from reaktoro_pse.examples.reaktoro_pse_to_phreeqc_comparison import (
    water_removal_comparison,
    vapor_pressure_comparison,
    precipitation_comparison,
    acid_base_addition_comparison,
)


def test_water_vapor_comp():
    result = vapor_pressure_comparison.main(True, False)
    print(result)

    expected_result = {"Vapor pressure": 9.243111447147829}
    for key in result:
        assert pytest.approx(result[key], 1e-3) == expected_result[key]


def test_water_removal_comp():
    result = water_removal_comparison.main(True, False)
    expected_result = {
        "Calcite": 0.2697638824647954,
        "pH": 0.08756483809897087,
        "Osmotic pressure": 0.045095730902661786,
    }
    for key in result:
        assert pytest.approx(result[key], 1e-3) == expected_result[key]
    print(result)


def test_precip_comp():
    result = precipitation_comparison.main(True, False)
    print(result)
    expected_result = {
        "Calcite": 4.9384302203136643e-08,
        "pH": 0.03764096550592823,
        "formed phase Calcite": 0.12859845852994584,
    }
    for key in result:
        assert pytest.approx(result[key], abs=1e-3) == expected_result[key]


def test_acid_base_comp():
    result = acid_base_addition_comparison.main(True, False)
    hcl_result = {
        "Calcite": 0.6170240328965659,
        "pH": 0.007244657764840176,
        "Osmotic pressure": 0.05902701717931409,
    }
    naoh_result = {
        "Calcite": 0.7957449866927796,
        "pH": 0.016994006109543942,
        "Osmotic pressure": 0.05628415918558991,
    }
    for key in result[0]:
        assert pytest.approx(result[0][key], 1e-3) == hcl_result[key]
    for key in result[1]:
        assert pytest.approx(result[1][key], 1e-3) == naoh_result[key]
