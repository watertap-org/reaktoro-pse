#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/reaktoro-pse/"
#################################################################################
import numpy as np
import pytest
from pyomo.environ import ConcreteModel, Var
import idaes.core.util.scaling as iscale

from reaktoro_pse.core.reaktoro_gray_box import (
    ReaktoroGrayBox,
    SpeciesAmountOffsetTypes,
    floor_species_amount,
    floor_species_amount_derivative,
)


def _gray_box_with_floors(offset_type, sf=1e3, offset=1e-4):
    """gray box with 3 outputs, output 0 and 2 are floored"""
    m = ConcreteModel()
    m.v = Var(range(3), initialize=1)
    for i in range(3):
        iscale.set_scaling_factor(m.v[i], sf)
    gb = ReaktoroGrayBox.__new__(ReaktoroGrayBox)
    gb.outputs = ["a", "b", "c"]
    gb._output_floors = {}
    gb.register_output_floor("a", m.v[0], offset, offset_type)
    gb.register_output_floor("c", m.v[2], offset, offset_type)
    return m, gb


@pytest.mark.parametrize(
    "offset_type",
    [SpeciesAmountOffsetTypes.additive, SpeciesAmountOffsetTypes.smooth_floor],
)
def test_floor_derivative_matches_finite_difference(offset_type):
    a = 1e-7
    for n in [1e-40, 1e-9, 1e-7, 1e-5, 1.0]:
        h = max(abs(n), a) * 1e-6
        fd = (
            floor_species_amount(n + h, a, offset_type)
            - floor_species_amount(n - h, a, offset_type)
        ) / (2 * h)
        assert floor_species_amount_derivative(n, a, offset_type) == pytest.approx(
            fd, rel=1e-6, abs=1e-12
        )


@pytest.mark.parametrize(
    "offset_type",
    [SpeciesAmountOffsetTypes.additive, SpeciesAmountOffsetTypes.smooth_floor],
)
def test_apply_output_floors(offset_type):
    m, gb = _gray_box_with_floors(offset_type)
    a = 1e-4 / 1e3
    outputs = np.array([1e-40, 2.0, 0.5])
    jac = np.array([[0.0, 0.0], [1.0, 2.0], [3.0, 4.0]])
    raw_out, raw_jac = outputs.copy(), jac.copy()
    new_jac, new_out = gb.apply_output_floors(jac, outputs)
    # raw arrays are not modified
    np.testing.assert_array_equal(outputs, raw_out)
    np.testing.assert_array_equal(jac, raw_jac)
    # absent species sits at a
    assert new_out[0] == pytest.approx(a, rel=1e-12)
    # unfloored output untouched
    assert new_out[1] == 2.0
    np.testing.assert_array_equal(new_jac[1], jac[1])
    expected = floor_species_amount(0.5, a, offset_type)
    assert new_out[2] == pytest.approx(expected, rel=1e-14)
    d = floor_species_amount_derivative(0.5, a, offset_type)
    np.testing.assert_allclose(new_jac[2], d * jac[2], rtol=1e-14)
    if offset_type == SpeciesAmountOffsetTypes.smooth_floor:
        # present species shifted by ~a**2/(2n)
        assert new_out[2] - 0.5 == pytest.approx(a**2 / (2 * 0.5), rel=1e-6)


def test_output_floor_follows_scaling_changes():
    m, gb = _gray_box_with_floors(SpeciesAmountOffsetTypes.smooth_floor)
    outputs = np.array([1e-40, 2.0, 0.5])
    jac = np.zeros((3, 2))
    _, out = gb.apply_output_floors(jac, outputs)
    assert out[0] == pytest.approx(1e-7, rel=1e-12)
    iscale.set_scaling_factor(m.v[0], 1e5)
    _, out = gb.apply_output_floors(jac, outputs)
    assert out[0] == pytest.approx(1e-9, rel=1e-12)


def test_no_floors_returns_inputs():
    gb = ReaktoroGrayBox.__new__(ReaktoroGrayBox)
    gb._output_floors = {}
    jac, out = np.ones((2, 2)), [1.0, 2.0]
    assert gb.apply_output_floors(jac, out) == (jac, out)


def test_register_output_floor_rejects_none_type():
    m, gb = _gray_box_with_floors(SpeciesAmountOffsetTypes.additive)
    with pytest.raises(ValueError):
        gb.register_output_floor("b", m.v[1], 1e-4, None)
