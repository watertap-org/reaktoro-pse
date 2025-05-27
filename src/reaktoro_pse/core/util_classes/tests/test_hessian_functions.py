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
from reaktoro_pse.core.util_classes.hessian_functions import (
    HessianApproximation,
    HessTypes,
)
import numpy as np


def get_jac_and_hessian(func, input_values):
    outputs = np.polynomial.polynomial.polyval2d(input_values[0], input_values[1], func)
    jac, hess = [], []
    for i in range(len(input_values)):
        first_order_derivative = np.polynomial.polynomial.polyder(func, m=1, axis=i)
        second_order_derivative = np.polynomial.polynomial.polyder(func, m=2, axis=i)

        jacobian = np.polynomial.polynomial.polyval2d(
            input_values[0], input_values[1], first_order_derivative
        )
        hessian = np.polynomial.polynomial.polyval2d(
            input_values[0], input_values[1], second_order_derivative
        )
        jac.append([jacobian])
        hess.append([hessian])
    jacobian = np.array(jac)
    hessian = np.array(hess)
    return outputs, jacobian, hessian


@pytest.fixture
def build_inputs():
    test_function = ((2, 3, 1), (1, 2, 3))
    input_values = np.array([1, 2])
    outputs, jacobian, hessian = get_jac_and_hessian(test_function, input_values)

    return input_values, outputs, jacobian, hessian


def test_gauss_newton(build_inputs):
    input_values, outputs, jacobian, hessian = build_inputs
    print(input_values, outputs, jacobian, hessian)
    hess = HessianApproximation(hessian_type=HessTypes.BFGS)
    approx_hessian = hess.get_hessian(input_values, outputs, jacobian)
    print("hessian", hessian)
    print("approx_hessian", approx_hessian)
