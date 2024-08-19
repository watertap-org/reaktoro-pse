import pytest
import reaktoro_pse.core.util_classes.rktInputs as rktInputs
from pyomo.environ import Var, units as pyunits

__author__ = "Alexander V. Dudchenko (SLAC)"


@pytest.fixture
def build_inputs():
    test_vars = {}
    test_vars["temperature"] = Var(initialize=293.15, units=pyunits.K)
    test_vars["temperature"].construct()

    test_vars["H2O"] = Var(initialize=1, units=pyunits.kg / pyunits.s)
    test_vars["H2O"].construct()

    return test_vars


def test_rktInput(build_inputs):
    inputs = build_inputs

    rkt_inputs = rktInputs.rktInputs()
    for key, var in inputs.items():
        rkt_inputs[key] = var

    assert rkt_inputs["temperature"].value == 293.15
    assert rkt_inputs["temperature"].timeUnit == None
    assert rkt_inputs["temperature"].mainUnit == "K"

    assert rkt_inputs["H2O"].value == 1
    assert rkt_inputs["H2O"].timeUnit == "s"
    assert rkt_inputs["H2O"].mainUnit == "kg"

    with pytest.raises(TypeError) as a:
        rkt_inputs["bad_var"] = 1
