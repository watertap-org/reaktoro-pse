import pytest
import pyomo.environ as pyo


@pytest.mark.parametrize(
    "solver",
    [
        "cyipopt",
    ],
)
def test_solver_available(solver: str):
    assert pyo.SolverFactory(solver).available()
