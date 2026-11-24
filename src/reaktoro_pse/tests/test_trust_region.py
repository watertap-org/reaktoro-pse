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
"""tests for trust region backend (ReaktoroSurrogateBlock + trust region solver)"""
import numpy as np
import pytest
import cyipopt
from pyomo.environ import Objective, assert_optimal_termination
from idaes.core.util.model_statistics import degrees_of_freedom

from reaktoro_pse.examples import simple_desalination, thermal_precipitation
from reaktoro_pse.parallel_tools.reaktoro_block_manager import ReaktoroBlockManager
from reaktoro_pse.core.reaktoro_surrogate import ReaktoroSurrogateBlockData
from reaktoro_pse.core.util_classes.cyipopt_solver import (
    get_cyipopt_watertap_solver,
    get_reaktoro_solver,
    ReaktoroSolverWrapper,
)
from reaktoro_pse.core.util_classes.trust_region_solver import (
    get_trust_region_stats,
    find_surrogate_blocks,
)
from reaktoro_pse.reaktoro_block_config.trust_region_options import SurrogateTypes


def check_desal_results(m):
    assert (
        pytest.approx(m.desal_properties[("scalingTendency", "Gypsum")].value, 1e-3)
        == 0.604051223942643
    )
    assert (
        pytest.approx(m.desal_properties[("osmoticPressure", "H2O")].value, 1e-1)
        == 1548396.415543
    )
    assert pytest.approx(m.desal_properties[("pH", None)].value, 1e-3) == 6.284055
    assert pytest.approx(m.water_recovery.value, 1e-3) == 0.899999
    assert pytest.approx(m.acid_addition.value, 1e-3) == 0.003043


def build_desal(
    parallel_mode=False, use_parallel_mode=True, hess_type=None, monkeypatch=None
):
    if parallel_mode:
        monkeypatch.setattr(
            simple_desalination,
            "ReaktoroBlockManager",
            lambda **kw: ReaktoroBlockManager(
                use_parallel_mode=use_parallel_mode, **kw
            ),
        )
    m = simple_desalination.build_simple_desal(
        hess_type, parallel_mode=parallel_mode, backend="trust_region"
    )
    simple_desalination.initialize(m)
    simple_desalination.setup_optimization(m)
    return m


def surrogate_block(m):
    if m.parallel_block_manager is not None:
        return m.parallel_block_manager.reaktoro_model
    return m.eq_desal_properties.reaktoro_model


@pytest.mark.parametrize(
    "hess_type, surrogate",
    [
        (None, SurrogateTypes.quadratic),  # SR1 default
        ("LBFGS", SurrogateTypes.quadratic),
        ("BFGS", SurrogateTypes.quadratic),
        ("ZeroHessian", SurrogateTypes.linear),
    ],
)
def test_desal_trust_region(hess_type, surrogate):
    m = simple_desalination.main(hess_type=hess_type, backend="trust_region")
    blk = m.eq_desal_properties.reaktoro_model
    assert isinstance(blk, ReaktoroSurrogateBlockData)
    assert blk.surrogate_type == surrogate
    check_desal_results(m)


def test_desal_trust_region_stats_and_cleanup():
    m = build_desal()
    blk = surrogate_block(m)
    solver = get_cyipopt_watertap_solver()
    assert isinstance(solver, ReaktoroSolverWrapper)
    result = solver.solve(m)
    assert_optimal_termination(result)
    check_desal_results(m)

    stats = get_trust_region_stats(result)
    assert stats is solver.trust_region_stats
    assert stats["termination_condition"] == "optimal"
    assert stats["outer_iterations"] == len(stats["iterations"]) - 1
    assert (
        stats["accepted_steps"] + stats["rejected_steps"] == stats["outer_iterations"]
    )
    assert stats["inner_solves"] == stats["outer_iterations"]
    assert stats["inner_iterations_total"] == sum(stats["inner_iterations"])
    # one reaktoro solve at start and one per trial point
    assert stats["reaktoro_solves"] == stats["outer_iterations"] + 1
    assert stats["reaktoro_failures"] == 0
    for key in ("feasibility", "trust_radius", "step_norm", "rejected"):
        assert key in stats["iterations"][0]
    assert "inner_status" in stats["iterations"][1]

    # trust region components are removed, user objective restored
    assert not blk.elastic_mode
    assert all(
        "tr_" not in c.local_name for c in blk.component_objects(descend_into=False)
    )
    objs = list(m.component_data_objects(Objective, active=True))
    assert objs == [m.objective]
    assert find_surrogate_blocks(m) == [blk]


def test_graybox_has_no_trust_region_stats():
    m = simple_desalination.build_simple_desal(None)
    simple_desalination.initialize(m)
    solver = get_reaktoro_solver()
    result = solver.solve(m)
    assert_optimal_termination(result)
    assert get_trust_region_stats(result) is None
    assert solver.trust_region_stats is None


def test_degrees_of_freedom_same_as_graybox():
    m_gb = simple_desalination.build_simple_desal(None, backend="graybox")
    m_tr = simple_desalination.build_simple_desal(None, backend="trust_region")
    assert degrees_of_freedom(m_tr) == degrees_of_freedom(m_gb)
    simple_desalination.setup_optimization(m_gb)
    simple_desalination.setup_optimization(m_tr)
    assert degrees_of_freedom(m_tr) == degrees_of_freedom(m_gb)


def test_surrogate_matches_truth_at_expansion_point():
    m = build_desal()
    blk = surrogate_block(m)
    get_cyipopt_watertap_solver().solve(m)
    w0 = blk.expansion_point
    y, J = blk.evaluate_truth(w0)
    np.testing.assert_allclose(blk.surrogate_outputs(w0), y, rtol=1e-10)
    # gradient of the surrogate at w0 is the reaktoro jacobian
    for i in range(len(w0)):
        h = 1e-6 * max(abs(w0[i]), 1e-12)
        wp, wm = w0.copy(), w0.copy()
        wp[i] += h
        wm[i] -= h
        fd = (blk.surrogate_outputs(wp) - blk.surrogate_outputs(wm)) / (2 * h)
        np.testing.assert_allclose(
            fd, J[:, i], rtol=1e-5, atol=1e-8 * np.max(np.abs(J))
        )


def test_reaktoro_failure_at_trial_point_rejects_step():
    m = build_desal()
    blk = surrogate_block(m)
    original = blk.evaluate_truth
    calls = {"n": 0}

    def failing(w=None):
        calls["n"] += 1
        if calls["n"] == 2:  # first trial point
            raise cyipopt.CyIpoptEvaluationError()
        return original(w)

    blk.evaluate_truth = failing
    try:
        result = get_cyipopt_watertap_solver().solve(m)
    finally:
        del blk.evaluate_truth
    assert_optimal_termination(result)
    stats = get_trust_region_stats(result)
    assert stats["reaktoro_failures"] == 1
    assert stats["rejected_steps"] >= 1
    assert stats["iterations"][1]["rejected"]
    check_desal_results(m)


@pytest.mark.parametrize("mode", ["fixed", "adaptive"])
def test_elastic_penalty_modes(mode):
    m = build_desal()
    surrogate_block(m).trust_region_options.elastic_penalty_mode = mode
    result = get_cyipopt_watertap_solver().solve(m)
    assert_optimal_termination(result)
    check_desal_results(m)


@pytest.mark.parametrize("keep", [False, True])
def test_keep_surrogate_curvature(keep):
    m = build_desal()
    blk = surrogate_block(m)
    blk.trust_region_options.keep_surrogate_curvature = keep
    # only count reset at start of solve, not resets after large steps
    blk.trust_region_options.curvature_reset_step = 1e10
    calls = {"n": 0}
    original = blk.reset_curvature

    def counting():
        calls["n"] += 1
        return original()

    blk.reset_curvature = counting
    try:
        result = get_cyipopt_watertap_solver().solve(m)
    finally:
        del blk.reset_curvature
    assert_optimal_termination(result)
    assert calls["n"] == (0 if keep else 1)


def test_backend_mismatch_with_manager_raises(monkeypatch):
    monkeypatch.setattr(
        simple_desalination,
        "ReaktoroBlockManager",
        lambda **kw: ReaktoroBlockManager(use_parallel_mode=False),
    )
    with pytest.raises(ValueError, match="differs from ReaktoroBlockManager backend"):
        simple_desalination.build_simple_desal(
            None, parallel_mode=True, backend="trust_region"
        )


@pytest.mark.parametrize("use_parallel_mode", [False, True])
def test_desal_trust_region_manager(monkeypatch, use_parallel_mode):
    m = build_desal(
        parallel_mode=True, use_parallel_mode=use_parallel_mode, monkeypatch=monkeypatch
    )
    try:
        blk = surrogate_block(m)
        assert isinstance(blk, ReaktoroSurrogateBlockData)
        result = get_cyipopt_watertap_solver().solve(m)
        assert_optimal_termination(result)
        check_desal_results(m)
        assert get_trust_region_stats(result)["number_of_surrogate_blocks"] == 1
    finally:
        m.parallel_block_manager.terminate_workers()


@pytest.mark.parametrize("parallel_mode", [False, True])
def test_thermal_precipitation_trust_region(parallel_mode):
    m = thermal_precipitation.main(
        hess_type="BFGS_mod", backend="trust_region", parallel_mode=parallel_mode
    )
    try:
        assert (
            pytest.approx(
                m.precipitation_properties[("speciesAmount", "Calcite")].value, 1e-2
            )
            == 8.38746210798396e-05
        )
        assert (
            pytest.approx(
                m.precipitation_properties[("vaporPressure", "H2O(g)")].value, 1e-3
            )
            == 12166.985206776764
        )
        assert (
            pytest.approx(m.precipitation_properties[("pH", None)].value, 1e-3)
            == 6.965357182573612
        )
        assert pytest.approx(m.Q_heating.value, abs=4e4) == 123.99859118830717 * 1000
    finally:
        if m.parallel_block_manager is not None:
            m.parallel_block_manager.terminate_workers()
