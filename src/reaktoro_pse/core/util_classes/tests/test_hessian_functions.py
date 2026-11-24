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
from scipy.sparse import coo_matrix

from reaktoro_pse.core.util_classes.hessian_functions import (
    HessianApproximation,
    HessTypes,
    QuasiNewtonInitializationTypes,
    _nonzero,
)

# single output y = 0.5 x^T A x, so jacobian row is (A x)^T and true hessian is A
A = np.array([[4.0, 1.0, 0.5], [1.0, 3.0, 0.2], [0.5, 0.2, 2.0]])
POINTS = [
    np.array([1.0, 2.0, 3.0]),
    np.array([1.3, 1.7, 3.2]),
    np.array([0.9, 2.4, 2.6]),
    np.array([1.5, 2.1, 3.4]),
    np.array([1.1, 1.6, 2.9]),
]


def jac(x):
    return np.array([A @ x])


def full(hess):
    """rebuild symmetric matrix from lower triangular coo"""
    h = hess.toarray()
    return h + h.T - np.diag(np.diag(h))


def run_points(hess_calc, points, duals=np.array([1.0])):
    hess = None
    for x in points:
        hess = hess_calc.get_hessian(x, jac(x), duals)
    return hess


def test_nonzero():
    assert _nonzero([1.0, -1.0])
    assert not _nonzero([0.0, 0.0])
    assert _nonzero(np.array([[0.0], [1e-300]]))


def test_check_step_detects_zero_sum_step():
    hc = HessianApproximation(hessian_type=HessTypes.BFGS)
    hc.hessian_memory.memorize(np.array([1.0, 1.0]), np.zeros((1, 2)))
    hc.hessian_memory.memorize(np.array([2.0, 0.0]), np.zeros((1, 2)))
    # s = [-1, 1] sums to zero but is a valid step
    assert hc.check_step()


@pytest.mark.parametrize(
    "init_type",
    [
        QuasiNewtonInitializationTypes.scalar1,
        QuasiNewtonInitializationTypes.scalar2,
        QuasiNewtonInitializationTypes.scalar3,
        QuasiNewtonInitializationTypes.scalar4,
        QuasiNewtonInitializationTypes.GaussNewton,
        QuasiNewtonInitializationTypes.constant,
    ],
)
def test_initializations_return_hessian(init_type):
    hc = HessianApproximation(
        hessian_type=HessTypes.LBFGS, quasi_newton_initialization_type=init_type
    )
    hess = run_points(hc, POINTS[:3])
    assert isinstance(hess, coo_matrix)
    assert np.all(np.isfinite(hess.toarray()))


@pytest.mark.parametrize(
    "hess_type",
    [
        HessTypes.BFGS,
        HessTypes.CBFGS,
        HessTypes.BFGS_mod,
        HessTypes.BFGS_damp,
        HessTypes.BFGS_ipopt,
        HessTypes.LBFGS,
    ],
)
def test_repeated_point_does_not_change_hessian(hess_type):
    hc = HessianApproximation(hessian_type=hess_type, lbfgs_memory=6)
    first = run_points(hc, POINTS[:3]).toarray()
    n_mem = len(hc.hessian_memory.inputs)
    # re-evaluate at the same point, e.g. with ipopt asking again
    second = hc.get_hessian(POINTS[2], jac(POINTS[2]), np.array([1.0])).toarray()
    assert len(hc.hessian_memory.inputs) == n_mem
    np.testing.assert_allclose(second, first)


def test_repeated_point_rescales_with_new_duals():
    hc = HessianApproximation(hessian_type=HessTypes.BFGS)
    first = run_points(hc, POINTS[:3]).toarray()
    second = hc.get_hessian(POINTS[2], jac(POINTS[2]), np.array([2.0])).toarray()
    np.testing.assert_allclose(second, 2 * first)


@pytest.mark.parametrize(
    "init_type",
    [
        QuasiNewtonInitializationTypes.scalar1,
        QuasiNewtonInitializationTypes.GaussNewton,
    ],
)
def test_lbfgs_satisfies_latest_secant(init_type):
    """with a consistent B0 the unrolled update must give B s_k = y_k for the latest pair"""
    hc = HessianApproximation(
        hessian_type=HessTypes.LBFGS,
        lbfgs_memory=6,
        quasi_newton_initialization_type=init_type,
    )
    hess = full(run_points(hc, POINTS))
    s = POINTS[-1] - POINTS[-2]
    y = A @ s
    np.testing.assert_allclose(hess @ s, y, rtol=1e-8)


def test_failure_falls_back_to_zero_hessian():
    hc = HessianApproximation(hessian_type=HessTypes.BFGS)
    good = hc.get_hessian(POINTS[0], jac(POINTS[0]), np.array([1.0]))

    def fail():
        raise ValueError("boom")

    hc.quasi_newton_hessian = fail
    hess = hc.get_hessian(POINTS[1], jac(POINTS[1]), np.array([1.0]))
    assert isinstance(hess, coo_matrix)
    assert hess.shape == (3, 3)
    assert not np.any(hess.toarray())
    # pynumero requires identical structure on every call
    np.testing.assert_array_equal(hess.row, good.row)
    np.testing.assert_array_equal(hess.col, good.col)


def test_unknown_hessian_type_raises():
    hc = HessianApproximation(hessian_type="not_a_type")
    with pytest.raises(NotImplementedError):
        hc.get_hessian(POINTS[0], jac(POINTS[0]), np.array([1.0]))


################################################################################
# SR1 / PSB / support / scaling (milestone 1)

# 4 points -> 3 linearly independent steps in 3 dimensions
SR1_POINTS = POINTS[:4]
# indefinite matrix (one negative eigenvalue)
A_INDEF = np.array([[2.0, 1.0, 0.0], [1.0, -3.0, 0.5], [0.0, 0.5, 1.0]])


def run_quadratic(hc, points, mat, duals=np.array([1.0])):
    hess = None
    for x in points:
        hess = hc.get_hessian(x, np.array([mat @ x]), duals)
    return hess


@pytest.mark.parametrize("mat", [A, A_INDEF], ids=["definite", "indefinite"])
def test_sr1_recovers_quadratic_hessian(mat):
    """SR1 terminates finitely on quadratics: after n linearly independent steps
    without skipped updates the matrix equals the true hessian (Nocedal & Wright thm 6.1)
    """
    hc = HessianApproximation(
        hessian_type=HessTypes.SR1,
        quasi_newton_initialization_type=QuasiNewtonInitializationTypes.zero,
    )
    hess = full(run_quadratic(hc, SR1_POINTS, mat))
    np.testing.assert_allclose(hess, mat, rtol=1e-8, atol=1e-10)


def test_bfgs_cannot_represent_indefinite_curvature():
    """reason for SR1/PSB, BFGS keeps matrix positive definite and gets sign wrong"""
    hc = HessianApproximation(hessian_type=HessTypes.BFGS_ipopt)
    hess = full(run_quadratic(hc, SR1_POINTS, A_INDEF))
    assert np.min(np.linalg.eigvalsh(hess)) > -1e-12
    assert np.min(np.linalg.eigvalsh(A_INDEF)) < 0


@pytest.mark.parametrize("mat", [A, A_INDEF], ids=["definite", "indefinite"])
@pytest.mark.parametrize(
    "init_type",
    [QuasiNewtonInitializationTypes.zero, QuasiNewtonInitializationTypes.scalar1],
)
def test_psb_satisfies_secant(mat, init_type):
    hc = HessianApproximation(
        hessian_type=HessTypes.PSB, quasi_newton_initialization_type=init_type
    )
    hess = full(run_quadratic(hc, POINTS, mat))
    s = POINTS[-1] - POINTS[-2]
    np.testing.assert_allclose(hess @ s, mat @ s, rtol=1e-8)
    np.testing.assert_allclose(hess, hess.T)


def test_sr1_skip_rule():
    hc = HessianApproximation(hessian_type=HessTypes.SR1, sr1_skip_tolerance=1e-8)
    B = np.identity(2)
    s = np.array([1.0, 0.0])
    # r = y - B s = [0, 1] is orthogonal to s, update undefined -> skipped
    assert hc._rule_sr1(B, s, np.array([1.0, 1.0]), None) is None
    # r = 0, secant already satisfied -> nothing to do
    assert hc._rule_sr1(B, s, np.array([1.0, 0.0]), None) is None
    # regular update is applied and satisfies secant
    y = np.array([3.0, 1.0])
    new_B = hc._rule_sr1(B, s, y, None)
    np.testing.assert_allclose(new_B @ s, y)


def two_block_jac(x):
    """output 0 depends on inputs 0, 1 and output 1 on inputs 2, 3"""
    J = np.zeros((2, 4))
    J[0, :2] = A[:2, :2] @ x[:2]
    J[1, 2:] = A_INDEF[:2, :2] @ x[2:]
    return J


TWO_BLOCK_POINTS = [
    np.array([1.0, 2.0, 3.0, 4.0]),
    np.array([1.2, 1.9, 3.3, 3.8]),
    np.array([0.9, 2.3, 2.8, 4.1]),
    np.array([1.1, 2.2, 3.1, 4.3]),
]


@pytest.mark.parametrize(
    "hess_type",
    [
        HessTypes.BFGS,
        HessTypes.BFGS_ipopt,
        HessTypes.LBFGS,
        HessTypes.SR1,
        HessTypes.PSB,
        HessTypes.SR1_lagrangian,
        HessTypes.GaussNewton,
    ],
)
def test_block_support_has_no_cross_block_entries(hess_type):
    support = [[0, 1], [2, 3]]
    hc = HessianApproximation(hessian_type=hess_type, hessian_support=support)
    duals = np.array([1.5, -0.5])
    structures = []
    for x in TWO_BLOCK_POINTS:
        hess = hc.get_hessian(x, two_block_jac(x), duals)
        structures.append((hess.row.copy(), hess.col.copy()))
        dense = full(hess)
        assert not np.any(dense[:2, 2:])
        assert not np.any(dense[2:, :2])
    # pattern is lower triangle of each block and never changes
    for rows, cols in structures:
        np.testing.assert_array_equal(rows, structures[0][0])
        np.testing.assert_array_equal(cols, structures[0][1])
    assert np.all(structures[0][0] >= structures[0][1])
    assert len(structures[0][0]) == 6


def test_full_support_without_cross_block_restriction_couples_blocks():
    """without support restriction, BFGS adds curvature between independent blocks
    when initial matrix is sigma*I (GaussNewton initialization, diag(g^2), is zero
    outside of the block, so it does not couple blocks)"""
    hc = HessianApproximation(
        hessian_type=HessTypes.BFGS,
        quasi_newton_initialization_type=QuasiNewtonInitializationTypes.scalar1,
    )
    duals = np.array([1.5, -0.5])
    for x in TWO_BLOCK_POINTS:
        dense = full(hc.get_hessian(x, two_block_jac(x), duals))
    assert np.any(dense[:2, 2:])


def test_support_validation():
    with pytest.raises(ValueError):
        HessianApproximation(
            hessian_type=HessTypes.SR1, hessian_support=[[0, 1]]
        ).get_hessian(TWO_BLOCK_POINTS[0], two_block_jac(TWO_BLOCK_POINTS[0]), [1, 1])
    with pytest.raises(ValueError):
        HessianApproximation(
            hessian_type=HessTypes.SR1, hessian_support=[[0, 1], [2, 4]]
        ).get_hessian(TWO_BLOCK_POINTS[0], two_block_jac(TWO_BLOCK_POINTS[0]), [1, 1])


def test_scaling_invariance_of_sr1():
    """SR1 is invariant to linear change of variables when the initial matrix is,
    so result in unscaled space must not depend on the scaling factors
    (safeguards compare norms, which are not scaling invariant, so they are off)"""
    results = []
    for scaling in [None, [1e3, 1.0, 1e-3]]:
        hc = HessianApproximation(
            hessian_type=HessTypes.SR1,
            quasi_newton_initialization_type=QuasiNewtonInitializationTypes.zero,
            quasi_newton_update_cap=0,
            quasi_newton_norm_cap=0,
            input_scaling=scaling,
        )
        results.append(full(run_quadratic(hc, SR1_POINTS, A_INDEF)))
    np.testing.assert_allclose(results[0], results[1], rtol=1e-8, atol=1e-10)


def test_scaling_changes_initial_matrix():
    """initial matrix sigma*I is built in scaled space, so becomes sigma*D^2 unscaled"""
    scaling = np.array([10.0, 1.0, 0.1])
    hc = HessianApproximation(
        hessian_type=HessTypes.SR1,
        quasi_newton_initialization_type=QuasiNewtonInitializationTypes.constant,
        quasi_newton_init_const_value=1.0,
        input_scaling=scaling,
    )
    hess = full(hc.get_hessian(POINTS[0], jac(POINTS[0]), np.array([1.0])))
    np.testing.assert_allclose(hess, np.diag(scaling**2))


def test_input_scaling_cannot_change_after_start():
    hc = HessianApproximation(hessian_type=HessTypes.SR1)
    run_points(hc, POINTS[:2])
    with pytest.raises(RuntimeError):
        hc.set_input_scaling([1.0, 1.0, 1.0])


def test_non_positive_scaling_replaced():
    hc = HessianApproximation(hessian_type=HessTypes.SR1, input_scaling=[0, -1, 2])
    np.testing.assert_array_equal(hc.input_scaling, [1, 1, 2])


def test_sr1_lagrangian_recovers_lagrangian_hessian():
    """one matrix per block, approximating sum_i lambda_i A_i directly"""
    mats = [A, A_INDEF]
    duals = np.array([0.7, -1.3])
    hc = HessianApproximation(
        hessian_type=HessTypes.SR1_lagrangian,
        quasi_newton_initialization_type=QuasiNewtonInitializationTypes.zero,
    )
    assert hc.hessian_matrix_type == HessTypes.SR1_lagrangian
    for x in SR1_POINTS:
        hess = hc.get_hessian(x, np.array([m @ x for m in mats]), duals)
    assert len(hc.elements) == 1
    np.testing.assert_allclose(
        full(hess), duals[0] * A + duals[1] * A_INDEF, rtol=1e-8, atol=1e-10
    )


def test_sr1_lagrangian_groups_by_support():
    hc = HessianApproximation(
        hessian_type=HessTypes.SR1_lagrangian,
        hessian_support=[[0, 1], [0, 1], [2, 3]],
    )
    x = TWO_BLOCK_POINTS[0]
    J = np.vstack([two_block_jac(x)[0], two_block_jac(x)])
    hc.get_hessian(x, J, np.ones(3))
    assert [el.outputs for el in hc.elements] == [[0, 1], [2]]


def test_non_finite_matrix_is_reset():
    hc = HessianApproximation(hessian_type=HessTypes.SR1)
    hc._update_rules[HessTypes.SR1] = lambda B, s, y, g: B * np.nan
    hess = run_points(hc, POINTS[:3])
    assert np.all(np.isfinite(hess.toarray()))
    assert hc.bfgs_matrix_not_initialized


def test_defaults_are_consistent():
    """HessianOptions, ReaktoroSolver and ReaktoroGrayBox share HESSIAN_DEFAULTS"""
    import inspect
    from reaktoro_pse.core.util_classes.hessian_functions import HESSIAN_DEFAULTS
    from reaktoro_pse.reaktoro_block_config.hessian_options import HessianOptions
    from reaktoro_pse.core.reaktoro_solver import ReaktoroSolver
    from reaktoro_pse.core.reaktoro_gray_box import ReaktoroGrayBox

    config = HessianOptions().get_dict()
    solver_defaults = inspect.signature(ReaktoroSolver.set_solver_options).parameters
    gray_box_defaults = inspect.signature(ReaktoroGrayBox.configure).parameters
    for key, value in HESSIAN_DEFAULTS.items():
        assert config[key] == value
        assert solver_defaults[key].default == value
        if key != "hessian_type":
            assert gray_box_defaults[key].default == value


def test_aggregate_solver_hessian_support():
    from reaktoro_pse.parallel_tools.reaktoro_block_manager import (
        AggregateSolverState,
    )

    agg = AggregateSolverState(parallel_mode=False)

    class _Input:
        pass

    for blk, (n_in, n_out) in enumerate([(2, 1), (3, 2)]):
        for i in range(n_in):
            agg.register_input(blk, f"in{i}", _Input())
        for o in range(n_out):
            agg.register_output(blk, f"out{o}")
    support = agg.get_hessian_support()
    assert [list(s) for s in support] == [[0, 1], [2, 3, 4], [2, 3, 4]]


###############################################################################
# deprecated bfgs_* option names


def test_deprecated_option_names_hessian_approximation():
    from reaktoro_pse.core.util_classes.hessian_functions import (
        BFGSInitializationTypes,
        QuasiNewtonInitializationTypes,
    )

    assert BFGSInitializationTypes is QuasiNewtonInitializationTypes
    hess = HessianApproximation(
        hessian_type=HessTypes.SR1,
        bfgs_initialization_type="zero",
        bfgs_hessian_memory=3,
        bfgs_epsilon=1e-12,
    )
    assert hess.quasi_newton_initialization_type == "zero"
    assert hess.hessian_memory.memory == 3
    assert hess.epsilon == 1e-12
    with pytest.raises(ValueError):
        HessianApproximation(
            bfgs_initialization_type="zero",
            quasi_newton_initialization_type="scalar1",
        )


def test_deprecated_option_names_config_and_solver():
    from reaktoro_pse.reaktoro_block_config.hessian_options import (
        HessianOptions,
        resolve_deprecated_hessian_options,
    )
    from reaktoro_pse.core.reaktoro_solver import ReaktoroSolver

    cfg = HessianOptions().get_dict()(
        {"bfgs_initialization_type": "scalar1", "bfgs_init_const_hessian_value": 1e-16}
    )
    resolve_deprecated_hessian_options(cfg)
    assert cfg.quasi_newton_initialization_type == "scalar1"
    assert cfg.quasi_newton_init_const_value == 1e-16
    assert cfg.bfgs_initialization_type is None
    # resolving again is a no-op
    resolve_deprecated_hessian_options(cfg)
    assert cfg.quasi_newton_initialization_type == "scalar1"
    cfg = HessianOptions().get_dict()(
        {
            "bfgs_initialization_type": "scalar1",
            "quasi_newton_initialization_type": "zero",
        }
    )
    with pytest.raises(ValueError):
        resolve_deprecated_hessian_options(cfg)
    # new name explicitly set to its default value still conflicts
    cfg = HessianOptions().get_dict()({"lbfgs_memory": 6, "bfgs_hessian_memory": 3})
    with pytest.raises(ValueError):
        resolve_deprecated_hessian_options(cfg)

    # solver keyword alias, without building a reaktoro solver
    solver = ReaktoroSolver.__new__(ReaktoroSolver)

    class _Opt:
        class optima:
            maxiters = None

            class convergence:
                tolerance = None

    solver.solver_options = _Opt()
    solver.presolve_options = _Opt()

    class _Solver:
        def setOptions(self, opts):
            pass

    solver.solver = _Solver()
    solver.set_solver_options(bfgs_initialization_type="zero", bfgs_hessian_memory=4)
    assert solver.quasi_newton_initialization_type == "zero"
    assert solver.lbfgs_memory == 4


################################################################################
# safeguards (quasi_newton_update_cap, quasi_newton_norm_cap)


def test_update_cap_skips_near_singular_sr1_update():
    """r's is small relative to |r||s| but above sr1_skip_tolerance, plain SR1 adds
    curvature ~1/cos(r, s) times larger than observed along the step"""
    points = [np.array([0.0, 0.0]), np.array([1.0, 0.0])]
    jacs = [np.array([[0.0, 0.0]]), np.array([[1e-4, 1.0]])]
    norms = {}
    for cap in [0.0, 100.0]:
        hc = HessianApproximation(
            hessian_type=HessTypes.SR1,
            quasi_newton_initialization_type=QuasiNewtonInitializationTypes.zero,
            quasi_newton_update_cap=cap,
            quasi_newton_norm_cap=0,
        )
        for x, J in zip(points, jacs):
            hess = hc.get_hessian(x, J, np.array([1.0]))
        norms[cap] = np.linalg.norm(full(hess))
        if cap > 0:
            assert hc.safeguard_stats["updates_capped"] == 1
    assert norms[0.0] > 1e3
    assert norms[100.0] == 0


@pytest.mark.parametrize(
    "hess_type", [HessTypes.SR1, HessTypes.BFGS_ipopt, HessTypes.PSB]
)
def test_safeguards_do_not_change_well_conditioned_updates(hess_type):
    """on a smooth quadratic the safeguards are inactive"""
    hessians = []
    for cap in [0.0, 100.0]:
        hc = HessianApproximation(
            hessian_type=hess_type,
            quasi_newton_initialization_type=QuasiNewtonInitializationTypes.scalar1,
            quasi_newton_update_cap=cap,
            quasi_newton_norm_cap=cap,
        )
        hessians.append(full(run_quadratic(hc, POINTS, A)))
        if cap > 0:
            assert hc.safeguard_stats == {"updates_capped": 0, "norms_capped": 0}
    np.testing.assert_allclose(hessians[0], hessians[1])


@pytest.mark.parametrize("hess_type", [HessTypes.BFGS, HessTypes.LBFGS])
def test_norm_cap_limits_gauss_newton_initial_matrix(hess_type):
    """GaussNewton initial matrix diag(g^2) is huge for outputs with large gradient but
    small curvature (e.g. y = b'x + 0.5 x'Ax with |b| >> |A|), norm cap bounds it
    (checked after first step, later updates correct B along the steps)"""
    b = np.array([100.0, -50.0, 80.0])
    points = POINTS[:2]
    cap = 10.0
    norms = {}
    for c in [0.0, cap]:
        hc = HessianApproximation(
            hessian_type=hess_type,
            quasi_newton_initialization_type=QuasiNewtonInitializationTypes.GaussNewton,
            quasi_newton_norm_cap=c,
        )
        for x in points:
            hess = hc.get_hessian(x, np.array([b + A @ x]), np.array([1.0]))
        norms[c] = np.linalg.norm(full(hess))
    step = points[1] - points[0]
    max_curv = np.linalg.norm(A @ step) / np.linalg.norm(step)
    assert norms[0.0] > 100 * max_curv
    assert norms[cap] <= cap * max_curv * (1 + 1e-12)


@pytest.mark.parametrize(
    "hess_type, expected",
    [
        (HessTypes.SR1, 10.0),
        (HessTypes.SR1_lagrangian, 10.0),
        (HessTypes.BFGS, 0.0),
        (HessTypes.LBFGS, 0.0),
        (HessTypes.PSB, 0.0),
    ],
)
def test_default_safeguards_only_for_sr1(hess_type, expected):
    hc = HessianApproximation(hessian_type=hess_type)
    assert hc.update_cap == expected
    assert hc.norm_cap == expected
    # explicit value always wins
    hc = HessianApproximation(
        hessian_type=hess_type, quasi_newton_update_cap=0, quasi_newton_norm_cap=5
    )
    assert (hc.update_cap, hc.norm_cap) == (0.0, 5.0)
