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

__author__ = "Alexander V. Dudchenko"


from watertap_solvers import get_solver
import warnings


class LinearSolverTypes:
    ma27 = "ma27"
    mumps = "mumps"


def get_cyipopt_watertap_solver(
    max_iter=500,
    linear_solver=LinearSolverTypes.mumps,
    limited_memory=False,
    solver_options=None,
    scalar_type="scalar1",
    dual_inf_tol=1e-1,
    constr_viol_tol=1e-8,
    tol=1e-8,
    pivtol=None,
    pivtolmax=None,
):
    """
    Helper function to get cyipopt-watertap solver with options commonly used in reaktoro_pse examples
    Args:
        max_iter: maximum number of iterations
        linear_solver: linear solver to use, see LinearSolverTypes class for options or solver supported by cyipopt configuration
        limited_memory: use limited memory BFGS hessian approximation
        solver_options: dictionary of additional solver options to set
        scalar_type: limited memory initialization type, see cyipopt documentation for options
        dual_inf_tol: dual infeasibility tolerance
        constr_viol_tol: constraint violation tolerance
        tol: overall convergence tolerance
        pivtol: pivot tolerance for linear solver (solver must accept option as <linear_solver>_pivtol)
        pivtolmax: maximum pivot tolerance for linear solver (solver must accept option as <linear_solver>_pivtolmax)
    Returns:
        cyipopt-watertap solver with specified options
    """

    cy_solver = get_solver(solver="cyipopt-watertap")
    cy_solver.options["max_iter"] = max_iter
    cy_solver.options["print_user_options"] = "yes"
    # helps handle property packages that have very small values requiring large steps
    cy_solver.options["diverging_iterates_tol"] = 1e30
    cy_solver.options["linear_solver"] = linear_solver
    if pivtol is not None:
        cy_solver.options[f"{linear_solver}_pivtol"] = float(pivtol)
    if pivtolmax is not None:
        cy_solver.options[f"{linear_solver}_pivtolmax"] = float(pivtolmax)
    if limited_memory:
        cy_solver.options["hessian_approximation"] = "limited-memory"
        cy_solver.options["limited_memory_initialization"] = scalar_type
    cy_solver.options["dual_inf_tol"] = dual_inf_tol
    cy_solver.options["constr_viol_tol"] = constr_viol_tol
    cy_solver.options["tol"] = tol
    # Ensure we never accept a "acceptable solution", which will be treated as infeasible"
    cy_solver.options["acceptable_dual_inf_tol"] = dual_inf_tol / 10
    cy_solver.options["acceptable_tol"] = tol / 10
    cy_solver.options["acceptable_constr_viol_tol"] = constr_viol_tol / 10
    cy_solver.options["acceptable_dual_inf_tol"] = dual_inf_tol / 10
    if solver_options is not None:
        for opt, value in solver_options.items():
            cy_solver.options[opt] = value
    return ReaktoroSolverWrapper(cy_solver)


class ReaktoroSolverWrapper:
    """Wraps a pyomo solver, models with reaktoro trust region surrogates
    (ReaktoroBlock(backend="trust_region") or ReaktoroBlockManager(backend="trust_region"))
    are solved with the trust region method using the wrapped solver for subproblems,
    all other models are solved directly with the wrapped solver.
    All attributes (options, available, ...) are those of the wrapped solver."""

    def __init__(self, solver):
        self.__dict__["_solver"] = solver
        self.__dict__["trust_region_stats"] = None

    @property
    def wrapped_solver(self):
        return self._solver

    def solve(self, model, *args, **kwds):
        # imported here, trust region solver imports reaktoro block modules
        from reaktoro_pse.core.util_classes.trust_region_solver import (
            ReaktoroTrustRegionSolver,
            find_surrogate_blocks,
        )

        if not find_surrogate_blocks(model):
            return self._solver.solve(model, *args, **kwds)
        tr_solver = ReaktoroTrustRegionSolver(self._solver)
        result = tr_solver.solve(model, *args, **kwds)
        self.__dict__["trust_region_stats"] = tr_solver.trust_region_stats
        return result

    def __getattr__(self, attr):
        return getattr(self.__dict__["_solver"], attr)

    def __setattr__(self, attr, value):
        setattr(self._solver, attr, value)


def get_reaktoro_solver(**kwargs):
    """same as get_cyipopt_watertap_solver, solves models with reaktoro blocks on
    graybox or trust_region backend"""
    return get_cyipopt_watertap_solver(**kwargs)
