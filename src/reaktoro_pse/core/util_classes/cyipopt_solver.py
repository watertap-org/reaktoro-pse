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
    return cy_solver
