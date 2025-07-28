from watertap_solvers import get_solver


def get_cyipopt_watertap_solver(
    max_iter=500,
    ma27=False,
    limited_memory=False,
    solver_options=None,
    scalar_type="scalar1",
    dual_inf_tol=1e-1,
    constr_viol_tol=1e-8,
    tol=1e-8,
):
    """general config for cyipopt solver"""
    cy_solver = get_solver(solver="cyipopt-watertap")
    cy_solver.options["max_iter"] = max_iter
    # only enable if avaialbe !
    cy_solver.options["print_user_options"] = "yes"
    # helps handle property packages that have very small values requiring large steps

    cy_solver.options["diverging_iterates_tol"] = 1e30
    if ma27:
        cy_solver.options["linear_solver"] = "ma27"
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
