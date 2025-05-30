from watertap_solvers import get_solver


def get_cyipopt_watertap_solver(
    max_iter=500, ma27=False, limited_memory=False, solver_args=None
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
        cy_solver.options["hessian_approximation_space"] = "all-variables"
        cy_solver.options["limited_memory_max_history"] = 5
        cy_solver.options["limited_memory_aug_solver"] = "extended"
        cy_solver.options["limited_memory_initialization"] = "scalar1"
    if solver_args is not None:
        for arg, value in solver_args.items():
            cy_solver.options[arg] = value
    return cy_solver
