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

from pyomo.common.config import ConfigValue, ConfigDict, In, IsInstance
from pyomo.contrib.trustregion.TRF import _trf_config

__author__ = "Alexander V. Dudchenko"


class BackendTypes:
    """solve backends for reaktoro blocks"""

    graybox = "graybox"
    trust_region = "trust_region"
    all = (graybox, trust_region)


class SurrogateTypes:
    """surrogate models used by the trust region backend"""

    auto = "auto"  # chosen from hessian_options.hessian_type
    linear = "linear"
    quadratic = "quadratic"
    all = (auto, linear, quadratic)


class TrustRegionFormulations:
    """how the trust region enters the subproblem"""

    soft = "soft"  # surrogate at actual inputs, step beyond radius is penalized
    trusted_point = "trusted_point"  # surrogate at w0 + h*s with |s| <= radius
    all = (soft, trusted_point)


class ElasticPenaltyModes:
    fixed = "fixed"
    adaptive = "adaptive"
    all = (fixed, adaptive)


# defaults that differ from pyomo.contrib.trustregion defaults
# (pyomo radius is unscaled on decision variables, ours is a relative step of the
# reaktoro inputs, see ReaktoroSurrogateBlock)
_TRF_DEFAULT_OVERRIDES = {
    "trust_radius": 0.5,
    "minimum_radius": 1e-8,
    "maximum_radius": 10.0,
    "maximum_iterations": 100,
    "feasibility_termination": 1e-7,
    "step_size_termination": 1e-6,
    "minimum_feasibility": 1e-3,
    "maximum_feasibility": 1e12,
    "verbose": True,
}


def backend_config_value():
    """ConfigValue for backend option (ReaktoroBlock and ReaktoroBlockManager)"""
    return ConfigValue(
        default=BackendTypes.graybox,
        domain=In(BackendTypes.all),
        description="Solve backend for reaktoro blocks",
        doc="""Defines how reaktoro is coupled to the pyomo model
        - graybox - reaktoro is a gray box model (ExternalGreyBoxBlock), ipopt calls reaktoro
          on every iteration and uses a quasi-Newton hessian approximation
        - trust_region - reaktoro is replaced by an algebraic surrogate (linear or quadratic
          taylor model) that is updated by a trust region filter method, ipopt solves
          the flowsheet with the surrogate, reaktoro is only evaluated at subproblem solutions.
          Solve with the solver returned by get_cyipopt_watertap_solver (or get_reaktoro_solver),
          which detects the surrogates and runs the trust region method.
        For blocks registered with ReaktoroBlockManager the manager's backend is used.
        """,
    )


class TrustRegionOptions:
    def __init__(self):
        pass

    def get_dict(self):
        # filter, radius and termination options of pyomo.contrib.trustregion
        CONFIG = _trf_config()
        for key in ("solver", "keepfiles", "tee"):
            # inner solver and tee are taken from solver calling the trust region method
            if key in CONFIG:
                del CONFIG[key]
        for key, value in _TRF_DEFAULT_OVERRIDES.items():
            CONFIG.get(key).set_default_value(value)
            CONFIG[key] = value
        # NOTE: trust_radius is the maximum relative step of reaktoro inputs
        # (amounts relative to their current value, pH, pE type inputs relative
        # to log_input_step)
        CONFIG.declare(
            "surrogate",
            ConfigValue(
                default=SurrogateTypes.auto,
                domain=In(SurrogateTypes.all),
                description="Surrogate model type",
                doc="""Surrogate of reaktoro outputs used in trust region subproblems
                - auto - quadratic for quasi-Newton hessian types (SR1, LBFGS, BFGS, ..),
                  linear for ZeroHessian, GaussNewton and no_hessian_estimation
                - linear - first order taylor model
                - quadratic - second order taylor model, curvature from hessian_options
                """,
            ),
        )
        CONFIG.declare(
            "log_input_step",
            ConfigValue(
                default=1.0,
                domain=float,
                description="Step scale for pH, pOH, pE, Eh inputs",
                doc="""Inputs on log scale (pH, pOH, pE, Eh) can move by trust_radius*log_input_step
                per trust region iteration, all other inputs can move by trust_radius*|value|
                (or trust_radius/scaling factor for values near zero)""",
            ),
        )
        CONFIG.declare(
            "trust_region_formulation",
            ConfigValue(
                default=TrustRegionFormulations.soft,
                domain=In(TrustRegionFormulations.all),
                description="Trust region formulation in the subproblem",
                doc="""In both cases (inputs - w0)/h == s + ep - en, |s| <= radius, ep, en >= 0
                with penalty on ep + en
                - soft - surrogate is evaluated at the actual inputs, so the trust region is
                  an exact l1 penalty on the step beyond the radius. Surrogate outputs can
                  always follow the flowsheet (e.g. large forced input changes)
                - trusted_point - surrogate is evaluated at w0 + h*s, so it is never
                  extrapolated beyond the radius, but outputs it can reach shrink with
                  the radius, which can make the subproblem infeasible
                """,
            ),
        )
        CONFIG.declare(
            "elastic_penalty_mode",
            ConfigValue(
                default=ElasticPenaltyModes.adaptive,
                domain=In(ElasticPenaltyModes.all),
                description="Elastic penalty mode",
                doc="""Inputs are allowed to move outside of the trust region at a penalty
                (elastic trust region), so the subproblem is never infeasible because of the trust region
                - fixed - penalty stays at elastic_penalty
                - adaptive - penalty is increased by elastic_penalty_increase when
                  elastic gap after an accepted step did not decrease compared to the
                  previous accepted step, up to elastic_penalty_max
                """,
            ),
        )
        CONFIG.declare(
            "elastic_penalty",
            ConfigValue(
                default=1e2,
                domain=float,
                description="Elastic penalty (initial value for adaptive mode)",
                doc="Penalty on scaled (relative) step outside of trust region, in scaled objective units",
            ),
        )
        CONFIG.declare(
            "elastic_penalty_increase",
            ConfigValue(
                default=10.0,
                domain=float,
                description="Elastic penalty increase factor (adaptive mode)",
            ),
        )
        CONFIG.declare(
            "elastic_penalty_max",
            ConfigValue(
                default=1e6,
                domain=float,
                description="Maximum elastic penalty (adaptive mode)",
            ),
        )
        CONFIG.declare(
            "keep_surrogate_curvature",
            ConfigValue(
                default=False,
                domain=bool,
                description="Keep quadratic surrogate curvature between solves",
                doc="""If False, curvature of quadratic surrogate (and quasi-Newton history) is reset
                at start of every solve. If True it is kept as a warm start for next solve""",
            ),
        )
        CONFIG.declare(
            "curvature_reset_step",
            ConfigValue(
                default=1.0,
                domain=float,
                description="Reset quadratic surrogate curvature after large steps",
                doc="""If an accepted step is larger than this (max relative input step, same
                units as trust_radius), or is a restoration step, curvature of the quadratic
                surrogate is reset (surrogate becomes linear) instead of a quasi-Newton update""",
            ),
        )
        CONFIG.declare(
            "warm_start_options",
            ConfigValue(
                default={"bound_push": 1e-8, "bound_frac": 1e-8},
                domain=dict,
                description="Ipopt options for inner solves after the first one",
                doc="""Applied to every inner solve after first one in a trust region solve,
                only if user did not set that option on the solver""",
            ),
        )
        CONFIG.declare(
            "maximum_failed_inner_solves",
            ConfigValue(
                default=10,
                domain=int,
                description="Maximum number of consecutive failed inner solves",
            ),
        )
        CONFIG.declare(
            "model_slack",
            ConfigValue(
                default=True,
                domain=bool,
                description="Elastic slack on surrogate equations after failed inner solve",
                doc="""If an inner solve fails (subproblem infeasible for any reaktoro inputs,
                e.g. a fixed spec far from the expansion point while an input is at its
                bound), it is retried with l1 penalized slack on surrogate equations (in
                units relative to output magnitude, same penalty as elastic trust region gap). Steps with
                slack are restoration steps, slack is fixed at zero again once it is not
                needed, and trust region can only terminate optimal with zero slack""",
            ),
        )
        return CONFIG
