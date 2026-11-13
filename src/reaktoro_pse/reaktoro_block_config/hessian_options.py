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

from pyomo.common.config import ConfigValue, IsInstance, ConfigDict
from reaktoro_pse.core.util_classes.hessian_functions import (
    HessTypes,
    QuasiNewtonInitializationTypes,
    HESSIAN_DEFAULTS,
    DEPRECATED_HESSIAN_OPTIONS,
    _warn_deprecated_option,
)


def resolve_deprecated_hessian_options(config):
    """copy values given with deprecated (bfgs_*) names in a hessian_options config
    block to their current names, the deprecated entries are reset to None"""
    for old, new in DEPRECATED_HESSIAN_OPTIONS.items():
        value = config[old]
        if value is None:
            continue
        _warn_deprecated_option(old, new)
        if config[new] != HESSIAN_DEFAULTS[new] and config[new] != value:
            raise ValueError(
                f"Both '{old}' (deprecated) and '{new}' were given with different values"
            )
        config[new] = value
        config[old] = None
    return config


def _optional_float(val):
    """config domain, float or None"""
    return None if val is None else float(val)


class HessianOptions:
    def __init__(self):
        pass

    def get_dict(self):
        CONFIG = ConfigDict()

        CONFIG.declare(
            "hessian_type",
            ConfigValue(
                default=HESSIAN_DEFAULTS["hessian_type"],
                domain=IsInstance((str, HessTypes)),
                description="Hessian type to use for reaktor gray box",
                doc="""Hessian type to use, some might provide better stability
                options:
                - ZeroHessian - no hessian
                - GaussNewton - Naive Gauss-Newton Hessian approximation (Jacobian^T * Jacobian)
                - LBFGS - Limited Memory Broyden-Fletcher-Goldfarb-Shanno
                - BFGS - Broyden-Fletcher-Goldfarb-Shanno
                - CBFGS - conditional BFGS
                - BFGS_mod - modified BFGS
                - BFGS_damp - damped BFGS
                - BFGS_ipopt - BFGS with ipopt update step
                - SR1 - symmetric rank one, per output, allows indefinite curvature
                - PSB - Powell symmetric Broyden, per output, allows indefinite curvature
                - SR1_lagrangian - symmetric rank one on the Lagrangian of each block
                    """,
            ),
        )
        CONFIG.declare(
            "quasi_newton_initialization_type",
            ConfigValue(
                default=HESSIAN_DEFAULTS["quasi_newton_initialization_type"],
                domain=IsInstance((str, QuasiNewtonInitializationTypes)),
                description="Initial matrix for quasi-Newton hessians",
                doc="""Initial matrix for quasi-Newton hessian types (BFGS family, LBFGS, SR1, PSB)
                - GaussNewton - Gauss-Newton Hessian approximation (Jacobian^T * Jacobian)
                - scalar1 - sTy/sTs
                - scalar2 - yTy/sTy
                - scalar3 - average of scalar1 and scalar2
                - scalar4 - geometric mean of scalar1 and scalar2
                - constant - quasi_newton_init_const_value
                - zero - start from zero matrix (SR1, PSB), adds no curvature in
                  directions that were not explored
                where s is change in step, and y is change in gradient
                    """,
            ),
        )
        CONFIG.declare(
            "quasi_newton_init_min_value",
            ConfigValue(
                default=HESSIAN_DEFAULTS["quasi_newton_init_min_value"],
                domain=float,
                description="Minimum value of quasi-Newton initial matrix",
                doc="""Minimum absolute value of quasi-Newton initial matrix diagonal (scalar initializations)
                    """,
            ),
        )
        CONFIG.declare(
            "quasi_newton_init_max_value",
            ConfigValue(
                default=HESSIAN_DEFAULTS["quasi_newton_init_max_value"],
                domain=float,
                description="Maximum value of quasi-Newton initial matrix",
                doc="""Maximum absolute value of quasi-Newton initial matrix diagonal (scalar initializations)
                    """,
            ),
        )
        CONFIG.declare(
            "quasi_newton_init_const_value",
            ConfigValue(
                default=HESSIAN_DEFAULTS["quasi_newton_init_const_value"],
                domain=float,
                description="Diagonal value of constant quasi-Newton initial matrix",
                doc="""Diagonal value of initial matrix for the constant initialization, also used
                before the first valid step
                    """,
            ),
        )
        CONFIG.declare(
            "lbfgs_memory",
            ConfigValue(
                default=HESSIAN_DEFAULTS["lbfgs_memory"],
                domain=int,
                description="Number of stored steps for LBFGS",
                doc="""Number of stored steps used to rebuild the LBFGS hessian (only used by LBFGS)
                    """,
            ),
        )
        CONFIG.declare(
            "quasi_newton_epsilon",
            ConfigValue(
                default=HESSIAN_DEFAULTS["quasi_newton_epsilon"],
                domain=float,
                description="Epsilon for BFGS family and LBFGS update tests",
                doc="""This is used to define when an update for BFGS hessian is accepted,
                in general should be an order of magnitude above then machine precision.
                    """,
            ),
        )
        CONFIG.declare(
            "sr1_skip_tolerance",
            ConfigValue(
                default=HESSIAN_DEFAULTS["sr1_skip_tolerance"],
                domain=float,
                description="Skip tolerance for SR1 updates",
                doc="""SR1 update is skipped if |r's| < tolerance*|s|*|r|, where
                r = y - B s (Nocedal & Wright eq. 6.26)
                    """,
            ),
        )
        CONFIG.declare(
            "quasi_newton_update_cap",
            ConfigValue(
                default=HESSIAN_DEFAULTS["quasi_newton_update_cap"],
                domain=_optional_float,
                description="Skip quasi-Newton updates that are large relative to observed curvature",
                doc="""Safeguard for cumulative quasi-Newton types (BFGS family, SR1, PSB,
                SR1_lagrangian). An update B -> B + dB is skipped if
                |dB| > cap * max(|B|, |y|/|s|) (Frobenius norms), where |y|/|s| is the
                curvature observed along the step. Catches SR1 updates with a small
                denominator (adds curvature far from the step direction, often with the
                wrong sign) and BFGS updates with tiny y's, both usually driven by noise in
                Reaktoro derivatives. 0 disables. None (default) uses 10 for SR1 and
                SR1_lagrangian and disables it for other types (on softening RO benchmarks it
                removed bad SR1 updates but made BFGS and LBFGS slower).
                    """,
            ),
        )
        CONFIG.declare(
            "quasi_newton_norm_cap",
            ConfigValue(
                default=HESSIAN_DEFAULTS["quasi_newton_norm_cap"],
                domain=_optional_float,
                description="Cap quasi-Newton matrix norm relative to observed curvature",
                doc="""Safeguard for all quasi-Newton types (including LBFGS). If the matrix of an
                element grows beyond cap * max |y|/|s| seen over the last lbfgs_memory steps
                (Frobenius norm), it is scaled down to that norm. Catches slow accumulation of
                updates that individually pass quasi_newton_update_cap, and initial matrices
                (e.g. GaussNewton) far above the observed curvature. 0 disables. None (default)
                uses 10 for SR1 and SR1_lagrangian and disables it for other types.
                Both caps compare norms in the space updates are done in, so they depend on
                hessian_input_scaling.
                    """,
            ),
        )
        CONFIG.declare(
            "hessian_input_scaling",
            ConfigValue(
                default=HESSIAN_DEFAULTS["hessian_input_scaling"],
                domain=bool,
                description="Do quasi-Newton updates in scaled inputs",
                doc="""If True, quasi-Newton updates are done in scaled inputs
                (x_scaled = sf*x, using scaling factors on gray box input variables),
                this makes the initial matrix and update tests consistent across inputs
                with very different magnitudes. Hessian returned to solver is unscaled.
                    """,
            ),
        )
        for old, new in DEPRECATED_HESSIAN_OPTIONS.items():
            CONFIG.declare(
                old,
                ConfigValue(
                    default=None,
                    description=f"Deprecated, use {new}",
                    doc=f"""Deprecated name of {new}, kept for backwards compatibility""",
                ),
            )
        return CONFIG
