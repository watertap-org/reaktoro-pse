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
    BFGSInitializationTypes,
    HESSIAN_DEFAULTS,
)


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
            "bfgs_initialization_type",
            ConfigValue(
                default=HESSIAN_DEFAULTS["bfgs_initialization_type"],
                domain=IsInstance((str, BFGSInitializationTypes)),
                description="Hessian initialzation type for BFGS",
                doc="""Scalar type for initialization BFGS hessian
                - GaussNewton - Gauss-Newton Hessian approximation (Jacobian^T * Jacobian)
                - scalar1 - sTy/sTs
                - scalar2 - yTy/sTy
                - scalar3 - average of scalar1 and scalar2
                - scalar4 - geometric mean of scalar1 and scalar2
                - constant - bfgs_init_const_hessian_value
                - zero - start from zero matrix (SR1, PSB), adds no curvature in
                  directions that were not explored
                where s is change in step, and y is change in gradient
                    """,
            ),
        )
        CONFIG.declare(
            "bfgs_init_min_hessian_value",
            ConfigValue(
                default=HESSIAN_DEFAULTS["bfgs_init_min_hessian_value"],
                domain=float,
                description="Minimum Hessian value for BFGS initialization",
                doc="""Minimum Hessian value for BFGS initialization
                    """,
            ),
        )
        CONFIG.declare(
            "bfgs_init_max_hessian_value",
            ConfigValue(
                default=HESSIAN_DEFAULTS["bfgs_init_max_hessian_value"],
                domain=float,
                description="Maximum Hessian value for BFGS initialization",
                doc="""Maximum Hessian value for BFGS initialization
                    """,
            ),
        )
        CONFIG.declare(
            "bfgs_init_const_hessian_value",
            ConfigValue(
                default=HESSIAN_DEFAULTS["bfgs_init_const_hessian_value"],
                domain=float,
                description="Constant Hessian value for BFGS initialization",
                doc="""Constant Hessian value for BFGS initialization
                    """,
            ),
        )
        CONFIG.declare(
            "bfgs_hessian_memory",
            ConfigValue(
                default=HESSIAN_DEFAULTS["bfgs_hessian_memory"],
                domain=int,
                description="Memory size for BFGS Hessian approximation",
                doc="""Memory size for BFGS Hessian approximation
                    """,
            ),
        )
        CONFIG.declare(
            "bfgs_epsilon",
            ConfigValue(
                default=HESSIAN_DEFAULTS["bfgs_epsilon"],
                domain=float,
                description="Epsilon value for BFGS Hessian approximation updates",
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
        return CONFIG
