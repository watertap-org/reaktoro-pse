from pyomo.common.config import ConfigValue, IsInstance, ConfigDict
from reaktoro_pse.core.reaktoro_jacobian import JacType
from reaktoro_pse.core.reaktoro_block_builder import JacScalingTypes
from reaktoro_pse.core.util_classes.hessian_functions import HessTypes


class HessianOptions:
    def __init__(self):
        pass

    def get_dict(self):
        CONFIG = ConfigDict()

        CONFIG.declare(
            "hessian_type",
            ConfigValue(
                default=HessTypes.LBFGS,
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
                    """,
            ),
        )
        CONFIG.declare(
            "bfgs_initialization_type",
            ConfigValue(
                default=HessTypes.GaussNewton,
                domain=IsInstance((str, HessTypes)),
                description="Hessian initialzation type for BFGS",
                doc="""Scalar type for initialization BFGS hessian
                - GaussNewton - Gauss-Newton Hessian approximation (Jacobian^T * Jacobian)
                - scalar1 - sTy/sTs
                - scalar2 - yTy/sTy
                where s is change in step, and y is change in gradient
                    """,
            ),
        )
        CONFIG.declare(
            "bfgs_init_min_hessian_value",
            ConfigValue(
                default=1e-64,
                domain=float,
                description="Minimum Hessian value for BFGS initialization",
                doc="""Minimum Hessian value for BFGS initialization
                    """,
            ),
        )
        CONFIG.declare(
            "bfgs_init_max_hessian_value",
            ConfigValue(
                default=1e8,
                domain=float,
                description="Maximum Hessian value for BFGS initialization",
                doc="""Maximum Hessian value for BFGS initialization
                    """,
            ),
        )
        CONFIG.declare(
            "bfgs_init_const_hessian_value",
            ConfigValue(
                default=1e-32,
                domain=float,
                description="Constant Hessian value for BFGS initialization",
                doc="""Constant Hessian value for BFGS initialization
                    """,
            ),
        )
        CONFIG.declare(
            "bfgs_hessian_memory",
            ConfigValue(
                default=6,
                domain=int,
                description="Memory size for BFGS Hessian approximation",
                doc="""Memory size for BFGS Hessian approximation
                    """,
            ),
        )
        CONFIG.declare(
            "bfgs_epsilon",
            ConfigValue(
                default=1e-12,
                domain=float,
                description="Epsilon value for BFGS Hessian approximation updates",
                doc="""This is used to define when an update for hesisan is accepted, 
                in general should be several orders of magnitude larger then machine precision.
                IPOPT uses 1e-12, so we use the same value here.
                    """,
            ),
        )
        return CONFIG
