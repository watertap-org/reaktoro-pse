from pyomo.common.config import ConfigValue, IsInstance, ConfigDict
from reaktoro_pse.core.reaktoro_jacobian import JacType
from reaktoro_pse.core.reaktoro_block_builder import JacScalingTypes
from reaktoro_pse.core.util_classes.hessian_functions import HessTypes


class JacobianOptions:
    def __init__(self):
        pass

    def get_dict(self):
        CONFIG = ConfigDict()
        CONFIG.declare(
            "numerical_type",
            ConfigValue(
                default=JacType.average,
                domain=IsInstance((str, JacType)),
                description="Defines method for numerical jacobian approximations",
                doc="""
                Derivatives for many of the properties in Reaktro are not directly available, 
                thus we numerically propagate derivatives from chemical state to methods for estimation of these properties. 
                Two methods are available, average and center_difference
                    - average methods takes defined number of derivatives by numerical_jacobian_order from center points and gets the average of them
                    - center_difference methods applies classical taylor difference approximation methods 
                In theory the two should yield same result- but due to round off errors the average method might provide better error dampening. 

                """,
            ),
        )
        CONFIG.declare(
            "numerical_order",
            ConfigValue(
                default=10,
                domain=int,
                description="Defines order of numerical jacobian (should be an even number)",
                doc="""
                This will define how many points to discretize the derivate over 
                - for numerical_jacobian_type==average - order can be any even number
                - for numerical_jacobian_type==center_difference - order can be 2, 4, 6, 8, 10
                """,
            ),
        )
        CONFIG.declare(
            "numerical_step",
            ConfigValue(
                default=1e-4,
                domain=float,
                description="Defines the step to use for numerical descritiazaiton",
                doc="""This will define how small of a step to use for numerical derivative propagation which takes
                the absolute chemical property and multiplies it by chemical property derivative multiplied by step 
                    chemical_property_step=chemical_property_absolute_value*chemical_property_derivative*step
                """,
            ),
        )
        CONFIG.declare(
            "scaling_type",
            ConfigValue(
                default=JacScalingTypes.jacobian_matrix_inverse_sum,
                domain=IsInstance((str, JacScalingTypes)),
                description="Defines how to scale Jacobian matrix",
                doc="""
                Defines methods for jacobian scaling:
                - if option is no_scaling, jacobian scale will == 1 for all outputs
                - if option is 'variable_scaling' will use output variable scaling factors
                - if option is 'inverse_variable_scaling' will use inverse of output variable scaling factors
                - if option is 'variable_io_scaling' will sums squared of input scales and output scales
                - if option is jacobian_matrix_inverse_sum will use inverse of sum of absolute values of jacobian matrix
                - if option is jacobian_matrix_square_sum will use squared sum of absolute values of jacobian matrix
                - if user_scaling is not None then uses user provided scaling
                """,
            ),
        )
        CONFIG.declare(
            "jacobian_scale_bounds",
            ConfigValue(
                default=(1e-10, 1e2),
                domain=IsInstance(tuple),
                description="Defines lower and upper bounds for jacobian scaling factors",
                doc="""
                This will clip jacobian scale by defined upper and lower bound (min, max).   
                Passing in None instead of a value will disable clipping for min or max (e.g. (None, 1e2) will disable lower bound clipping).             
                """,
            ),
        )
        CONFIG.declare(
            "jacobian_scaling_bounds_output_based",
            ConfigValue(
                default=True,
                domain=bool,
                description="Defines if lower and upper bounds for jacobian scaling factors should be baseded on output scale",
                doc="""
                If True, the jacobian is clipped based on jacbian_scale_bounds multiplied by output variable scaling factors. 
                If False, the jacobian is clipped based on jacbian_scale_bounds only.            
                """,
            ),
        )
        CONFIG.declare(
            "update_jacobian_scale_every_solve",
            ConfigValue(
                default=False,
                domain=bool,
                description="Defines if jacobian scale should be updated every solve",
                doc="""
                This will recalculate jacobian scale every time a new solve is started. 
                This only works if user updates output/input variable scaling between solves 
                or if user uses any of the jacobian_matrix scaling methods, otherwise the jacobian scale factors will not change            
                """,
            ),
        )
        CONFIG.declare(
            "user_scaling",
            ConfigValue(
                default=None,
                domain=IsInstance((float, list, dict)),
                description="Manual scaling factors for jacobian",
                doc="""
                Applies user provided jacobian scaling values:
                - either single value that will be applied to all outputs in jacobian
                - array applied across jacobian
                - dict that specifics output and scaling factor to which apply scaling, (variable_scaling will be applied to non specified outputs)
                    e.g. {output_name: scaling_factor} applies to specific jac output 
                """,
            ),
        )
        CONFIG.declare(
            "hessian_type",
            ConfigValue(
                default=HessTypes.LBFGS,
                domain=IsInstance((str, HessTypes)),
                description="Hessian type to use for reaktor gray box",
                doc="""Hessian type to use, some might provide better stability
                options:                
                - ZeroHessian - no hessian
                - GaussNewton - default
                - BFGS - Broyden-Fletcher-Goldfarb-Shanno   
                - CBFGS - conditional BFGS
                - BFGS_mod - modified BFGS
                - BFGS_damp - damped BFGS   
                - BFGS_ipopt - BFGS with ipopt update step
                    """,
            ),
        )
        return CONFIG
