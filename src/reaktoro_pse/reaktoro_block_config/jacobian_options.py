from pyomo.common.config import ConfigValue, IsInstance, ConfigDict
from reaktoro_pse.core.reaktoro_jacobian import JacType
from reaktoro_pse.core.reaktoro_block_builder import JacScalingTypes
from reaktoro_pse.core.reaktoro_gray_box import HessTypes


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
                Derivatives for many of the properties in reaktro are not directly available, 
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
                default=8,
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
                default=1e-5,
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
                default=JacScalingTypes.variable_scaling,
                domain=IsInstance((str, JacScalingTypes)),
                description="Defines how to scale Jacobian matrix",
                doc="""
                Defines methods for jacobian scaling:
                - if option is no_scaling, jacobian scale will == 1 for all outputs
                - if option is 'variable_scaling' will use output variable scaling factors
                - if option is jacobian_matrix will use actual jac matrix to calculate scaling factors
                - if user_scaling is not None then uses user provided scaling
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
                    e.g. {output_name:scaling_factor} applies to specific jac output 
                """,
            ),
        )
        CONFIG.declare(
            "hessian_type",
            ConfigValue(
                default="BFGS",
                domain=IsInstance((str, HessTypes)),
                description="Hessian type to use for reaktor gray box",
                doc="""Hessian type to use, some might provide better stability
                options (Jt.J, BFGS, BFGS-mod,BFGS-damp, BFGS-ipopt""",
            ),
        )
        return CONFIG
