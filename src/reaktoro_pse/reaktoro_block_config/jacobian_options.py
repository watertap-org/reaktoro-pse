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

from pyomo.common.config import ConfigValue, IsInstance, ConfigDict
from reaktoro_pse.core.reaktoro_jacobian import JacType
from reaktoro_pse.core.reaktoro_block_builder import JacScalingTypes
from reaktoro_pse.core.util_classes.rkt_inputs import RktInputTypes


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
                default={
                    RktInputTypes.pH: 1e-4,
                    RktInputTypes.temperature: 1e-4,
                    RktInputTypes.pressure: 1e-4,
                    RktInputTypes.enthalpy: 1e-4,
                    RktInputTypes.species: 1e-4,
                },
                domain=IsInstance((dict, float)),
                description="Defines the step to use for numerical descritiazaiton based on input type, if None, automatically found",
                doc="""This will define how small of a step to use for numerical derivative propagation which takes
                the absolute chemical property and multiplies it by chemical property derivative multiplied by step size. 
                    chemical_property_step=chemical_input*chemical_property_derivative*step
                """,
            ),
        )

        CONFIG.declare(
            "scaling_type",
            ConfigValue(
                default=JacScalingTypes.variable_oi_scaling_square_sum,
                domain=IsInstance((str, JacScalingTypes)),
                description="Defines how to scale Jacobian matrix",
                doc="""
                Defines methods for jacobian scaling:
                - no_scaling -- jacobian scale == 1 for all outputs
                - variable_output_scaling -- use output variable scaling factors (output_scale_i)                
                - variable_oi_scaling_inverse_sum -- sum squared of output/input variable scaling factors output_scale_i/((sum(input_scales_i)**-1)**-1)
                - variable_oi_scaling_square_sum --  (default) use inverse of sum squared of output/input variable scaling factors output_scale_i/((sum(input_scales_i)**2)**0.5)
                - jacobian_matrix_inverse_sum -- use inverse of sum of absolute values of jacobian matrix
                - jacobian_matrix_square_sum -- use squared sum of absolute values of jacobian matrix
                - user_scaling -- Use user provided scaling
                """,
            ),
        )
        CONFIG.declare(
            "jacobian_scale_bounds",
            ConfigValue(
                default=(1e-8, 1e2),
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
                default=False,
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
        return CONFIG
