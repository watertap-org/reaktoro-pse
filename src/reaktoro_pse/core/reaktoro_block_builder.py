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
from pyomo.contrib.pynumero.interfaces.external_grey_box import (
    ExternalGreyBoxBlock,
)
from pyomo.environ import Var

import numpy as np

from reaktoro_pse.core.reaktoro_outputs import PropTypes
from reaktoro_pse.core.reaktoro_solver import (
    ReaktoroSolver,
)
from reaktoro_pse.core.reaktoro_coupled_solver import (
    ReaktoroCoupledSolver,
)
from reaktoro_pse.core.reaktoro_gray_box import (
    ReaktoroGrayBox,
)

from pyomo.util.calc_var_value import calculate_variable_from_constraint

import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
import math

__author__ = "Alexander V. Dudchenko"


# class to build reaktoro block and constraints on provided block
_log = idaeslog.getLogger(__name__)


class JacScalingTypes:
    no_scaling = "no_scaling"
    variable_output_scaling = "variable_output_scaling"
    jacobian_matrix_square_sum = "jacobian_matrix_square_sum"
    jacobian_matrix_inverse_sum = "jacobian_matrix_inverse_sum"
    manual_scaling = "manual_scaling"
    variable_oi_scaling_square_sum = "variable_oi_scaling_square_sum"
    variable_oi_scaling_inverse_sum = "variable_oi_scaling_inverse_sum"


class ReaktoroBlockBuilder:
    def __init__(self, block, reaktoro_solver, build_on_init=True):
        """build reaktoro block builder
        Args:
            block: pyomo block to build reaktoro model on
            reaktoro_solver: ReaktoroSolver or ReaktoroCoupledSolver object
            build_on_init: if True will build reaktoro block right away, otherwise will require explicit call to build_reaktoro_block
        """
        self.block = block
        # TODO: add check to make sure block is a pyomo block, or model
        self.solver = reaktoro_solver
        if (
            isinstance(self.solver, ReaktoroSolver) == False
            and isinstance(self.solver, ReaktoroCoupledSolver) == False
        ):
            raise TypeError("Reaktoro block builder requires a ReaktoroSolver class")
        self.configure_jacobian_scaling()
        self.reaktoro_initialize_function = None  # used to provide external solve call
        self.get_jacobian_matrix_function = None
        self.display_reaktoro_state_function = (
            None  # used to specifying external function to display rkt state
        )
        self.build_output_vars()
        if build_on_init:  # option to support legacy implementation
            self.build_reaktoro_block()

    def build_reaktoro_block(
        self,
        gray_box_model=None,
        reaktoro_initialize_function=None,
        display_reaktoro_state_function=None,
        get_jacobian_matrix_function=None,
    ):
        """build reaktoro model"""
        if gray_box_model is None:
            external_model = ReaktoroGrayBox()
            external_model.configure(self.solver)
            self.block.reaktoro_model = ExternalGreyBoxBlock(
                external_model=external_model
            )
        else:
            self.block.reaktoro_model = gray_box_model
        if reaktoro_initialize_function is not None:
            self.reaktoro_initialize_function = reaktoro_initialize_function
        if display_reaktoro_state_function is not None:
            self.display_reaktoro_state_function = display_reaktoro_state_function
        if get_jacobian_matrix_function is not None:
            self.get_jacobian_matrix_function = get_jacobian_matrix_function
        self.build_input_constraints()
        self.build_output_constraints()
        self.solver.get_jacobian_scaling = self.get_jacobian_scaling
        self.solver.get_input_scaling = self.get_input_scaling

    def configure_jacobian_scaling(
        self,
        jacobian_scaling_type=None,
        user_scaling=None,
        jacobian_scaling_bounds=(1e-16, 1e2),
        jacobian_scaling_bounds_output_based=True,
        update_jacobian_scale_every_solve=False,
    ):
        """define scaling for jacobian, defaults to useing variable scaling

        Keyword:
        jacobianScalingType -- defines type of scaling to use (default: variable_scaling)
            - if option is 'variable_scaling' will use output scaling factors
            - if option is jacobian_matrix will use actual jac matrix
            - if user_scaling is not None then uses user provided scaling
        user_scaling -- either a single value or array with length of rkt outputs defining scaling
        """
        if jacobian_scaling_type is None:
            self.jacobian_scaling_type = JacScalingTypes.no_scaling
        else:
            self.jacobian_scaling_type = jacobian_scaling_type
        if isinstance(user_scaling, float):
            self.solver.jacobian_scaling_values = (
                np.ones(len(self.solver.output_specs.rkt_outputs)) + user_scaling
            )
            self.jacobian_scaling_type = JacScalingTypes.manual_scaling
        elif isinstance(user_scaling, list):
            self.solver.jacobian_scaling_values = user_scaling
            self.jacobian_scaling_type = JacScalingTypes.manual_scaling
        else:
            self.solver.jacobian_scaling_values = np.ones(
                len(self.solver.output_specs.rkt_outputs.keys())
            )
        if isinstance(user_scaling, dict):
            self.user_scaling = user_scaling
        else:
            self.user_scaling = {}
        self.jacobian_scaling_bounds = jacobian_scaling_bounds
        self.update_jacobian_scale_every_solve = update_jacobian_scale_every_solve
        self.jacobian_scaling_bounds_output_based = jacobian_scaling_bounds_output_based

    def build_input_constraints(self):
        """build input constraints for reaktoro model"""
        if self.solver.input_specs.dissolve_species_in_rkt:

            @self.block.Constraint(self.solver.input_specs.rkt_inputs.rkt_input_list)
            def input_constraints(fs, *kwargs):
                if len(kwargs) == 1:
                    kwargs = kwargs[0]

                return (
                    self.block.reaktoro_model.inputs[kwargs]
                    == self.solver.input_specs.rkt_inputs[
                        kwargs
                    ].get_pyomo_with_required_units()
                )

        else:
            # only build these if we are summing species to elements in pyomo
            constraint_dict = self.solver.input_specs.constraint_dict
            self._input_constraint_scaling = {}
            #  connect rektor model vars to our inputs
            for element in constraint_dict:
                self.solver.input_specs.rkt_inputs[element].set_pyomo_var(
                    self.block.reaktoro_model.inputs[element]
                )

            @self.block.Expression(constraint_dict)
            def inputs(fs, element):
                sum_species = []
                for mol, specie in constraint_dict[element]:
                    sum_species.append(mol * self.get_specie_object(specie))
                return sum(sum_species)

            for element in constraint_dict:
                self.solver.input_specs.rkt_inputs[element].set_pyomo_var(
                    self.block.reaktoro_model.inputs[element]
                )

            @self.block.Constraint(self.solver.input_specs.rkt_inputs.rkt_input_list)
            def input_constraints(fs, key):
                if key in constraint_dict:
                    return (
                        self.block.reaktoro_model.inputs[key] == self.block.inputs[key]
                    )
                else:
                    return (
                        self.block.reaktoro_model.inputs[key]
                        == self.solver.input_specs.user_inputs[
                            key
                        ].get_pyomo_with_required_units()
                    )

    def get_specie_object(self, specie):
        """get specie object from input dicts"""
        if specie in self.solver.input_specs.user_inputs:
            pyo_obj = self.solver.input_specs.user_inputs[
                specie
            ].get_pyomo_with_required_units()
        elif specie in self.solver.input_specs.rkt_chemical_inputs:
            pyo_obj = self.solver.input_specs.rkt_chemical_inputs[
                specie
            ].get_pyomo_with_required_units()
        else:
            raise KeyError(f"specie {specie} not found in input dicts")
        return pyo_obj

    def build_output_vars(self):
        """build output variables for reaktoro model"""
        new_output_vars = {}

        for key, obj in self.solver.output_specs.user_outputs.items():
            # NOTE: We do not set rkt_outputs to reaktoro_model outputs as they
            # same as user inputs - we want RKt model to update "user provided vars"
            # rather then pyomo vars in reaktoro model (e.g. reaktor_block.outputs)
            if obj.get_pyomo_var() is None:
                new_output_vars[key] = obj
        if new_output_vars != {}:
            self.block.outputs = Var(new_output_vars.keys(), initialize=1)
            for key, obj in new_output_vars.items():
                obj.set_pyomo_var(self.block.outputs[key])
        self.new_output_vars = new_output_vars

    def build_output_constraints(self):
        """first update rktOuptutObjects for pyomoBuildProperties with reaktoro pyomo variables as
        they will be used in construction of constraints
        The is will also check if user provided an output pyomo var and if not will
        add them to new_output_var dict, which will be used to create new output variables on the block
        """
        for key, obj in self.solver.output_specs.user_outputs.items():
            if PropTypes.pyomo_built_prop == obj.property_type:
                for (
                    pyoPropKey,
                    pyoPropObj,
                ) in obj.pyomo_build_options.properties.items():
                    if pyoPropObj.get_pyomo_var() is None:
                        pyoPropObj.set_pyomo_var(
                            self.block.reaktoro_model.outputs[pyoPropKey]
                        )

        @self.block.Constraint(self.solver.output_specs.user_outputs)
        def output_constraints(fs, prop, prop_index):
            prop_object = self.solver.output_specs.user_outputs[(prop, prop_index)]
            if prop_object.property_type == PropTypes.pyomo_built_prop:
                return prop_object.pyomo_build_options.build_constraint_function(
                    prop_object
                )
            else:
                return (
                    prop_object.get_pyomo_var()
                    == self.block.reaktoro_model.outputs[(prop, prop_index)]
                )

    def initialize(self, presolve_during_initialization=False):
        """initialize reaktoro block
        Args:
            presolve_during_initialization: if True will call presolve before initialization
        """
        self.initialize_input_variables_and_constraints()
        if self.reaktoro_initialize_function is None:
            self.solver.equilibrate_state()
            self.solver.solve_reaktoro_block(presolve=presolve_during_initialization)
        else:
            self.reaktoro_initialize_function(presolve=presolve_during_initialization)

        self.initialize_output_variables_and_constraints()

        self.set_jacobian_scaling()
        self.set_user_jacobian_scaling()
        _log.info(f"Initialized rkt block")

    def get_rkt_scale(self, rkt_input_output, use_default_scaling=True):
        """util function to get scaling for reaktoro vars"""
        sf = self.get_sf(rkt_input_output.get_pyomo_var(), use_default_scaling)
        return sf

    def get_sf(self, pyo_var, use_default_scaling, return_none=1):
        """get scaling factor for pyomo variable"""
        dsf = iscale.get_scaling_factor(pyo_var)
        # only return default scaling factor if we request to use default scaling and its not None,
        # otherwise use the pyomo variable value to calculate scaling factor
        if dsf is not None and use_default_scaling:
            return dsf
        else:
            if pyo_var.value == 0:
                if return_none is None:
                    return None
                _log.warning(f"Var {pyo_var} value is 0")
                return 1

            sf = 1 / (abs(pyo_var.value))
            # Magic Numbers! -  generally for species amounts.
            max_scale = 1e32
            min_scale = 1e-32
            if sf > max_scale:
                _log.warning(
                    f"Var {pyo_var} scale {sf:e}>{max_scale:e}, applied max scale of {max_scale:e}"
                )
                sf = max_scale
            if sf < min_scale:
                _log.warning(
                    f"Var {pyo_var} scale {sf:e}<{min_scale:e}, applied min scale of {min_scale:e}"
                )
                sf = min_scale
            return sf

    def set_output_vars_and_scale(self, use_default_scaling=True):
        for key, obj in self.solver.output_specs.user_outputs.items():
            """update vars scaling in pyomo build constraints
            these are updated to actual value when we call solve_reaktoro_block"""
            if PropTypes.pyomo_built_prop == obj.property_type:
                for (
                    pyoPropKey,
                    pyoPropObj,
                ) in obj.pyomo_build_options.properties.items():
                    val = pyoPropObj.value
                    pyoPropObj.set_pyomo_var_value(val)
                    iscale.set_scaling_factor(
                        pyoPropObj.get_pyomo_var(),
                        self.get_rkt_scale(pyoPropObj, use_default_scaling),
                    )
                output_constraint = self.block.output_constraints[key]
                calculate_variable_from_constraint(
                    obj.get_pyomo_var(), output_constraint
                )
                sf = self.get_rkt_scale(obj, use_default_scaling)
                iscale.constraint_scaling_transform(
                    output_constraint,
                    sf,
                )
            else:
                obj.set_pyomo_var_value(obj.value)
                rkt_var = self.block.reaktoro_model.outputs[key]
                output_constraint = self.block.output_constraints[key]
                calculate_variable_from_constraint(rkt_var, output_constraint)
                sf = self.get_rkt_scale(obj, use_default_scaling)
                iscale.constraint_scaling_transform(
                    output_constraint,
                    sf,
                )
                iscale.set_scaling_factor(rkt_var, sf)
            iscale.set_scaling_factor(
                obj.get_pyomo_var(),
                sf,
            )

    def initialize_output_variables_and_constraints(self):
        self.set_output_vars_and_scale(True)

    def set_jacobian_scaling(self, use_default_scaling=True):
        """Function to calculate jacobian scaling values
        Args:
            use_default_scaling: if True will use default scaling factors for jacobian scaling when
            using variable scaling methods
        """
        output_scales = [
            1 / iscale.get_scaling_factor(obj.get_pyomo_var(), default=1)
            for _, obj in self.solver.output_specs.rkt_outputs.items()
        ]
        output_keys = list(self.solver.output_specs.rkt_outputs.keys())
        if self.jacobian_scaling_type == JacScalingTypes.no_scaling:
            for i, (key, obj) in enumerate(
                self.solver.output_specs.rkt_outputs.items()
            ):
                self.solver.jacobian_scaling_values[i] = 1
        elif self.jacobian_scaling_type == JacScalingTypes.variable_output_scaling:
            for i, (key, obj) in enumerate(
                self.solver.output_specs.rkt_outputs.items()
            ):
                out_sf = self.get_rkt_scale(
                    obj, use_default_scaling=use_default_scaling
                )
                sf = out_sf
                self.solver.jacobian_scaling_values[i] = sf
        elif (
            self.jacobian_scaling_type
            == JacScalingTypes.variable_oi_scaling_inverse_sum
        ):
            for i, (key, obj) in enumerate(
                self.solver.output_specs.rkt_outputs.items()
            ):
                input_scales = []
                out_sf = iscale.get_scaling_factor(obj.get_pyomo_var(), default=1)
                for input_key, input_obj in self.solver.input_specs.rkt_inputs.items():
                    sf = self.get_rkt_scale(
                        input_obj, use_default_scaling=use_default_scaling
                    )
                    input_scales.append(sf)

                sf = np.sum(np.array(input_scales) ** -1) ** -1
                self.solver.jacobian_scaling_values[i] = out_sf / sf
        elif (
            self.jacobian_scaling_type == JacScalingTypes.variable_oi_scaling_square_sum
        ):
            for i, (key, obj) in enumerate(
                self.solver.output_specs.rkt_outputs.items()
            ):
                input_scales = []
                out_sf = iscale.get_scaling_factor(obj.get_pyomo_var(), default=1)
                for input_key, input_obj in self.solver.input_specs.rkt_inputs.items():
                    sf = self.get_rkt_scale(
                        input_obj, use_default_scaling=use_default_scaling
                    )
                    input_scales.append(sf)

                sf = np.sum(np.array(input_scales) ** 2) ** 0.5
                self.solver.jacobian_scaling_values[i] = out_sf / sf
        elif self.jacobian_scaling_type == JacScalingTypes.jacobian_matrix_square_sum:
            jac_matrix = self.get_jacobian_matrix().copy()
            scale_factors = np.sum(np.abs(jac_matrix) ** 2, axis=1) ** 0.5
            scale_factors[scale_factors != 0] = scale_factors[scale_factors != 0] ** -1
            self.solver.jacobian_scaling_values = scale_factors
        elif self.jacobian_scaling_type == JacScalingTypes.jacobian_matrix_inverse_sum:

            jac_matrix = self.get_jacobian_matrix().copy()
            jac_matrix[jac_matrix != 0] = jac_matrix[jac_matrix != 0] ** -1
            scale_factors = np.sum(np.abs(jac_matrix), axis=1)
            scale_factors[scale_factors != 0] = scale_factors[scale_factors != 0] ** -1

            self.solver.jacobian_scaling_values = scale_factors

        max_scale = self.jacobian_scaling_bounds[1]
        min_scale = self.jacobian_scaling_bounds[0]

        for i, scale in enumerate(self.solver.jacobian_scaling_values):
            if self.jacobian_scaling_bounds_output_based:
                mx_multiplier = output_scales[i]
            else:
                mx_multiplier = 1
            if min_scale is not None and scale < mx_multiplier * min_scale:
                self.solver.jacobian_scaling_values[i] = mx_multiplier * min_scale
                _log.warning(
                    f"Jacobian scale for {output_keys[i]} below {min_scale*mx_multiplier }, set to {mx_multiplier * min_scale}"
                )
            if max_scale is not None and scale > mx_multiplier * max_scale:
                self.solver.jacobian_scaling_values[i] = mx_multiplier * max_scale
                _log.warning(
                    f"Jacobian scale for {output_keys[i]} above {max_scale*mx_multiplier }, set to {mx_multiplier * max_scale}"
                )

    def get_jacobian_matrix(self):
        """get jacobian matrix from reaktoro solver"""
        if self.get_jacobian_matrix_function is not None:
            return self.get_jacobian_matrix_function()
        else:
            return self.solver.jacobian_matrix

    def get_jacobian_scaling(self):
        """get jacobian scaling values from reaktoro solver
        generally used parallel manager"""
        if self.update_jacobian_scale_every_solve:
            self.set_jacobian_scaling()
        return self.solver.jacobian_scaling_values

    def get_input_scaling(self):
        """utility function for getting input scaling values from solver
        generally used parallel manager"""
        return self.solver.input_scaling_values

    def set_user_jacobian_scaling(self, user_scaling=None):
        """apply user scaling to jacobian scaling values
        Args:
            user_scaling: dict with keys as output names and values as scaling factors
        """

        if user_scaling is None:
            user_scaling = self.user_scaling
        for i, (key, obj) in enumerate(self.solver.output_specs.rkt_outputs.items()):
            if user_scaling.get(key) != None:
                scale = user_scaling[key]
                self.solver.jacobian_scaling_values[i] = scale

    def display_jacobian_scaling(self):
        """display jacobian scaling values for each output variable"""
        jac_scale = {}
        for i, (key, obj) in enumerate(self.solver.output_specs.rkt_outputs.items()):
            scale = self.solver.jacobian_scaling_values[i]
            _log.info(
                f"Jacobian scale for {key} : {self.solver.jacobian_scaling_values[i]}, IDX: {i}"
            )
            jac_scale[key] = scale
        return jac_scale

    def initialize_input_variables_and_constraints(self, use_default_scaling=True):
        """initialize input variables and constraints"""
        self.solver.input_scaling_values = []
        for key in self.solver.input_specs.rkt_inputs.rkt_input_list:

            if key in self.block.input_constraints:
                calculate_variable_from_constraint(
                    self.block.reaktoro_model.inputs[key],
                    self.block.input_constraints[key],
                )

                sf = self.get_rkt_scale(
                    self.solver.input_specs.rkt_inputs[key],
                    use_default_scaling,
                )
                if self.block.reaktoro_model.inputs[key].value == 0:
                    self.block.reaktoro_model.inputs[key].value = (
                        self.solver.input_specs.rkt_inputs[key].get_value(
                            apply_conversion=True
                        )
                    )

                iscale.set_scaling_factor(self.block.reaktoro_model.inputs[key], sf)
                iscale.constraint_scaling_transform(
                    self.block.input_constraints[key], sf
                )
                self.solver.input_scaling_values.append(sf)

    def display_state(self):
        """display reaktoro state"""
        if self.display_reaktoro_state_function is None:
            self.solver.display_state()
        else:
            self.display_reaktoro_state_function()
