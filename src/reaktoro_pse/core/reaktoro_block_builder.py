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
from pyomo.environ import Var, value, units as pyunits
from pyomo.contrib.pynumero.interfaces.external_grey_box import (
    ExternalGreyBoxBlock,
)
import numpy as np
from reaktoro_pse.core.reaktoro_outputs import PropTypes

from reaktoro_pse.core.reaktoro_solver import (
    ReaktoroSolver,
)

from reaktoro_pse.core.reaktoro_gray_box import (
    ReaktoroGrayBox,
)

from pyomo.util.calc_var_value import calculate_variable_from_constraint

import idaes.core.util.scaling as iscale
import cyipopt
import idaes.logger as idaeslog


__author__ = "Alexander Dudchenko"


"""class to build reaktoro block and constraints on provided block"""
_log = idaeslog.getLogger(__name__)


class JacScalingTypes:
    no_scaling = "no_scaling"
    variable_scaling = "variable_scaling"
    jacobian_matrix = "jacobian_matrix"
    manual_scaling = "manual_scaling"


class ReaktoroBlockBuilder:
    def __init__(self, block, reaktoro_solver, build_on_init=True):
        self.block = block
        # TODO: add check to make sure block is a pyomo block, or model
        self.solver = reaktoro_solver
        if isinstance(self.solver, ReaktoroSolver) == False:
            raise TypeError("Reaktoro block builder requires a ReaktoroSolver class")
        self.configure_jacobian_scaling()
        if build_on_init:  # option to support legacy implementation
            self.build_reaktoro_block()

    def build_reaktoro_block(self):
        """build reaktoro model"""
        external_model = ReaktoroGrayBox()
        external_model.configure(self.solver)
        self.block.reaktoro_model = ExternalGreyBoxBlock(external_model=external_model)
        self.build_input_constraints()
        self.build_output_constraints()

    def configure_jacobian_scaling(self, jacobian_scaling_type=None, user_scaling=None):
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

    def build_input_constraints(self):

        if self.solver.input_specs.dissolve_species_in_rkt:

            @self.block.Constraint(self.solver.input_specs.rkt_inputs.rkt_input_list)
            def input_constraints(fs, key):
                return (
                    self.block.reaktoro_model.inputs[key]
                    == self.solver.input_specs.rkt_inputs[
                        key
                    ].get_pyomo_with_required_units()
                )

        else:
            """only build these if we are summing species to elements in pyomo"""
            constraint_dict = self.solver.input_specs.constraint_dict
            self._input_constraint_scaling = {}
            """ connect rektor model vars to our inputs"""
            for element in constraint_dict:
                self.solver.input_specs.rkt_inputs[element].set_pyomo_var(
                    self.block.reaktoro_model.inputs[element]
                )

            @self.block.Expression(constraint_dict)
            def inputs(fs, element):
                sum_species = []
                for mol, specie in constraint_dict[element]:
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
                    sum_species.append(mol * pyo_obj)
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

    def build_output_constraints(self):
        """first update rktOuptutObjects for pyomoBuildProperties with reaktoro pyomo variables as
        they will be used in construction of constraints
        The is will also check if user provided an output pyomo var and if not will
        add them to new_output_var dict, which will be used to create new output variables on the block
        """
        new_output_vars = {}
        for key, obj in self.solver.output_specs.user_outputs.items():
            if PropTypes.pyomo_built_prop == obj.property_type:
                for (
                    pyoPropKey,
                    pyoPropObj,
                ) in obj.pyomo_build_options.properties.items():
                    pyoPropObj.set_pyomo_var(
                        self.block.reaktoro_model.outputs[pyoPropKey]
                    )
            # NOTE: We do not set rkt_outputs to reaktoro_model outputs as they
            # same as user inputs - we want RKt model to update "user provided vars"
            # rather then pyomo vars in reaktoro model (e.g. reaktor_block.outputs)
            if obj.get_pyomo_var() is None:
                new_output_vars[key] = obj
        if new_output_vars != {}:
            self.block.outputs = Var(new_output_vars.keys(), initialize=1)
            for key, obj in new_output_vars.items():
                obj.set_pyomo_var(self.block.outputs[key])

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

    def initialize(self, presolveDuringInitialization=False):
        self.initialize_input_variables_and_constraints()
        self.solver.state.equilibrate_state()
        self.solver.solve_reaktoro_block(presolve=presolveDuringInitialization)
        self.initialize_output_variables_and_constraints()
        _log.info(f"Initialized rkt block")

    def get_sf(self, pyo_var):
        if pyo_var.value == 0:
            return 1
        else:
            if iscale.get_scaling_factor(pyo_var) is not None:
                return iscale.get_scaling_factor(pyo_var)

            sf = 1 / abs(pyo_var.value)
            if sf > 1e16:
                _log.warning(f"Var {pyo_var} scale >1e16")
                # sf = 1e16
            if sf < 1e-16:
                _log.warning(f"Var {pyo_var} scale <1e-16")
                # sf = 1e-16
            return sf

    def initialize_output_variables_and_constraints(self):
        for key, obj in self.solver.output_specs.user_outputs.items():
            """update vars scaling in pyomo build ocnstraints
            these are updated to actual value when we call solve_rektoro_block"""
            if PropTypes.pyomo_built_prop == obj.property_type:
                for (
                    pyoPropKey,
                    pyoPropObj,
                ) in obj.pyomo_build_options.properties.items():
                    val = pyoPropObj.value
                    pyoPropObj.set_pyomo_var_value(val)
                    if iscale.get_scaling_factor(pyoPropObj.get_pyomo_var()) is None:
                        iscale.set_scaling_factor(
                            pyoPropObj.get_pyomo_var(),
                            self.get_sf(pyoPropObj.get_pyomo_var()),
                        )
                output_constraint = self.block.output_constraints[key]
                calculate_variable_from_constraint(
                    obj.get_pyomo_var(), output_constraint
                )
                iscale.constraint_scaling_transform(
                    output_constraint, self.get_sf(obj.get_pyomo_var())
                )
            else:
                obj.set_pyomo_var_value(obj.value)
                rkt_var = self.block.reaktoro_model.outputs[key]
                output_constraint = self.block.output_constraints[key]
                calculate_variable_from_constraint(rkt_var, output_constraint)
                iscale.constraint_scaling_transform(
                    output_constraint, self.get_sf(obj.get_pyomo_var())
                )

                if iscale.get_scaling_factor(rkt_var) is None:
                    iscale.set_scaling_factor(rkt_var, self.get_sf(obj.get_pyomo_var()))

            """ scale user provided vars if they are not scaled"""
            if iscale.get_scaling_factor(obj.get_pyomo_var()) is None:
                iscale.set_scaling_factor(
                    obj.get_pyomo_var(), self.get_sf(obj.get_pyomo_var())
                )

        """ update jacobian scaling """
        self.get_jacobian_scaling()

    def get_jacobian_scaling(self):
        if self.jacobian_scaling_type == JacScalingTypes.no_scaling:
            for i, (key, obj) in enumerate(
                self.solver.output_specs.rkt_outputs.items()
            ):
                self.solver.jacobian_scaling_values[i] = 1
        elif self.jacobian_scaling_type == JacScalingTypes.variable_scaling:
            for i, (key, obj) in enumerate(
                self.solver.output_specs.rkt_outputs.items()
            ):
                out_sf = iscale.get_scaling_factor(obj.get_pyomo_var(), default=1)
                self.solver.jacobian_scaling_values[i] = 1 / out_sf
        elif self.jacobian_scaling_type == JacScalingTypes.jacobian_matrix:
            self.solver.jacobian_scaling_values = (
                np.sum(np.abs(self.solver.jacobian_matrix) ** 2, axis=1) ** 0.5
            )
        self.set_user_jacobian_scaling()

    def set_user_jacobian_scaling(self, user_scaling=None):
        if user_scaling is None:
            user_scaling = self.user_scaling

        for i, (key, obj) in enumerate(self.solver.output_specs.rkt_outputs.items()):
            if key in user_scaling:
                scale = user_scaling[key]
                self.solver.jacobian_scaling_values[i] = scale

    def display_jacobian_scaling(self):
        jac_scale = {}
        for i, (key, obj) in enumerate(self.solver.output_specs.rkt_outputs.items()):
            scale = self.solver.jacobian_scaling_values[i]
            _log.info(
                f"Jacobian scale for {key} : {self.solver.jacobian_scaling_values[i]}, IDX: {i}"
            )
            jac_scale[key] = scale
        return jac_scale

    def initialize_input_variables_and_constraints(self):
        """initialize input variables and constraints"""
        # if self.solver.input_specs.dissolve_species_in_rkt == False:
        #     for element in self.solver.input_specs.constraint_dict:
        #         calculate_variable_from_constraint(
        #             self.block.reaktoro_model.inputs[element],
        #             self.block.input_constraints[element],
        #         )
        #         val = value(self.block.reaktoro_model.inputs[element])
        #         self.solver.input_specs.rkt_inputs[element].set_pyomo_var_value(val)
        for key in self.solver.input_specs.rkt_inputs.rkt_input_list:
            pyo_var = self.solver.input_specs.rkt_inputs[key].get_pyomo_var()
            calculate_variable_from_constraint(
                self.block.reaktoro_model.inputs[key],
                self.block.input_constraints[key],
            )
            # self.solver.input_specs.rkt_inputs[key].set_pyomo_var_value(
            #     self.block.reaktoro_model.inputs[key].value
            # )
            sf = iscale.get_scaling_factor(pyo_var)
            if sf == None:
                sf = self.get_sf(self.block.reaktoro_model.inputs[key])
            iscale.set_scaling_factor(self.block.reaktoro_model.inputs[key], sf)
            iscale.constraint_scaling_transform(self.block.input_constraints[key], sf)
