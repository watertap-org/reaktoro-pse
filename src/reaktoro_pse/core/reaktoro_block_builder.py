from pyomo.environ import Var, value
from pyomo.contrib.pynumero.interfaces.external_grey_box import (
    ExternalGreyBoxBlock,
)
import numpy as np
from reaktoro_pse.core.reaktoro_outputs import propTypes

from reaktoro_pse.core.reaktoro_solver import (
    reaktoroSolver,
)

from reaktoro_pse.core.reaktoro_gray_box import (
    reaktoroGrayBox,
)

from pyomo.util.calc_var_value import calculate_variable_from_constraint

import idaes.core.util.scaling as iscale
import cyipopt
import idaes.logger as idaeslog


__author__ = "Alexander Dudchenko"


"""class to build reaktoro block and constraints on provided block"""
_log = idaeslog.getLogger(__name__)


class jacScalingTypes:
    variable_scaling = "variable_scaling"
    jacobian_matrix = "jacobian_matrix"
    manual_scaling = "manual_scaling"


class reaktoroBlockBuilder:
    def __init__(self, block, rktSolver, build_on_init=True):
        self.block = block
        # TODO: add check to make sure block is a pyomo block, or model
        self.rktSolver = rktSolver
        if isinstance(self.rktSolver, reaktoroSolver) == False:
            raise TypeError("Reaktoro block builder requires a reaktoroSolver class")
        self.configure_jacobian_scaling()
        if build_on_init:  # option to support legacy implementation
            self.build_reaktoro_block()

    def build_reaktoro_block(self):
        """build reaktoro model"""
        external_model = reaktoroGrayBox()
        external_model.configure(self.rktSolver)
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
            self.jacobianScalingType = jacScalingTypes.variable_scaling
        else:
            self.jacobianScalingType = jacobian_scaling_type
        if isinstance(user_scaling, float):
            self.rktSolver.jacobianScalingValues = (
                np.ones(len(self.rktSolver.rktOutputSpec.rktOutputs)) + user_scaling
            )
            self.jacobianScalingType = jacScalingTypes.manual_scaling
        elif isinstance(user_scaling, list):
            self.rktSolver.jacobianScalingValues = user_scaling
            self.jacobianScalingType = jacScalingTypes.manual_scaling
        else:
            self.rktSolver.jacobianScalingValues = np.ones(
                len(self.rktSolver.rktOutputSpec.rktOutputs.keys())
            )
        if isinstance(user_scaling, dict):
            self.userScaling = user_scaling
        else:
            self.userScaling = {}

    def build_input_constraints(self):

        if self.rktSolver.rktInputSpec.dissolveSpeciesInRkt:

            @self.block.Constraint(self.rktSolver.rktInputSpec.rktInputs.rktInputList)
            def input_constraints(fs, key):
                return (
                    self.block.reaktoro_model.inputs[key]
                    == self.rktSolver.rktInputSpec.rktInputs[key].pyomoVar
                )

        else:
            """only build these if we are summing species to elements in pyomo"""
            constraint_dict = self.rktSolver.rktInputSpec.constraintDict
            # self.block.inputs = Var(
            #     list(constraint_dict.keys()), initialize=1, bounds=(0, None)
            # )
            self._input_constraint_scaling = {}
            """ connect rektor model vars to our inputs"""
            for element in constraint_dict:
                self.rktSolver.rktInputSpec.rktInputs[element].set_pyomoVar(
                    self.block.reaktoro_model.inputs[element]
                )

            @self.block.Expression(constraint_dict)
            def inputs(fs, element):
                sum_species = []
                # print(self.rktInputSpec.rktInputs.rktInputs.keys())
                for mol, specie in constraint_dict[element]:
                    if specie in self.rktSolver.rktInputSpec.userInputs:
                        pyo_obj = self.rktSolver.rktInputSpec.userInputs[
                            specie
                        ].pyomoVar
                    elif specie in self.rktSolver.rktInputSpec.rktChemicalInputs:
                        pyo_obj = self.rktSolver.rktInputSpec.rktChemicalInputs[
                            specie
                        ].pyomoVar
                    else:
                        raise KeyError(f"specie {specie} not found in input dicts")
                    sum_species.append(mol * pyo_obj)
                return sum(sum_species)

            for element in constraint_dict:
                self.rktSolver.rktInputSpec.rktInputs[element].set_pyomoVar(
                    self.block.reaktoro_model.inputs[element]
                )

            @self.block.Constraint(self.rktSolver.rktInputSpec.rktInputs.rktInputList)
            def input_constraints(fs, key):
                if key in constraint_dict:
                    return (
                        self.block.reaktoro_model.inputs[key] == self.block.inputs[key]
                    )
                else:
                    return (
                        self.block.reaktoro_model.inputs[key]
                        == self.rktSolver.rktInputSpec.userInputs[key].pyomoVar
                    )

    def build_output_constraints(self):
        """first update rktOuptutObjects for pyomoBuildProperties with reaktoro pyomo variables as
        they will be used in construction of constraints
        The is will also check if user provided an output pyomo var and if not will
        add them to new_output_var dict, which will be used to create new output variables on the block
        """
        new_output_vars = {}
        for key, obj in self.rktSolver.rktOutputSpec.userOutputs.items():
            if propTypes.pyomoBuiltProperties == obj.propertyType:
                for (
                    pyoPropKey,
                    pyoPropObj,
                ) in obj.pyomoBuildOptions.properties.items():
                    pyoPropObj.set_pyomo_var(
                        self.block.reaktoro_model.outputs[pyoPropKey]
                    )
            # NOTE: We do not set rktOutputs to reaktoro_model outputs as they
            # same as user inputs - we want RKt model to update "user provided vars"
            # rather then pyomo vars in reaktoro model (e.g. reaktor_block.outputs)
            if obj.pyomoVar is None:
                new_output_vars[key] = obj
        if new_output_vars != {}:
            self.block.outputs = Var(new_output_vars.keys(), initialize=1)
            for key, obj in new_output_vars.items():
                obj.set_pyomo_var(self.block.outputs[key])

        @self.block.Constraint(self.rktSolver.rktOutputSpec.userOutputs)
        def output_constraints(fs, prop, prop_index):
            prop_object = self.rktSolver.rktOutputSpec.userOutputs[(prop, prop_index)]
            if prop_object.propertyType == propTypes.pyomoBuiltProperties:
                return prop_object.pyomoBuildOptions.build_constraint_function(
                    prop_object
                )
            else:
                return (
                    prop_object.pyomoVar
                    == self.block.reaktoro_model.outputs[(prop, prop_index)]
                )

    def initialize(self, presolveDuringInitialization=False):
        self.initialize_input_variables_and_constraints()
        self.rktSolver.rktBase.equilibrate_state()
        self.rktSolver.solve_reaktoro_block(presolve=presolveDuringInitialization)
        self.initialize_output_variables_and_constraints()
        _log.info(f"Initialized rkt block")

    def initialize_output_variables_and_constraints(self):
        def get_sf(val):
            if val == 0:
                return 1
            else:
                return abs(val)

        for key, obj in self.rktSolver.rktOutputSpec.userOutputs.items():
            """update vars scaling in pyomo build ocnstraints
            these are updated to actual value when we call solve_rektoro_block"""
            if propTypes.pyomoBuiltProperties == obj.propertyType:
                for (
                    pyoPropKey,
                    pyoPropObj,
                ) in obj.pyomoBuildOptions.properties.items():
                    val = pyoPropObj.value
                    pyoPropObj.set_pyomo_var_value(val)
                    if iscale.get_scaling_factor(pyoPropObj.pyomoVar) is None:
                        iscale.set_scaling_factor(pyoPropObj.pyomoVar, 1 / get_sf(val))
                output_constraint = self.block.output_constraints[key]
                calculate_variable_from_constraint(obj.pyomoVar, output_constraint)
                iscale.constraint_scaling_transform(
                    output_constraint, 1 / get_sf(obj.pyomoVar.value)
                )
            else:
                obj.set_pyomo_var_value(obj.value)
                rkt_var = self.block.reaktoro_model.outputs[key]
                output_constraint = self.block.output_constraints[key]
                # obj.pyomo_var.display()
                # rkt_var.display()
                calculate_variable_from_constraint(rkt_var, output_constraint)

                if iscale.get_scaling_factor(rkt_var) is None:
                    iscale.set_scaling_factor(rkt_var, 1 / get_sf(obj.value))
            val = obj.pyomoVar.value
            """ scale user provided vars if they are not scaled"""
            if iscale.get_scaling_factor(obj.pyomoVar) is None:
                iscale.set_scaling_factor(obj.pyomoVar, 1 / get_sf(val))

        """ update jacobian scaling """
        self.get_jacobian_scaling()

    def get_jacobian_scaling(self):
        if self.jacobianScalingType == jacScalingTypes.variable_scaling:
            for i, (key, obj) in enumerate(
                self.rktSolver.rktOutputSpec.rktOutputs.items()
            ):
                out_sf = iscale.get_scaling_factor(obj.pyomoVar, default=1)
                self.rktSolver.jacobianScalingValues[i] = 1 / out_sf
        elif self.jacobianScalingType == jacScalingTypes.jacobian_matrix:
            self.rktSolver.jacobianScalingValues = (
                np.sum(np.abs(self.jacobian_matrix) ** 2, axis=1) ** 0.5
            )
        for i, (key, obj) in enumerate(self.rktSolver.rktOutputSpec.rktOutputs.items()):
            if key in self.userScaling:
                scale = self.userScaling[key]
                self.rktSolver.jacobianScalingValues[i] = scale

        _log.info(
            f"Jacobian output order: {list(self.rktSolver.rktOutputSpec.rktOutputs.keys())}"
        )
        _log.info(f"Jacobian scaleing: {self.rktSolver.jacobianScalingValues}")

    def initialize_input_variables_and_constraints(self):
        """intialize input variables and constraints"""
        if self.rktSolver.rktInputSpec.dissolveSpeciesInRkt == False:
            for element in self.rktSolver.rktInputSpec.constraintDict:
                val = value(self.block.inputs[element])
                self.rktSolver.rktInputSpec.rktInputs[element].set_pyomo_var_value(val)
        for key in self.rktSolver.rktInputSpec.rktInputs.rktInputList:
            pyo_var = self.rktSolver.rktInputSpec.rktInputs[key].pyomoVar
            self.block.reaktoro_model.inputs[key] = pyo_var.value
            sf = iscale.get_scaling_factor(pyo_var)
            if sf == None:
                sf = abs(1 / pyo_var.value)

            iscale.set_scaling_factor(self.block.reaktoro_model.inputs[key], sf)
            iscale.constraint_scaling_transform(self.block.input_constraints[key], sf)
