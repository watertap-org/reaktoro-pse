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
import numpy as np
from reaktoro_pse.core.util_classes.rkt_inputs import (
    RktInputs,
    RktInputTypes,
    DummyPyomoVar,
)
from reaktoro_pse.core.reaktoro_outputs import (
    ReaktoroOutputSpec,
)
from reaktoro_pse.core.reaktoro_inputs import (
    ReaktoroInputSpec,
)
import idaes.logger as idaeslog

__author__ = "Alexander V. Dudchenko, Ben Knueven"

_log = idaeslog.getLogger(__name__)


class ReaktoroCoupledSolver:
    """Couples speciation and property reaktoro solvers into a single solver.
    This in general mimics regular Reaktoro solver but couples the
    speciation solver with a property solver."""

    def __init__(self, speciation_solvers):
        """
        Args:
            speciation_solvers: single or list of speciation solvers to couple with property solver
        """
        if not isinstance(speciation_solvers, list):
            speciation_solvers = [speciation_solvers]
        self.speciation_solvers = speciation_solvers
        self.get_master_inputs()
        self.get_speciation_outputs()

    def register_property_solver(self, property_solver):
        """
        Register a property solver to the coupled solver.
        Args:
            property_solver: a ReaktoroPropertySolver instance to couple with speciation solvers
        """
        self.property_solver = property_solver
        self.prop_inputs_idx = []
        self.prop_jac_idx = []
        # used to store property indexes
        self.prop_jac_propagation_idx = np.zeros(len(self.output_key_order), dtype=int)
        # used to store property keys for propagation
        self.prop_jac_propagation_keys = np.empty(
            len(self.output_key_order), dtype=object
        )
        # get user inputs into property block, we will track these manually in the
        # coupled solver
        for key, obj in self.property_solver.input_specs.user_inputs.items():
            if (
                obj.io_type != RktInputTypes.specie
                and obj.io_type != RktInputTypes.element
            ):
                new_key = self.modify_key("prop", key)
                self.input_specs.user_inputs[new_key] = obj
                self.master_mapping[new_key] = key
        # get all inputs.
        for idx, key in enumerate(
            self.property_solver.input_specs.rkt_inputs.rkt_input_list
        ):
            if key in self.property_solver.input_specs.rkt_inputs:
                obj = self.property_solver.input_specs.rkt_inputs[key]
                if (
                    obj.io_type != RktInputTypes.specie
                    and obj.io_type != RktInputTypes.element
                    or obj.dummy_var_key == None
                ):
                    new_key = self.modify_key("prop", key)
                    self.input_specs.rkt_inputs[new_key] = obj
                    self.master_mapping[new_key] = key
                    self.update_input_list(new_key)
                    self.prop_jac_idx.append(idx)
                else:
                    self.prop_jac_propagation_idx[
                        self.output_key_order[obj.dummy_var_key]
                    ] = idx
                    self.prop_jac_propagation_keys[
                        self.output_key_order[obj.dummy_var_key]
                    ] = key
            elif key in self.property_solver.input_specs.rkt_chemical_inputs:
                obj = self.property_solver.input_specs.rkt_chemical_inputs[key]
                if (
                    obj.io_type != RktInputTypes.specie
                    and obj.io_type != RktInputTypes.element
                    or obj.dummy_var_key == None
                ):
                    new_key = self.modify_key("prop", key)
                    self.input_specs.rkt_chemical_inputs[new_key] = obj
                    self.master_mapping[new_key] = key
                    self.update_input_list(new_key)
                    self.prop_jac_idx.append(idx)
                else:
                    self.prop_jac_propagation_idx[
                        self.output_key_order[obj.dummy_var_key]
                    ] = idx
                    self.prop_jac_propagation_keys[
                        self.output_key_order[obj.dummy_var_key]
                    ] = key
        # rkt_inputs are not used in property solver, so we do not register them
        self.get_master_outputs()
        for key in self.prop_jac_propagation_keys:
            new_key = self.modify_key("prop", key)
            self.input_specs.user_inputs[new_key] = (
                self.property_solver.input_specs.user_inputs[key]
            )
            self.master_mapping[new_key] = key
        # register hessian options for proerty solver, these will be used
        # as master options for our graybox
        self.hessian_type = self.property_solver.hessian_type
        self.bfgs_initialization_type = self.property_solver.bfgs_initialization_type
        self.bfgs_init_min_hessian_value = (
            self.property_solver.bfgs_init_min_hessian_value
        )
        self.bfgs_init_max_hessian_value = (
            self.property_solver.bfgs_init_max_hessian_value
        )
        self.bfgs_init_const_hessian_value = (
            self.property_solver.bfgs_init_const_hessian_value
        )
        self.bfgs_hessian_memory = self.property_solver.bfgs_hessian_memory
        self.bfgs_epsilon = self.property_solver.bfgs_epsilon

    def modify_key(self, index, key):
        """utility for generating single tuple key"""
        new_index = [index]
        if isinstance(key, str):
            new_index.append(key)
        elif isinstance(key, (list, tuple)):
            for k in key:
                new_index.append(k)
        return tuple(new_index)

    def update_input_list(self, new_key):
        # update input list for master inputs
        if new_key not in self.input_specs.rkt_inputs.rkt_input_list:
            self.input_specs.rkt_inputs.rkt_input_list.append(new_key)

    def get_master_inputs(self):
        # create master inputs for all speciation solvers
        self.input_specs = ReaktoroInputSpec()
        self.input_specs.user_inputs = RktInputs()
        self.input_specs.rkt_chemical_inputs = RktInputs()
        self.input_specs.rkt_inputs = RktInputs()
        self.master_mapping = {}
        self.speciation_jac_idx = {}

        def update_input_list(new_key):
            if new_key not in self.input_specs.rkt_inputs.rkt_input_list:
                self.input_specs.rkt_inputs.rkt_input_list.append(new_key)

        for i, solver in enumerate(self.speciation_solvers):
            self.speciation_jac_idx[i] = []
            for key, obj in solver.input_specs.user_inputs.items():
                new_key = self.modify_key(f"s_{i}", key)
                self.master_mapping[new_key] = key
                self.input_specs.user_inputs[new_key] = obj
            for idx, key in enumerate(solver.input_specs.rkt_inputs.rkt_input_list):
                if key in solver.input_specs.rkt_inputs:
                    new_key = self.modify_key(f"s_{i}", key)
                    self.master_mapping[new_key] = key
                    self.input_specs.rkt_inputs[new_key] = (
                        solver.input_specs.rkt_inputs[key]
                    )
                    self.update_input_list(new_key)
                elif key in solver.input_specs.rkt_chemical_inputs:
                    new_key = self.modify_key(f"s_{i}", key)
                    self.master_mapping[new_key] = key
                    self.input_specs.rkt_chemical_inputs[new_key] = (
                        solver.input_specs.rkt_chemical_inputs[key]
                    )
                    self.update_input_list(new_key)
                self.speciation_jac_idx[i].append(idx)
        self.input_specs.dissolve_species_in_rkt = self.speciation_solvers[
            0
        ].input_specs.dissolve_species_in_rkt
        self.input_specs.exact_speciation = self.speciation_solvers[
            0
        ].input_specs.exact_speciation

    def get_speciation_outputs(self):
        self.outputs = {}
        self.output_key_order = {}
        self.output_specs = ReaktoroOutputSpec()
        for idx, output in enumerate(
            self.speciation_solvers[0].output_specs.rkt_outputs
        ):
            self.outputs[output] = DummyPyomoVar()
            self.outputs[output].original_key = output
            self.output_key_order[output] = idx

    def get_master_outputs(self):
        """create master inputs/outputs for property solver"""
        self.output_specs = ReaktoroOutputSpec()
        self.output_specs.user_outputs = self.property_solver.output_specs.user_outputs
        self.output_specs.rkt_outputs = self.property_solver.output_specs.rkt_outputs

    def equilibrate_state(
        self,
    ):
        """Initialize all reaktoro states"""
        for solver in self.speciation_solvers:
            solver.equilibrate_state()
        self.prop_block_not_equilibrated = True

    def display_state(self):
        """Display the current state of the Reaktoro coupled solver."""
        _log.info("Reaktoro Coupled Solver State:")
        for i, solver in enumerate(self.speciation_solvers):
            _log.info(f"Speciation Solver {i}:")
            solver.display_state()
        _log.info("Property Solver State:")
        self.property_solver.display_state()

    def equilibrate_property_state(self):
        """Equilibrate the property solver state."""
        if self.prop_block_not_equilibrated:
            self.property_solver.equilibrate_state()
        self.prop_block_not_equilibrated = False

    def propagate_speciation_outputs(self):
        """Propagate outputs from speciation solvers to the inputs for main prop block."""
        for output, obj in self.outputs.items():
            sum_element = 0
            for i in self.speciation_solvers:
                sum_element += i.output_specs.rkt_outputs[output].value
            self.outputs[output].set_value(sum_element)

    def compute_combined_jacobian(self, speciation_jacs, property_jac):
        """manually propagate speciation jacobians to property jacobian"""
        self.jacobian_matrix = np.zeros(
            (
                len(self.output_specs.rkt_outputs),
                len(self.input_specs.rkt_inputs.rkt_input_list),
            )
        )

        prop_prop_jac = property_jac.T[[self.prop_jac_propagation_idx]][0].T
        end_idx = 0
        start_idx = 0
        for i, spc in enumerate(self.speciation_solvers):
            spec_prop_jack = prop_prop_jac @ speciation_jacs[i]
            end_idx += speciation_jacs[i].shape[1]
            self.jacobian_matrix[:, start_idx:end_idx] = spec_prop_jack
            start_idx = end_idx
        sub_prop_jac = property_jac.T[[self.prop_jac_idx]][0].T
        self.jacobian_matrix[:, end_idx:] = sub_prop_jac

    def solve_reaktoro_block(self, params=None, presolve=False):
        """
        Solve the coupled Reaktoro block.
        Args:
            params: dictionary of parameters to update the inputs, if None, use current inputs
            presolve: boolean flag to indicate if presolve should be used"""
        if params is None:
            use_temp = False
        else:
            use_temp = True
        self.update_inputs(params)
        outputs, jacs = {}, {}
        for i, solver in enumerate(self.speciation_solvers):
            jac, out = solver.solve_reaktoro_block(presolve=presolve, use_temp=use_temp)
            outputs[i] = out
            jacs[i] = jac
        self.propagate_speciation_outputs()
        self.equilibrate_property_state()
        jac, outputs = self.property_solver.solve_reaktoro_block(
            presolve=presolve, use_temp=use_temp
        )

        self.compute_combined_jacobian(jacs, jac)
        return self.jacobian_matrix, outputs

    def update_inputs(self, params):
        """Propagate inputs to all speciation solvers."""
        if params != None:
            for key, value in params.items():
                self.input_specs.rkt_inputs[key].set_temp_value(value)

    def get_jacobian_scaling(self):
        raise NotImplementedError(
            "This method gets updated by ReaktoroBlockBuilder, did you build the builder?"
        )

    def get_input_scaling(self):
        raise NotImplementedError(
            "This method gets updated by ReaktoroBlockBuilder, did you build the builder?"
        )
