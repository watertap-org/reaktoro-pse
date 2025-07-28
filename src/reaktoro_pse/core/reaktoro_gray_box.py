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

import pyomo.environ as pyo
from pyomo.contrib.pynumero.interfaces.external_grey_box import (
    ExternalGreyBoxModel,
)
import numpy as np
from scipy.sparse import coo_matrix
import copy
import idaes.logger as idaeslog
from reaktoro_pse.core.util_classes.hessian_functions import (
    HessianApproximation,
    HessTypes,
    BFGSInitializationTypes,
)

__author__ = "Ilayda Akkor, Alexander V. Dudchenko, Paul Vecchiarelli, Ben Knueven"
_log = idaeslog.getLogger(__name__)


import idaes.core.util.scaling as iscale


class ReaktoroGrayBox(ExternalGreyBoxModel):
    ########################################################################################
    # custom Grey Box functions
    def configure(
        self,
        reaktoro_solver=None,
        inputs=None,
        input_dict=None,
        outputs=None,
        hessian_type=None,
        bfgs_initialization_type=BFGSInitializationTypes.GaussNewton,
        bfgs_init_min_hessian_value=1e-32,
        bfgs_init_max_hessian_value=1e8,
        bfgs_init_const_hessian_value=1e-16,
        bfgs_hessian_memory=3,
        bfgs_epsilon=1e-12,
    ):
        # assign a Reaktoro state object to instance
        self.reaktoro_solver = reaktoro_solver
        if hessian_type is None:
            hess_type = reaktoro_solver.hessian_type
            bfgs_hessian_memory = reaktoro_solver.bfgs_hessian_memory
            bfgs_init_const_hessian_value = (
                reaktoro_solver.bfgs_init_const_hessian_value
            )
            bfgs_init_min_hessian_value = reaktoro_solver.bfgs_init_min_hessian_value
            bfgs_init_max_hessian_value = reaktoro_solver.bfgs_init_max_hessian_value
            bfgs_initialization_type = reaktoro_solver.bfgs_initialization_type
            bfgs_epsilon = reaktoro_solver.bfgs_epsilon
        else:
            hess_type = hessian_type
        if inputs is None:
            self.inputs = reaktoro_solver.input_specs.rkt_inputs.rkt_input_list
        else:
            self.inputs = inputs
        if input_dict is None:
            self.input_dict = reaktoro_solver.input_specs.rkt_inputs
        else:
            self.input_dict = input_dict
        if outputs is None:
            self.outputs = list(reaktoro_solver.output_specs.rkt_outputs.keys())
        else:
            self.outputs = outputs
        self._outputs_dual_multipliers = np.ones(len(self.outputs))
        self.header_saved = False
        self.step = 0
        self.old_params = None

        _log.info(f"RKT gray box using {hess_type} hessian type")
        if hess_type != HessTypes.no_hessian_estimation:
            self.hessian_calculator = HessianApproximation(
                hessian_type=hess_type,
                bfgs_hessian_memory=bfgs_hessian_memory,
                bfgs_init_const_hessian_value=bfgs_init_const_hessian_value,
                bfgs_init_min_hessian_value=bfgs_init_min_hessian_value,
                bfgs_init_max_hessian_value=bfgs_init_max_hessian_value,
                bfgs_initialization_type=bfgs_initialization_type,
                bfgs_epsilon=bfgs_epsilon,
            )
            setattr(self, "evaluate_hessian_outputs", self._evaluate_hessian_outputs)

    ########################################################################################
    # standard Grey Box functions
    def input_names(self):
        # get input names (required by Grey Box)
        return self.inputs

    def output_names(self):
        # get output names (not required, but helpful)
        return self.outputs

    def set_input_values(self, input_values):
        # set input values from Pyomo as inputs to External Model (required by Grey Box)
        self._input_values = list(input_values)

    def finalize_block_construction(self, pyomo_block):
        # initialize Pyomo block for External Model
        block_components = [obj for obj in pyomo_block.component_objects(pyo.Var)]
        for block in block_components:
            if "inputs" in block.name:
                for var in self.inputs:
                    block[var].value = 1  # self.input_dict[var].get_pyomo_var_value()
                    block[var].setlb(self.input_dict[var].get_lower_bound())
                    block[var].setub(self.input_dict[var].get_upper_bound())
            elif "outputs" in block.name:
                for prop in self.outputs:
                    block[prop].setlb(None)
                    block[prop].setub(None)
                    block[prop].value = 0.1

    def evaluate_outputs(self):
        # update Reaktoro state with current inputs (this function runs repeatedly)
        self.params = dict(zip(self.inputs, self._input_values))
        self.get_last_output(self.params)
        return np.array(self.rkt_result, dtype=np.float64)

    def get_last_output(self, new_params):
        """only eval reaktoro if params changed!"""
        if self.old_params is None or any(
            new_params[key] != self.old_params[key] for key in new_params
        ):
            self.jacobian_matrix, self.rkt_result = (
                self.reaktoro_solver.solve_reaktoro_block(params=new_params)
            )
            self.step += 1
        self.old_params = copy.deepcopy(new_params)

    def evaluate_jacobian_outputs(self):
        self.evaluate_outputs()
        jm = np.array(self.jacobian_matrix)
        i = np.array([i for i in range(jm.shape[0]) for j in range(jm.shape[1])])
        j = np.array([j for i in range(jm.shape[0]) for j in range(jm.shape[1])])

        cm = coo_matrix((jm.flatten(), (i, j)))
        return cm

    def set_output_constraint_multipliers(self, _outputs_dual_multipliers):
        np.copyto(self._outputs_dual_multipliers, _outputs_dual_multipliers)

    def get_output_constraint_scaling_factors(self):

        return self.reaktoro_solver.get_jacobian_scaling()

    def _evaluate_hessian_outputs(self):
        """Evaluate the Hessian matrix of the outputs with respect to the inputs."""

        low_triangular_hessian = self.hessian_calculator.get_hessian(
            self._input_values,
            self.jacobian_matrix,
            self._outputs_dual_multipliers,
        )

        return low_triangular_hessian
