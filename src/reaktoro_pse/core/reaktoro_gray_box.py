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
import reaktoro as rkt
import pyomo.environ as pyo
from pyomo.contrib.pynumero.interfaces.external_grey_box import (
    ExternalGreyBoxModel,
)
import numpy as np
from scipy.sparse import coo_matrix, tril

import copy
import idaes.logger as idaeslog

__author__ = "Ilayda Akkor, Alexander Dudchenko, Paul Vecchiarelli, Ben Knueven"
_log = idaeslog.getLogger(__name__)


class HessTypes:
    JtJ = "J.tJ"
    BFGS = "BFGS"
    BFGS_mod = "BFGS-mod"
    BFGS_damp = "BFGS-damp"
    BFGS_ipopt = "BFGS-ipopt"


class ReaktoroGrayBox(ExternalGreyBoxModel):
    ########################################################################################
    # custom Grey Box functions
    def configure(self, reaktoro_solver):
        # assign a Reaktoro state object to instance
        self.reaktoro_solver = reaktoro_solver
        self.hess_type = reaktoro_solver.hessian_type
        self.inputs = reaktoro_solver.input_specs.rkt_inputs.rkt_input_list
        self.input_dict = reaktoro_solver.input_specs.rkt_inputs
        self.outputs = list(reaktoro_solver.output_specs.rkt_outputs.keys())
        self._outputs_dual_multipliers = np.ones(len(self.outputs))
        self._hess = np.zeros((len(self.inputs), len(self.inputs)))
        self.header_saved = False
        self.step = 0
        self.old_params = None

        _log.info(f"RKT gray box using {self.hess_type} hessian type")

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
            # 1e-16 is Reaktoro's epsilon value
            if "inputs" in block.name:
                for var in self.inputs:
                    block[var].value = 1  # self.input_dict[var].get_pyomo_var_value()
                    block[var].setlb(self.input_dict[var].get_lower_bound())
                    block[var].setub(None)
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
        return self.reaktoro_solver.jacobian_scaling_values

    def hessian_gauss_newton_version(self, sparse_jac, threshold=7):

        s = np.zeros((len(self.inputs), len(self.inputs)))
        for i in range(self.jacobian_matrix.shape[0]):
            row = self.jacobian_matrix[i, :]
            if sparse_jac:
                row = np.round(row, decimals=threshold)
            s += self._outputs_dual_multipliers[i] * np.outer(row.T, row)
        hess = s
        return hess

    def hessian_bfgs(self):
        # Cautious BFGS update implementation (Li and Fukushima)
        if not hasattr(self, "x"):
            self.x = None
        if not hasattr(self, "H"):
            self.H = []
            H_i = np.identity(len(self.inputs))
            for i in range(self.jacobian_matrix.shape[0]):
                self.H.append(H_i.copy())
        if self.x is not None:
            s_k = (np.array([self._input_values]) - self.x).T
            alpha = 1
            eps = 1e-6
            for i in range(self.jacobian_matrix.shape[0]):
                y_k = np.array([self.jacobian_matrix[i, :] - self.del_f[i, :]]).T
                y_s = y_k.T @ s_k
                H_s = self.H[i] @ s_k
                w = (np.linalg.norm(s_k) ** 2) * (
                    np.linalg.norm(self.del_f[i, :]) ** alpha
                )

                if y_s > eps * w:
                    self.H[i] = (
                        self.H[i]
                        + (y_k @ y_k.T) / (y_s)
                        - (H_s @ H_s.T) / (s_k.T @ H_s)
                    )

        h_sum = np.zeros((len(self.inputs), len(self.inputs)))
        for i in range(self.jacobian_matrix.shape[0]):
            h_sum += self._outputs_dual_multipliers[i] * self.H[i]

        # print(np.linalg.cond(h_sum))
        self.x = self._input_values
        self.del_f = self.jacobian_matrix
        self.hessian = h_sum
        return self.hessian

    def hessian_modified_bfgs(self):
        # Modified BFGS update implementation (Li and Fukushima)
        if not hasattr(self, "x"):
            self.x = None
        if not hasattr(self, "H"):
            self.H = []
            H_i = np.identity(len(self.inputs))
            for i in range(self.jacobian_matrix.shape[0]):
                self.H.append(H_i.copy())
        if self.x is not None:
            s_k = (np.array([self._input_values]) - self.x).T
            for i in range(self.jacobian_matrix.shape[0]):
                y_k = np.array([self.jacobian_matrix[i, :] - self.del_f[i, :]]).T
                y_s = y_k.T @ s_k
                H_s = self.H[i] @ s_k
                if s_k.any():
                    t_k = 1 + max(0, -y_s / (np.linalg.norm(s_k) ** 2))
                    z_k = y_k + t_k * np.linalg.norm(self.del_f[i, :]) * s_k
                    self.H[i] = (
                        self.H[i]
                        + (z_k @ z_k.T) / (z_k.T @ s_k)
                        - (H_s @ H_s.T) / (s_k.T @ H_s)
                    )

        h_sum = np.zeros((len(self.inputs), len(self.inputs)))
        for i in range(self.jacobian_matrix.shape[0]):
            h_sum += self._outputs_dual_multipliers[i] * self.H[i]

        self.x = self._input_values
        self.del_f = self.jacobian_matrix
        self.hessian = h_sum
        return self.hessian

    def hessian_damped_bfgs(self):
        # apply Powell's damping on the BFGS update
        if not hasattr(self, "x"):
            self.x = None
        if not hasattr(self, "H"):
            self.H = []
            H_i = np.identity(len(self.inputs))
            for i in range(self.jacobian_matrix.shape[0]):
                self.H.append(H_i.copy())
        if self.x is not None:
            s_k = (np.array([self._input_values]) - self.x).T
            phi = 0.9
            for i in range(self.jacobian_matrix.shape[0]):
                y_k = np.array([self.jacobian_matrix[i, :] - self.del_f[i, :]]).T
                y_s = y_k.T @ s_k
                H_s = self.H[i] @ s_k

                # new
                s_H_s = s_k.T @ H_s
                if y_s >= phi * s_H_s:
                    delta_k = 1
                else:
                    delta_k = (1 - phi) * s_H_s / (s_H_s - y_s)
                z_k = delta_k * y_k + (1 - delta_k) * H_s
                z_s = z_k.T @ s_k
                if z_k.shape != y_k.shape:
                    raise RuntimeError()
                if s_k.any():  # extra
                    self.H[i] = (
                        self.H[i] + (z_k @ z_k.T) / (z_s) - (H_s @ H_s.T) / (s_H_s)
                    )
                ###########################################
        h_sum = np.zeros((len(self.inputs), len(self.inputs)))
        for i in range(self.jacobian_matrix.shape[0]):
            h_sum += self._outputs_dual_multipliers[i] * self.H[i]

        self.x = self._input_values
        self.del_f = self.jacobian_matrix
        self.hessian = h_sum
        # print(self.hessian)
        return self.hessian

    def hessian_ipopt_bfgs_modification(self):
        # BFGS update is only done on certain conditions (taken from IPOPT's implementation)
        if not hasattr(self, "x"):
            self.x = None
        if not hasattr(self, "H"):
            self.H = []
            H_i = np.identity(len(self.inputs))
            for i in range(self.jacobian_matrix.shape[0]):
                self.H.append(H_i.copy())
        if self.x is not None:
            s_k = (np.array([self._input_values]) - self.x).T
            for i in range(self.jacobian_matrix.shape[0]):
                y_k = np.array([self.jacobian_matrix[i, :] - self.del_f[i, :]]).T
                y_s = y_k.T @ s_k
                H_s = self.H[i] @ s_k
                mach_eps = np.finfo(float).eps
                if (
                    y_s.T
                    > np.sqrt(mach_eps) * np.linalg.norm(s_k) * np.linalg.norm(y_k)
                ) and (np.linalg.norm(s_k, np.inf) >= 100 * mach_eps):
                    self.H[i] = (
                        self.H[i]
                        + (y_k @ y_k.T) / (y_s)
                        - (H_s @ H_s.T) / (s_k.T @ H_s)
                    )
        h_sum = np.zeros((len(self.inputs), len(self.inputs)))
        for i in range(self.jacobian_matrix.shape[0]):
            h_sum += self._outputs_dual_multipliers[i] * self.H[i]

        self.x = self._input_values
        self.del_f = self.jacobian_matrix
        self.hessian = h_sum
        return self.hessian

    def evaluate_hessian_outputs(self):
        if self.hess_type == HessTypes.JtJ:
            self._hess = self.hessian_gauss_newton_version(sparse_jac=False)
        if self.hess_type == HessTypes.BFGS:
            self._hess = self.hessian_bfgs()
        if self.hess_type == HessTypes.BFGS_mod:
            self._hess = self.hessian_modified_bfgs()
        if self.hess_type == HessTypes.BFGS_damp:
            self._hess = self.hessian_damped_bfgs()
        if self.hess_type == HessTypes.BFGS_ipopt:
            self._hess = self.hessian_ipopt_bfgs_modification()
        jm = np.array(self._hess)
        cm = tril(jm)
        return cm
