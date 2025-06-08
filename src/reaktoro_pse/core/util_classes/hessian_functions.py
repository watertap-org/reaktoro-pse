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

from scipy.sparse import coo_matrix, tril

__author__ = "Ilayda Akkor, Alexander V. Dudchenko, Paul Vecchiarelli, Ben Knueven"


class HessTypes:
    GaussNewton = "GaussNewton"

    LBFGS = "LBFGS"
    BFGS = "BFGS"
    CBFGS = "CBFGS"
    BFGS_mod = "BFGS_mod"
    BFGS_damp = "BFGS_damp"
    BFGS_ipopt = "BFGS_ipopt"
    no_hessian_estimation = "no_hessian_estimation"
    ZeroHessian = "ZeroHessian"
    sparse_16 = "sparse_16"
    diag_inv = "diag_inv"


class HessianMemory:
    def __init__(self, memory=6):
        self.memory = memory
        self.inputs = []
        self.jacobian = []
        self.hessian = []
        self.memory_limited_arrays = [
            self.inputs,
            self.jacobian,
            self.hessian,
        ]

    def trim_memory(self, force_trim=False):
        if len(self.inputs) > self.memory or force_trim:
            for mla in self.memory_limited_arrays:
                mla.pop(0)

    def reset_memory(self):
        print("hess memory resdet")
        for mla in self.memory_limited_arrays:
            mla.clear()

    def memorize(self, inputs, jacobian, hessian):
        self.inputs.append(inputs.copy())
        self.jacobian.append(jacobian.copy())
        self.hessian.append(hessian.copy())
        # print(self.hessian)  #
        self.trim_memory()  #

    def get_last_input(self):
        if len(self.inputs) > 0:
            return self.inputs[-1]
        else:
            return None

    def get_last_jacobian(self):
        if len(self.jacobian) > 0:
            return self.jacobian[-1]
        else:
            return None

    def get_last_hessian(self):
        if len(self.hessian) > 0:
            return self.hessian[-1]
        else:
            return None

    def get_range(self, start=0, end_offset=0):
        return range(start, len(self.inputs) + end_offset)

    def set_initial_matrix(self, hessian):
        self.initial_hessian = hessian.copy()


class HessianApproximation:
    def __init__(self, hessian_type=None):
        self.hessian_memory = HessianMemory(memory=4)
        if hessian_type is None:
            self.hessian_matrix_type = HessTypes.ZeroHessian
        else:
            self.hessian_matrix_type = hessian_type
        self.iters = 0
        self.search_started = False
        self.bfgs_hessian = None
        self.old_inputs = None
        self.s = None
        self.bfgs_matrix_not_initialized = True
        self.reset_counts = 0

    def hessian_gauss_newton_version(self, sparse_jac, threshold=1e-16):
        """standard gauss newton hessian approximation"""
        hess = np.zeros((len(self.inputs), len(self.inputs)))
        for i in range(self.jacobian_matrix.shape[0]):
            row = self.jacobian_matrix[i, :]
            if sparse_jac:
                row[np.abs(row) < threshold] = 0
            hess += self._outputs_dual_multipliers[i] * np.outer(row.T, row)

        self.hessian_matrix = hess

    def get_initial_hessian(self, old_step, new_step, old_jacobian, new_jacobian):
        bfgs_hessian = []
        for i in range(self.jacobian_matrix.shape[0]):
            bfgs_hessian.append(np.identity(len(self.inputs)) * 0)
        bfgs_hessian = np.array(bfgs_hessian)

        s = np.array([new_step]) - old_step
        for i in range(self.jacobian_matrix.shape[0]):
            y = np.array([new_jacobian[i, :] - old_jacobian[i, :]])
            sTs = s.T @ s
            yTy = y.T @ y
            sTy = s.T @ y
            if np.sum(y) != 0 and np.sum(s) != 0:
                sigma = np.zeros(bfgs_hessian[i].shape)
                sigma[sTs != 0] = sTy[sTs != 0] / sTs[sTs != 0]
                sigma[sigma == np.inf] = 0
                sigma[sigma == -np.inf] = 0
                bfgs_hessian[i] = np.identity(len(old_step)) * sigma
        return bfgs_hessian

    def create_bfgs_matrix(self):
        if self.bfgs_hessian is None:
            self.bfgs_hessian = []
            for i in range(self.jacobian_matrix.shape[0]):
                self.bfgs_hessian.append(np.identity(len(self.inputs)) * 0)
            self.bfgs_hessian = np.array(self.bfgs_hessian)

        if (
            self.check_step()
            and self.bfgs_matrix_not_initialized
            or self.reset_counts == 2
        ):
            self.bfgs_matrix_not_initialized = False
            self.bfgs_hessian = self.get_initial_hessian(
                self.hessian_memory.get_last_input(),
                self.inputs,
                self.hessian_memory.get_last_jacobian(),
                self.jacobian_matrix,
            )

            self.hessian_memory.set_initial_matrix(self.bfgs_hessian)
            self.hessian_memory.reset_memory()
            self.reset_counts = 0
        if np.sum(self._outputs_dual_multipliers) == 0:
            self.reset_counts += 1
        else:
            self.reset_counts = 0
        self.iters += 1

    def check_step(self):
        if self.hessian_memory.get_last_input() is None:
            return False
        self.s = np.array([self.inputs]) - self.hessian_memory.get_last_input()

        if np.sum(self.s) != 0:
            return True
        else:
            return False

    def update_bfgs_matrix(self):
        self.hessian_memory.memorize(
            self.inputs, self.jacobian_matrix, self.bfgs_hessian
        )
        h_sum = np.zeros((len(self.inputs), len(self.inputs)))
        for i in range(self.jacobian_matrix.shape[0]):
            h_sum += self._outputs_dual_multipliers[i] * self.bfgs_hessian[i]
        self.hessian_matrix = h_sum.copy()

    def hessian_lbfgs(self):
        """Direct LBFGS based on Byrd Representations of quasi-newton matrices and their use in limited memory methods"""
        self.create_bfgs_matrix()
        if len(self.hessian_memory.get_range()) > 1:

            self.bfgs_hessian = self.get_initial_hessian(
                self.hessian_memory.inputs[-2],
                self.hessian_memory.inputs[-1],
                self.hessian_memory.jacobian[-2],
                self.hessian_memory.jacobian[-1],
            )

            initial_hessians = [None]
            for r in self.hessian_memory.get_range(1, -1):
                initial_hessians.append(
                    self.get_initial_hessian(
                        self.hessian_memory.inputs[r - 1],
                        self.hessian_memory.inputs[r],
                        self.hessian_memory.jacobian[r - 1],
                        self.hessian_memory.jacobian[r],
                    ).copy()
                )
            for i in range(self.jacobian_matrix.shape[0]):
                bk = []
                ak = []
                for r in self.hessian_memory.get_range(1, -1):
                    sk = (
                        self.hessian_memory.inputs[r]
                        - self.hessian_memory.inputs[r - 1]
                    )
                    yk = (
                        self.hessian_memory.jacobian[r][i, :]
                        - self.hessian_memory.jacobian[r - 1][i, :]
                    )
                    y_s = yk.T @ sk
                    mach_eps = np.finfo(float).eps
                    if (
                        yk.T @ sk > 1e-6
                        and np.sum(yk) != 0
                        and np.sum(sk) != 0
                        and np.sum(initial_hessians[r]) != 0
                    ):
                        b = yk / (yk.T @ sk) ** 0.5
                        b[b != b] = 0
                        b[b == np.inf] = 0
                        b[b == -np.inf] = 0
                        bk.append(b.copy())
                        if len(ak) > 0:
                            s = np.zeros(ak[0].shape)
                            for k in range(1, r):
                                _bk = bk[k]
                                _ak = ak[k]
                                s += (_bk.T @ sk) @ _bk - (_ak.T @ sk) @ _ak
                            _ak = initial_hessians[r][i] @ sk + s
                        else:
                            _ak = initial_hessians[r][i] @ sk

                        ska = sk / (sk.T @ _ak) ** 0.5
                        ska[ska != ska] = 0
                        _ak[_ak == np.inf] = 0
                        _ak[_ak == -np.inf] = 0
                        ak.append(_ak.copy())
                for m in range(len(bk)):
                    sum_ak_bk = bk[m] @ bk[m].T - ak[m] @ ak[m].T

                    self.bfgs_hessian[i] = self.bfgs_hessian[i] + sum_ak_bk
                    if (
                        np.isnan(self.bfgs_hessian).any()
                        or np.isinf(self.bfgs_hessian).any()
                    ):
                        print(self.bfgs_hessian[i], sum_ak_bk)

            # print("bfgs_hessian", self.bfgs_hessian)
        self.update_bfgs_matrix()

    def hessian_bfgs(self):
        """Vanilla BFGS update implementation"""
        self.create_bfgs_matrix()
        if len(self.hessian_memory.get_range()) > 1:
            s_k = (np.array([self.inputs]) - self.hessian_memory.get_last_input()).T
            for i in range(self.jacobian_matrix.shape[0]):
                y_k = np.array(
                    [
                        self.jacobian_matrix[i, :]
                        - self.hessian_memory.get_last_jacobian()[i, :]
                    ]
                ).T
                y_s = y_k.T @ s_k
                H_s = self.bfgs_hessian[i] @ s_k

                if y_k.T @ s_k > 1e-32 and np.sum(H_s) != 0:
                    self.bfgs_hessian[i] = (
                        self.bfgs_hessian[i]
                        + (y_k @ y_k.T) / (y_s)
                        - (H_s @ H_s.T) / (s_k.T @ H_s)
                    )
        self.update_bfgs_matrix()

    def hessian_cbfgs(self):
        """Cautious BFGS update implementation (Li and Fukushima)"""
        self.create_bfgs_matrix()
        if len(self.hessian_memory.get_range()) > 1:
            s_k = (np.array([self.inputs]) - self.hessian_memory.get_last_input()).T
            eps = 100 * np.finfo(float).eps
            for i in range(self.jacobian_matrix.shape[0]):
                y_k = np.array(
                    [
                        self.jacobian_matrix[i, :]
                        - self.hessian_memory.get_last_jacobian()[i, :]
                    ]
                ).T
                y_s = y_k.T @ s_k
                H_s = self.bfgs_hessian[i] @ s_k
                if np.linalg.norm(self.jacobian_matrix[i, :]) >= 1:
                    alpha = 0.01
                if np.linalg.norm(self.jacobian_matrix[i, :]) < 1:
                    alpha = 3
                min_test = eps * np.linalg.norm(self.jacobian_matrix[i, :]) ** alpha
                max_test = y_s / np.linalg.norm(s_k) ** 2
                if max_test > min_test and np.sum(H_s) != 0:
                    self.bfgs_hessian[i] = (
                        self.bfgs_hessian[i]
                        + (y_k @ y_k.T) / (y_s)
                        - (H_s @ H_s.T) / (s_k.T @ H_s)
                    )
        self.update_bfgs_matrix()

    def hessian_modified_bfgs(self):
        """Modified BFGS update implementation (Li and Fukushima)"""
        self.create_bfgs_matrix()
        if len(self.hessian_memory.get_range()) > 1:
            s_k = (np.array([self.inputs]) - self.hessian_memory.get_last_input()).T
            for i in range(self.jacobian_matrix.shape[0]):
                y_k = np.array(
                    [
                        self.jacobian_matrix[i, :]
                        - self.hessian_memory.get_last_jacobian()[i, :]
                    ]
                ).T
                y_s = y_k.T @ s_k
                H_s = self.bfgs_hessian[i] @ s_k
                if s_k.any() and np.sum(H_s) != 0:
                    t_k = 1 + max(0, -y_s / (np.linalg.norm(s_k) ** 2))
                    z_k = (
                        y_k
                        + t_k
                        * np.linalg.norm(self.hessian_memory.get_last_jacobian()[i, :])
                        * s_k
                    )
                    self.bfgs_hessian[i] = (
                        self.bfgs_hessian[i]
                        + (z_k @ z_k.T) / (z_k.T @ s_k)
                        - (H_s @ H_s.T) / (s_k.T @ H_s)
                    )
        self.update_bfgs_matrix()

    def hessian_damped_bfgs(self):
        """apply Powell's damping on the BFGS update"""
        self.create_bfgs_matrix()
        if len(self.hessian_memory.get_range()) > 1:
            s_k = (np.array([self.inputs]) - self.hessian_memory.get_last_input()).T
            phi = 0.75
            for i in range(self.jacobian_matrix.shape[0]):
                y_k = np.array(
                    [
                        self.jacobian_matrix[i, :]
                        - self.hessian_memory.get_last_jacobian()[i, :]
                    ]
                ).T
                y_s = y_k.T @ s_k
                H_s = self.bfgs_hessian[i] @ s_k

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
                if s_k.any() and np.sum(H_s) != 0 and np.sum(z_s) != 0:  # extra
                    self.bfgs_hessian[i] = (
                        self.bfgs_hessian[i]
                        + (z_k @ z_k.T) / (z_s)
                        - (H_s @ H_s.T) / (s_H_s)
                    )
                ###########################################
        self.update_bfgs_matrix()

    def hessian_ipopt_bfgs_modification(self):
        """BFGS update is only done on certain conditions (taken from IPOPT's implementation)"""
        self.create_bfgs_matrix()
        if len(self.hessian_memory.get_range()) > 1:
            s_k = (np.array([self.inputs]) - self.hessian_memory.get_last_input()).T
            for i in range(self.jacobian_matrix.shape[0]):
                y_k = np.array(
                    [
                        self.jacobian_matrix[i, :]
                        - self.hessian_memory.get_last_jacobian()[i, :]
                    ]
                ).T
                y_s = y_k.T @ s_k
                H_s = self.bfgs_hessian[i] @ s_k
                mach_eps = np.finfo(float).eps
                if (
                    (
                        y_s.T
                        > np.sqrt(mach_eps) * np.linalg.norm(s_k) * np.linalg.norm(y_k)
                    )
                    and (np.linalg.norm(s_k, np.inf) >= 1 * mach_eps)
                    and np.sum(H_s) != 0
                ):
                    self.bfgs_hessian[i] = (
                        self.bfgs_hessian[i]
                        + (y_k @ y_k.T) / (y_s)
                        - (H_s @ H_s.T) / (s_k.T @ H_s)
                    )
        self.update_bfgs_matrix()

    def hessian_diag_inv_value(self):
        hessian = np.zeros((len(self.inputs), len(self.inputs)))
        for idx, v in enumerate(self.inputs):
            hessian[idx, idx] = 1.0 / v
        h_sum = np.zeros((len(self.inputs), len(self.inputs)))
        for i in range(self.jacobian_matrix.shape[0]):
            h_sum += self._outputs_dual_multipliers[i] * hessian[i]
        self.hessian_matrix = h_sum.copy()

    def sparse_diagonal(self, shape, value=1e-16):
        rows = []
        cols = []
        vals = []
        for i in range(shape):
            rows.append(i)
            cols.append(i)
            vals.append(value)

        self.hessian_matrix = coo_matrix((vals, (rows, cols)), shape=(shape, shape))

    def get_hessian(self, input_values, output_values, jacobian, dual_multipliers):
        self.inputs = np.array(input_values)
        self.outputs = np.array(output_values)
        self.jacobian_matrix = np.array(jacobian)
        self._outputs_dual_multipliers = dual_multipliers
        if self.hessian_matrix_type == HessTypes.ZeroHessian:
            self.sparse_diagonal(len(self.inputs), 0)
        elif self.hessian_matrix_type == HessTypes.sparse_16:
            self.sparse_diagonal(len(self.inputs), 1e-16)
        elif self.hessian_matrix_type == HessTypes.GaussNewton:
            self.hessian_gauss_newton_version(sparse_jac=False)
        elif self.hessian_matrix_type == HessTypes.LBFGS:
            self.hessian_lbfgs()
        elif self.hessian_matrix_type == HessTypes.BFGS:
            self.hessian_bfgs()
        elif self.hessian_matrix_type == HessTypes.CBFGS:
            self.hessian_cbfgs()
        elif self.hessian_matrix_type == HessTypes.BFGS_mod:
            self.hessian_modified_bfgs()
        elif self.hessian_matrix_type == HessTypes.BFGS_damp:
            self.hessian_damped_bfgs()
        elif self.hessian_matrix_type == HessTypes.BFGS_ipopt:
            self.hessian_ipopt_bfgs_modification()
        elif self.hessian_matrix_type == HessTypes.diag_inv:
            self.hessian_diag_inv_value()
        else:
            raise NotImplementedError(
                f"Hessian type {self.hessian_matrix_type} not implemented"
            )
        if isinstance(self.hessian_matrix, coo_matrix):
            return self.hessian_matrix
        else:
            low_triangular_hessian = _hand_tril(np.array(self.hessian_matrix))
            return low_triangular_hessian


def _hand_tril(jm):

    assert jm.shape[0] == jm.shape[1]
    shape = jm.shape[0]
    row = []
    col = []
    val = []

    for i in range(shape):
        for j in range(i + 1):
            row.append(i)
            col.append(j)
            v = jm[i, j]

            val.append(v)

    return coo_matrix((val, (row, col)), shape=(shape, shape))
