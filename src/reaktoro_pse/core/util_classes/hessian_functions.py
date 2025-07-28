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

from scipy.sparse import coo_matrix

__author__ = "Ilayda Akkor, Alexander V. Dudchenko, Paul Vecchiarelli, Ben Knueven"


class HessTypes:
    GaussNewton = "GaussNewton"
    BFGS = "BFGS"
    CBFGS = "CBFGS"
    BFGS_mod = "BFGS_mod"
    BFGS_damp = "BFGS_damp"
    BFGS_ipopt = "BFGS_ipopt"
    LBFGS = "LBFGS"
    no_hessian_estimation = "no_hessian_estimation"
    ZeroHessian = "ZeroHessian"
    sparse_16 = "sparse_16"


class BFGSInitializationTypes:
    scalar1 = "scalar1"
    scalar2 = "scalar2"
    scalar3 = "scalar3"
    scalar4 = "scalar4"
    GaussNewton = "GaussNewton"
    constant = "constant"


class HessianMemory:
    def __init__(self, memory=6):
        self.memory = memory
        self.inputs = []
        self.jacobian = []
        self.memory_limited_arrays = [
            self.inputs,
            self.jacobian,
        ]

    def trim_memory(self, force_trim=False):
        if len(self.inputs) > self.memory or force_trim:
            for mla in self.memory_limited_arrays:
                mla.pop(0)

    def reset_memory(self, leave_last=False):
        if leave_last:
            mem_to_reset = self.memory_limited_arrays[:-2]
        else:
            mem_to_reset = self.memory_limited_arrays
        for mla in mem_to_reset:
            mla.clear()

    def memorize(self, inputs, jacobian):
        self.inputs.append(inputs.copy())
        self.jacobian.append(jacobian.copy())
        self.trim_memory()

    def get_range(self, start=0, end_offset=0):
        return range(start, len(self.inputs) + end_offset)


class HessianApproximation:
    """general classs for hessian approximation methods"""

    def __init__(
        self,
        hessian_type=None,
        bfgs_init_min_hessian_value=1e-32,
        bfgs_init_max_hessian_value=1e8,
        bfgs_init_const_hessian_value=1e-32,
        bfgs_initialization_type=BFGSInitializationTypes.GaussNewton,
        bfgs_hessian_memory=3,
        bfgs_epsilon=1e-12,  # same as ipopt!
    ):
        """initialize hessian approximation class
        Args:
            hessian_type: type of hessian approximation to use
            bfgs_init_min_hessian_value: minimum value for hessian diagonal elements
            bfgs_init_max_hessian_value: maximum value for hessian diagonal elements
            bfgs_init_const_hessian_value: constant value for hessian diagonal elements
            bfgs_initialization_type: type of initialization for BFGS hessian
            bfgs_hessian_memory: memory size for BFGS hessian approximation
            bfgs_epsilon: epsilon value for numerical stability in BFGS updates
        """
        self.hessian_memory = HessianMemory(memory=bfgs_hessian_memory)
        if hessian_type is None:
            self.hessian_matrix_type = HessTypes.ZeroHessian
        else:
            self.hessian_matrix_type = hessian_type
        self.iters = 0
        self.search_started = False
        self.bfgs_hessian = None
        self.old_inputs = None
        self.s = None
        self.epsilon = bfgs_epsilon  # same as ipopt!   # np.finfo(float).eps
        self.init_min_hessian_value = bfgs_init_min_hessian_value
        self.init_max_hessian_value = bfgs_init_max_hessian_value
        self.init_const_hessian_value = bfgs_init_const_hessian_value
        self.bfgs_matrix_not_initialized = True
        self.bfgs_initialization_type = bfgs_initialization_type

    def apply_sigma_bounds(self, sigma, abs_test=True):
        """methods for bounding the sigma values in BFGS initialization"""
        sigma_signs = np.sign(sigma)
        zeros = sigma == 0
        if self.init_min_hessian_value > 0:
            if abs_test:
                test = np.abs(sigma) < self.init_min_hessian_value
            else:
                test = sigma < self.init_min_hessian_value
            sigma[test] = self.init_min_hessian_value
            if abs_test:
                sigma[test] *= sigma_signs[test]

        if self.init_max_hessian_value > 0:
            if abs_test:
                test = np.abs(sigma) > self.init_max_hessian_value
            else:
                test = sigma > self.init_max_hessian_value
            sigma[test] = self.init_max_hessian_value
            if abs_test:
                sigma[test] *= sigma_signs[test]
        sigma[zeros] = 0
        return sigma

    def get_initial_hessian(self, old_step, new_step, old_jacobian, new_jacobian):
        """grab initial hessian matrix based on the old and new steps and jacobian matrices for BFGS initialization"""
        bfgs_hessian = []
        for i in range(new_jacobian.shape[0]):
            bfgs_hessian.append(np.identity(len(new_step)))
        bfgs_hessian = np.array(bfgs_hessian)

        s = (np.array([new_step]) - old_step).T
        for i in range(new_jacobian.shape[0]):
            y = np.array([new_jacobian[i, :] - old_jacobian[i, :]]).T
            if (
                np.sum(y) != 0
                and np.sum(s) != 0
                and self.bfgs_initialization_type != BFGSInitializationTypes.constant
            ):
                if self.bfgs_initialization_type == BFGSInitializationTypes.scalar1:
                    sTy = s.T @ y
                    sTs = s.T @ s
                    if sTs == 0:
                        sigma = self.init_const_hessian_value
                    else:
                        sigma = sTy / sTs
                elif self.bfgs_initialization_type == BFGSInitializationTypes.scalar2:
                    sTy = s.T @ y
                    yTy = y.T @ y
                    if sTy == 0:
                        sigma = self.init_const_hessian_value
                    else:
                        sigma = yTy / sTy
                elif self.bfgs_initialization_type == BFGSInitializationTypes.scalar3:
                    sTy = s.T @ y
                    yTy = y.T @ y
                    sTs = s.T @ s
                    if sTs == 0:
                        sigma = self.init_const_hessian_value
                    else:
                        sigma = sTy / sTs / 2
                        if sTy != 0:
                            sigma += (yTy / sTy) / 2

                elif self.bfgs_initialization_type == BFGSInitializationTypes.scalar4:
                    yTy = y.T @ y
                    sTs = s.T @ s
                    if sTs == 0:
                        sigma = self.init_const_hessian_value
                    else:
                        sigma = sTy / sTs
                    if sTy != 0:
                        sigma *= yTy / sTy
                    sigma_sign = np.sign(sigma)
                    sigma = np.sqrt(np.abs(sigma))
                    sigma *= sigma_sign
                elif (
                    self.bfgs_initialization_type == BFGSInitializationTypes.GaussNewton
                ):
                    sigma = np.outer(new_jacobian[i, :].T, new_jacobian[i, :])

            else:
                sigma = self.init_const_hessian_value
            bfgs_hessian[i] *= sigma
        bfgs_hessian = self.apply_sigma_bounds(bfgs_hessian)
        return bfgs_hessian

    def create_bfgs_matrix(self):
        """create and amanage BFGS matrix
        - initialize matrix if detected that all duals are zero or no matrix exists"""
        self.hessian_memory.memorize(self.inputs, self.jacobian_matrix)

        if self.bfgs_hessian is None:
            self.bfgs_hessian = []
            for i in range(self.hessian_memory.jacobian[-1].shape[0]):
                self.bfgs_hessian.append(
                    np.identity(len(self.inputs)) * self.init_const_hessian_value
                )
            self.bfgs_hessian = np.array(self.bfgs_hessian)

        if self.check_step() and self.bfgs_matrix_not_initialized:
            self.bfgs_matrix_not_initialized = False
            self.bfgs_hessian = self.get_initial_hessian(
                self.hessian_memory.inputs[-2],
                self.hessian_memory.inputs[-1],
                self.hessian_memory.jacobian[-2],
                self.hessian_memory.jacobian[-1],
            )
            self.hessian_memory.reset_memory(leave_last=True)
            self.reset_counts = 0
        if np.sum(self._outputs_dual_multipliers) == 0:
            self.bfgs_matrix_not_initialized = True
        self.iters += 1

    def check_step(self):
        """check if the step is valid for BFGS update or if hessian memory is too small
        - if the step is zero, return False
        """
        if len(self.hessian_memory.inputs) < 2:
            return False
        self.s = self.hessian_memory.inputs[-2] - self.hessian_memory.inputs[-1]
        if np.sum(self.s) != 0:
            return True
        else:
            return False

    def update_bfgs_matrix(self):
        """update matrix with dual multipliers"""
        h_sum = np.zeros((len(self.inputs), len(self.inputs)))
        for i in range(self.hessian_memory.jacobian[-1].shape[0]):
            h_sum += self._outputs_dual_multipliers[i] * self.bfgs_hessian[i]
        self.hessian_matrix = h_sum.copy()

    def hessian_gauss_newton_version(self, sparse_jac, threshold=1e-8):
        """standard gauss newton hessian approximation"""
        hess = np.zeros((len(self.inputs), len(self.inputs)))

        for i in range(self.jacobian_matrix.shape[0]):
            row = self.jacobian_matrix[i, :]
            if sparse_jac:
                row[np.abs(row) < threshold] = 0
            hess_row = np.outer(row.T, row)
            hess += self._outputs_dual_multipliers[i] * hess_row

        self.hessian_matrix = hess

    def hessian_bfgs(self):
        """Vanilla BFGS update implementation"""
        self.create_bfgs_matrix()
        if len(self.hessian_memory.get_range()) > 1:
            s_k = (
                np.array([self.hessian_memory.inputs[-1]])
                - self.hessian_memory.inputs[-2]
            ).T
            for i in range(self.hessian_memory.jacobian[-1].shape[0]):
                y_k = np.array(
                    [
                        self.hessian_memory.jacobian[-1][i, :]
                        - self.hessian_memory.jacobian[-2][i, :]
                    ]
                ).T
                B_k = self.bfgs_hessian[i]
                H_s = self.bfgs_hessian[i] @ s_k
                if y_k.T @ s_k > self.epsilon and np.sum(H_s) != 0:
                    update_pos = (y_k @ y_k.T) / (y_k.T @ s_k)
                    update_neg = (B_k @ s_k @ s_k.T @ B_k) / (s_k.T @ B_k @ s_k)

                    self.bfgs_hessian[i] = (
                        self.bfgs_hessian[i] + update_pos - update_neg
                    )
        self.update_bfgs_matrix()

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
            for r in self.hessian_memory.get_range(1, 0):
                initial_hessians.append(
                    self.get_initial_hessian(
                        self.hessian_memory.inputs[r - 1],
                        self.hessian_memory.inputs[r],
                        self.hessian_memory.jacobian[r - 1],
                        self.hessian_memory.jacobian[r],
                    ).copy()
                )

            def ipopt_update_test(sn, so, jn, jo, hessian):
                sk = (np.array([sn]) - so).T
                yk = np.array([jn - jo]).T
                y_s = yk.T @ sk
                if (
                    (
                        y_s
                        > np.sqrt(self.epsilon)
                        * np.linalg.norm(sk)
                        * np.linalg.norm(yk)
                    )
                    and (np.linalg.norm(sk, np.inf) >= self.epsilon)
                    and np.sum(yk) != 0
                    and np.sum(sk) != 0
                    and np.sum(hessian) != 0
                ):
                    return True
                else:
                    return False

            for i in range(self.hessian_memory.jacobian[-1].shape[0]):
                # only update if current step is good
                if ipopt_update_test(
                    self.hessian_memory.inputs[-1],
                    self.hessian_memory.inputs[-2],
                    self.hessian_memory.jacobian[-1][i, :],
                    self.hessian_memory.jacobian[-2][i, :],
                    initial_hessians[-1][i],
                ):
                    bk = []
                    ak = []
                    for r in self.hessian_memory.get_range(1, 0):
                        sk = (
                            np.array([self.hessian_memory.inputs[r]])
                            - self.hessian_memory.inputs[r - 1]
                        ).T
                        yk = np.array(
                            [
                                self.hessian_memory.jacobian[r][i, :]
                                - self.hessian_memory.jacobian[r - 1][i, :]
                            ]
                        ).T
                        # only include update if sub step is good
                        if ipopt_update_test(
                            self.hessian_memory.inputs[r],
                            self.hessian_memory.inputs[r - 1],
                            self.hessian_memory.jacobian[r][i, :],
                            self.hessian_memory.jacobian[r - 1][i, :],
                            initial_hessians[r][i],
                        ):
                            b = yk / np.sqrt(yk.T @ sk)
                            _ak = initial_hessians[r][i] @ sk
                            for k in range(len(ak)):
                                _ak += (bk[k].T @ sk) * bk[k] - (ak[k].T @ sk) * ak[k]
                            if sk.T @ _ak > 0:
                                _ak = _ak / np.sqrt(sk.T @ _ak)
                                ak.append(_ak.copy())
                                bk.append(b.copy())
                    sum_ak_bk = np.zeros(self.bfgs_hessian[i].shape)

                    for m in range(len(bk)):
                        sum_ak_bk += bk[m] @ bk[m].T - ak[m] @ ak[m].T
                    self.bfgs_hessian[i] = self.bfgs_hessian[i] + sum_ak_bk
        self.update_bfgs_matrix()

    def hessian_cbfgs(self):
        """Cautious BFGS update implementation (Li and Fukushima)"""
        self.create_bfgs_matrix()
        if len(self.hessian_memory.get_range()) > 1:
            s_k = (
                np.array([self.hessian_memory.inputs[-1]])
                - self.hessian_memory.inputs[-2]
            ).T
            for i in range(self.hessian_memory.jacobian[-1].shape[0]):
                y_k = np.array(
                    [
                        self.hessian_memory.jacobian[-1][i, :]
                        - self.hessian_memory.jacobian[-2][i, :]
                    ]
                ).T
                y_s = y_k.T @ s_k
                H_s = self.bfgs_hessian[i] @ s_k
                if np.linalg.norm(self.hessian_memory.jacobian[-1][i, :]) >= 1:
                    alpha = 0.01
                if np.linalg.norm(self.hessian_memory.jacobian[-1][i, :]) < 1:
                    alpha = 3
                min_test = (
                    self.epsilon
                    * np.linalg.norm(self.hessian_memory.jacobian[-1][i, :]) ** alpha
                )
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
            s_k = (
                np.array([self.hessian_memory.inputs[-1]])
                - self.hessian_memory.inputs[-2]
            ).T
            for i in range(self.hessian_memory.jacobian[-1].shape[0]):
                y_k = np.array(
                    [
                        self.hessian_memory.jacobian[-1][i, :]
                        - self.hessian_memory.jacobian[-2][i, :]
                    ]
                ).T
                y_s = y_k.T @ s_k
                H_s = self.bfgs_hessian[i] @ s_k
                if s_k.any() and np.sum(H_s) != 0:
                    t_k = 1 + max(0, -y_s / (np.linalg.norm(s_k) ** 2))
                    z_k = (
                        y_k
                        + t_k
                        * np.linalg.norm(self.hessian_memory.jacobian[-1][i, :])
                        * s_k
                    )
                    if z_k.T @ s_k != 0 and s_k.T @ H_s != 0:
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
            s_k = (
                np.array([self.hessian_memory.inputs[-1]])
                - self.hessian_memory.inputs[-2]
            ).T
            phi = 0.9
            for i in range(self.hessian_memory.jacobian[-1].shape[0]):
                y_k = np.array(
                    [
                        self.hessian_memory.jacobian[-1][i, :]
                        - self.hessian_memory.jacobian[-2][i, :]
                    ]
                ).T
                y_s = y_k.T @ s_k

                H_s = self.bfgs_hessian[i] @ s_k
                if np.sum(H_s) != 0 and np.sum(s_k) != 0:
                    s_H_s = s_k.T @ H_s
                if np.sum(H_s) != 0 and np.sum(s_k) != 0 and s_H_s - y_s > 0:
                    if y_s >= phi * s_H_s:
                        delta_k = 1
                    else:
                        delta_k = (1 - phi) * s_H_s / (s_H_s - y_s)
                    z_k = delta_k * y_k + (1 - delta_k) * H_s
                    z_s = z_k.T @ s_k
                    if z_k.shape != y_k.shape:
                        raise RuntimeError()
                    if s_k.any() and np.sum(z_s) != 0:
                        self.bfgs_hessian[i] = (
                            self.bfgs_hessian[i]
                            + (z_k @ z_k.T) / (z_s)
                            - (H_s @ H_s.T) / (s_H_s)
                        )
        self.update_bfgs_matrix()

    def hessian_ipopt_bfgs_modification(self):
        """BFGS update is only done on certain conditions (taken from IPOPT's implementation)"""
        self.create_bfgs_matrix()
        if len(self.hessian_memory.get_range()) > 1:
            s_k = (
                np.array([self.hessian_memory.inputs[-1]])
                - self.hessian_memory.inputs[-2]
            ).T
            for i in range(self.hessian_memory.jacobian[-1].shape[0]):
                y_k = np.array(
                    [
                        self.hessian_memory.jacobian[-1][i, :]
                        - self.hessian_memory.jacobian[-2][i, :]
                    ]
                ).T
                y_s = y_k.T @ s_k
                H_s = self.bfgs_hessian[i] @ s_k
                mach_eps = self.epsilon
                if (
                    (
                        y_s.T
                        > np.sqrt(mach_eps) * np.linalg.norm(s_k) * np.linalg.norm(y_k)
                    )
                    and (np.linalg.norm(s_k, np.inf) >= self.epsilon)
                    and np.sum(H_s) != 0
                ):
                    self.bfgs_hessian[i] = (
                        self.bfgs_hessian[i]
                        + (y_k @ y_k.T) / (y_s)
                        - (H_s @ H_s.T) / (s_k.T @ H_s)
                    )
        self.update_bfgs_matrix()

    def sparse_diagonal(self, shape, value=1e-16):
        rows = []
        cols = []
        vals = []
        for i in range(shape):
            rows.append(i)
            cols.append(i)
            vals.append(value)

        self.hessian_matrix = coo_matrix((vals, (rows, cols)), shape=(shape, shape))

    def get_hessian(self, input_values, jacobian, dual_multipliers):
        self.inputs = np.array(input_values)
        self.jacobian_matrix = np.array(jacobian)
        self._outputs_dual_multipliers = dual_multipliers
        try:
            if self.hessian_matrix_type == HessTypes.ZeroHessian:
                self.sparse_diagonal(len(self.inputs), 0)
            elif self.hessian_matrix_type == HessTypes.sparse_16:
                self.sparse_diagonal(len(self.inputs), 1e-16)
            elif self.hessian_matrix_type == HessTypes.GaussNewton:
                self.hessian_gauss_newton_version(sparse_jac=False)
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
            elif self.hessian_matrix_type == HessTypes.LBFGS:
                self.hessian_lbfgs()
            else:
                raise NotImplementedError(
                    f"Hessian type {self.hessian_matrix_type} not implemented"
                )

            if isinstance(self.hessian_matrix, coo_matrix):
                return self.hessian_matrix
            else:

                low_triangular_hessian = _hand_tril(np.array(self.hessian_matrix))
                return low_triangular_hessian
        except Exception as e:
            print(
                f"Error in Hessian approximation: {e}. "
                f"Hessian type: {self.hessian_matrix_type}, "
                f"Inputs: {self.inputs}, "
                f"Jacobian: {self.jacobian_matrix}, "
                f"Dual multipliers: {self._outputs_dual_multipliers}"
            )


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
