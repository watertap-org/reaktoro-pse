#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/reaktoro-pse/"
#################################################################################
"""Hessian approximations for the Reaktoro gray box.

The gray box has outputs y_i(x) and PyNumero asks for sum_i lambda_i * d2y_i/dx2.
Reaktoro only provides first derivatives, so the curvature is approximated from
changes in the Jacobian between iterates.

Quasi-Newton types keep one matrix per *element*:
 - per-output types (BFGS family, SR1, PSB, LBFGS): one matrix per output,
   assembled as sum_i lambda_i * B_i (partitioned quasi-Newton, Griewank & Toint)
 - SR1_lagrangian: one matrix per group of outputs that share an input support,
   updated with the change in the Lagrangian gradient, assembled as sum_g B_g

Each element only covers the inputs its outputs depend on (``hessian_support``),
which keeps blocks in an aggregated gray box decoupled. Updates can optionally be
done in scaled inputs x_s = D x (``input_scaling``) which makes the initial matrix
and the skip tests consistent across inputs of very different magnitudes.

The returned matrix always has the same lower triangular sparsity pattern,
as required by PyNumero.
"""

import functools
import numpy as np

from scipy.sparse import coo_matrix
from pyomo.common.deprecation import deprecation_warning
import idaes.logger as idaeslog

__author__ = "Ilayda Akkor, Alexander V. Dudchenko, Paul Vecchiarelli, Ben Knueven"
_log = idaeslog.getLogger(__name__)


def _nonzero(vec):
    """True if any element of vec is non-zero (unlike np.sum(vec) != 0, which
    is False for vectors such as [1, -1])"""
    return bool(np.any(np.asarray(vec) != 0))


class HessTypes:
    GaussNewton = "GaussNewton"
    BFGS = "BFGS"
    CBFGS = "CBFGS"
    BFGS_mod = "BFGS_mod"
    BFGS_damp = "BFGS_damp"
    BFGS_ipopt = "BFGS_ipopt"
    LBFGS = "LBFGS"
    SR1 = "SR1"
    PSB = "PSB"
    SR1_lagrangian = "SR1_lagrangian"
    no_hessian_estimation = "no_hessian_estimation"
    ZeroHessian = "ZeroHessian"
    sparse_16 = "sparse_16"


class QuasiNewtonInitializationTypes:
    scalar1 = "scalar1"
    scalar2 = "scalar2"
    scalar3 = "scalar3"
    scalar4 = "scalar4"
    GaussNewton = "GaussNewton"
    constant = "constant"
    zero = "zero"


# single source of default values, used by HessianOptions, ReaktoroSolver,
# ReaktoroGrayBox and HessianApproximation
HESSIAN_DEFAULTS = {
    "hessian_type": HessTypes.LBFGS,
    "quasi_newton_initialization_type": QuasiNewtonInitializationTypes.GaussNewton,
    "quasi_newton_init_min_value": 1e-32,
    "quasi_newton_init_max_value": 1e8,
    "quasi_newton_init_const_value": 1e-64,
    "lbfgs_memory": 6,
    "quasi_newton_epsilon": 1e-16,
    "sr1_skip_tolerance": 1e-8,
    # None: SR1_SAFEGUARD_CAP for SR1 types, off for others (benchmark on softening RO
    # flowsheets: removes bad SR1 updates, makes BFGS/LBFGS slower)
    "quasi_newton_update_cap": None,
    "quasi_newton_norm_cap": None,
    "hessian_input_scaling": False,
}

# options were named bfgs_* before non BFGS quasi-Newton types (SR1, PSB) were added,
# the old names are still accepted (with a deprecation warning)
DEPRECATED_HESSIAN_OPTIONS = {
    "bfgs_initialization_type": "quasi_newton_initialization_type",
    "bfgs_init_min_hessian_value": "quasi_newton_init_min_value",
    "bfgs_init_max_hessian_value": "quasi_newton_init_max_value",
    "bfgs_init_const_hessian_value": "quasi_newton_init_const_value",
    "bfgs_hessian_memory": "lbfgs_memory",
    "bfgs_epsilon": "quasi_newton_epsilon",
}
BFGSInitializationTypes = QuasiNewtonInitializationTypes  # deprecated name


def _warn_deprecated_option(old, new):
    deprecation_warning(
        f"Hessian option '{old}' is deprecated, use '{new}' instead",
        version="0.4.1",
    )


def accepts_deprecated_hessian_options(func):
    """decorator, maps deprecated bfgs_* keyword arguments to current names"""

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        for old, new in DEPRECATED_HESSIAN_OPTIONS.items():
            if old in kwargs:
                value = kwargs.pop(old)
                _warn_deprecated_option(old, new)
                if new in kwargs and kwargs[new] != value:
                    raise ValueError(
                        f"Both '{old}' (deprecated) and '{new}' were given with different values"
                    )
                kwargs[new] = value
        return func(*args, **kwargs)

    return wrapper


_DIAGONAL_TYPES = (HessTypes.ZeroHessian, HessTypes.sparse_16)
_LAGRANGIAN_TYPES = (HessTypes.SR1_lagrangian,)
_SR1_TYPES = (HessTypes.SR1, HessTypes.SR1_lagrangian)
# types whose update needs B s != 0, they never move away from a zero matrix
_BFGS_FAMILY_TYPES = (
    HessTypes.BFGS,
    HessTypes.CBFGS,
    HessTypes.BFGS_mod,
    HessTypes.BFGS_damp,
    HessTypes.BFGS_ipopt,
    HessTypes.LBFGS,
)
_INITIALIZATION_TYPES = tuple(
    v for k, v in vars(QuasiNewtonInitializationTypes).items() if not k.startswith("_")
)
SR1_SAFEGUARD_CAP = 10.0


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

    def reset_memory(self):
        for mla in self.memory_limited_arrays:
            mla.clear()

    def memorize(self, inputs, jacobian):
        self.inputs.append(inputs.copy())
        self.jacobian.append(jacobian.copy())
        self.trim_memory()


class _Element:
    """one quasi-Newton matrix: the inputs it covers and the outputs it represents"""

    def __init__(self, support, outputs):
        self.support = np.asarray(support, dtype=int)
        self.outputs = list(outputs)
        self.matrix = None
        self.curvatures = []  # recent |y|/|s|, used by quasi_newton_norm_cap

    @property
    def size(self):
        return len(self.support)


class HessianApproximation:
    """general classs for hessian approximation methods"""

    @accepts_deprecated_hessian_options
    def __init__(
        self,
        hessian_type=None,
        quasi_newton_init_min_value=HESSIAN_DEFAULTS["quasi_newton_init_min_value"],
        quasi_newton_init_max_value=HESSIAN_DEFAULTS["quasi_newton_init_max_value"],
        quasi_newton_init_const_value=HESSIAN_DEFAULTS["quasi_newton_init_const_value"],
        quasi_newton_initialization_type=HESSIAN_DEFAULTS[
            "quasi_newton_initialization_type"
        ],
        lbfgs_memory=HESSIAN_DEFAULTS["lbfgs_memory"],
        quasi_newton_epsilon=HESSIAN_DEFAULTS["quasi_newton_epsilon"],
        sr1_skip_tolerance=HESSIAN_DEFAULTS["sr1_skip_tolerance"],
        quasi_newton_update_cap=HESSIAN_DEFAULTS["quasi_newton_update_cap"],
        quasi_newton_norm_cap=HESSIAN_DEFAULTS["quasi_newton_norm_cap"],
        hessian_support=None,
        input_scaling=None,
    ):
        """initialize hessian approximation class
        Args:
            hessian_type: type of hessian approximation to use (see HessTypes)
            quasi_newton_init_min_value: minimum absolute value for initial matrix elements
            quasi_newton_init_max_value: maximum absolute value for initial matrix elements
            quasi_newton_init_const_value: constant value for initial matrix diagonal
            quasi_newton_initialization_type: type of initialization for quasi-Newton matrix
            lbfgs_memory: number of iterates stored (used by LBFGS)
            quasi_newton_epsilon: epsilon value for numerical stability in BFGS updates
            sr1_skip_tolerance: SR1 update is skipped if |r's| < tol*|s|*|r|
            quasi_newton_update_cap: skip updates with |dB| > cap*max(|B|, |y|/|s|),
                cumulative types only, 0 disables, None uses SR1_SAFEGUARD_CAP for SR1
                types and disables for others
            quasi_newton_norm_cap: scale matrices down to cap*max recent |y|/|s|,
                all quasi-Newton types, 0 disables, None same as update cap
            hessian_support: list (one entry per output) of input indices the output
                depends on, None means every output depends on every input
            input_scaling: array of input scaling factors D (x_scaled = D*x), None
                means no scaling, can also be set later with set_input_scaling
        """
        self.hessian_memory = HessianMemory(memory=lbfgs_memory)
        if hessian_type is None:
            self.hessian_matrix_type = HESSIAN_DEFAULTS["hessian_type"]
        else:
            self.hessian_matrix_type = hessian_type
        self.iters = 0
        self.new_point = True
        self.epsilon = quasi_newton_epsilon
        self.sr1_skip_tolerance = sr1_skip_tolerance
        self.update_cap = self._resolve_cap(quasi_newton_update_cap)
        self.norm_cap = self._resolve_cap(quasi_newton_norm_cap)
        # safeguard counters, for diagnostics
        self.safeguard_stats = {"updates_capped": 0, "norms_capped": 0}
        self.init_min_hessian_value = quasi_newton_init_min_value
        self.init_max_hessian_value = quasi_newton_init_max_value
        self.init_const_hessian_value = quasi_newton_init_const_value
        self.bfgs_matrix_not_initialized = True
        if quasi_newton_initialization_type not in _INITIALIZATION_TYPES:
            raise ValueError(
                f"Unknown quasi-Newton initialization type {quasi_newton_initialization_type}, "
                f"expected one of {_INITIALIZATION_TYPES}"
            )
        if (
            quasi_newton_initialization_type == QuasiNewtonInitializationTypes.zero
            and self.hessian_matrix_type in _BFGS_FAMILY_TYPES
        ):
            _log.warning(
                f"Initialization type 'zero' with {self.hessian_matrix_type}: BFGS updates "
                "of a zero matrix are always skipped, the hessian will stay zero"
            )
        self.quasi_newton_initialization_type = quasi_newton_initialization_type
        self.hessian_support = hessian_support
        self.elements = None
        self.hessian_matrix = None
        self.input_scaling = None
        if input_scaling is not None:
            self.set_input_scaling(input_scaling)
        self._update_rules = {
            HessTypes.BFGS: self._rule_bfgs,
            HessTypes.CBFGS: self._rule_cbfgs,
            HessTypes.BFGS_mod: self._rule_modified_bfgs,
            HessTypes.BFGS_damp: self._rule_damped_bfgs,
            HessTypes.BFGS_ipopt: self._rule_ipopt_bfgs,
            HessTypes.SR1: self._rule_sr1,
            HessTypes.PSB: self._rule_psb,
            HessTypes.SR1_lagrangian: self._rule_sr1,
        }

    ########################################################################################
    # set up
    def _resolve_cap(self, cap):
        """None selects the default safeguard cap for the hessian type"""
        if cap is None:
            return SR1_SAFEGUARD_CAP if self.hessian_matrix_type in _SR1_TYPES else 0.0
        return float(cap)

    def set_input_scaling(self, input_scaling):
        """set scaling factors D for inputs (x_scaled = D*x), must be called before
        quasi-Newton matrices are created, otherwise stored matrices would be inconsistent
        """
        if self.elements is not None and self.elements[0].matrix is not None:
            raise RuntimeError(
                "Input scaling must be set before hessian approximation is started"
            )
        scaling = np.array(input_scaling, dtype=float)
        bad = ~np.isfinite(scaling) | (scaling <= 0)
        if np.any(bad):
            _log.warning(
                f"{np.sum(bad)} input scaling factors are not positive and finite, using 1 for those"
            )
            scaling[bad] = 1.0
        self.input_scaling = scaling

    def _setup(self, n_inputs, n_outputs):
        if self.hessian_support is None:
            supports = [np.arange(n_inputs)] * n_outputs
        else:
            if len(self.hessian_support) != n_outputs:
                raise ValueError(
                    f"hessian_support has {len(self.hessian_support)} entries, expected {n_outputs} (one per output)"
                )
            supports = [
                np.unique(np.asarray(s, dtype=int)) for s in self.hessian_support
            ]
            for s in supports:
                if s.size and (s.min() < 0 or s.max() >= n_inputs):
                    raise ValueError("hessian_support has input index out of range")
        if self.input_scaling is None:
            self.input_scaling = np.ones(n_inputs)
        elif len(self.input_scaling) != n_inputs:
            raise ValueError(
                f"input_scaling has {len(self.input_scaling)} entries, expected {n_inputs}"
            )
        if self.hessian_matrix_type in _LAGRANGIAN_TYPES:
            groups = {}
            for i, s in enumerate(supports):
                groups.setdefault(tuple(s), []).append(i)
            self.elements = [_Element(s, outs) for s, outs in groups.items()]
        else:
            self.elements = [_Element(s, [i]) for i, s in enumerate(supports)]

        # fixed sparsity pattern, PyNumero requires same structure on every call
        if self.hessian_matrix_type in _DIAGONAL_TYPES:
            rows = np.arange(n_inputs)
            cols = np.arange(n_inputs)
        else:
            mask = np.zeros((n_inputs, n_inputs), dtype=bool)
            for s in supports:
                mask[np.ix_(s, s)] = True
            rows, cols = np.nonzero(np.tril(mask))
        self._pattern = (rows, cols)
        self._shape = (n_inputs, n_inputs)

    def _pattern_matrix(self, values):
        rows, cols = self._pattern
        return coo_matrix((values, (rows, cols)), shape=self._shape)

    def zero_hessian(self):
        """zero valued hessian with the fixed sparsity pattern"""
        return self._pattern_matrix(np.zeros(len(self._pattern[0])))

    ########################################################################################
    # element vectors (in scaled space)
    def _element_vectors(self, element, idx_old, idx_new):
        """return scaled step s, gradient change y and current gradient g of element
        between memory entries idx_old and idx_new"""
        sup = element.support
        d = self.input_scaling[sup]
        x_old = self.hessian_memory.inputs[idx_old][sup]
        x_new = self.hessian_memory.inputs[idx_new][sup]
        j_old = self.hessian_memory.jacobian[idx_old][:, sup]
        j_new = self.hessian_memory.jacobian[idx_new][:, sup]
        if self.hessian_matrix_type in _LAGRANGIAN_TYPES:
            # gradient of the element's part of the Lagrangian, using current multipliers
            # for both points (same as ipopt)
            lam = self._outputs_dual_multipliers[element.outputs]
            g_old = lam @ j_old[element.outputs, :]
            g_new = lam @ j_new[element.outputs, :]
        else:
            g_old = j_old[element.outputs[0], :]
            g_new = j_new[element.outputs[0], :]
        s = (x_new - x_old) * d
        y = (g_new - g_old) / d
        g = g_new / d
        return s, y, g

    ########################################################################################
    # initial matrix
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

    def get_initial_matrix(self, s, y, g):
        """initial matrix of an element from step s, gradient change y and current gradient g"""
        n = len(s)
        init_type = self.quasi_newton_initialization_type
        if init_type == QuasiNewtonInitializationTypes.zero:
            return np.zeros((n, n))
        sigma = self.init_const_hessian_value
        if (
            _nonzero(y)
            and _nonzero(s)
            and init_type != QuasiNewtonInitializationTypes.constant
        ):
            sTy = s @ y
            sTs = s @ s
            yTy = y @ y
            if init_type == QuasiNewtonInitializationTypes.scalar1:
                if sTs != 0:
                    sigma = sTy / sTs
            elif init_type == QuasiNewtonInitializationTypes.scalar2:
                if sTy != 0:
                    sigma = yTy / sTy
            elif init_type == QuasiNewtonInitializationTypes.scalar3:
                if sTs != 0:
                    sigma = sTy / sTs / 2
                    if sTy != 0:
                        sigma += (yTy / sTy) / 2
            elif init_type == QuasiNewtonInitializationTypes.scalar4:
                if sTs != 0:
                    sigma = sTy / sTs
                if sTy != 0:
                    sigma *= yTy / sTy
                sigma = np.sign(sigma) * np.sqrt(np.abs(sigma))
            elif init_type == QuasiNewtonInitializationTypes.GaussNewton:
                # diagonal of g g^T
                return self.apply_sigma_bounds(np.diag(g**2))
            else:
                raise NotImplementedError(
                    f"BFGS initialization type {init_type} not implemented"
                )
        return self.apply_sigma_bounds(np.identity(n) * sigma)

    ########################################################################################
    # update rules, all operate on 1D vectors s, y, g of one element and return new matrix
    # or None if update is skipped
    def _rule_bfgs(self, B, s, y, g):
        """Vanilla BFGS update"""
        Hs = B @ s
        sHs = s @ Hs
        if y @ s > self.epsilon and _nonzero(Hs) and sHs != 0:
            return B + np.outer(y, y) / (y @ s) - np.outer(Hs, Hs) / sHs
        return None

    def _rule_cbfgs(self, B, s, y, g):
        """Cautious BFGS update (Li and Fukushima)"""
        Hs = B @ s
        sHs = s @ Hs
        y_s = y @ s
        g_norm = np.linalg.norm(g)
        alpha = 0.01 if g_norm >= 1 else 3
        min_test = self.epsilon * g_norm**alpha
        max_test = y_s / np.linalg.norm(s) ** 2
        if max_test > min_test and _nonzero(Hs) and sHs != 0:
            return B + np.outer(y, y) / y_s - np.outer(Hs, Hs) / sHs
        return None

    def _rule_modified_bfgs(self, B, s, y, g):
        """Modified BFGS update (Li and Fukushima)"""
        Hs = B @ s
        sHs = s @ Hs
        if _nonzero(s) and _nonzero(Hs):
            y_s = y @ s
            t_k = 1 + max(0, -y_s / (np.linalg.norm(s) ** 2))
            z = y + t_k * np.linalg.norm(g) * s
            z_s = z @ s
            if z_s != 0 and sHs != 0:
                return B + np.outer(z, z) / z_s - np.outer(Hs, Hs) / sHs
        return None

    def _rule_damped_bfgs(self, B, s, y, g):
        """Powell's damping on the BFGS update"""
        phi = 0.9
        Hs = B @ s
        if _nonzero(Hs) and _nonzero(s):
            y_s = y @ s
            sHs = s @ Hs
            if sHs > 0:
                # full BFGS update when curvature condition holds (delta = 1), damped
                # otherwise (y_s < phi*sHs guarantees sHs - y_s > 0)
                if y_s >= phi * sHs:
                    delta = 1
                else:
                    delta = (1 - phi) * sHs / (sHs - y_s)
                z = delta * y + (1 - delta) * Hs
                z_s = z @ s
                if z_s != 0:
                    return B + np.outer(z, z) / z_s - np.outer(Hs, Hs) / sHs
        return None

    def _rule_ipopt_bfgs(self, B, s, y, g):
        """BFGS update only done on certain conditions (taken from IPOPT's implementation)"""
        Hs = B @ s
        sHs = s @ Hs
        y_s = y @ s
        if (
            y_s > np.sqrt(self.epsilon) * np.linalg.norm(s) * np.linalg.norm(y)
            and np.linalg.norm(s, np.inf) >= self.epsilon
            and _nonzero(Hs)
            and sHs != 0
        ):
            return B + np.outer(y, y) / y_s - np.outer(Hs, Hs) / sHs
        return None

    def _rule_sr1(self, B, s, y, g):
        """Symmetric rank one update, allows indefinite matrices.
        Skipped if |r's| < tol*|s|*|r| (Nocedal & Wright, eq. 6.26)"""
        r = y - B @ s
        r_s = r @ s
        if (
            _nonzero(r)
            and r_s != 0
            and abs(r_s)
            >= self.sr1_skip_tolerance * np.linalg.norm(s) * np.linalg.norm(r)
        ):
            return B + np.outer(r, r) / r_s
        return None

    def _rule_psb(self, B, s, y, g):
        """Powell symmetric Broyden update, always defined for s != 0 and does
        not require positive definiteness"""
        s_s = s @ s
        if s_s > 0:
            r = y - B @ s
            if not _nonzero(r):
                return None
            r_s = r @ s
            return (
                B
                + (np.outer(r, s) + np.outer(s, r)) / s_s
                - r_s * np.outer(s, s) / s_s**2
            )
        return None

    ########################################################################################
    # quasi-Newton driver
    def _memorize_point(self):
        """only store a new (s, y) pair if inputs changed since last call, re-evaluating
        at the same point (e.g. with new multipliers) must not consume history or
        re-apply the last update"""
        self.new_point = len(self.hessian_memory.inputs) == 0 or not np.array_equal(
            self.hessian_memory.inputs[-1], self.inputs
        )
        if self.new_point:
            self.hessian_memory.memorize(self.inputs, self.jacobian_matrix)

    def check_step(self):
        """check if last stored step is valid for update (non zero)"""
        if len(self.hessian_memory.inputs) < 2:
            return False
        s = self.hessian_memory.inputs[-1] - self.hessian_memory.inputs[-2]
        return _nonzero(s)

    def _initialize_elements(self):
        """create matrices for elements
        - before first valid step, use constant diagonal
        - on first valid step (and after all duals were zero), scale initial matrix from
          that step, that same step is then also used as first update"""
        for el in self.elements:
            if el.matrix is None:
                if (
                    self.quasi_newton_initialization_type
                    == QuasiNewtonInitializationTypes.zero
                ):
                    el.matrix = np.zeros((el.size, el.size))
                else:
                    el.matrix = np.identity(el.size) * self.init_const_hessian_value
        if self.hessian_matrix_type == HessTypes.LBFGS:
            # LBFGS rebuilds matrices from its own initial matrix on every update
            return
        if self.check_step() and self.bfgs_matrix_not_initialized:
            self.bfgs_matrix_not_initialized = False
            n = len(self.hessian_memory.inputs)
            for el in self.elements:
                el.matrix = self.get_initial_matrix(
                    *self._element_vectors(el, n - 2, n - 1)
                )
        if not _nonzero(self._outputs_dual_multipliers):
            self.bfgs_matrix_not_initialized = True

    ########################################################################################
    # safeguards
    @staticmethod
    def _curvature(s, y):
        """magnitude of curvature observed along step, |y|/|s|"""
        s_norm = np.linalg.norm(s)
        return np.linalg.norm(y) / s_norm if s_norm > 0 else 0.0

    def _record_curvature(self, el, s, y):
        el.curvatures.append(self._curvature(s, y))
        if len(el.curvatures) > self.hessian_memory.memory:
            el.curvatures.pop(0)

    def _update_too_large(self, B, new_matrix, s, y):
        """quasi_newton_update_cap test, True if update should be skipped"""
        if self.update_cap <= 0:
            return False
        ref = max(np.linalg.norm(B), self._curvature(s, y))
        if ref == 0:
            return True
        return np.linalg.norm(new_matrix - B) > self.update_cap * ref

    def _apply_norm_cap(self, el, curvatures):
        """quasi_newton_norm_cap, scale matrix down to cap * max recent curvature"""
        if self.norm_cap <= 0 or not curvatures:
            return
        limit = self.norm_cap * max(curvatures)
        norm = np.linalg.norm(el.matrix)
        if limit > 0 and norm > limit:
            el.matrix = el.matrix * (limit / norm)
            self.safeguard_stats["norms_capped"] += 1

    def _cumulative_update(self, rule):
        """update every element matrix with most recent stored step"""
        if not (self.new_point and len(self.hessian_memory.inputs) > 1):
            return
        n = len(self.hessian_memory.inputs)
        for el in self.elements:
            s, y, g = self._element_vectors(el, n - 2, n - 1)
            if not _nonzero(s):
                continue
            self._record_curvature(el, s, y)
            new_matrix = rule(el.matrix, s, y, g)
            if new_matrix is not None:
                if self._update_too_large(el.matrix, new_matrix, s, y):
                    self.safeguard_stats["updates_capped"] += 1
                else:
                    el.matrix = new_matrix
            self._apply_norm_cap(el, el.curvatures)

    def _ipopt_update_test(self, s, y, initial_matrix):
        y_s = y @ s
        return bool(
            y_s > np.sqrt(self.epsilon) * np.linalg.norm(s) * np.linalg.norm(y)
            and np.linalg.norm(s, np.inf) >= self.epsilon
            and _nonzero(y)
            and _nonzero(s)
            and _nonzero(initial_matrix)
        )

    def _lbfgs_update(self):
        """Direct LBFGS based on Byrd, Nocedal & Schnabel, Representations of quasi-newton
        matrices and their use in limited memory methods. Matrix is rebuilt from the
        stored pairs on every call, using single initial matrix scaled from most recent pair
        """
        n = len(self.hessian_memory.inputs)
        # matrices only depend on stored pairs, which do not change on a repeated point
        if n < 2 or not self.new_point:
            return
        for el in self.elements:
            pairs = [self._element_vectors(el, r - 1, r) for r in range(1, n)]
            s_last, y_last, g_last = pairs[-1]
            initial_matrix = self.get_initial_matrix(s_last, y_last, g_last)
            el.matrix = initial_matrix.copy()
            # only update if current step is good
            if not self._ipopt_update_test(s_last, y_last, initial_matrix):
                continue
            ak = []
            bk = []
            for s, y, _ in pairs:
                # only include update if sub step is good
                if self._ipopt_update_test(s, y, initial_matrix):
                    b = y / np.sqrt(y @ s)
                    a = initial_matrix @ s
                    for k in range(len(ak)):
                        a = a + (bk[k] @ s) * bk[k] - (ak[k] @ s) * ak[k]
                    s_a = s @ a
                    if s_a > 0:
                        ak.append(a / np.sqrt(s_a))
                        bk.append(b)
            for a, b in zip(ak, bk):
                el.matrix = el.matrix + np.outer(b, b) - np.outer(a, a)
            self._apply_norm_cap(el, [self._curvature(s, y) for s, y, _ in pairs])

    def _reset_non_finite(self):
        """reset elements whose matrix became non finite"""
        bad = [el for el in self.elements if not np.all(np.isfinite(el.matrix))]
        if bad:
            _log.warning(
                f"{len(bad)} hessian approximation matrices were not finite, resetting them"
            )
            for el in bad:
                if (
                    self.quasi_newton_initialization_type
                    == QuasiNewtonInitializationTypes.zero
                ):
                    el.matrix = np.zeros((el.size, el.size))
                else:
                    el.matrix = np.identity(el.size) * self.init_const_hessian_value
            self.bfgs_matrix_not_initialized = True

    def _assemble(self):
        """assemble hessian in unscaled space: H = sum lambda_i D B_i D"""
        H = np.zeros(self._shape)
        for el in self.elements:
            d = self.input_scaling[el.support]
            unscaled = el.matrix * np.outer(d, d)
            if self.hessian_matrix_type in _LAGRANGIAN_TYPES:
                weight = 1.0
            else:
                weight = self._outputs_dual_multipliers[el.outputs[0]]
            H[np.ix_(el.support, el.support)] += weight * unscaled
        return H

    def quasi_newton_hessian(self):
        self._memorize_point()
        self._initialize_elements()
        if self.hessian_matrix_type == HessTypes.LBFGS:
            self._lbfgs_update()
        else:
            self._cumulative_update(self._update_rules[self.hessian_matrix_type])
        self._reset_non_finite()
        self.iters += 1
        self.hessian_matrix = self._assemble()

    def hessian_gauss_newton_version(self):
        """naive gauss newton hessian approximation, sum lambda_i J_i^T J_i"""
        H = np.zeros(self._shape)
        for i, el in enumerate(self.elements):
            row = self.jacobian_matrix[i, el.support]
            H[np.ix_(el.support, el.support)] += self._outputs_dual_multipliers[
                i
            ] * np.outer(row, row)
        self.hessian_matrix = H

    ########################################################################################
    # main entry
    def get_hessian(self, input_values, jacobian, dual_multipliers):
        """return lower triangular hessian of sum_i lambda_i y_i(x) as coo matrix
        with a fixed sparsity pattern"""
        self.inputs = np.array(input_values, dtype=float)
        self.jacobian_matrix = np.array(jacobian, dtype=float)
        self._outputs_dual_multipliers = np.asarray(dual_multipliers, dtype=float)
        if self.elements is None:
            self._setup(len(self.inputs), self.jacobian_matrix.shape[0])
        rows, cols = self._pattern
        try:
            if self.hessian_matrix_type == HessTypes.ZeroHessian:
                return self._pattern_matrix(np.zeros(len(rows)))
            elif self.hessian_matrix_type == HessTypes.sparse_16:
                return self._pattern_matrix(np.full(len(rows), 1e-16))
            elif self.hessian_matrix_type == HessTypes.GaussNewton:
                self.hessian_gauss_newton_version()
            elif self.hessian_matrix_type in self._update_rules or (
                self.hessian_matrix_type == HessTypes.LBFGS
            ):
                self.quasi_newton_hessian()
            else:
                raise NotImplementedError(
                    f"Hessian type {self.hessian_matrix_type} not implemented"
                )
            return self._pattern_matrix(self.hessian_matrix[rows, cols])
        except NotImplementedError:
            raise
        except Exception as e:
            # never hand None back to PyNumero, fall back to a zero hessian (same
            # sparsity pattern) so the solve can continue, but make the failure visible
            _log.warning(
                f"Error in Hessian approximation: {e!r}. "
                f"Hessian type: {self.hessian_matrix_type}, "
                f"Inputs: {self.inputs}, "
                f"Jacobian: {self.jacobian_matrix}, "
                f"Dual multipliers: {self._outputs_dual_multipliers}. "
                f"Falling back to zero hessian for this evaluation.",
                exc_info=True,
            )
            return self.zero_hessian()
