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
"""Algebraic surrogate of a reaktoro gray box, used by the trust region backend.

The block has the same `inputs` and `outputs` variables as the ExternalGreyBoxBlock it
replaces, so all linking constraints and scaling built by ReaktoroBlockBuilder work
unchanged. Outputs are given by a taylor model around an expansion point w0

    outputs[o] == y0[o] + sum_i J[o,i] d_i + 0.5 sum_ij B[o,i,j] d_i d_j

where d = inputs - w0. During a trust region solve (see trust_region_solver.py)

    (inputs_i - w0_i) / h_i == s_i + ep_i - en_i,  |s_i| <= radius, ep, en >= 0 (penalized)

is added (elastic trust region). With the soft formulation (default) the model is
evaluated at the actual inputs, so the step beyond the radius is an exact l1 penalty
and the subproblem is never infeasible because of the trust region. With the
trusted_point formulation the model is evaluated at d_i = h_i * s_i instead.
h_i is the step scale of input i, its current magnitude (or log_input_step for pH like inputs).

The truth model (reaktoro), output floors and quasi-Newton curvature are provided
by the ReaktoroGrayBox external model, which is not added to the pyomo model.
"""

import numpy as np
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from pyomo.core.base.block import BlockData, declare_custom_block
from pyomo.environ import (
    Set,
    Var,
    Param,
    Expression,
    Constraint,
    NonNegativeReals,
    value,
)
from pyomo.common.modeling import unique_component_name

from reaktoro_pse.core.util_classes.rkt_inputs import RktInputTypes
from reaktoro_pse.core.util_classes.hessian_functions import (
    HessTypes,
    _LAGRANGIAN_TYPES,
)
from reaktoro_pse.reaktoro_block_config.trust_region_options import (
    SurrogateTypes,
    TrustRegionOptions,
    TrustRegionFormulations,
)

__author__ = "Alexander V. Dudchenko"
_log = idaeslog.getLogger(__name__)

# inputs on a log scale, trust region on these is absolute (log_input_step)
LOG_INPUT_TYPES = (RktInputTypes.pH, RktInputTypes.pOH, RktInputTypes.pE, "Eh")
# hessian types that do not provide output curvature for quadratic surrogate
LINEAR_HESSIAN_TYPES = (
    HessTypes.ZeroHessian,
    HessTypes.GaussNewton,
    HessTypes.no_hessian_estimation,
    HessTypes.sparse_16,
)


def resolve_surrogate_type(surrogate, hessian_type):
    """returns linear or quadratic for configured surrogate and hessian type"""
    if surrogate == SurrogateTypes.auto:
        if hessian_type in LINEAR_HESSIAN_TYPES or hessian_type in _LAGRANGIAN_TYPES:
            return SurrogateTypes.linear
        return SurrogateTypes.quadratic
    if surrogate == SurrogateTypes.quadratic and (
        hessian_type in LINEAR_HESSIAN_TYPES or hessian_type in _LAGRANGIAN_TYPES
    ):
        raise ValueError(
            f"Quadratic surrogate requires a per output quasi-Newton hessian_type "
            f"(e.g. SR1, LBFGS, BFGS), got {hessian_type}"
        )
    return surrogate


def _input_key(name):
    """original reaktoro input key, aggregated (manager) inputs are (block index, key)"""
    if isinstance(name, tuple) and len(name) >= 2 and isinstance(name[0], int):
        return name[1] if len(name) == 2 else name[1:]
    return name


@declare_custom_block(name="ReaktoroSurrogateBlock")
class ReaktoroSurrogateBlockData(BlockData):
    """pyomo block with algebraic surrogate of reaktoro outputs (see module doc)"""

    def build_surrogate(
        self,
        external_model,
        trust_region_options=None,
        hessian_type=None,
    ):
        """build surrogate variables and constraints

        Args:
            external_model: configured ReaktoroGrayBox, provides truth model, output floors,
                scaling and quasi-Newton curvature
            trust_region_options: trust region ConfigDict (TrustRegionOptions)
            hessian_type: hessian type from hessian_options, sets surrogate type for 'auto'
        """
        self.external_model = external_model
        self.block_index = None
        if trust_region_options is None:
            trust_region_options = TrustRegionOptions().get_dict()()
        self.trust_region_options = trust_region_options
        if hessian_type is None:
            hessian_type = getattr(
                getattr(external_model, "hessian_calculator", None),
                "hessian_matrix_type",
                HessTypes.no_hessian_estimation,
            )
        self.hessian_type = hessian_type
        self.surrogate_type = resolve_surrogate_type(
            trust_region_options.surrogate, hessian_type
        )
        self._input_names = list(external_model.input_names())
        self._output_names = list(external_model.output_names())
        n_in = len(self._input_names)
        n_out = len(self._output_names)
        self._input_names_set = Set(
            initialize=self._input_names, ordered=True, dimen=None
        )
        self.inputs = Var(self._input_names_set)
        self._output_names_set = Set(
            initialize=self._output_names, ordered=True, dimen=None
        )
        self.outputs = Var(self._output_names_set)
        # sets bounds and initial values the same way as for the gray box
        external_model.finalize_block_construction(self)
        self._input_vars = [self.inputs[n] for n in self._input_names]
        self._output_vars = [self.outputs[n] for n in self._output_names]
        self._log_inputs = np.array(
            [self._is_log_input(n) for n in self._input_names], dtype=bool
        )

        calc = getattr(external_model, "hessian_calculator", None)
        support = None if calc is None else calc.hessian_support
        if support is None:
            self._support = [np.arange(n_in) for _ in range(n_out)]
        else:
            self._support = [np.unique(np.asarray(s, dtype=int)) for s in support]

        self._in_idx = Set(initialize=range(n_in), ordered=True)
        self._out_idx = Set(initialize=range(n_out), ordered=True)
        self.w0 = Param(self._in_idx, mutable=True, initialize=0)
        self.h = Param(self._in_idx, mutable=True, initialize=1)
        self.y0 = Param(self._out_idx, mutable=True, initialize=0)
        self._jac_index = Set(
            initialize=[(o, int(i)) for o in range(n_out) for i in self._support[o]],
            dimen=2,
            ordered=True,
        )
        self.J = Param(self._jac_index, mutable=True, initialize=0)
        if self.surrogate_type == SurrogateTypes.quadratic:
            self._hess_index = Set(
                initialize=[
                    (o, int(i), int(j))
                    for o in range(n_out)
                    for a, i in enumerate(self._support[o])
                    for j in self._support[o][: a + 1]
                ],
                dimen=3,
                ordered=True,
            )
            self.B = Param(self._hess_index, mutable=True, initialize=0)

        @self.Expression(self._in_idx)
        def delta(b, i):
            return b._input_vars[i] - b.w0[i]

        # elastic slack on surrogate equations, zero unless model slack is enabled
        # during trust region solves (see enable_model_slack)
        @self.Expression(self._out_idx)
        def out_slack(b, o):
            return 0

        @self.Constraint(self._out_idx)
        def surrogate_eq(b, o):
            sup = b._support[o]
            expr = b.y0[o] + sum(b.J[o, int(i)] * b.delta[int(i)] for i in sup)
            if b.surrogate_type == SurrogateTypes.quadratic:
                quad = 0
                for a, i in enumerate(sup):
                    i = int(i)
                    quad += 0.5 * b.B[o, i, i] * b.delta[i] ** 2
                    for j in sup[:a]:
                        j = int(j)
                        quad += b.B[o, i, j] * b.delta[i] * b.delta[j]
                expr = expr + quad
            return b._output_vars[o] == expr + b.out_slack[o]

        self._elastic_names = None
        self.model_slack_enabled = False
        self.expansion_point = None  # np array of inputs at expansion point
        self.reset_curvature()

    ####################################################################################
    # compatibility with ExternalGreyBoxBlock
    def get_external_model(self):
        return self.external_model

    ####################################################################################
    # helpers
    def _is_log_input(self, name):
        key = _input_key(name)
        if key in LOG_INPUT_TYPES:
            return True
        # keys of speciation/multi system inputs are (system, input type)
        if isinstance(key, tuple) and key[-1] in LOG_INPUT_TYPES:
            return True
        obj = self.external_model.input_dict[name]
        get_type = getattr(obj, "get_input_type", None)
        return get_type is not None and get_type() in LOG_INPUT_TYPES

    def input_values(self):
        return np.array([value(v) for v in self._input_vars], dtype=float)

    def output_values(self):
        return np.array([value(v) for v in self._output_vars], dtype=float)

    def input_scaling(self):
        return np.array(
            [iscale.get_scaling_factor(v, default=1) for v in self._input_vars],
            dtype=float,
        )

    def output_scaling(self):
        return np.array(
            [iscale.get_scaling_factor(v, default=1) for v in self._output_vars],
            dtype=float,
        )

    def step_scale(self, w):
        """per input step scale h, trust region is |w - w0| <= radius * h"""
        sf = self.input_scaling()
        sf[sf <= 0] = 1
        h = np.maximum(np.abs(w), 1 / sf)
        h[self._log_inputs] = self.trust_region_options.log_input_step
        return h

    ####################################################################################
    # truth model
    def evaluate_truth(self, w=None):
        """evaluate reaktoro (with output floors) at inputs w (default: current inputs)
        returns (outputs, jacobian), raises on reaktoro failure"""
        if w is None:
            w = self.input_values()
        em = self.external_model
        em.set_input_values(w)
        outputs = np.array(em.evaluate_outputs(), dtype=float)
        jacobian = np.array(em.jacobian_matrix, dtype=float)
        return outputs, jacobian

    def set_expansion_point(self, w, outputs, jacobian, update_curvature=True):
        """set taylor model around w"""
        w = np.asarray(w, dtype=float)
        h = self.step_scale(w)
        for i in range(len(w)):
            self.w0[i].set_value(float(w[i]))
            self.h[i].set_value(float(h[i]))
        for o in range(len(outputs)):
            self.y0[o].set_value(float(outputs[o]))
        for o, i in self._jac_index:
            self.J[o, i].set_value(float(jacobian[o, i]))
        if self.surrogate_type == SurrogateTypes.quadratic and update_curvature:
            self._update_curvature(w, jacobian)
        self.expansion_point = w.copy()

    def _update_curvature(self, w, jacobian):
        calc = self.external_model.hessian_calculator
        em = self.external_model
        if em.hessian_input_scaling and calc.input_scaling is None:
            calc.set_input_scaling(em.get_input_scaling_factors())
        # unit multipliers, all outputs keep their own curvature
        calc.get_hessian(w, jacobian, np.ones(len(self._output_names)))
        matrices = calc.get_output_matrices()
        for o, (sup, mat) in matrices.items():
            mat = np.where(np.isfinite(mat), mat, 0.0)
            for a, i in enumerate(sup):
                for b in range(a + 1):
                    self.B[o, int(i), int(sup[b])].set_value(float(mat[a, b]))

    def reset_curvature(self):
        """reset quadratic surrogate curvature and quasi-Newton history"""
        if self.surrogate_type != SurrogateTypes.quadratic:
            return
        calc = self.external_model.hessian_calculator
        calc.reset()
        for idx in self._hess_index:
            self.B[idx].set_value(0)

    def surrogate_outputs(self, w):
        """value of surrogate at inputs w (outside of elastic mode)"""
        d = np.asarray(w, dtype=float) - self.expansion_point
        y = np.array([value(self.y0[o]) for o in self._out_idx])
        for o, i in self._jac_index:
            y[o] += value(self.J[o, i]) * d[i]
        if self.surrogate_type == SurrogateTypes.quadratic:
            for o, i, j in self._hess_index:
                f = 0.5 if i == j else 1.0
                y[o] += f * value(self.B[o, i, j]) * d[i] * d[j]
        return y

    def set_constraint_scaling(self):
        """scale surrogate constraints the same way gray box output constraints are"""
        try:
            sf = self.external_model.get_output_constraint_scaling_factors()
        except Exception:  # scaling not set up (block not initialized)
            sf = None
        if sf is None:
            sf = self.output_scaling()
        for o in self._out_idx:
            iscale.set_scaling_factor(self.surrogate_eq[o], float(sf[o]))

    ####################################################################################
    # elastic trust region (only active during trust region solves)
    @property
    def elastic_mode(self):
        return self._elastic_names is not None

    def enable_elastic(self):
        """add trust region step and elastic gap variables, returns penalty expression
        (sum of elastic gaps in units of the step scale h)"""
        if self.elastic_mode:
            return self._penalty_expr
        names = {}
        for base in (
            "tr_step",
            "tr_gap_pos",
            "tr_gap_neg",
            "tr_elastic_eq",
            "tr_model_slack_pos",
            "tr_model_slack_neg",
        ):
            names[base] = unique_component_name(self, base)
        self.add_component(names["tr_step"], Var(self._in_idx, initialize=0))
        self.add_component(
            names["tr_gap_pos"],
            Var(self._in_idx, initialize=0, domain=NonNegativeReals),
        )
        self.add_component(
            names["tr_gap_neg"],
            Var(self._in_idx, initialize=0, domain=NonNegativeReals),
        )
        s = self.component(names["tr_step"])
        ep = self.component(names["tr_gap_pos"])
        en = self.component(names["tr_gap_neg"])
        w = self.input_values()
        for i in self._in_idx:
            # start at current point, all movement is inside the trust region
            s[i].set_value((w[i] - value(self.w0[i])) / value(self.h[i]))
            iscale.set_scaling_factor(s[i], 1)
            iscale.set_scaling_factor(ep[i], 1)
            iscale.set_scaling_factor(en[i], 1)

        def _elastic_rule(b, i):
            return (b._input_vars[i] - b.w0[i]) / b.h[i] == s[i] + ep[i] - en[i]

        self.add_component(
            names["tr_elastic_eq"], Constraint(self._in_idx, rule=_elastic_rule)
        )
        con = self.component(names["tr_elastic_eq"])
        for i in self._in_idx:
            iscale.set_scaling_factor(con[i], 1)
        if (
            self.trust_region_options.trust_region_formulation
            == TrustRegionFormulations.trusted_point
        ):
            for i in self._in_idx:
                self.delta[i].set_value(self.h[i] * s[i])
        # model slack, in units of scaled surrogate equation residual, fixed at zero
        # until enabled (when the subproblem is infeasible for any input, e.g. a
        # fixed spec far from the expansion point with an input at its bound)
        for base in ("tr_model_slack_pos", "tr_model_slack_neg"):
            self.add_component(
                names[base],
                Var(self._out_idx, initialize=0, domain=NonNegativeReals),
            )
        mp = self.component(names["tr_model_slack_pos"])
        mn = self.component(names["tr_model_slack_neg"])
        for o in self._out_idx:
            mp[o].fix(0)
            mn[o].fix(0)
            iscale.set_scaling_factor(mp[o], 1)
            iscale.set_scaling_factor(mn[o], 1)
            sf = iscale.get_scaling_factor(self.surrogate_eq[o], default=1)
            sf = sf if sf is not None and sf > 0 else 1
            self.out_slack[o].set_value((mp[o] - mn[o]) / sf)
        self.model_slack_enabled = False
        self._elastic_names = names
        self._penalty_expr = sum(ep[i] + en[i] for i in self._in_idx) + sum(
            mp[o] + mn[o] for o in self._out_idx
        )
        return self._penalty_expr

    def enable_model_slack(self, enable=True):
        """free (or fix at zero) elastic slack on surrogate equations"""
        mp = self.component(self._elastic_names["tr_model_slack_pos"])
        mn = self.component(self._elastic_names["tr_model_slack_neg"])
        for o in self._out_idx:
            if enable:
                mp[o].unfix()
                mn[o].unfix()
                mp[o].set_value(0)
                mn[o].set_value(0)
            else:
                mp[o].fix(0)
                mn[o].fix(0)
        self.model_slack_enabled = enable

    def model_slack(self):
        """sum of surrogate equation slack (scaled residual units)"""
        if not self.elastic_mode:
            return 0.0
        mp = self.component(self._elastic_names["tr_model_slack_pos"])
        mn = self.component(self._elastic_names["tr_model_slack_neg"])
        return float(sum(abs(value(mp[o]) - value(mn[o])) for o in self._out_idx))

    def set_trust_radius(self, radius):
        s = self.component(self._elastic_names["tr_step"])
        for i in self._in_idx:
            s[i].setlb(-radius)
            s[i].setub(radius)
            if s[i].value > radius:
                s[i].set_value(radius)
            elif s[i].value < -radius:
                s[i].set_value(-radius)

    def reset_elastic_values(self):
        """put step and gaps at values consistent with current inputs and radius"""
        s = self.component(self._elastic_names["tr_step"])
        ep = self.component(self._elastic_names["tr_gap_pos"])
        en = self.component(self._elastic_names["tr_gap_neg"])
        w = self.input_values()
        for i in self._in_idx:
            d = (w[i] - value(self.w0[i])) / value(self.h[i])
            lb, ub = s[i].lb, s[i].ub
            sd = min(max(d, lb if lb is not None else d), ub if ub is not None else d)
            s[i].set_value(sd)
            ep[i].set_value(max(d - sd, 0))
            en[i].set_value(max(sd - d, 0))

    def elastic_gap(self):
        """sum of elastic gaps, computed from current inputs"""
        if not self.elastic_mode:
            return 0.0
        s = self.component(self._elastic_names["tr_step"])
        w = self.input_values()
        gap = 0.0
        for i in self._in_idx:
            d = (w[i] - value(self.w0[i])) / value(self.h[i])
            gap += abs(d - value(s[i]))
        return gap + self.model_slack()

    def trusted_point(self):
        """inputs at which surrogate is evaluated in elastic mode"""
        s = self.component(self._elastic_names["tr_step"])
        return np.array(
            [value(self.w0[i]) + value(self.h[i]) * value(s[i]) for i in self._in_idx]
        )

    def disable_elastic(self):
        if not self.elastic_mode:
            return
        for i in self._in_idx:
            self.delta[i].set_value(self._input_vars[i] - self.w0[i])
        for o in self._out_idx:
            self.out_slack[o].set_value(0)
        self.model_slack_enabled = False
        for name in self._elastic_names.values():
            self.del_component(name)
        self._elastic_names = None
        self._penalty_expr = None

    def relative_step(self, w_new, w_old):
        """max relative step of inputs (in units of step scale at expansion point)"""
        h = np.array([value(self.h[i]) for i in self._in_idx])
        if len(h) == 0:
            return 0.0
        return float(np.max(np.abs(np.asarray(w_new) - np.asarray(w_old)) / h))
