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
"""Trust region filter method for models with ReaktoroSurrogateBlocks.

The algorithm is the trust region filter method of pyomo.contrib.trustregion
(Eason & Biegler 2016, 2018, Yoshio & Biegler 2021). Filter, iteration logger and options
are imported from pyomo, the main loop (trust_region_method) is a copy of
pyomo.contrib.trustregion.TRF.trust_region_method with these changes:
    - takes an interface instance (pyomo creates its own TRFInterface, which clones the
      model, counts fixed variables as degrees of freedom, puts the trust region on
      decision variables and does not support square problems)
    - surrogate is rebuilt before trust region bounds are set (bounds depend on the
      expansion point)
    - first step norm is inf, so at least one subproblem is solved
    - trial point has to be acceptable to the current iterate as well as the filter
      (Fletcher-Leyffer), pyomo only checks the filter
    - restoration step: a rejected step that needed elastic slack (subproblem is
      infeasible within the trust region) is accepted while the elastic gap decreases,
      the filter is reset and quadratic curvature is reset. Rejection shrinks the
      radius from min(step, radius)
    - no model.display() on exit, returns termination status and statistics
    - adaptive elastic penalty update after accepted steps
ReaktoroTRInterface replaces pyomo's TRFInterface for reaktoro surrogate blocks.
"""

import time
import numpy as np
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale
from pyomo.contrib.trustregion.filter import Filter, FilterElement
from pyomo.contrib.trustregion.util import IterationLogger
from pyomo.common.modeling import unique_component_name
from pyomo.environ import (
    Block,
    Objective,
    Param,
    Var,
    value,
    minimize,
    maximize,
    check_optimal_termination,
)
from pyomo.opt import SolverResults, SolverStatus, TerminationCondition

from reaktoro_pse.core.reaktoro_surrogate import ReaktoroSurrogateBlockData
from reaktoro_pse.reaktoro_block_config.trust_region_options import (
    TrustRegionOptions,
    ElasticPenaltyModes,
)

__author__ = "Alexander V. Dudchenko"
_log = idaeslog.getLogger(__name__)

try:
    import cyipopt

    _TRUTH_ERRORS = (cyipopt.CyIpoptEvaluationError, RuntimeError, ValueError)
except ImportError:  # pragma: no cover
    _TRUTH_ERRORS = (RuntimeError, ValueError)


def find_surrogate_blocks(model):
    """active ReaktoroSurrogateBlocks in model (including model itself)"""
    blocks = []
    if isinstance(model, ReaktoroSurrogateBlockData) and model.active:
        blocks.append(model)
    for blk in model.component_data_objects(Block, active=True, descend_into=True):
        if isinstance(blk, ReaktoroSurrogateBlockData):
            blocks.append(blk)
    return blocks


def get_trust_region_stats(results):
    """trust region statistics of a solve, None for solves without trust region"""
    return getattr(getattr(results, "solver", None), "trust_region", None)


class ReaktoroTRInterface:
    """Trust region interface for models with reaktoro surrogate blocks, implements the
    methods of pyomo.contrib.trustregion.interface.TRFInterface used by the main loop"""

    def __init__(
        self, model, surrogate_blocks, inner_solver, config, tee=False, solve_kwds=None
    ):
        self.model = model
        self.blocks = surrogate_blocks
        self.inner_solver = inner_solver
        self.config = config
        self.tee = tee
        self.solve_kwds = {} if solve_kwds is None else dict(solve_kwds)
        self.inner_solves = 0
        self.consecutive_failed_inner_solves = 0
        self._last_accepted_gap = None
        self.stats = {
            "inner_iterations": [],
            "inner_status": [],
            "inner_time": [],
            "truth_time": [],
            "reaktoro_failures": [],
            "theta_per_block": [],
            "elastic_gap": [],
            "elastic_penalty": [],
            "model_slack": [],
            "model_slack_enabled": 0,
        }
        self.total_reaktoro_solves = 0
        self.total_reaktoro_failures = 0
        self._last_iteration_count = None
        self._objective_component = None
        self._rho = None
        self._user_objective = None
        self._sense = 1
        self._variables = None
        self.previous_model_state = None
        self._current_radius = config.trust_radius
        # accepted truth (inputs, outputs, jacobian) per block, and trial truth
        self.accepted_truth = {}
        self.trial_truth = {}

    ####################################################################################
    # setup and clean up
    def _setup_objective(self):
        objs = list(
            self.model.component_data_objects(Objective, active=True, descend_into=True)
        )
        if len(objs) > 1:
            raise ValueError(
                "Trust region solve supports models with at most one active objective"
            )
        rho = self.config.elastic_penalty
        self._rho = Param(mutable=True, initialize=rho)
        self._rho_name = unique_component_name(
            self.model, "_reaktoro_trf_elastic_penalty"
        )
        self.model.add_component(self._rho_name, self._rho)
        penalty = sum(blk.enable_elastic() for blk in self.blocks)
        if objs:
            self._user_objective = objs[0]
            self._sense = 1 if self._user_objective.sense == minimize else -1
            obj_sf = iscale.get_scaling_factor(self._user_objective, default=1)
            expr = (
                self._sense * self._user_objective.expr + self._rho / obj_sf * penalty
            )
            self._user_objective.deactivate()
        else:
            obj_sf = None
            expr = self._rho * penalty
        self._obj_name = unique_component_name(self.model, "_reaktoro_trf_objective")
        self._objective_component = Objective(expr=expr, sense=minimize)
        self.model.add_component(self._obj_name, self._objective_component)
        if obj_sf is not None:
            iscale.set_scaling_factor(self._objective_component, obj_sf)

    def restore(self):
        """remove trust region components, reactivate user objective"""
        for blk in self.blocks:
            blk.disable_elastic()
        if self._objective_component is not None:
            self.model.del_component(self._objective_component)
            self._objective_component = None
        if self._rho is not None:
            self.model.del_component(self._rho)
            self._rho = None
        if self._user_objective is not None:
            self._user_objective.activate()

    def objective_value(self):
        if self._user_objective is None:
            return 0.0
        return self._sense * value(self._user_objective.expr)

    ####################################################################################
    # truth model
    def _evaluate_truth(self, blk, w):
        self.total_reaktoro_solves += 1
        return blk.evaluate_truth(w)

    def _theta(self, blk, truth_outputs):
        """surrogate error, relative to output magnitude (or 1/scaling factor for
        outputs near zero), so it does not depend on output units or on outputs
        without scaling factors (same idea as input step scale)"""
        sf = blk.output_scaling()
        sf[sf <= 0] = 1
        scale = np.maximum(np.abs(truth_outputs), 1 / sf)
        return float(np.sum(np.abs(blk.output_values() - truth_outputs) / scale))

    ####################################################################################
    # methods used by trust region main loop (same as TRFInterface)
    def initializeProblem(self):
        """evaluate reaktoro at current inputs, build surrogates and elastic trust region,
        returns (objective, feasibility)"""
        theta = 0.0
        for blk in self.blocks:
            w = blk.input_values()
            outputs, jacobian = self._evaluate_truth(blk, w)
            if not self.config.keep_surrogate_curvature:
                blk.reset_curvature()
            blk.set_expansion_point(w, outputs, jacobian)
            blk.set_constraint_scaling()
            self.accepted_truth[blk] = (w, outputs, jacobian)
            theta += self._theta(blk, outputs)
        self._setup_objective()
        self._variables = [
            v
            for v in self.model.component_data_objects(Var, descend_into=True)
            if not v.fixed
        ]
        return self.objective_value(), theta

    def updateDecisionVariableBounds(self, radius):
        for blk in self.blocks:
            blk.set_trust_radius(radius)
            blk.reset_elastic_values()

    def updateSurrogateModel(self, restoration=False):
        """move expansion point to accepted trial point, after a restoration step
        or a step larger than curvature_reset_step quadratic curvature is reset
        instead of updated, as a quasi-Newton update over such a step is not
        representative of local curvature"""
        for blk in self.blocks:
            w, outputs, jacobian = self.trial_truth[blk]
            reset = restoration or (
                blk.relative_step(w, blk.expansion_point)
                > self.config.curvature_reset_step
            )
            if reset:
                blk.reset_curvature()
            blk.set_expansion_point(w, outputs, jacobian, update_curvature=not reset)
            self.accepted_truth[blk] = (w, outputs, jacobian)

    def getCurrentModelState(self):
        return [v.value for v in self._variables]

    def rejectStep(self):
        for v, val in zip(self._variables, self.previous_model_state):
            v.set_value(val, skip_validation=True)
        if self.model_slack_enabled:
            self.set_model_slack(True)  # back to zero slack

    ####################################################################################
    # model slack (elastic surrogate equations)
    @property
    def model_slack_enabled(self):
        return any(blk.model_slack_enabled for blk in self.blocks)

    def set_model_slack(self, enable):
        for blk in self.blocks:
            blk.enable_model_slack(enable)

    def model_slack(self):
        return float(sum(blk.model_slack() for blk in self.blocks))

    def _inner_options(self):
        """warm start options for inner solves after first one, only if not set by user"""
        if self.inner_solves == 0:
            return {}
        user = getattr(self.inner_solver, "options", {}) or {}
        user_kwds = self.solve_kwds.get("options", {}) or {}
        return {
            k: v
            for k, v in self.config.warm_start_options.items()
            if k not in user and k not in user_kwds
        }

    def _inner_solve(self):
        kwds = dict(self.solve_kwds)
        options = dict(kwds.pop("options", {}) or {})
        options.update(self._inner_options())
        if options:
            kwds["options"] = options
        self._last_iteration_count = None
        if "intermediate_callback" not in kwds:

            def _count(
                nlp,
                alg_mod,
                iter_count,
                obj_value,
                inf_pr,
                inf_du,
                mu,
                d_norm,
                regularization_size,
                alpha_du,
                alpha_pr,
                ls_trials,
            ):
                self._last_iteration_count = iter_count
                return True

            kwds["intermediate_callback"] = _count
        self.inner_solves += 1
        t0 = time.time()
        try:
            results = self.inner_solver.solve(self.model, tee=self.tee, **kwds)
            ok = check_optimal_termination(results)
            status = str(results.solver.termination_condition)
        except Exception as e:  # solver errors count as failed inner solve
            _log.warning(f"Trust region inner solve failed: {type(e).__name__}: {e}")
            ok = False
            status = f"error: {type(e).__name__}"
        self.stats["inner_time"].append(time.time() - t0)
        self.stats["inner_iterations"].append(self._last_iteration_count)
        self.stats["inner_status"].append(status)
        return ok

    def solveModel(self):
        """solve trust region subproblem, returns (objective, step norm, feasibility),
        feasibility is inf if the inner solve or reaktoro failed"""
        self.previous_model_state = self.getCurrentModelState()
        radius = self._current_radius
        ok = self._inner_solve()
        if not ok and self.config.model_slack and not self.model_slack_enabled:
            # subproblem is infeasible for any inputs (e.g. surrogate can not reach a
            # fixed spec far from the expansion point because an input is at its
            # bound), retry with elastic slack on surrogate equations, the step is
            # then a restoration step towards the spec and reaktoro re-linearizes
            _log.info("Trust region inner solve failed, enabling model slack")
            self.rejectStep()
            self.set_model_slack(True)
            self.stats["model_slack_enabled"] += 1
            ok = self._inner_solve()
            # one outer iteration, merge stats of both inner solves
            it = self.stats["inner_iterations"]
            first, second = it[-2], it.pop()
            it[-1] = None if first is None and second is None else (first or 0) + (second or 0)
            t = self.stats["inner_time"]
            second_time = t.pop()
            t[-1] += second_time
            st = self.stats["inner_status"]
            second_status = st.pop()
            st[-1] = f"{second_status} (model slack)"
        failures = 0
        thetas = []
        truth_time = 0.0
        if not ok:
            self.consecutive_failed_inner_solves += 1
            self.stats["truth_time"].append(truth_time)
            self.stats["reaktoro_failures"].append(failures)
            self.stats["theta_per_block"].append([np.inf] * len(self.blocks))
            self.stats["elastic_gap"].append(np.nan)
            self.stats["elastic_penalty"].append(value(self._rho))
            self.stats["model_slack"].append(np.nan)
            return self.objective_value(), radius, np.inf
        self.consecutive_failed_inner_solves = 0
        step_norm = 0.0
        t0 = time.time()
        for blk in self.blocks:
            w = blk.input_values()
            step_norm = max(step_norm, blk.relative_step(w, blk.expansion_point))
            try:
                outputs, jacobian = self._evaluate_truth(blk, w)
                self.trial_truth[blk] = (w, outputs, jacobian)
                thetas.append(self._theta(blk, outputs))
            except _TRUTH_ERRORS as e:
                failures += 1
                self.total_reaktoro_failures += 1
                thetas.append(np.inf)
                _log.warning(
                    f"Reaktoro failed at trust region trial point ({type(e).__name__}), step will be rejected"
                )
                # put reaktoro back at accepted state, so next solve starts from good state
                w_ok = self.accepted_truth[blk][0]
                try:
                    self._evaluate_truth(blk, w_ok)
                except _TRUTH_ERRORS:
                    _log.warning("Reaktoro failed to re-solve at last accepted point")
        truth_time = time.time() - t0
        self.stats["truth_time"].append(truth_time)
        self.stats["reaktoro_failures"].append(failures)
        self.stats["theta_per_block"].append(thetas)
        gap = float(sum(blk.elastic_gap() for blk in self.blocks))
        self.stats["elastic_gap"].append(gap)
        self.stats["elastic_penalty"].append(value(self._rho))
        self.stats["model_slack"].append(self.model_slack())
        theta = float(np.sum(thetas))
        return self.objective_value(), step_norm, theta

    def after_accepted_step(self):
        """adaptive elastic penalty, model slack is fixed at zero again once it is not
        needed"""
        if (
            self.model_slack_enabled
            and self.model_slack() <= self.config.step_size_termination
        ):
            _log.info("Trust region model slack no longer needed, disabling")
            self.set_model_slack(False)
        if self.config.elastic_penalty_mode != ElasticPenaltyModes.adaptive:
            return
        gap = self.stats["elastic_gap"][-1]
        if gap is None or not np.isfinite(gap):
            return
        previous = self._last_accepted_gap
        self._last_accepted_gap = gap
        # forced moves (e.g. a fixed variable changed between solves) need elastic
        # steps for a few iterations, only increase penalty if gap is not decreasing
        if gap > self.config.step_size_termination and (
            previous is None or gap >= previous
        ):
            rho = min(
                value(self._rho) * self.config.elastic_penalty_increase,
                self.config.elastic_penalty_max,
            )
            self._rho.set_value(rho)


def trust_region_method(interface, config):
    """main trust region filter loop, copy of pyomo.contrib.trustregion.TRF.trust_region_method
    (see module doc for changes), returns (termination condition, message, logger)"""
    TRFLogger = IterationLogger()
    TRFilter = Filter()

    rebuildSM = False
    obj_val, feasibility = interface.initializeProblem()
    feasibility_k = feasibility
    obj_val_k = obj_val
    # changed: inf, so that at least one subproblem is solved
    step_norm_k = np.inf
    trust_radius = config.trust_radius
    subopt_flag = False
    restoration_gap = np.inf
    last_step_restoration = False
    iteration = 0
    termination = TerminationCondition.maxIterations
    message = "maximum iterations reached"

    TRFLogger.newIteration(iteration, feasibility_k, obj_val_k, trust_radius, 0)
    TRFLogger.logIteration()
    _print_header(config)
    _print_iteration(config, TRFLogger.iterrecord, interface, initial=True)
    while iteration < config.maximum_iterations:
        iteration += 1

        # Check termination conditions
        # changed: model slack has to be zero (surrogate equations exact)
        if (
            (feasibility_k <= config.feasibility_termination)
            and (step_norm_k <= config.step_size_termination)
            and (interface.model_slack() <= config.step_size_termination)
        ):
            termination = TerminationCondition.optimal
            message = "optimal solution found"
            break

        if (trust_radius <= config.minimum_radius) and (
            abs(feasibility_k - feasibility) < config.feasibility_termination
        ):
            if subopt_flag:
                _log.warning("Trust region: insufficient progress")
                termination = TerminationCondition.other
                message = "insufficient progress, trust radius at minimum"
                break
            else:
                subopt_flag = True
        else:
            subopt_flag = False

        if (
            interface.consecutive_failed_inner_solves
            > config.maximum_failed_inner_solves
        ):
            termination = TerminationCondition.other
            message = "maximum number of failed inner solves exceeded"
            break

        # changed: surrogate is rebuilt before trust region bounds (bounds depend on it)
        if rebuildSM:
            interface.updateSurrogateModel(restoration=last_step_restoration)
        interface._current_radius = trust_radius
        interface.updateDecisionVariableBounds(trust_radius)

        # Solve the Trust Region Subproblem (TRSP)
        obj_val_k, step_norm_k, feasibility_k = interface.solveModel()

        TRFLogger.newIteration(
            iteration, feasibility_k, obj_val_k, trust_radius, step_norm_k
        )

        # Check filter acceptance
        filterElement = FilterElement(obj_val_k, feasibility_k)
        # changed: trial point must also be acceptable to the current iterate
        # (Fletcher-Leyffer), otherwise an objective decrease starting from a feasible
        # point is accepted as an f-step no matter how much feasibility degrades
        # trial points that are nearly feasible (theta <= minimum_feasibility, where
        # f-steps are allowed) are not checked, otherwise steps near the optimum,
        # where objective decrease is small compared to theta, are all rejected
        acceptable_to_current = (
            (feasibility_k <= config.minimum_feasibility)
            or (feasibility_k <= (1 - config.param_filter_gamma_theta) * feasibility)
            or (obj_val_k <= obj_val - config.param_filter_gamma_f * feasibility_k)
        )
        TRFLogger.iterrecord.restoration = False
        gap_k = interface.stats["elastic_gap"][-1]
        if (
            not np.isfinite(feasibility_k)
            or not acceptable_to_current
            or not TRFilter.isAcceptable(filterElement, config.maximum_feasibility)
        ):
            # changed: restoration step. If the subproblem needed elastic slack, the
            # flowsheet forces reaktoro inputs outside of the trust region (e.g. a
            # fixed variable changed between solves). Rejecting such a step repeats the
            # same subproblem forever, instead accept it as long as the elastic gap
            # decreases, move the surrogate to the trial point and reset the filter.
            if (
                np.isfinite(feasibility_k)
                and gap_k is not None
                and np.isfinite(gap_k)
                and gap_k > config.step_size_termination
                and gap_k < restoration_gap
            ):
                restoration_gap = gap_k
                last_step_restoration = True
                TRFLogger.iterrecord.restoration = True
                TRFLogger.iterrecord.thetaStep = True
                TRFilter = Filter()
                TRFLogger.updateIteration(trustRadius=trust_radius)
                rebuildSM = True
                feasibility = feasibility_k
                obj_val = obj_val_k
                interface.after_accepted_step()
                TRFLogger.logIteration()
                _print_iteration(config, TRFLogger.iterrecord, interface)
                continue
            # Reject the step
            TRFLogger.iterrecord.rejected = True
            # changed: with elastic slack the step can be larger than the radius
            trust_radius = max(
                config.minimum_radius,
                min(step_norm_k, trust_radius) * config.radius_update_param_gamma_c,
            )
            rebuildSM = False
            interface.rejectStep()
            # changed: rejected step keeps last accepted values for termination test
            feasibility_k = feasibility
            obj_val_k = obj_val
            step_norm_k = np.inf
            TRFLogger.logIteration()
            _print_iteration(config, TRFLogger.iterrecord, interface)
            continue

        if (obj_val - obj_val_k) >= (
            config.switch_condition_kappa_theta
            * pow(feasibility, config.switch_condition_gamma_s)
        ) and (feasibility <= config.minimum_feasibility):
            TRFLogger.iterrecord.fStep = True
            trust_radius = min(
                max(step_norm_k * config.radius_update_param_gamma_e, trust_radius),
                config.maximum_radius,
            )
        else:
            TRFLogger.iterrecord.thetaStep = True
            filterElement = FilterElement(
                obj_val_k - config.param_filter_gamma_f * feasibility_k,
                (1 - config.param_filter_gamma_theta) * feasibility_k,
            )
            TRFilter.addToFilter(filterElement)
            rho_k = (
                feasibility - feasibility_k + config.feasibility_termination
            ) / max(feasibility, config.feasibility_termination)
            if (rho_k < config.ratio_test_param_eta_1) or (
                feasibility > config.minimum_feasibility
            ):
                trust_radius = max(
                    config.minimum_radius,
                    (config.radius_update_param_gamma_c * step_norm_k),
                )
            elif rho_k >= config.ratio_test_param_eta_2:
                trust_radius = min(
                    config.maximum_radius,
                    max(
                        trust_radius, (config.radius_update_param_gamma_e * step_norm_k)
                    ),
                )

        TRFLogger.updateIteration(trustRadius=trust_radius)
        rebuildSM = True
        last_step_restoration = False
        feasibility = feasibility_k
        obj_val = obj_val_k
        interface.after_accepted_step()
        TRFLogger.logIteration()
        _print_iteration(config, TRFLogger.iterrecord, interface)

    if termination == TerminationCondition.maxIterations:
        _log.warning(
            "Trust region: maximum iterations reached: {}".format(
                config.maximum_iterations
            )
        )
    if config.verbose:
        print(f"Trust region: {message}")
    return termination, message, TRFLogger


def _print_header(config):
    if config.verbose:
        print(
            f"{'iter':>5} {'theta':>11} {'objective':>13} {'radius':>9} {'step':>9} "
            f"{'gap':>9} {'penalty':>8} {'inner':>6} {'type':>5}"
        )


def _print_iteration(config, record, interface, initial=False):
    if not config.verbose:
        return
    if initial:
        gap, rho, inner = "-", "-", "-"
    else:
        gap = interface.stats["elastic_gap"][-1]
        gap = "nan" if gap is None or not np.isfinite(gap) else f"{gap:.2e}"
        rho = f"{interface.stats['elastic_penalty'][-1]:.1e}"
        inner = interface.stats["inner_iterations"][-1]
        inner = "-" if inner is None else str(inner)
    if record.rejected:
        kind = "rej"
    elif getattr(record, "restoration", False):
        kind = "rest"
    elif record.fStep:
        kind = "f"
    elif record.thetaStep:
        kind = "theta"
    else:
        kind = ""
    print(
        f"{record.iteration:>5} {record.feasibility:>11.4e} {record.objectiveValue:>13.6e} "
        f"{record.trustRadius:>9.2e} {record.stepNorm:>9.2e} {gap:>9} {rho:>8} {inner:>6} {kind:>5}"
    )


class ReaktoroTrustRegionSolver:
    """Solves models with ReaktoroSurrogateBlocks with the trust region filter method,
    inner_solver (e.g. cyipopt-watertap) solves the subproblems"""

    def __init__(self, inner_solver, trust_region_options=None):
        self.inner_solver = inner_solver
        self.trust_region_options = trust_region_options
        self.trust_region_stats = None

    def _get_config(self, blocks):
        if self.trust_region_options is not None:
            config = TrustRegionOptions().get_dict()()
            config.set_value(self.trust_region_options)
            return config
        configs = [blk.trust_region_options for blk in blocks]
        config = configs[0]
        for other in configs[1:]:
            if other.value() != config.value():
                _log.warning(
                    "Reaktoro blocks have different trust region options, using options of first block"
                )
                break
        return config

    def solve(self, model, tee=False, **kwds):
        blocks = find_surrogate_blocks(model)
        if not blocks:
            return self.inner_solver.solve(model, tee=tee, **kwds)
        config = self._get_config(blocks)
        interface = ReaktoroTRInterface(
            model, blocks, self.inner_solver, config, tee=tee, solve_kwds=kwds
        )
        t0 = time.time()
        logger = None
        try:
            termination, message, logger = trust_region_method(interface, config)
        finally:
            interface.restore()
        stats = self._build_stats(
            interface, logger, termination, message, time.time() - t0
        )
        self.trust_region_stats = stats
        results = SolverResults()
        results.solver.name = "reaktoro-trust-region"
        results.solver.termination_condition = termination
        results.solver.status = (
            SolverStatus.ok
            if termination == TerminationCondition.optimal
            else SolverStatus.warning
        )
        results.solver.message = message
        results.solver.wallclock_time = stats["wall_time"]
        results.solver.trust_region = stats
        return results

    @staticmethod
    def _build_stats(interface, logger, termination, message, wall_time):
        iterations = []
        records = [] if logger is None else logger.iterations
        s = interface.stats
        for k, rec in enumerate(records):
            row = {
                "iteration": rec.iteration,
                "feasibility": rec.feasibility,
                "objective": rec.objectiveValue,
                "trust_radius": rec.trustRadius,
                "step_norm": rec.stepNorm,
                "f_step": rec.fStep,
                "theta_step": rec.thetaStep,
                "rejected": rec.rejected,
                "restoration": getattr(rec, "restoration", False),
            }
            if k > 0 and k - 1 < len(s["inner_status"]):
                j = k - 1
                row.update(
                    {
                        "feasibility_per_block": s["theta_per_block"][j],
                        "elastic_gap": s["elastic_gap"][j],
                        "elastic_penalty": s["elastic_penalty"][j],
                        "model_slack": s["model_slack"][j],
                        "inner_iterations": s["inner_iterations"][j],
                        "inner_status": s["inner_status"][j],
                        "inner_time": s["inner_time"][j],
                        "truth_time": s["truth_time"][j],
                        "reaktoro_failures": s["reaktoro_failures"][j],
                    }
                )
            iterations.append(row)
        inner_iters = [i for i in s["inner_iterations"] if i is not None]
        return {
            "termination_condition": str(termination),
            "message": message,
            "outer_iterations": max(len(records) - 1, 0),
            "accepted_steps": sum(1 for r in records[1:] if not r.rejected),
            "rejected_steps": sum(1 for r in records[1:] if r.rejected),
            "inner_solves": interface.inner_solves,
            "inner_iterations_total": int(np.sum(inner_iters)) if inner_iters else None,
            "inner_iterations": list(s["inner_iterations"]),
            "reaktoro_solves": interface.total_reaktoro_solves,
            "reaktoro_failures": interface.total_reaktoro_failures,
            "model_slack_enabled": s["model_slack_enabled"],
            "wall_time": wall_time,
            "inner_time": float(np.sum(s["inner_time"])),
            "truth_time": float(np.sum(s["truth_time"])),
            "number_of_surrogate_blocks": len(interface.blocks),
            "iterations": iterations,
        }
