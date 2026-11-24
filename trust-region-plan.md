# Trust-region backend for ReaktoroBlock: plan

Status: phases 0–2 implemented on branch `try_trust_region` (see section 10); validated on reaktoro-pse examples only.
Date: 2026-11-23

## 1. Goal

Add a second solve backend for `ReaktoroBlock` and `ReaktoroBlockManager`: a trust-region method (TRF) with algebraic surrogates instead of a gray box. Switching should take **one setting and nothing else**:

```python
ReaktoroBlockManager(backend="trust_region", hessian_options=...)  # managed blocks
ReaktoroBlock(..., backend="trust_region")                          # blocks without a manager
rkt_options["backend"] = "trust_region"                             # same, via REW reaktoro_options
```

Nothing else changes: inputs, outputs, scaling, speciation blocks, the parallel manager, Hessian options, solver helpers, results handling. The REW reference case is `flowsheets/softening_acid_ro`, which does square solves (initialization, fixed recovery) and optimization (RO design and chemical dosing at fixed recovery).

## 2. Summary

| Question | Answer |
|---|---|
| Can we use `pyomo.contrib.trustregion` (TRF) as is? | **No.** Six problems block it for reaktoro-pse/REW models (section 4). The worst: it clones the model, counts fixed variables as degrees of freedom, crashes on square solves, and cannot pass Ipopt options. |
| Can we use the same algorithm? | **Yes, mostly reused (option B, chosen).** All the blockers are in Pyomo's model plumbing (`TRFInterface`), not in the algorithm. We import Pyomo's `Filter`, `IterationLogger` and `_trf_config()` and replace only `TRFInterface` with our own adapter (about 150 lines). The one duplication is a copy of the main loop `trust_region_method` (about 110 lines), because it creates `TRFInterface` itself. Follow-up option C: an upstream Pyomo PR that makes the interface class injectable, so the copy can be deleted (section 6.4). |
| How much of reaktoro-pse changes? | Little. `ReaktoroSolver.solve_reaktoro_block(params) -> (jacobian, outputs)` is already the "truth model" interface the algorithm needs. This holds for the single, coupled and parallel-manager paths. Only the gray-box block and the Ipopt solve call get replaced. |
| Does `softening_acid_ro` work with one switch? | **Yes**, with three requirements on reaktoro-pse (section 7): automatic solver selection inside `get_cyipopt_watertap_solver`; the same degrees of freedom outside a solve as with the gray box; and the surrogate type taken from the existing `hessian_options`. Trust-region statistics are available from the solver and the results object (section 6.6). |
| What does it buy? | Reaktoro is only evaluated at subproblem solutions, and failures become rejected steps instead of Ipopt evaluation errors. Steps in dosing, pH and concentration are bounded directly. The inner Ipopt solves see an exact Hessian (of the surrogate). No gray-box Hessian callback, no zero-dual reset. |
| What doesn't it fix? | Ill-posed chemistry coupling. The surrogate has the same Jacobian as the gray box, so a numerically singular Jacobian (as diagnosed in the Kenya flowsheet) remains singular. |

Recommendation: build the phase 0 prototype on reaktoro-pse examples and compare it with the gray box before integrating.

## 3. How the algorithm works

For each Reaktoro block, the gray box today enforces `y = d(w)`. Here `w` are the inputs (`reaktoro_model.inputs`), `y` the outputs (`reaktoro_model.outputs`), and `d` the Reaktoro equilibrium. Ipopt calls `d`, `dd/dw` and a quasi-Newton Hessian on every iteration and line-search trial.

The trust-region method replaces `y = d(w)` with an algebraic surrogate. That surrogate is exact in value and gradient at the current point `w_k`:

```
y = r_k(w) = d(w_k) + J_k (w - w_k)              (linear, "zero basis")
y = d(w_k) + J_k (w - w_k) + 0.5 (w - w_k)' B_k (w - w_k)   (quadratic, B_k from quasi-Newton)
```

Outer iteration k:

1. Evaluate the truth `d(w_k), J_k` for all blocks (one Reaktoro solve per block; can run in parallel).
2. Build the subproblem: the full flowsheet with surrogates in place of the Reaktoro blocks, plus a trust region on the chemistry inputs (section 6.3). Solve it with Ipopt. The model is fully algebraic and needs no gray box.
3. At the subproblem solution `w*`, evaluate the truth again. Feasibility measure: `θ = ||D_c (y* - d(w*))||`.
4. Use a filter on (objective, θ) and a ratio test to accept or reject the step and update `Δ`.
   - If Reaktoro fails at `w*`, reject the step and shrink `Δ`. This is the built-in recovery path.
5. Stop when `θ <= tol` and the step is small.

Convergence: with a surrogate that matches value and gradient (first-order consistent), the method converges to a KKT point of the original problem. With the linear surrogate the subproblem has no chemistry curvature, so convergence near an optimum is linear; the quadratic surrogate is meant to fix that. Pure feasibility (square) solves behave like a damped Newton method and converge fast.

## 4. Why the stock Pyomo TRF does not fit

These were checked in the installed source (`pyomo/contrib/trustregion`, TRF version 0.2.0, `recycled_softening_analysis` env):

| # | Limitation (where) | Effect on reaktoro-pse/REW |
|---|---|---|
| 1 | Truth models must be `ExternalFunction` expressions inside constraints (`EFReplacement.exitNode`). Value and gradient come from `value()` and `differentiate()` on each one separately (`updateSurrogateModel`). | We use `ExternalGreyBoxBlock` with many outputs. We would need one `ExternalFunction` per output with a gradient callback and an input cache, so we don't run one Reaktoro solve per output. |
| 2 | `TRFInterface.__init__` calls `model.clone()` and the solve returns the **clone** (`trust_region_method` returns `interface.model`). | Deep-copying Reaktoro objects and the parallel manager (worker processes) is unsupported or very expensive. Results never reach the user's model, and REW reads results from `m`. |
| 3 | The degrees-of-freedom check counts **every** variable met in the expression walk, including fixed ones (`beforeChild` adds any `is_variable_type()`). It then requires `DOF == len(decision_variables)` (`replaceExternalFunctionsWithVariables`). | REW flowsheets have many fixed variables, so this raises `ValueError`. |
| 4 | The trust region is applied only to user-given decision variables, in **unscaled** units (`updateDecisionVariableBounds`; its docstring says "we assume that the users have correctly scaled their variables"). The step norm is `max()` over those variables (`calculateStepSizeInfNorm`). | Square solves (`softening_acid_ro` initialization and fixed-recovery solves) have no decision variables, so `max([])` crashes. Even with degrees of freedom, the region doesn't limit what breaks Reaktoro (pH, concentrations at the RO interface, dosing near equivalence points). |
| 5 | The subsolver is `SolverFactory(config.solver)` with no options, and `solve(model, keepfiles, tee)` only. | We can't pass `cyipopt-watertap` settings (user scaling, linear solver, `mumps_pivtol`, tolerances, `output_file`). |
| 6 | Exactly one active objective is required. Any subproblem failure raises `ArithmeticError`. | No recovery or fallback. (`softening_acid_ro` always has `lcow_objective` active, so the objective requirement itself is fine there.) |

All six are in `TRFInterface` (`interface.py`), the model-plumbing layer. The algorithm layer (`trust_region_method` in `TRF.py`, `Filter` in `filter.py`, `IterationLogger` in `util.py`, `_trf_config()`) has no blockers. It calls only these interface methods:
`initializeProblem`, `updateDecisionVariableBounds`, `updateSurrogateModel`, `solveModel`, `rejectStep`, and the attribute `model`.

So we keep the algorithm and replace the interface (option B, section 6.4). Two things in `trust_region_method` stop us from calling it directly: it creates `TRFInterface(...)` itself, and on success it calls `interface.model.display()` (which would print the whole flowsheet) and returns no termination status. That is why we copy the loop.

## 5. What already fits in reaktoro-pse

- **Truth interface.** `ReaktoroSolver.solve_reaktoro_block(params)` ([core/reaktoro_solver.py](src/reaktoro_pse/core/reaktoro_solver.py)), `ReaktoroCoupledSolver.solve_reaktoro_block` ([core/reaktoro_coupled_solver.py](src/reaktoro_pse/core/reaktoro_coupled_solver.py)) and `AggregateSolverState.solve_reaktoro_block` ([parallel_tools/reaktoro_block_manager.py](src/reaktoro_pse/parallel_tools/reaktoro_block_manager.py)) all return `(jacobian, outputs)` for a dict of input values. The parallel path already solves every block at once, which is what a TR outer iteration needs.
- **Pyomo-side structure.** `ReaktoroBlockBuilder` ([core/reaktoro_block_builder.py](src/reaktoro_pse/core/reaktoro_block_builder.py)) creates `block.reaktoro_model` with `.inputs[...]` and `.outputs[...]` variables. It adds `input_constraints` that link user variables to `reaktoro_model.inputs` and `output_constraints` that link `reaktoro_model.outputs` to user variables. Scaling is set on those variables (`initialize_input_variables_and_constraints`, `set_output_vars_and_scale`). If the TR backend builds a `reaktoro_model` **with the same component names**, all of this works unchanged.
- **Output scaling.** The gray box's output-constraint scaling comes from `get_jacobian_scaling()`, and the same values can scale the surrogate constraints.
- **Output floors.** `ReaktoroGrayBox.apply_output_floors` (species amount offset) is a pure function of `(jacobian, outputs)`. It can be moved to a shared helper and applied to truth values before they go into the surrogate.
- **Hessian machinery.** `HessianApproximation` ([core/util_classes/hessian_functions.py](src/reaktoro_pse/core/util_classes/hessian_functions.py)) supplies the per-output matrices `B_o` for the quadratic surrogate (any of its types: SR1, LBFGS, BFGS variants), updated once per *accepted outer step* instead of per Ipopt iteration. That gives much cleaner (s, y) pairs.

## 6. What needs to be added or changed

### 6.1 Configuration: one switch
- `backend`: `"graybox"` (default) or `"trust_region"`. Resolution order:
  1. `ReaktoroBlockManager(backend=...)`, for blocks registered with a manager. A block-level setting that conflicts raises an error.
  2. `ReaktoroBlock(backend=...)` / `rkt_options["backend"]`, for blocks without a manager.
  - There is no global default (decided): the backend is an option of `ReaktoroBlock` or `ReaktoroBlockManager` only.
  - Note for REW: the feed unit's charge-balance sub-model forces `reaktoro_block_manager=None` and copies the user's `reaktoro_options`. It therefore picks up the block-level option, never the manager's backend. Either backend is correct for that small square solve. To run everything on one backend, set `backend` in `rkt_options` as well as on the manager.
- **Surrogate type comes from the existing `hessian_options`, so no new option is needed to switch:**

  | `hessian_options.hessian_type` | Surrogate |
  |---|---|
  | `SR1` (reaktoro-pse default), `SR1_lagrangian`, `PSB` | quadratic, `B_o` from that update rule |
  | `LBFGS` (what `softening_acid_ro` uses), `BFGS*` | quadratic, `B_o` from that rule (positive definite) |
  | `GaussNewton` | quadratic, `B_o = J_o' J_o` per output |
  | `ZeroHessian` (also set by `softening_acid_ro` for `rkt_hessian_type="limited-memory"`) | linear |

  The initialization type, caps and input scaling from `hessian_options` all apply. An explicit `trust_region_options.surrogate` can override the choice.
- New `trust_region_options` ConfigDict (all optional, defaults chosen for REW flowsheets):
  - Radius settings: initial, minimum and maximum radius, and update factors.
  - Filter and ratio-test parameters, feasibility tolerance, and maximum outer iterations.
  - Trust-region shape (section 6.3).
  - Elastic penalty: `elastic_penalty_mode` = `"fixed"` or `"adaptive"`, with `elastic_penalty` (the fixed value, or the starting value for adaptive) and `elastic_penalty_increase` (factor applied when the gap `e` stays non-zero over consecutive accepted steps). Both modes are implemented; the default is picked in phase 0.
  - `keep_surrogate_curvature` (default `False`): `B_o` and the quasi-Newton history are reset at the start of every solve, like the gray box today. With `True`, they carry over as a warm start into the next solve, e.g. along a recovery sweep. The first-order terms `y0`, `J` are always refreshed from a new truth evaluation.
  - It lives on the manager and on the block, alongside `hessian_options`.

### 6.2 Surrogate block (replaces `ExternalGreyBoxBlock`)
New `core/reaktoro_surrogate.py`, used by `ReaktoroBlockBuilder.build_reaktoro_block` when the backend is `trust_region`:
- `reaktoro_model = Block()` with `inputs = Var(input_names)` and `outputs = Var(output_names)`, the same index sets and names as the gray box.
- Mutable `Param`s: `y0[o]`, `w0[i]`, `J[o, i]` (pattern from the Hessian support, fixed at build time, so the Ipopt structure stays constant) and `B[o, i, j]` (quadratic only, lower triangle over the output's support).
- `surrogate_eq[o]: outputs[o] == y0[o] + J[o,:] s + 0.5 s' B_o s`, with `s = inputs - w0`. It is scaled with the gray box's output-constraint scaling (`get_jacobian_scaling()`). For the linear surrogate, `B` is omitted.
  - `B_o` comes from the existing `HessianApproximation`, updated once per accepted outer step from consecutive truth points `(w_k, J_k)`. We need a small accessor that returns the per-output element matrices (`el.matrix` already holds them); the dual-weighted assembly isn't used.
  - Because the surrogate is an algebraic Pyomo expression, Ipopt gets its **exact** Hessian through automatic differentiation, combined with Ipopt's own multipliers.
  - SR1 can be indefinite. That is acceptable inside a trust region, which is the setting SR1 is designed for.
  - First outer iteration: no step pair exists yet, so `B_o` starts from the configured initialization (GaussNewton by default) or zero.
- `update_surrogate(jacobian, outputs, w)`: sets the Params after applying the output floors.
- **Equation count must match the gray box.** Outside a TR solve the block must have the same number of variables and equality constraints as the gray-box version: one `surrogate_eq` per output. Otherwise `degrees_of_freedom(m)` changes and `softening_acid_ro.fix_and_scale` (`assert degrees_of_freedom(m) == 0`) and `report_global_state` break. The elastic variables (section 6.3) therefore exist only during a TR solve: they are fixed at 0 with their constraints deactivated outside it. A test checks that `degrees_of_freedom` is the same on both backends for every example.
- The block's `activate`/`deactivate` (currently iterating `ExternalGreyBoxModel`, [reaktoro_block.py:1354-1372](src/reaktoro_pse/reaktoro_block.py#L1354-L1372)) must also handle `surrogate_eq`.
- **Manager: one aggregate surrogate (decided).** This is the least code and gives the most parallelism:
  - `ReaktoroBlockManager.build_reaktoro_blocks` already builds one aggregate gray box from `AggregateSolverState`, and `PseudoGrayBox` already maps each block into the aggregate input and output windows. The TR backend builds one aggregate surrogate in the same place, with the same windows, and the per-block linking code is unchanged.
  - One truth evaluation per outer iteration = one `AggregateSolverState.solve_reaktoro_block(params)` call. That call already sends every block to the workers in parallel (`parallel_solver`), so the adapter needs no parallel code of its own.
  - `J` and `B_o` are block-diagonal through the existing windows and `get_hessian_support()`.
  - One surrogate per registered block would need the adapter to collect results from N solvers and run them in parallel itself, which duplicates `parallel_solver`.
- Blocks without a manager: one surrogate per `ReaktoroBlock`, truth evaluated serially. The adapter treats both cases the same way, as a list of (surrogate block, `solve_reaktoro_block`).

### 6.3 Trust region: on chemistry inputs, elastic, mixed relative/absolute
- **Variables.** The trust region applies to the chemistry inputs only (the variables Reaktoro sees), not to the flowsheet decision variables. Units, costing, recycles and pumps stay exact algebra and need no trust region. The size of the region depends only on the number of Reaktoro inputs, not on the flowsheet size.
  - For `softening_acid_ro`, the decision variables (reagent doses, RO pressure and area) act on the chemistry through the flowsheet. A box on the decisions would not bound what Reaktoro sees, for example the RO interface concentrations at high recovery, or the pH near an equivalence point as the dose changes.
  - Square solves (initialization, fixed recovery before `set_optimization`) have no decisions at all.
- **Elastic, so the subproblem is never infeasible because of the trust region.** In a square or tightly specified flowsheet the equations fix the chemistry inputs, so a hard box could make the subproblem infeasible. Instead:

  ```
  y = r_k(w_k + s)            surrogate evaluated at the trusted point
  |s| <= Δ (per input, see shape below)
  w - (w_k + s) = e,  e free  elastic gap, penalized
  minimize  f + ρ·||D_c e||_1  (f = the user's objective, e.g. LCOW)
  ```

  - The surrogate is never evaluated outside the region. `e` measures how far the flowsheet wants to move past it, and is driven to zero by the outer loop (it enters θ).
  - In square solves this acts like Levenberg–Marquardt.
  - Cost: 2·n_inputs variables and n_inputs constraints, added only during the TR solve (section 6.2).
- **Shape (mixed).** Amounts in REW flowsheets span many orders of magnitude, so the box is:
  - Relative for species and element amounts: `|s_i| <= Δ·max(|w_k,i|, floor_i)`. With Δ = 0.5, each amount may change by 50% per step.
  - Absolute for pH, pE, temperature and pressure: `|s_i| <= Δ_abs/sf_i`.
  - Both overridable per input through `trust_region_options`.
- **One global radius first.** Per-block mismatch `θ_b` is logged every outer iteration. If one block (e.g. the RO scaling block) is clearly holding back all the others, add per-block radii in phase 3, adjusted by each block's own mismatch.

### 6.4 Trust-region driver (option B: reuse Pyomo's algorithm, replace its interface)

What is reused, what is new, and what is copied:

| Piece | Source | Lines (approx.) |
|---|---|---|
| `Filter`, `FilterElement` | imported from `pyomo.contrib.trustregion.filter` | 0 |
| `IterationLogger` | imported from `pyomo.contrib.trustregion.util` | 0 |
| Config (radius, filter, ratio-test, termination options) | `_trf_config()` imported, then our own options declared on top (section 6.1) | about 20 |
| Main loop | **copied** from `trust_region_method`, with three changes: it accepts the interface object as an argument, it drops `interface.model.display()`, and it returns a status (optimal / feasible-not-optimal / max iterations / failed) | about 110 (the only duplication) |
| `ReaktoroTRInterface` adapter | new; implements the same methods as `TRFInterface` for our surrogate blocks (table below) | about 150, plus about 40 for the elastic setup |

Adapter methods, in the same shape as `TRFInterface`:

| Method | Pyomo version does | Ours does |
|---|---|---|
| `__init__` | clones the model | keeps the model; finds active surrogate blocks inside the block being solved; stores the inner solver |
| `initializeProblem` | replaces `ExternalFunction`s with holder variables, builds `basis_constraint`/`sm_constraint_basis`, solves once | truth evaluation at the current inputs, `update_surrogate`, activates the elastic part, returns (objective, θ) without a solve |
| `updateSurrogateModel` | `value()` and `differentiate()` on each `ExternalFunction` | one `solve_reaktoro_block` per block (or one manager call), then output floors, then Params and the `B_o` update |
| `calculateFeasibility` | sum of `|y - d(w)|` | scaled surrogate mismatch plus the elastic gap `e`; a Reaktoro failure gives +inf |
| `updateDecisionVariableBounds(radius)` | bounds on decision variables, unscaled | bounds on the trust-region step `s` (section 6.3) |
| `solveModel` | `SolverFactory(name).solve`; raises `ArithmeticError` on failure | inner solver with the user's options plus the driver's warm-start defaults; on failure, restores the state and returns θ = +inf, so the filter rejects the step and the radius shrinks |
| `rejectStep`, `getCurrentModelState` | over `all_variables` from the expression walk | over all unfixed variables of the solved block |
| `calculateStepSizeInfNorm` | `max()` over decision variables (crashes for square problems) | infinity norm of the scaled `s` |
| `restore` (new) | – | deactivates the elastic part, fixes `e`/`s` at 0; called in `finally` |

Follow-up option C (after B works): a Pyomo PR that lets `trust_region_method` take an interface class or instance, makes the final `display()` optional and returns a status. Then the copied loop is removed and we call Pyomo's loop with our adapter.

Solver: new `core/util_classes/trust_region_solver.py` with `ReaktoroTrustRegionSolver`:
1. Collect the active surrogate blocks inside the block being solved (so sub-block solves, like the REW feed reconciliation, work).
2. Truth at the current inputs; update the surrogates. A Reaktoro failure at the start point raises, as with the gray box.
3. Run the loop. Each subproblem is solved with the user's solver and options (e.g. `cyipopt-watertap` with `linear_solver`, `mumps_pivtol`, tolerances).
4. Return a `SolverResults` with the termination condition (optimal, infeasible or maxIterations), so `assert_optimal_termination(result)` works. Statistics are attached as described in section 6.6.
5. Logging: a one-line outer-iteration table (k, θ, θ per block, objective, Δ, inner iterations and status, Reaktoro failures). `tee` passes through to the inner Ipopt solves.

Inner solves:
- **Warm starts (decided: driver defaults, user selections always win).** Every outer iteration is a new Ipopt solve, so Ipopt's default `bound_push`/`bound_frac` (1e-2, scaled) would move near-bound variables each time.
  - For inner solves after the first, the driver adds `warm_start_init_point=yes` (primal and dual values from the previous accepted subproblem), small `warm_start_bound_push`/`warm_start_bound_frac`/`warm_start_mult_bound_push`, and `mu_init` near the last barrier value.
  - Each default is applied only if the user didn't set that option (`setdefault` on a copy of the user's options); the user's options are never changed. The first inner solve uses the user's options unchanged.
  - Phase 0 check: how `cyipopt-watertap` takes initial multipliers (Pyomo `ipopt_zL_in`/`ipopt_zU_in`/`dual` suffixes vs PyNumero `init_duals`). Fallback: primal warm start with a small bound push only.
- A subproblem failure (restoration failure, max iterations) counts as a rejected step with a smaller `Δ`, up to N times, then the solver stops with a non-optimal status.
- User choices that would be odd but are respected: `hessian_approximation=limited-memory` (from `softening_acid_ro`'s `rkt_hessian_type="limited-memory"`). The inner Ipopt then uses L-BFGS instead of the surrogate's exact Hessian. The plan leaves it as the user's choice and logs a note.

### 6.5 Solver selection (automatic)
- `get_cyipopt_watertap_solver` ([cyipopt_solver.py](src/reaktoro_pse/core/util_classes/cyipopt_solver.py)) returns a thin wrapper around the current `cyipopt-watertap` solver.
  - Its `.solve(model, ...)` checks for active surrogate blocks in `model`: if there are any, it runs `ReaktoroTrustRegionSolver` with the wrapped solver as the inner solver; otherwise it solves directly, as now.
  - The wrapper passes `.options`, `.available()` and `.solve(m, tee=...)` straight through, because callers set `solver.options[...]` after creating it (`softening_acid_ro.solve_model` sets `output_file`).
- `get_reaktoro_solver(...)`: same thing, under a clearer name for new code.
- **Guard.** A plain Ipopt solve of a model containing surrogates would treat the frozen surrogate as truth and give a wrong "optimal" answer. Surrogate blocks register with a module-level registry, and the wrapper is the only path that activates the elastic part and updates the surrogates. For solves outside the wrapper (direct `watertap_solvers.get_solver()` calls, section 7.2), the surrogate block checks at the start of the next TR solve whether its inputs moved since the last truth evaluation. If they did, it logs a warning that the previous solve used a stale surrogate.

### 6.6 Trust-region solve statistics
The stock Pyomo TRF records per-iteration data in an `IterationLogger` (`IterationRecord`: iteration, feasibility θ, objective, trust radius, step norm, and whether the step was an f-step, a θ-step or rejected). That logger is local to `trust_region_method` and is never returned, because the function returns only `interface.model`. So these statistics are **not available** from the default TRF solver. Our copied loop keeps the logger and returns it, along with what the adapter records:

| Level | Statistics |
|---|---|
| Per outer iteration | Pyomo's `IterationRecord` fields (θ, objective, radius, step norm, step type, rejected), plus: θ per block, elastic gap `‖e‖`, elastic penalty `ρ`, inner Ipopt iterations, status and time, Reaktoro solves and failures per block, truth-evaluation time |
| Per solve (totals) | outer iterations, accepted/rejected steps, inner Ipopt iterations (total and per inner solve), Reaktoro solves and failures (total and per block), wall time split into inner solves / truth evaluations / overhead, termination reason |

Access:
- `results.solver.trust_region`: a dict with the totals and the per-iteration list (list of dicts, easy to turn into a DataFrame).
- `solver.trust_region_stats`: the same object for the last solve, for callers that don't keep the results.
- `reaktoro_pse.core.util_classes.trust_region_solver.get_trust_region_stats(results)`: a helper that returns it, or `None` for a gray-box solve.

The user's Ipopt options, including `output_file`, are passed unchanged to every inner solve, with no special handling. With `output_file` set, each inner solve overwrites it, so the file holds the last inner solve's log. The trust-region statistics above are the supported way to get totals.

### 6.7 What the TR backend does not use
- The gray-box Hessian callback and the per-solve zero-dual reset: not needed, because Ipopt differentiates the surrogate itself.
- `CyIpoptEvaluationError` handling and `max_reaktoro_failed_solves` change meaning: a Reaktoro failure becomes a rejected step, never an Ipopt evaluation error.

## 7. `softening_acid_ro` (REW) with the TR backend

### 7.1 Walk through the flowsheet's steps

| Step in [softening_acid_ro.py](../reaktoro_enabled_watertap/src/reaktoro_enabled_watertap/flowsheets/softening_acid_ro/softening_acid_ro.py) | With the gray box today | With the TR backend | Change in REW |
|---|---|---|---|
| `build_model(... rkt_hessian_type="LBFGS")` creates `ReaktoroBlockManager(hessian_options=...)` or per-block `rkt_options` | aggregate gray box | aggregate surrogate; LBFGS `hessian_options` → quadratic surrogate with LBFGS `B_o` | one argument (`backend=` on the manager, and in `rkt_options` for the feed sub-model) |
| Units: `PrecipitationUnit` (softening), `ChemicalAdditionUnit` (acid), `MultiCompROUnit` (scaling block), `hpro_unit` (Seawater) | Reaktoro blocks in the manager | same blocks, same inputs and outputs | none |
| `MultiCompFeed(reconcile_using_reaktoro=True)`: charge-balance sub-model, own `ReaktoroBlock` (manager forced off), solved with `get_cyipopt_watertap_solver()` | small gray-box solve | small TR solve on the sub-model (wrapper finds its surrogate), if the backend reaches the block (section 6.1) | none |
| `m.reaktoro_manager.build_reaktoro_blocks()` | builds the aggregate gray box | builds the aggregate surrogate | none |
| `fix_and_scale` → `assert degrees_of_freedom(m) == 0` | passes | passes, provided the equation count matches (section 6.2) | none |
| Output bounds: alkalinity ≥ 10 ppm, softening and RO feed pH 6–11 | bounds on gray-box output variables | the same bounds on the surrogate output variables; the truth may violate them slightly at a trial point, which shows up in θ and is corrected by the outer loop | none |
| `initialize`: `unit.initialize()` for every unit, then `costing.initialize()` | Reaktoro solves and IDAES initialization; no Ipopt solve with a gray box except in the feed | unchanged; these Reaktoro solves also set the surrogates' first points | none |
| `initialize`: first `solve_model` (square, LCOW objective active, DOF 0) | one Ipopt solve | TR solve, square: elastic region, damped Newton; typically a few outer iterations | none |
| `set_optimization` → units unfix dosing, pH, RO pressure and area | – | – | none |
| `solve_model` at recovery 0.5, then 0.6/0.7/0.8 in `main` (optimization: RO design and dosing for minimum LCOW) | one Ipopt solve per recovery; Hessian reset each solve | one TR solve per recovery; fresh truth evaluation at the start, `B_o` reset each solve (option to keep it across solves) | none |
| Active RO scaling constraints (Calcite, Gypsum `scaling_tendency ≤ 1`) at the optimum | gray-box outputs | surrogate outputs; this is where the quadratic surrogate matters most (curvature of ST in dose and recovery) | none |
| `solve_model`: `get_cyipopt_watertap_solver(linear_solver, max_iter=1000, limited_memory, pivtol, pivtolmax)` | cyipopt | wrapper → TR with these as the inner options; `max_iter` applies per inner solve | none |
| `solve_model`: `solver.options["output_file"]` + `ipopt_perf_utils.get_ipopt_performance_data` | iterations and final errors of the one solve | `output_file` passed through unchanged, so it holds the last inner solve's log; trust-region totals come from the statistics (section 6.6) | none required; optional to read the statistics |
| `print_variables_close_to_bounds`, `print_constraints_close_to_bounds` | – | elastic part deactivated and fixed before returning, so the reports see the same model | none |
| `assert_optimal_termination(result)` | Ipopt status | TR status (optimal only if θ ≤ tol and the last inner solve was optimal) | none |
| `m.reaktoro_manager.terminate_workers()` | – | unchanged | none |
| Sweeps (`analysis_scripts/softening_acid_ro/data_generation/*_sweep.py`, `loopTool`) | call `solve_model` | unchanged | none |

### 7.2 Direct `watertap_solvers.get_solver()` calls in REW
These bypass `get_cyipopt_watertap_solver`: `multi_comp_feed_unit.py:310` and `:504`, `multi_comp_mvc.py:1205`, `utils/translator_utils.py:58`.
- In `softening_acid_ro`, the feed calls (`initialize_unit`, and the solve before reconciliation) solve the feed block itself, whose Reaktoro block lives in a separate sub-model (section 6.1). They don't touch a surrogate.
- `multi_comp_mvc` and the translator are not used in `softening_acid_ro`.
- Phase 2 still adds the stale-surrogate guard (section 6.5) as a safety net. A small REW change, using `get_cyipopt_watertap_solver` in those four places, would remove the issue for good, but it isn't needed for `softening_acid_ro`.

### 7.3 Performance tracking
`softening_acid_ro.solve_model` parses Ipopt's `output_file`. The driver doesn't handle that file specially (section 6.6), so the parse reflects the last inner solve. For trust-region totals, REW can read `results.solver.trust_region` (optional; not needed for the flowsheet to work).

### 7.4 Expected behaviour on this flowsheet

| Aspect | Gray box | Trust region |
|---|---|---|
| Reaktoro fails inside an Ipopt line search (dosing near an equivalence point, high-recovery RO interface) | Ipopt evaluation error; solve aborted after `max_reaktoro_failed_solves` | trial point rejected, radius shrunk; Reaktoro only at subproblem solutions |
| Large steps in dose, pH or RO concentrations | fraction-to-boundary only | bounded directly by the input trust region |
| Hessian | LBFGS approximation per Ipopt iteration, multiplier-weighted | exact Hessian of the surrogate inside Ipopt; LBFGS/SR1 `B_o` updated per outer step |
| Reaktoro solves per flowsheet solve | one per Ipopt iteration and line-search trial, all blocks | one per outer iteration, all blocks in parallel; outer iterations expected 5–30 (to be measured) |
| Recovery sweeps (60 → 70 → 80%) | warm start from the previous recovery | same, and the trust region protects the first step after a recovery change |
| Ill-posed chemistry (singular Jacobian, trace species) | fails | same Jacobian, so it still fails; needs a formulation or scaling fix |

## 8. Phases

Validation in all phases uses **reaktoro-pse only**: examples and the test suite. `softening_acid_ro` (section 7) is the design target for interfaces and behaviour, and is validated after these phases.

| Phase | Content | Exit criterion |
|---|---|---|
| 0. Prototype (scratch, no package changes) | The copied loop, the adapter, the elastic mixed trust region, and hand-built linear **and** quadratic surrogates (SR1, LBFGS) for one block, reusing `solve_reaktoro_block`. Warm-start defaults under the user's options. Run on `examples/simple_desalination` (single block) and `examples/thermal_precipitation` (manager, parallel), each with a square solve followed by an optimization. Compare with the gray box: Reaktoro calls, outer and inner iterations, wall time, robustness to perturbed starts. | Same solution as the gray box. Call counts known. Linear vs quadratic compared. Warm-dual mechanism settled. |
| 1. Core backend | The `backend` option with its resolution order. Surrogate type from `hessian_options`. The surrogate block (linear and quadratic) with the same equation count as the gray box. The adapter, the copied loop (option B), the elastic region, warm starts, the results object, activate/deactivate, and output floors moved to a shared helper. | `examples/tests/test_examples.py` passes on both backends (parametrized). `degrees_of_freedom` is the same on both backends. |
| 2. Manager and solver selection | The aggregate surrogate in `build_reaktoro_blocks`, using the existing windows and `parallel_solver`. The `get_cyipopt_watertap_solver` wrapper, the statistics (section 6.6), the stale-surrogate guard. | Manager examples (`simple_desalination`, `thermal_precipitation`, `redox_tracking`) pass with TR in parallel and serial mode, and parallel matches serial. A solve through `get_cyipopt_watertap_solver` of a model with TR blocks uses TR. |
| 3. Robustness and speed | Criticality and termination tuning, tuning of the `B_o` update safeguards for once-per-step updates, per-block radius if the logs call for it, failure handling (Reaktoro failure at a trial point, failed inner solve). | Full reaktoro-pse test suite green on both backends. TR results match the gray box within tolerance. |

Tests to add:
- The surrogate matches truth value and Jacobian at `w_k`.
- A rejected step restores the model state.
- The elastic part is deactivated and fixed after success and after an exception.
- `degrees_of_freedom` is the same on both backends.
- Square solves work.
- A Reaktoro failure at a trial point is handled.
- A results object from a non-converged run fails `assert_optimal_termination`.
- Backend resolution (manager > block), and the surrogate type from `hessian_options`.
- Solver selection: a model with surrogates goes through the TR driver.
- Statistics: present in `results.solver.trust_region` and `solver.trust_region_stats`, consistent totals, `None` for a gray-box solve.
- Elastic penalty: fixed and adaptive both converge on the examples.
- `keep_surrogate_curvature`: `B_o` is reset by default and kept when set.
- Parallel and serial manager give the same answers.

## 9. Decisions

Decided (2026-11-23):
- **Option B:** reuse Pyomo's algorithm pieces with our adapter and a copied main loop. Option C (upstream PR to remove the copy) comes after B is shown to work.
- **One switch:** `backend` on the manager or block (no global default); everything else unchanged (section 6.1).
- **Surrogates:** linear and quadratic, chosen from the existing `hessian_options` (SR1 default; LBFGS for `softening_acid_ro`).
- **Inner warm starts:** driver defaults under the user's options; user selections never overridden (section 6.4).
- **Parallelism:** one aggregate surrogate for the manager; per-block surrogates only without a manager (section 6.2).
- **Trust-region variables:** chemistry inputs, elastic, mixed relative/absolute box (section 6.3).
- **Solver selection:** automatic inside `get_cyipopt_watertap_solver` (section 6.5).
- **Validation:** reaktoro-pse only, first (section 8).
- **Elastic penalty:** both fixed and adaptive modes are implemented; the default is chosen from phase 0 results (section 6.1).
- **Curvature across solves:** reset each solve by default; `keep_surrogate_curvature=True` keeps it as a warm start (section 6.1).
- **Statistics:** no special handling of the user's `output_file`. Trust-region statistics, which the stock TRF doesn't return, are exposed through the results object, the solver and a helper (section 6.6).

Still open: none at the planning level. Defaults (penalty mode and values, radii) are set in phase 0.

## 10. Implementation status (2026-11-23, branch `try_trust_region`)

### New and changed files

| File | Content |
|---|---|
| `reaktoro_block_config/trust_region_options.py` (new) | `BackendTypes`, `SurrogateTypes`, `ElasticPenaltyModes`, `backend` ConfigValue, `TrustRegionOptions` (Pyomo's `_trf_config()` with our defaults, plus surrogate, elastic penalty (fixed/adaptive), `keep_surrogate_curvature`, `warm_start_options`, `log_input_step`, `maximum_failed_inner_solves`) |
| `core/reaktoro_surrogate.py` (new) | `ReaktoroSurrogateBlock`: same `inputs`/`outputs` as the gray box, linear or quadratic surrogate, elastic trust region added only during a solve, `get_external_model()` for compatibility. Truth, output floors, scaling and curvature come from the existing `ReaktoroGrayBox` object (not added to the model), so floors were not moved to a shared helper. |
| `core/util_classes/trust_region_solver.py` (new) | `ReaktoroTRInterface` (adapter), copied main loop `trust_region_method`, `ReaktoroTrustRegionSolver`, statistics, `get_trust_region_stats`, `find_surrogate_blocks` |
| `core/util_classes/cyipopt_solver.py` | `get_cyipopt_watertap_solver` returns `ReaktoroSolverWrapper` (TR when the model has surrogates, plain solve otherwise); `get_reaktoro_solver` alias |
| `core/util_classes/hessian_functions.py` | `HessianApproximation.reset()`, `get_output_matrices()` |
| `core/reaktoro_block_builder.py`, `reaktoro_block.py` | `backend` and `trust_region_options` options; mismatch with the manager's backend raises |
| `parallel_tools/reaktoro_block_manager.py` | `backend`/`trust_region_options` on the manager; aggregate surrogate; `parallel_solver` collects every worker before re-raising an error; serial mode uses the builder's coupled solver for direct-coupling blocks |
| `parallel_tools/parallel_manager.py` | `LocalWorker.display_state` reads the worker's reply (bug fix, see below) |
| `examples/simple_desalination.py`, `examples/thermal_precipitation.py` | `backend` (and `parallel_mode` for `main`) arguments |
| tests | `tests/test_trust_region.py` (18 tests), plus unit tests in `test_hessian_functions.py` and `test_manager.py` |

### Results on reaktoro-pse examples

| Case | Gray box | Trust region |
|---|---|---|
| `simple_desalination`, no manager (square solve + optimization) | 4.5 s | 4.3 s; square: 1 outer iteration; optimization: 6 outer iterations, 7 Reaktoro solves |
| `simple_desalination`, serial manager | 9.5 s | 1.9 s |
| `simple_desalination`, parallel manager | – | same iterations and result as no manager |
| `simple_desalination` with SR1, LBFGS, BFGS (quadratic) and ZeroHessian (linear) | – | all converge to the reference values |
| `thermal_precipitation`, no manager | 10.2 s | 10.5 s, 3 outer iterations |
| `thermal_precipitation`, parallel manager | fails at build (see below) | converges, 3 outer iterations |

### Differences from the plan
- Warm starts: only primal warm start with `bound_push`/`bound_frac` = 1e-8 for inner solves after the first (driver defaults, user options win). Dual warm start is not implemented.
- The stale-surrogate guard (section 6.5) and per-block radii (section 6.3) are not implemented.
- Inner Ipopt iteration counts come from an `intermediate_callback`; if the caller passes their own callback, counts are `None`.

### Existing problems found while testing
- **Fixed:** `LocalWorker.display_state` never read the worker's reply. Every later solve read that stale reply and returned shared-memory results before the worker finished, so results after any `display_reaktoro_state()` call in parallel mode could be one evaluation behind (both backends).
- **Fixed:** the serial manager (`use_parallel_mode=False`) failed for direct-coupling blocks (`block.solver` is `None`).
- **Not fixed (gray box only):** `ExternalGreyBoxBlock` cannot build the aggregate gray box when managed blocks have input keys of different tuple lengths (e.g. `thermal_precipitation` with the manager). The surrogate uses `dimen=None` and is not affected. When this fails after workers are started, the process hangs until the workers time out.
- **Not fixed:** `thermal_precipitation` with the serial manager fails in Reaktoro during block initialization (this path was not reachable before).

## 11. REW `softening_acid_ro` benchmark (2026-11-24)

Run with the REW test setup (`multi_process_reaktoro=True`, LBFGS + GaussNewton, `watertap_default` costing, 5000 m3/d, no HPRO): `initialize`, then optimization at water recovery 0.7 and 0.8 (0.6 only for Seawater). REW was not modified, the backend was switched by setting the `ReaktoroBlockManager` `backend` default to `trust_region` (same as passing `backend="trust_region"`). Script: `scratchpad/trf/sar_bench.py`.

| Water | Gray box total (init / opt 0.7 / opt 0.8) | Trust region total (init / opt 0.7 / opt 0.8), outer iterations | Result |
|---|---|---|---|
| USDA_brackish | 49.0 s (16.7 / 15.5 / 14.4) | 73.3 s (24.5 / 36.6 / 9.8), 13 / 32 / 8 | same LCOW 0.476189 and doses |
| sample_500_hardness | 44.4 s (29.5 / 6.8 / 5.8) | fails in opt 0.7 (every inner solve fails) | - |
| sample_1500_hardness | 30.9 s (17.7 / 5.3 / 5.6) | 100.5 s (84.2 / 7.2 / 6.9), 66 / 6 / 6 | same LCOW 1.183421 and doses |
| Seawater (opt 0.6) | 22.0 s (16.1 / 3.5) | 23.2 s (14.0 / 6.8), 6 outer in opt | same LCOW 0.587003 and doses |

### Driver changes needed for this flowsheet
1. **Current iterate in acceptance test.** Pyomo only checks the filter, so from a feasible start any objective decrease was accepted as an f-step even when theta went from 1e-10 to 1e2. Trial points must now also be acceptable to the current iterate (Fletcher-Leyffer), unless theta <= `minimum_feasibility`.
2. **Restoration step.** `initialize` fixes water recovery from about 0.92 to 0.5 between solves, so the flowsheet forces the RO inputs outside of the trust region. Rejecting these steps repeated the same subproblem forever. A rejected step that needed elastic slack is accepted while the elastic gap decreases, and the filter is reset.
3. **Soft trust region (`trust_region_formulation="soft"`, new default).** With the surrogate evaluated at the trusted point `w0 + h*s`, the outputs it can reach shrink with the radius, and the subproblem became infeasible after large forced moves. Now the surrogate is evaluated at the actual inputs and the step beyond the radius is an exact l1 penalty. `trusted_point` is kept as an option.
4. **Curvature reset after large steps (`curvature_reset_step`, default 1).** A quasi-Newton update over a restoration step or a large forced step (relative step 1.9) gave curvature that made every later step fail. Curvature is reset (linear surrogate) instead.
5. **Relative theta.** theta was `sum(sf*|y_surr - y_truth|)`, outputs without scaling factors (e.g. `alkalinityAsCaCO3` in mg/L, sf=1) dominated it, so steps were limited to about 1e-3. theta is now `sum(|y_surr - y_truth| / max(|y_truth|, 1/sf))`.
6. **Adaptive penalty** only increases when the elastic gap did not decrease, `elastic_penalty_max` is now 1e6 (1e8 made inner solves fail).
7. **Bug:** pH inputs of managed speciation blocks (names like `(2, 's_0', 'pH')`) were not recognised as log inputs. `log_input_step` default changed from 10 to 1 (one radius unit is 1 pH unit).

### Remaining problem
- `sample_500_hardness`, optimization at 0.7 from the 0.5 solution: all inner solves fail at every radius (0.5 down to 1e-3), so no step is ever accepted. The gray box evaluates reaktoro along the Ipopt path, the trust region has to solve the jump with chemistry linearized at the 0.5 point.
- Both hardness samples need 66 outer iterations in `initialize` (gray box: one Ipopt solve).
