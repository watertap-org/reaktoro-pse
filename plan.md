# Grey-box Hessian improvement plan

Two milestones, plus a small shared step 0:

- **Milestone 1:** partitioned SR1/PSB quasi-Newton (fixes the positive-definite updates on non-convex outputs).
- **Milestone 2:** Hessian from finite differences of Reaktoro's exact sensitivity Jacobian.

Items 1–6 of the original review (scalar4 crash, `None` return, `np.sum` zero checks,
L-BFGS starting matrix, repeated-point history, stale Jacobian) are already fixed on
`stability-revisit`, with tests in
[test_hessian_functions.py](src/reaktoro_pse/core/util_classes/tests/test_hessian_functions.py).

## Constraints found in the code

- **Cross-block coupling in the aggregated grey box.** `ReaktoroBlockManager` puts every
  block under one `ReaktoroGrayBox`
  ([reaktoro_block_manager.py:338-352](src/reaktoro_pse/parallel_tools/reaktoro_block_manager.py#L338-L352)),
  so each output's BFGS matrix covers *all* inputs of *all* blocks. Output yᵢ only
  depends on its own block's inputs, but the update term B s sᵀB adds curvature between
  blocks that has no physical basis. It also costs n_out × n_total² memory.
- **PyNumero requires a fixed Hessian sparsity pattern.** It asserts that `row`/`col` are
  identical on every call (`pyomo_grey_box_nlp.py`, `evaluate_hessian_lag`). Any new
  method has to pick its sparsity pattern once, at configure time.
- **Multipliers arrive unscaled.** PyNumero passes Ipopt's multipliers through unchanged,
  and the Hessian is expected in the original, unscaled variables. Any scaling we add has
  to be internal and mapped back before returning.

---

## Step 0 (shared): benchmark harness

Both milestones need the same yardstick, so this comes first.

- [x] A dev script, not shipped, that runs the examples and tutorials against a matrix of
      `hessian_type` values:
  - current LBFGS, BFGS_ipopt, ZeroHessian;
  - Ipopt's own `limited-memory` with `no_hessian_estimation`;
  - new types as they are added.
  - Done as `dev/bench_softening_ro.py` on the reaktoro_enabled_watertap softening-acid-RO
    flowsheet (5 waters, recovery sweep 0.5-0.9, warm started); `dev/summarize_bench.py`
    builds the comparison table. Examples/tutorials are covered by the pytest suite instead.
- [ ] Record per run: solve status, Ipopt iterations, wall time, number of Reaktoro solves,
      restoration-phase entries.
  - Status, iterations, time and LCOW are recorded. Reaktoro solve counts and restoration
    entries are not yet.
- [x] Write the results as a CSV table, so every later change is compared against the
      same baseline.

---

## Milestone 1: partitioned SR1/PSB quasi-Newton

**Goal:** per-output curvature estimates that are allowed to be indefinite, only cover
the inputs each output depends on, and are computed in scaled variables.

### M1.1 Refactor the update loop

- [x] Pull the shared machinery in `HessianApproximation` (history storage, `new_point`,
      starting matrix, assembly with λ) into one place.
- [x] Each method then only supplies `update(B, s, y) -> B`.
- [x] Port BFGS, CBFGS, BFGS_mod, BFGS_damp and BFGS_ipopt onto it unchanged, so the old
      behaviour stays available for comparison. Check the port with the existing tests.
- [x] Decide on `reset_memory(leave_last=True)`. It currently does nothing, because it
      clears `memory_limited_arrays[:-2]` and that slice is always empty. Fixing it
      changes BFGS right after initialization.

  - Resolved: the no-op call was removed. Behaviour is unchanged (the first pair both
    scales B0 and is applied as the first update).
### M1.2 SR1 and PSB update rules (new `HessTypes`)

- [x] **SR1:** with r = y − Bs, apply B += r rᵀ / (rᵀs) only if |rᵀs| ≥ τ‖s‖‖r‖,
      with τ ≈ 1e-8 (Nocedal & Wright). The update can be indefinite.
- [x] **PSB:** B += (r sᵀ + s rᵀ)/(sᵀs) − (rᵀs) s sᵀ/(sᵀs)². Always defined when s ≠ 0,
      symmetric, no positive-definite requirement.
- [x] **Starting matrix:** allow a signed scalar sᵀy/sᵀs, plus a new "zero" start so no
      curvature is invented in unexplored directions. Make it configurable; step 0 picks
      the default.
- [ ] **Safeguards:** clip entries at `init_max`, reset any matrix that becomes
      non-finite, keep the existing reset when all multipliers are zero.
  - Non-finite reset and zero-multiplier reset done. Clipping at `init_max` was **not**
    added to SR1/PSB: real unscaled curvature can legitimately exceed 1e8.

### M1.3 Per-output input support and a fixed sparsity pattern

- [x] Add an optional `hessian_support` argument to `ReaktoroGrayBox.configure`: for each
      output, the list of input indices it depends on.
  - Single block: all inputs.
  - Aggregated manager: built from the existing `input_windows` / `output_windows`.
- [x] Updates for each output work on s[support] and y[support] and only touch
      B[support, support].
- [x] At configure time, compute the lower-triangular pattern as the union over outputs of
      support × support. Always return exactly that pattern (this replaces `_hand_tril`).
- [x] Apply this to every quasi-Newton type, the old BFGS ones included. On its own it
      removes the spurious cross-block coupling and the memory blow-up.

### M1.4 Scaled variables

- [x] Work in scaled variables x̃ = Dx, using the `get_input_scaling()` values that are
      already wired up: s̃ = Ds, ỹ = D⁻¹y, and the returned Hessian is B = D B̃ D.
- [ ] This makes τ and ε mean the same thing across inputs that range from about 300 K to
      1e-5 mol.

### M1.5 (stretch) Lagrangian-level SR1

- [x] One matrix per block, updated with y_L = (J₊ − J)ᵀλ₊ (new multipliers on both sides,
      as Ipopt does). Cheaper, and a useful comparison against the per-output approach.

### M1.6 Config plumbing

- [x] Route the new options (`sr1_skip_tol`, the zero-start choice, the scaling switch)
      through `HessianOptions`, `set_solver_options`, `export_config` /
      `load_from_export_object`, and the block-manager `configure` call.
- [x] Unify the three conflicting default sets:

  | Setting | `hessian_options.py` | `reaktoro_solver.py` | `HessianApproximation` |
  |---|---|---|---|
  | type | LBFGS | ZeroHessian | — |
  | const | 1e-64 | 1e-16 | 1e-32 |
  | memory | 6 | 3 | — |
  | ε | 1e-16 | 1e-12 | — |

### M1 tests

- [x] On a quadratic with linearly independent steps, SR1 recovers A exactly in n steps
      (its finite-termination property).
- [x] Indefinite (saddle) A: SR1 and PSB match it, while BFGS stays positive definite and
      gets the sign wrong.
- [x] PSB satisfies the secant condition.
- [ ] The SR1 skip rule fires when it should.
- [x] Two blocks: no entries between blocks, and an identical COO pattern on every call.
- [x] Rescaling the inputs gives the same Hessian once mapped back.
- [x] Integration: the existing examples run with SR1 and PSB.

### M1 done when

On the step-0 benchmark, SR1 or PSB with support restriction and scaling is at least as
robust as the current LBFGS default and needs fewer iterations on most cases. If it is,
switch the default.

**Status (M1 implemented, benchmark run, default not yet switched).** Softening-acid-RO,
5 waters × recovery sweep = 38 solve points, MUMPS (MA27 unavailable in this env).
Iterations summed over all 38 points, relative to the baseline (`HEAD`) LBFGS + GaussNewton:

| case | solves ok | iters | vs baseline LBFGS |
|---|---|---|---|
| SR1 + GaussNewton | 38/38 | 1339 | 0.55 |
| SR1 + scalar1, scaled | 38/38 | 1381 | 0.57 |
| SR1 + zero | 38/38 | 1407 | 0.58 |
| BFGS_ipopt + GaussNewton | 38/38 | 1457 | 0.60 |
| SR1_lagrangian + zero | 38/38 | 1535 | 0.63 |
| LBFGS + GaussNewton (new code) | 38/38 | 1908 | 0.78 |
| LBFGS + GaussNewton (baseline) | 38/38 | 2435 | 1.00 |
| SR1 + scalar1, unscaled | 29/38 | – | init fails on BGW_500 |
| PSB (any start) | ≤ 29/38 | – | not robust |

- With SR1, the maximum iterations over all points is 51, against 394 for baseline LBFGS.
- SR1 reaches the same LCOW as BFGS_ipopt at every point. LBFGS (old and new code) lands at
  a worse point at BGW 0.85 (+17%) and SW_HPRO 0.8 (+28%).
- Input scaling matters only for the scalar1 start. GaussNewton and zero starts are robust
  without it.
- PSB is not robust here and should not be a default.
- Recommended default: `hessian_type = SR1`, `bfgs_initialization_type = GaussNewton`,
  `hessian_input_scaling = False`.

---

## Milestone 2: finite-differenced exact Jacobian

**Goal:** a near-true Hessian of the Lagrangian, computed by finite-differencing
Reaktoro's analytic sensitivity Jacobian (`sensitivity.dudw()` plus the existing
property chain).

### M2.1 Standalone finite-difference utility

For a solver at a converged point, compute the perturbed Jacobians {J(w + hₖeₖ)} for
k = 1…n, and from them H(λ) = ½(G + Gᵀ), where column k of G is (Jₖ − J₀)ᵀλ / hₖ.

- [ ] **Separate "shadow" solver.** `solve_reaktoro_block` changes many things on the main
      solver: the warm-start `ChemicalState`, input temp values, `rkt_scaling_factor`s,
      `_old_input_params`, output values, the Jacobian spec's property states and the
      failure counter. So perturbations run on a separate solver built from
      `export_config()`.
  - [ ] Move `RktModelData` from `parallel_manager.py` into `core` so both paths share it.
  - [ ] The coupled solver uses the same dict-of-configs construction as
        `RemoteWorker.process_config_data`. Its combined Jacobian is then differentiated
        directly, which picks up the speciation × property second-order terms that are
        currently dropped.
  - [ ] Before each column, copy the main solver's converged state into the shadow solver
        as the warm start.
- [ ] **Step size:** hₖ = δ·max(|wₖ|, typical value), with δ ≈ √(Reaktoro tolerance)
      (about 1e-4 at the default tolerance of 1e-8).
  - Per-input-type overrides through a dict, mirroring `numerical_step`.
  - Forward differences by default (n solves), central as an option (2n).
  - Species near their lower bound only step in the positive direction.
- [ ] **Diagnostic:** log the asymmetry ‖G − Gᵀ‖/‖G‖. A large value points to noise or a
      phase boundary being crossed.
- [ ] Use this utility as the **ground truth for evaluating M1 offline**: compare the
      SR1/PSB/BFGS matrices against it at saved iterates. Build it alongside M1 for that
      reason.

### M2.2 Grey-box integration (new type, e.g. `HessTypes.FD_exact_jacobian`)

- [ ] **Caching:** cache the perturbed Jacobians per input vector. A call with only new λ
      just re-contracts and needs no new solves.
- [ ] **Failed solves:** if a perturbed solve fails, retry with the opposite sign. If that
      fails too, use the M1 quasi-Newton result for that evaluation, or zero that column.
  - A Hessian failure never increments the main solver's `_sequential_fails` and never
    raises `CyIpoptEvaluationError`.
- [ ] **Sparsity:** return the same fixed pattern defined in M1.3.

### M2.3 Parallel manager

- [ ] Compute each block's Hessian inside its worker, so blocks run in parallel and the
      result is block-diagonal.
- [ ] Add a worker message for the Hessian and shared-memory buffers for λ and the
      n_in × n_in result.
- [ ] Workers keep their perturbed Jacobians cached, so a call with only new λ costs one
      message and no solves.

### M2.4 Hybrid refresh

- [ ] Recompute the finite-difference Hessian every N iterations, or when the M1 skip
      tests or asymmetry checks fail, and use SR1 updates in between. This bounds the
      cost to roughly n/N extra solves per iteration.

### Accuracy caveat

For outputs whose Jacobian is itself numerical (`JacType.numeric`, step 1e-8), this is a
difference of differences, with an expected relative error of about 1e-8 / 1e-4 = 1e-4.
Report which outputs are affected using `display_jacobian_output_types()`.

### M2 tests

- [ ] Mock solver with J(w) linear in w: H equals Σλᵢ Aᵢ exactly, is symmetric, keeps a
      fixed pattern, causes no re-solves when only λ changes, and falls back correctly on
      failure.
- [ ] Real Reaktoro on the `test_reaktoro_solver` fixture: compare with second differences
      of the outputs (loose tolerance) and check the asymmetry is small.
- [ ] Real Reaktoro: the main solver's outputs and J are bit-identical before and after a
      Hessian evaluation (the warm start isn't disturbed).
- [ ] Coupled solver and parallel manager cases.

### M2 done when

On the step-0 benchmark it needs clearly fewer Ipopt iterations than M1's best option,
with no new failures. Report wall time separately, since this method trades iterations
for Reaktoro solves.

---

## Open questions

- **Order:** M1 first, with M2.1 built alongside it as the yardstick for M1?
- **`reset_memory(leave_last=True)`:** fix it as part of M1.1?
- **Benchmark script location:** an untracked `dev/` folder, or committed under
  `examples/benchmarks/`?
- **Longer term (upstream):** an analytic second-order sensitivity needs Reaktoro to expose
  its KKT factorization or second-order sensitivities. That's a feature request for the
  Reaktoro project, not part of this plan.
