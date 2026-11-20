How to troubleshoot solver issues
=================================

Dual infeasibility does not close
---------------------------------

Ipopt reduces ``inf_pr`` below 1e-8, but ``inf_du`` stops decreasing, and the solver does
not terminate or ends with *Solved to Acceptable level*. This is caused by the
approximated Hessian of the gray box, which affects the estimate of dual infeasibility
(see this `Ipopt mailing list discussion <https://list.coin-or.org/pipermail/ipopt/2007-February/000700.html>`_).

**A. Recompute multipliers with least squares** (`recalc_y <https://coin-or.github.io/Ipopt/OPTIONS.html#OPT_recalc_y>`_):

.. code-block:: python

   solver = get_cyipopt_watertap_solver()
   solver.options["recalc_y"] = "yes"
   solver.options["recalc_y_feas_tol"] = 1e-2

**B. Use outputs with exact derivatives.** Numerical derivatives carry errors that reduce
the accuracy of dual infeasibility estimates. Check which outputs are exact, calculated
(chain rule), or numerical with ``m.block.display_jacobian_outputs()``. Properties that
Reaktoro differentiates exactly can be used to build the desired property with exact
derivatives, see :doc:`/technical_reference/outputs`.

**C. Try a different Hessian approximation**, e.g.
``hessian_options={"hessian_type": "BFGS"}``. The default (SR1) is the most robust in
benchmarks, but individual problems can behave differently.

**D. Check for absent phases** pinning bounded variables, see
:ref:`how_to_guides/how_to_scale_and_initialize:Species amount offset`.

Iterates diverging
------------------

.. code-block:: text

   EXIT: Iterates diverging; problem might be unbounded.

This is generally caused by bad scaling. If the model solves without the ReaktoroBlock,
the likely cause is the automatic Jacobian (and possibly variable) scaling. Provide
Jacobian scaling for the problematic outputs:

.. code-block:: python

   jacobian_options={
       "user_scaling": {
           ("molarEnthalpy", None): 1,
           ("specificHeatCapacityConstP", None): 1,
       },
   }

Check the Jacobian scaling with ``m.block.display_jacobian_scaling()``. Alternatively
increase the divergence tolerance, ``solver.options["diverging_iterates_tol"] = 1e30``.

Reaktoro fails to solve during optimization
-------------------------------------------

When Reaktoro fails during an Ipopt solve, the block raises an evaluation error and Ipopt
cuts back its step. After ``max_reaktoro_failed_solves`` (default 2) consecutive failures
the solve stops. Common causes:

* inputs reach unphysical values (negative amounts, extreme pH), add bounds on the Pyomo
  variables that feed the block
* the redox state is not defined for a redox sensitive water, see
  :doc:`how_to_track_redox`
* species that never form cause trouble in the property block, remove them with
  ``speciation_block_exclude_species`` (no checks are made that they stay at 0) or from
  the system with ``exclude_species_list``

Inspecting a block
------------------

.. code-block:: python

   m.block.display_reaktoro_state()     # Reaktoro chemical state (all species, properties)
   m.block.display_jacobian_outputs()   # derivative type of each output
   m.block.display_jacobian_scaling()   # Jacobian scaling factors
