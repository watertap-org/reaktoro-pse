How to scale and initialize
===========================

Scaling
-------

Reaktoro-PSE uses IDAES scaling factors (``iscale.set_scaling_factor``) on input and
output variables to scale the gray box:

* **Inputs**: set scaling factors on composition, chemical additions, and other unfixed
  inputs, e.g. ``1 / expected value``. If not set, they are computed from current values.
* **Outputs**: scaling factors on output variables are computed from their values on
  initialization if not set. For ``speciesAmount`` of phases that can be absent, set them
  from the amount of the limiting element (e.g. ``1 / feed Ca`` for Calcite), as a value
  based factor of an absent phase is meaningless (see below).
* **Jacobian**: the gray box Jacobian is scaled with ``jacobian_options["scaling_type"]``
  (default combines output and input scaling factors, see
  :doc:`/technical_reference/reaktoro_block_options`). Check it with
  ``m.block.display_jacobian_scaling()``, and override individual outputs with
  ``jacobian_options={"user_scaling": {("molarEnthalpy", None): 1}}``.

When scaling factors change after the block is built (e.g. after initialization of the
rest of the flowsheet), update the block:

.. code-block:: python

   m.block.update_jacobian_scaling()   # or m.block.update_block_scaling()

Initialization
--------------

``m.block.initialize()`` solves the Reaktoro equilibrium problem at current input values
and sets all outputs (and internal variables) to the result. Initialize upstream units
first so inputs have realistic values, then initialize each ReaktoroBlock, then solve the
flowsheet:

.. code-block:: python

   calculate_variable_from_constraint(m.concentrate["Ca"], m.eq_concentrate["Ca"])
   m.block.initialize()
   solver.solve(m)

To speed up flowsheet initialization, a block can be temporarily removed from the
problem: ``m.block.deactivate()`` fixes its outputs at current values and deactivates
the gray box, ``m.block.activate()`` restores it.

If Reaktoro fails to converge during initialization, presolving to a lower tolerance can
help: ``reaktoro_presolve_options={"presolve_during_initialization": True}``.

Species amount offset
---------------------

When a mineral or other species is not present, Reaktoro returns its amount as epsilon
(~1e-40) with zero derivatives. If its output variable (or a variable linked to it, e.g.
precipitated solids in a mass balance) has a lower bound of 0, it sits at the bound,
which makes the problem degenerate and can prevent Ipopt from converging.

``speciesAmount`` outputs are therefore offset by ``a = species_amount_zero_offset /
scaling factor`` (default ``1e-4`` in scaled units) with the ``smooth_floor`` type
(``output = sqrt(amount**2 + a**2)``): absent species sit at ``a``, and present species
(amount ``n``) are changed by only ``~a**2/(2n)``. This relies on output scaling factors
set from the limiting element, as above. Set
``reaktoro_solve_options={"species_amount_zero_offset": 0}`` to disable it.

After solving, check the error the offset introduced:

.. code-block:: python

   status = m.block.check_species_amount_offsets(raise_error=True)

It flags present species changed by more than 1% and elements of which the offsets take
up more than 1% of the system amount. See ``examples/evaporative_crystallizer.py`` for a
case that only converges with the offset.
