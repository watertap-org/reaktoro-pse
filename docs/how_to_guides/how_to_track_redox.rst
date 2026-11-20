How to specify and track the redox state
========================================

With apparent species (``exact_speciation=False``), H and O are open and pH fixes only
one of the two degrees of freedom they leave. The other one is the **redox state**: how
oxidizing or reducing the water is, which sets how Fe, Mn, S, N, and C split between
oxidation states (e.g. Fe(II)/Fe(III), sulfide/sulfate, ammonium/nitrate). If it is not
specified, Reaktoro finds a state that is numerically valid but often physically wrong
(e.g. nitrate reduced to N2), and results can change between solves and platforms.

For waters with redox sensitive species, and a database that contains them
(see :doc:`/technical_reference/databases`), define the redox state with **one** of:

.. list-table::
   :header-rows: 1
   :widths: 20 30 50

   * - ``system_state`` key
     - Units
     - Use when
   * - ``pE``
     - dimensionless
     - pE (or ORP converted to pE) is known, or the water is anoxic
   * - ``Eh``
     - V
     - redox potential is known
   * - ``dissolved_O2``
     - same basis as composition (e.g. mol/s, kg/s)
     - dissolved oxygen is measured (most common for treatment trains), water is oxic
   * - ``O2_fugacity``
     - bar
     - water is in equilibrium with a gas, e.g. 0.21 for air (O2 saturation)

Only one can be given, more raises an error. ``dissolved_O2`` and ``O2_fugacity`` need a
database with O2(aq) (not ``pitzer.dat``). Concentrations (mg/L) are not supported and
raise an error, convert dissolved oxygen to a flow (e.g. mg/L times L/s).

.. code-block:: python

   ReaktoroBlock(
       aqueous_phase={"composition": m.feed, "convert_to_rkt_species": True,
                      "activity_model": "ActivityModelPhreeqc"},
       system_state={
           "temperature": m.temperature,
           "pressure": m.pressure,
           "pH": m.pH,
           "dissolved_O2": m.feed_dissolved_O2,   # mol/s
       },
       database_file="wateq4f.dat",
       ...
   )

Aeration and oxygen dosing
--------------------------

In a block with a speciation block, the property block can change the oxygen content:

* **Equilibrium with air** (open basin, aerated softening): give ``O2_fugacity`` in
  ``system_state_modifier``. O2 is taken up or stripped until the water is saturated.

  .. code-block:: python

     system_state_modifier={"O2_fugacity": m.O2_fugacity},  # 0.21 bar

  0.21 bar is a good default at sea level and moderate temperature. For warm water
  correct for water vapor, f(O2) ~ 0.2095 (P_air - p_H2O(T)) (about 7% lower at 40 C,
  20% lower at 60 C), and for elevation (P_air ~ 1 atm exp(-elevation / 8400 m)).
  Partial aeration (e.g. 80% of saturation) can be modeled as 0.8 times the fugacity of
  air. Salting out is only approximately captured by activity models.
* **Known oxygen transfer** (aerator, pure O2 dosing): add ``O2`` as a chemical,
  ``chemistry_modifier={"O2": m.O2_dose}``. All of it dissolves, and dissolved oxygen
  leaving the block is found from equilibrium (O2 consumed by oxidation of Fe, Mn, etc.).
  There is no gas phase, so O2 is not stripped if the water is supersaturated.
  Fixing the dissolved O2 output and unfixing the dose solves for the dose that gives a
  target dissolved oxygen.

Passing the redox state between blocks
--------------------------------------

Pass the redox state between units the same way as pH: request it as an output of one
block and connect it to the input of the next.

.. code-block:: python

   m.softening_outputs = Var([("pH", None), ("speciesAmount", "O2"), ...])
   ...
   m.neutralization = ReaktoroBlock(
       system_state={
           "pH": m.softening_outputs["pH", None],
           "dissolved_O2": m.softening_outputs["speciesAmount", "O2"],
           ...
       },
       ...
   )

Within a block, changes in the property block (acid, base, aeration) update pE and
dissolved oxygen consistently, as the property block conserves all elements and charge.
Between blocks, what to pass depends on the water:

.. list-table::
   :header-rows: 1
   :widths: 35 20 45

   * - Connection
     - Pass
     - Why
   * - ReaktoroBlock to ReaktoroBlock
     - either
     - both reproduce the upstream state, choose by the rows below
   * - Oxic water (DO above ~0.01-0.1 mg/L)
     - ``dissolved_O2``
     - conserved, mixes linearly, measurable
   * - Through mixers, recycles, or units without Reaktoro (e.g. MCAS RO)
     - ``dissolved_O2`` as a component
     - follows the mass balance, pE does not mix and goes stale in units that change
       composition
   * - Anoxic or reduced water (DO ~ 0)
     - ``pE``
     - DO near 0 is a poor specification (pE depends on log(DO)), and pE carries the
       oxidation state of Fe, Mn, S
   * - ``pitzer.dat``
     - neither
     - no redox chemistry in the database

Waters with nothing left to hold their redox state (poorly poised, e.g. dissolved Mn(II)
after all oxidants are consumed) have a pE set by trace species, and adding chemicals to
them can be numerically fragile regardless of what is passed.

See ``examples/oxygen_tracking.py`` (:doc:`/technical_reference/examples`) for a complete
treatment train, and :doc:`/background/redox_state` for the theory.
