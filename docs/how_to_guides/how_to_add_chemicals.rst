How to add chemicals and modify the state
=========================================

Most water treatment steps change an equilibrated water: acid or base dosing,
softening with lime, heating, evaporation, or aeration. Because the lab composition only
defines an apparent state at a measured pH, the block first has to find the true
speciation of the feed, and then the new equilibrium after the change.

Speciation and property blocks
------------------------------

With ``build_speciation_block=True`` the ReaktoroBlock contains two equilibrium problems:

#. **speciation block**: speciates the feed at the given pH (apparent species, H and O
   open, charge balanced on ``charge_neutrality_ion``), and outputs the amounts of all
   true species.
#. **property block**: takes the true species from the speciation block (exact
   speciation, all elements fixed), adds chemicals, applies state changes, finds the new
   equilibrium (pH is now a result), and computes the requested outputs.

.. image:: ../tutorials/notebooks/spec_block.png
   :width: 100%
   :alt: ReaktoroBlock with speciation block

By default (``direct_speciation_to_property_block_coupling=True``) both are solved in
one gray box and derivatives are combined by chain rule. With ``False`` they are two
gray boxes linked by Pyomo constraints, which is less stable, as the solver can take
steps that give an infeasible true composition.

Mineral and gas phases are only added to the property block, unless
``build_speciation_block_with_phases=True``.

Chemical additions
------------------

``chemistry_modifier`` takes a dict (or indexed variable) of chemical name to amount
added, on the same basis as composition:

.. code-block:: python

   ReaktoroBlock(
       ...,
       chemistry_modifier={"CaO": m.lime_dose, "Na2CO3": m.soda_ash_dose},
       build_speciation_block=True,
   )

Available by default: ``HCl``, ``H2SO4``, ``CaO``, ``Ca(OH)2``, ``Na2CO3``, ``CO2``,
``O2``, ``NaOH``, ``H``, ``OH``, and ``H2O_evaporation`` (removes water, amount in mol of
H2O). Register other chemicals by their element composition:

.. code-block:: python

   ReaktoroBlock(
       ...,
       chemistry_modifier={"FeCl3": m.coagulant_dose},
       register_new_chemistry_modifiers={"FeCl3": {"Fe": 1, "Cl": 3}},
       build_speciation_block=True,
   )

.. note::

   Chemicals are added as elements. In a block without a speciation block and without
   exact speciation, H and O are open, so the H and O of added chemicals are ignored
   (e.g. only Cl of HCl is added). Use ``build_speciation_block=True`` to add chemicals
   to an equilibrated water.

State changes
-------------

``system_state_modifier`` changes the state in the property block:
``temperature``, ``pressure``, ``enthalpy``, and ``O2_fugacity`` (aeration, see
:doc:`how_to_track_redox`). pH can not be changed directly, dose acid or base instead
(and solve for the dose that gives the target pH, as in
:doc:`/getting_started/first_reaktoro_block`).

.. code-block:: python

   ReaktoroBlock(
       ...,
       system_state={"temperature": m.feed_temperature, "pressure": m.pressure, "pH": m.pH},
       system_state_modifier={"temperature": m.heater_temperature},
       mineral_phase={"phase_components": ["Calcite"]},
       build_speciation_block=True,
   )

Precipitating phases
--------------------

Phases given in ``mineral_phase``, ``solid_phase``, or ``gas_phase`` can form in the
property block. Request their amounts with ``("speciesAmount", "Calcite")`` and use them
in mass balances. Absent phases are returned at a small offset from 0 to keep bounded
solids variables away from their bounds, see
:ref:`how_to_guides/how_to_scale_and_initialize:Species amount offset`.
