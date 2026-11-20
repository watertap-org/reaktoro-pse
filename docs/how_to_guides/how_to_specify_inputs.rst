How to specify inputs
=====================

A ``ReaktoroBlock`` needs a composition for at least one phase and a system state
(temperature, pressure, and usually pH). All inputs are Pyomo variables (fixed or
unfixed), so they can be decision variables or results of other unit models.

Composition
-----------

Composition is given per phase (``aqueous_phase``, ``liquid_phase``, ``gas_phase``,
``condensed_phase``, ``mineral_phase``, ``solid_phase``, ``ion_exchange_phase``) as an
indexed variable or a dict of variables, keyed by species name:

.. code-block:: python

   ReaktoroBlock(
       aqueous_phase={
           "composition": m.composition,       # e.g. m.composition["Ca"], ...
           "convert_to_rkt_species": True,     # Ca -> Ca+2, HCO3 -> HCO3-, ...
           "activity_model": "ActivityModelPitzer",
       },
       ...
   )

* **Units**: composition can be in mol or mass (e.g. ``mol/s``, ``kg/s``, ``mol``), on any
  consistent basis; all species are converted to mol for Reaktoro. Concentrations
  (``mg/L``, ``mol/m**3``) are not supported and raise an error, convert them to flows
  or amounts first.
* **Species names**: with ``convert_to_rkt_species`` the default converter maps common
  ions to PHREEQC database names. Provide ``species_to_rkt_species_dict``
  (``{"Ca": "Ca+2", ...}``) for other names or databases.
* **Elements**: set ``composition_is_elements=True`` if composition is element amounts.
* **Phases**: give only ``phase_components`` (e.g. a list of minerals) for phases that can
  form but are not in the inlet, e.g.
  ``mineral_phase={"phase_components": ["Calcite", "Gypsum"]}``.

Apparent or exact speciation
----------------------------

By default (``exact_speciation=False``) the composition is treated as **apparent
species** from a lab analysis: each species only fixes the total amount of its elements
(``Fe`` fixes total Fe, regardless of Fe(II)/Fe(III) split), the solvent (H2O) amount
is fixed, and H and O are left open. The equilibrium state is then defined by pH (and pE
or dissolved oxygen for redox sensitive waters, see :doc:`how_to_track_redox`).

Set ``exact_speciation=True`` when the composition is the true equilibrium speciation,
e.g. the ``speciesAmount`` outputs of another ReaktoroBlock. All elements, including H
and O, are then fixed, and pH is a result of the equilibrium (it should not be given).

System state
------------

``system_state`` takes ``temperature`` (K), ``pressure`` (Pa), ``pH``, and optionally a
redox specification (``pE``, ``Eh``, ``dissolved_O2``, or ``O2_fugacity``, only one) and
``enthalpy`` (to solve for temperature, e.g. in combustion or heating).

.. code-block:: python

   system_state={
       "temperature": m.temperature,
       "pressure": m.pressure,
       "pH": m.pH,
   }

Charge balance
--------------

Lab analyses are rarely charge balanced. With ``assert_charge_neutrality=True`` (default)
the block adjusts ``charge_neutrality_ion`` (default ``Cl``) until the solution is neutral:

* an **element** (e.g. ``Cl``): the amount of Cl (and all Cl species) is adjusted. This is
  the most stable choice when exact speciation is not provided.
* an **ion** (e.g. ``Cl-``): species ratios are shifted without changing element amounts.
  Use only with exact speciation.
* ``pH``: H+ is adjusted without changing total H.

To get the charge balanced composition, request ``("elementAmount", "Cl")`` as an output.
The property block built after a speciation block is charge balanced on pH, unless
``assert_charge_neutrality_on_property_block=False``.

Indexed blocks
--------------

An indexed ReaktoroBlock (e.g. ``ReaktoroBlock(m.fs.time, ...)`` or along the length of an
RO module) passes each index of indexed inputs to the matching block. Set
``<input>_indexed=False`` (e.g. ``composition_indexed``, ``temperature_indexed``) to pass the
same variable to all blocks.
