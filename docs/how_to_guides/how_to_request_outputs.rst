How to request outputs
======================

Outputs are given with the ``outputs`` option as an indexed variable or a dict, keyed
by ``(property, index)``, where index is the species, phase, or element the property
applies to, or ``None``:

.. code-block:: python

   m.properties = Var(
       [
           ("scalingTendency", "Calcite"),   # phase
           ("speciesAmount", "Calcite"),     # species
           ("elementAmount", "Cl"),          # element
           ("osmoticPressure", "H2O"),
           ("pH", None),
       ],
       initialize=1,
   )
   ReaktoroBlock(..., outputs=m.properties)

* Existing variables are used as outputs and are bound to the gray box with equality
  constraints, so outputs can be variables of other unit models.
* A dict can mix variables from different parts of the model,
  ``{("pH", None): m.unit.pH, ("osmoticPressure", "H2O"): m.unit.pressure_osm}``.
  Use ``None`` as a value to have the block create the output variable on
  ``m.block.outputs``.
* ``{"speciesAmount": True}`` returns amounts of all species in the system, useful to
  pass exact speciation to another block without knowing the species beforehand.

Available properties
--------------------

Any Reaktoro property that returns a single value can be used, from
`ChemicalProps <https://reaktoro.org/api/classReaktoro_1_1ChemicalProps.html>`_
(e.g. ``speciesAmount``, ``speciesActivityCoefficient``, ``molarEnthalpy``) and
`AqueousProps <https://reaktoro.org/api/classReaktoro_1_1AqueousProps.html>`_
(e.g. ``pH``, ``pE``, ``ionicStrength``, ``saturationIndex``). Some properties are
computed by Reaktoro-PSE with exact derivatives. See :doc:`/technical_reference/outputs`
for the full list and how derivatives are obtained.

Checking derivative types
-------------------------

Derivatives are exact for properties Reaktoro differentiates, computed by chain rule for
converted properties, and numerical otherwise. Numerical derivatives are less accurate
and can slow convergence. To see which outputs use which:

.. code-block:: python

   m.block.display_jacobian_outputs()
