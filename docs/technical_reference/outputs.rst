Outputs and derivatives
=======================

Outputs are requested as ``(property, index)`` keys (see
:doc:`/how_to_guides/how_to_request_outputs`). For each key, Reaktoro-PSE looks for the
property in this order, and uses the first match:

#. **Converted properties** (Reaktoro-PSE): computed from exact Reaktoro derivatives with
   the chain rule.
#. **Pyomo properties** (Reaktoro-PSE): built as Pyomo constraints from outputs with
   exact derivatives.
#. **Chemical properties**: any single valued method of Reaktoro
   `ChemicalProps <https://reaktoro.org/api/classReaktoro_1_1ChemicalProps.html>`_.
#. **Aqueous properties**: any single valued method of Reaktoro
   `AqueousProps <https://reaktoro.org/api/classReaktoro_1_1AqueousProps.html>`_
   (only when an aqueous phase is defined).

Only properties that return a single real value are supported, as each output is matched
to one Pyomo variable. Methods that return arrays or strings are not supported.

Converted properties
--------------------

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Property
     - Index
     - Description
   * - ``pH``
     - ``None``
     - solution pH
   * - ``scalingTendency``
     - phase (e.g. ``Calcite``)
     - scaling tendency (saturation ratio, ion activity product / K), computed from
       species activities and the phase equilibrium constant with exact derivatives
       (PHREEQC databases)
   * - ``scalingTendencySaturationIndex``
     - phase
     - scaling tendency as 10^saturation index from Reaktoro ``saturationIndex``
       (numeric derivatives), works with any database
   * - ``osmoticPressure``
     - solvent species (``H2O``)
     - osmotic pressure (Pa)
   * - ``vaporPressure``
     - species (e.g. ``H2O(g)``)
     - vapor pressure (Pa), requires the gas species in ``gas_phase``
   * - ``elementAmount``
     - element
     - amount of element in the system (mol)
   * - ``charge``
     - ``None``
     - solution charge
   * - ``alkalinityAsCaCO3``
     - ``None``
     - alkalinity as CaCO3 (mg/L)

Pyomo properties
----------------

Alternative implementations built as Pyomo constraints from exact outputs:
``scalingTendencyPyomo``, ``osmoticPressurePyomo``, ``vaporPressurePyomo``. They are
mostly useful for debugging and comparison with converted properties.

Derivative types
----------------

Each output gets one of these derivative types (shown by
``m.block.display_jacobian_outputs()``):

* **exact**: Reaktoro provides the derivative with respect to the equilibrium inputs.
* **calculated**: derivative computed by chain rule from exact derivatives (converted
  properties).
* **numeric**: derivative propagated numerically from the chemical state, controlled by
  ``jacobian_options`` (``numerical_type``, ``numerical_order``, ``numerical_step``).

Properties with exact derivatives, which can also be used to build custom properties
with exact derivatives:

* ``temperature``, ``pressure`` (system and per phase)
* ``speciesAmount``, ``speciesMoleFraction`` (per species)
* ``amount``, ``mass`` (per phase)
* ``speciesStandardGibbsEnergy``, ``speciesStandardEnthalpy``, ``speciesStandardVolume``,
  ``speciesStandardVolumeT``, ``speciesStandardVolumeP``,
  ``speciesStandardHeatCapacityConstP`` (per species)
* ``correctiveMolarVolume``, ``correctiveMolarVolumeT``, ``correctiveMolarVolumeP``,
  ``correctiveMolarGibbsEnergy``, ``correctiveMolarEnthalpy``,
  ``correctiveMolarHeatCapacityConstP`` (per phase)
* ``speciesActivityCoefficientLn``, ``speciesActivityLn``, ``speciesChemicalPotential``
  (per species)

See :class:`reaktoro_pse.core.reaktoro_outputs.ConvertedPropTypes` and
:class:`reaktoro_pse.core.reaktoro_outputs.PyomoProperties` for implementations.
