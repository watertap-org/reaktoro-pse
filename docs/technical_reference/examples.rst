Examples
========

Examples are Python scripts in
`src/reaktoro_pse/examples <https://github.com/watertap-org/reaktoro-pse/tree/main/src/reaktoro_pse/examples>`_,
each with a ``main()`` function, and are tested in ``examples/tests``. Run one with, e.g.:

.. code-block:: shell

   python -m reaktoro_pse.examples.simple_desalination

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Example
     - Demonstrates
   * - ``simple_desalination.py``
     - basic ReaktoroBlock, scaling tendency, pH, and osmotic pressure outputs, and
       optimizing acid dose for operation at a target scaling tendency
   * - ``thermal_precipitation.py``
     - using a different database, enthalpy and vapor pressure outputs, precipitation of
       calcite, and energy needed to remove calcite at different temperatures
   * - ``simple_ion_exchange.py``
     - charge balancing the feed, passing outputs of a speciation block into a second
       property block, ion exchange phases, and optimizing acid and base addition for
       calcium over magnesium selectivity
   * - ``biogas_combustion.py``
     - custom database, condensed phase, exhaust temperature from combustion, and the air
       to fuel ratio needed to reach 2000 K (based on the
       `Reaktoro biomass gasification example <https://reaktoro.org/applications/biomass-gasification/biomass-gasification.html>`_)
   * - ``redox_tracking.py``
     - why pE has to be passed between blocks for redox sensitive species in databases
       such as ``wateq4f.dat``
   * - ``oxygen_tracking.py``
     - defining the redox state with measured dissolved oxygen, equilibrating lime
       softening with air (``O2_fugacity``), and passing dissolved oxygen to the next
       step; shows how oxygen limits Mn removal
   * - ``evaporative_crystallizer.py``
     - why the species amount offset is needed when absent minerals feed bounded
       variables, optimizing evaporation of seawater with lime and soda ash dosing

Comparisons with PHREEQC
------------------------

Scripts in ``examples/reaktoro_pse_to_phreeqc_comparison`` compare Reaktoro-PSE to
PHREEQC (``phreeqcinwt``) for common calculations, using ``pitzer.dat`` and
``wateq4f.dat``. They also show how to set up each type of calculation. Run ``main()``
with ``save_fig=True`` to plot results.

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Script
     - Calculation
   * - ``water_removal_comparison.py``
     - removing water from a solution (evaporative processes)
   * - ``vapor_pressure_comparison.py``
     - vapor pressure at different temperatures
   * - ``precipitation_comparison.py``
     - precipitation of mineral phases
   * - ``acid_base_addition_comparison.py``
     - addition of HCl and NaOH
   * - ``solution_mixing_comparison.py``
     - mixing of two solutions
