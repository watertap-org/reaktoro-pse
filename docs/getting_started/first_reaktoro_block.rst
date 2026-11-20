Your first ReaktoroBlock
========================

This page builds a ``ReaktoroBlock`` that computes the scaling tendency of calcite and
gypsum, the osmotic pressure, and the pH of a brackish water, and then finds the acid
dose that keeps calcite below saturation.

.. testsetup::

   # silence logs (idaes logs to stdout, which doctest compares to expected output),
   # setting logger levels does not work as importing idaes resets them
   import logging
   logging.disable(logging.WARNING)

.. testcleanup::

   logging.disable(logging.NOTSET)

Define the inputs
-----------------

Inputs are ordinary Pyomo variables. Composition is given as amounts of apparent species
(as reported by a lab analysis), on any consistent basis, here mol/s for about 1 kg/s
of water.

.. testcode::

   from pyomo.environ import ConcreteModel, Var, Objective, value, units as pyunits
   import idaes.core.util.scaling as iscale
   from reaktoro_pse.reaktoro_block import ReaktoroBlock
   from reaktoro_pse.core.util_classes.cyipopt_solver import get_cyipopt_watertap_solver

   m = ConcreteModel()
   feed = {"H2O": 55.5, "Na": 0.025, "Cl": 0.025, "Ca": 0.004, "Mg": 0.002,
           "SO4": 0.01, "HCO3": 0.005}
   m.composition = Var(list(feed), initialize=feed, units=pyunits.mol / pyunits.s)
   m.composition.fix()
   m.temperature = Var(initialize=298.15, units=pyunits.K)
   m.temperature.fix()
   m.pressure = Var(initialize=1e5, units=pyunits.Pa)
   m.pressure.fix()
   m.pH = Var(initialize=7.5, units=pyunits.dimensionless)
   m.pH.fix()
   m.acid_dose = Var(initialize=1e-5, units=pyunits.mol / pyunits.s)
   m.acid_dose.fix()

Define the outputs
------------------

Outputs are declared as an indexed variable (or a dict of variables) keyed by
``(property, index)``. Any property Reaktoro reports can be requested, see
:doc:`/technical_reference/outputs`.

.. testcode::

   m.properties = Var(
       [
           ("scalingTendency", "Calcite"),
           ("scalingTendency", "Gypsum"),
           ("osmoticPressure", "H2O"),
           ("pH", None),
       ],
       initialize=1,
   )

Build the block
---------------

The feed is first speciated at the measured pH (speciation block), then HCl is added and
the new equilibrium state is found (property block). Species names are converted to
Reaktoro names (``Ca`` to ``Ca+2``, etc.) for the PHREEQC databases.

.. testcode::

   m.water_chemistry = ReaktoroBlock(
       aqueous_phase={
           "composition": m.composition,
           "convert_to_rkt_species": True,
           "activity_model": "ActivityModelPitzer",
       },
       system_state={
           "temperature": m.temperature,
           "pressure": m.pressure,
           "pH": m.pH,
       },
       database="PhreeqcDatabase",
       database_file="pitzer.dat",
       chemistry_modifier={"HCl": m.acid_dose},
       outputs=m.properties,
       build_speciation_block=True,
   )

Scale, initialize, and solve
----------------------------

Scaling factors on inputs are used to scale the gray box. ``initialize`` solves Reaktoro
at current input values and sets outputs to the results.

.. testcode::

   for key, var in m.composition.items():
       iscale.set_scaling_factor(var, 1 / var.value)
   iscale.set_scaling_factor(m.acid_dose, 1e4)
   m.water_chemistry.initialize()

   solver = get_cyipopt_watertap_solver()
   result = solver.solve(m)
   print(f"pH {value(m.properties['pH', None]):.2f}")
   print(f"Calcite scaling tendency {value(m.properties['scalingTendency', 'Calcite']):.2f}")

.. testoutput::

   pH 7.49
   Calcite scaling tendency 2.07

Optimize
--------

Find the smallest acid dose that keeps calcite at or below saturation (scaling
tendency of 1):

.. testcode::

   m.acid_dose.unfix()
   m.properties["scalingTendency", "Calcite"].setub(1)
   m.objective = Objective(expr=m.acid_dose * 1e4)
   result = solver.solve(m)
   print(f"HCl dose {value(m.acid_dose):.2e} mol/s, pH {value(m.properties['pH', None]):.2f}")

.. testoutput::

   HCl dose 2.72e-04 mol/s, pH 7.19

Next steps
----------

* :doc:`/tutorials/index` walk through these steps in more detail and couple
  ``ReaktoroBlock`` to WaterTAP unit models.
* :doc:`/how_to_guides/index` cover inputs, outputs, chemical additions, redox, scaling,
  and troubleshooting.
* :doc:`/technical_reference/reaktoro_block_options` lists all options.
