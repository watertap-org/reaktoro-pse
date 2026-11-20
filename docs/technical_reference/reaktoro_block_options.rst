ReaktoroBlock options
=====================

.. code-block:: python

   from reaktoro_pse.reaktoro_block import ReaktoroBlock

All options are passed as keyword arguments when building the block. Options on this
page are generated from the code (``ReaktoroBlock.CONFIG``).

Inputs marked with an *indexed option* (e.g. ``temperature_indexed``) define how the
input is used by an indexed ReaktoroBlock: if ``True`` (default) the input has the same
indexing as the block and each index is passed to the matching block, if ``False`` all
blocks get the same input.

.. contents:: On this page
   :local:
   :depth: 1

Phases
------

Each phase (``aqueous_phase``, ``liquid_phase``, ``gas_phase``, ``condensed_phase``,
``mineral_phase``, ``solid_phase``, ``ion_exchange_phase``) is configured with a dict of
the options below. ``aqueous_phase`` and ``liquid_phase`` also accept
``fixed_solvent_specie`` and ``free_element``, and ``gas_phase`` accepts
``speciate_phase_component``.

.. code-block:: python

   ReaktoroBlock(
       aqueous_phase={
           "composition": m.composition,
           "convert_to_rkt_species": True,
           "activity_model": "ActivityModelPitzer",
       },
       mineral_phase={"phase_components": ["Calcite", "Gypsum"]},
       ...
   )

.. rkt-config:: reaktoro_pse.reaktoro_block.ReaktoroBlockData.CONFIG
   :path: aqueous_phase

``gas_phase`` only:

.. rkt-config:: reaktoro_pse.reaktoro_block.ReaktoroBlockData.CONFIG
   :path: gas_phase
   :exclude: composition, composition_indexed, phase_components, convert_to_rkt_species, species_to_rkt_species_dict, composition_is_elements, activity_model

System state
------------

``system_state`` defines the state of the inlet (or the only) equilibrium problem.

.. code-block:: python

   system_state={"temperature": m.temperature, "pressure": m.pressure, "pH": m.pH}

.. rkt-config:: reaktoro_pse.reaktoro_block.ReaktoroBlockData.CONFIG
   :path: system_state

System state modifier
---------------------

``system_state_modifier`` changes the state in the property block, when it is built after
a speciation block (``build_speciation_block=True``). It accepts the same options as
``system_state`` except pH and the redox specifications (pE, Eh, dissolved_O2), which are
results of the property block equilibrium. ``O2_fugacity`` here models aeration.

Block options
-------------

.. rkt-config:: reaktoro_pse.reaktoro_block.ReaktoroBlockData.CONFIG
   :values-only:

Jacobian options
----------------

``jacobian_options`` configure numerical derivatives and scaling of the gray box
Jacobian.

.. rkt-config:: reaktoro_pse.reaktoro_block.ReaktoroBlockData.CONFIG
   :path: jacobian_options

Hessian options
---------------

``hessian_options`` configure the Hessian approximation of the gray box.

.. rkt-config:: reaktoro_pse.reaktoro_block.ReaktoroBlockData.CONFIG
   :path: hessian_options

Reaktoro solver options
-----------------------

``reaktoro_solve_options`` configure the Reaktoro equilibrium solver and the species
amount offset.

.. rkt-config:: reaktoro_pse.reaktoro_block.ReaktoroBlockData.CONFIG
   :path: reaktoro_solve_options

Reaktoro presolve options
-------------------------

``reaktoro_presolve_options`` configure an optional low tolerance solve before the main
Reaktoro solve.

.. rkt-config:: reaktoro_pse.reaktoro_block.ReaktoroBlockData.CONFIG
   :path: reaktoro_presolve_options

ReaktoroBlockManager options
----------------------------

.. code-block:: python

   from reaktoro_pse.parallel_tools.reaktoro_block_manager import ReaktoroBlockManager

See :doc:`/how_to_guides/how_to_run_in_parallel`. ``hessian_options`` are the same as
above.

.. rkt-config:: reaktoro_pse.parallel_tools.reaktoro_block_manager.ReaktoroBlockManagerData.CONFIG
   :values-only:

ReaktoroBlock methods
---------------------

.. autoclass:: reaktoro_pse.reaktoro_block.ReaktoroBlockData
   :members: initialize, update_block_scaling, update_jacobian_scaling, check_species_amount_offsets, get_species_amount_offset_report, get_species_amount_offset_element_report, display_jacobian_outputs, display_jacobian_scaling, display_reaktoro_state, fix_outputs, unfix_outputs, activate, deactivate
   :noindex:
