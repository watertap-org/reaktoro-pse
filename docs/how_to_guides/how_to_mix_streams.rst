How to mix streams
==================

Mixing waters at different pH can not be done by adding apparent compositions, as each
water has to be speciated at its own pH first. ReaktoroBlock does this by speciating each
stream separately and passing all true species into a single property block.

#. Build a ReaktoroBlock for each additional stream with ``build_graybox_model=False``
   and ``outputs={"speciesAmount": True}``. This configures the speciation problem
   without building a gray box.
#. Build the mixing ReaktoroBlock with ``build_speciation_block=True`` for its own feed,
   and pass the other streams with ``external_speciation_reaktoro_blocks``. All speciation
   problems and the mixed property block are solved in one gray box.

.. code-block:: python

   m.stream_2 = ReaktoroBlock(
       aqueous_phase={"composition": m.composition_2, "convert_to_rkt_species": True},
       system_state={"temperature": m.temperature, "pressure": m.pressure,
                     "pH": m.pH_2},
       outputs={"speciesAmount": True},
       build_graybox_model=False,
   )
   m.mixed = ReaktoroBlock(
       aqueous_phase={"composition": m.composition_1, "convert_to_rkt_species": True},
       system_state={"temperature": m.temperature, "pressure": m.pressure,
                     "pH": m.pH_1},
       outputs=m.mixed_properties,
       build_speciation_block=True,
       external_speciation_reaktoro_blocks=[m.stream_2],
   )

Use the same database and activity models for all blocks. For redox sensitive waters
give each stream its own redox specification (see :doc:`how_to_track_redox`).

Examples: ``examples/reaktoro_pse_to_phreeqc_comparison/solution_mixing_comparison.py``
(compared to PHREEQC) and ``examples/simple_ion_exchange.py``.
