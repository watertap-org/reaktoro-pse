Tutorials
=========

The tutorials are Jupyter notebooks in
`src/reaktoro_pse/tutorials <https://github.com/watertap-org/reaktoro-pse/tree/main/src/reaktoro_pse/tutorials>`_.
They are shown here with their stored outputs, download them from GitHub to run them.
Tutorials 2 to 4 need WaterTAP (``pip install watertap``).

#. **Interacting with ReaktoroBlock**: getting properties, charge balancing the feed,
   providing exact speciation, chemical addition with a speciation block, and the
   effect of concentrating a feed on scaling tendency.
#. **Integration with reverse osmosis**: indexed ReaktoroBlocks on the
   :ref:`WaterTAP RO model <watertap:technical_reference/unit_models/reverse_osmosis_1D:Reverse Osmosis (1D)>`
   for osmotic pressure and scaling tendencies along the membrane.
#. **Integration with nanofiltration**: ReaktoroBlocks on a WaterTAP NF zero order model.
#. **Softening, acid addition, and desalination**: a treatment train with lime and soda
   ash softening, acid addition, and a concentrator, optimized to control scaling at
   minimum chemical cost.

.. toctree::
   :maxdepth: 1
   :caption: Notebooks

   notebooks/basic_reaktoro_block_interaction
   notebooks/integration_with_ro
   notebooks/integration_with_nf
   notebooks/softening_acid_ro_example
