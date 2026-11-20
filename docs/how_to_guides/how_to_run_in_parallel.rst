How to solve ReaktoroBlocks in parallel
=======================================

Flowsheets with many ReaktoroBlocks (e.g. indexed blocks along an RO module, or multiple
treatment steps) spend most of their solve time in Reaktoro. ``ReaktoroBlockManager``
collects all blocks into a single gray box and solves their Reaktoro problems in
parallel worker processes.

.. code-block:: python

   from reaktoro_pse.reaktoro_block import ReaktoroBlock
   from reaktoro_pse.parallel_tools.reaktoro_block_manager import ReaktoroBlockManager

   m.reaktoro_manager = ReaktoroBlockManager()
   m.softening = ReaktoroBlock(..., reaktoro_block_manager=m.reaktoro_manager)
   m.desalination = ReaktoroBlock(..., reaktoro_block_manager=m.reaktoro_manager)

   # after all ReaktoroBlocks are built, before initializing or using them
   m.reaktoro_manager.build_reaktoro_blocks()

   m.softening.initialize()
   ...
   solver.solve(m)

   # shut down worker processes when done
   m.reaktoro_manager.terminate_workers()

Options (see :doc:`/technical_reference/reaktoro_block_options`):

* ``use_parallel_mode`` (default ``True``): set to ``False`` to solve blocks in series in
  the main process, with the same single gray box.
* ``maximum_number_of_parallel_solves``: limits how many blocks are solved at once.
* ``worker_timeout``: workers shut down after this many seconds without commands.
* ``hessian_options``: Hessian approximation for the combined gray box, which replaces
  the options of individual blocks.

On Windows and macOS, workers are started with ``spawn``, so scripts that use the manager
must protect their entry point with ``if __name__ == "__main__":``.
