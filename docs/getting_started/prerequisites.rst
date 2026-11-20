Choosing databases and activity models
======================================

.. important::

   Familiarize yourself with Reaktoro options, especially database, database file, and
   activity models for aqueous, solid, and gas phases. Reaktoro-PSE does not select
   these for you, and uses ideal activity models by default, which are not accurate for
   most real waters.

All Reaktoro databases and activity models are supported, see the Reaktoro tutorials on
`loading databases <https://reaktoro.org/tutorials/basics/loading-databases.html>`_ and
`specifying activity models <https://reaktoro.org/tutorials/basics/specifying-activity-models.html>`_.
Activity models can be given as a name, an initialized Reaktoro activity model, a list,
or a chain (``reaktoro.chain(...)``).

Defaults used by ``ReaktoroBlock``:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Option
     - Default
   * - Database
     - ``PhreeqcDatabase``
   * - Database file
     - ``pitzer.dat``
   * - Aqueous activity model
     - ``ActivityModelIdealAqueous``
   * - Gas activity model
     - ``ActivityModelIdealGas``
   * - Solid activity model
     - ``ActivityModelIdealSolution``
   * - Ion exchange activity model
     - ``ActivityModelIonExchange``

Common choices for water treatment:

* ``pitzer.dat`` with ``ActivityModelPitzer`` for brines and high salinity waters.
  ``pitzer.dat`` has no redox chemistry (see :doc:`/technical_reference/databases`).
* ``phreeqc.dat`` or ``wateq4f.dat`` with ``ActivityModelPhreeqc`` for fresh and brackish
  waters, including redox sensitive species (Fe, Mn, S, N, O2).

See :doc:`/technical_reference/databases` for what each PHREEQC database supports.
