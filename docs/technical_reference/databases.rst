Databases
=========

Reaktoro-PSE works with any Reaktoro database (``database`` option, as a name with
``database_file``, or an initialized Reaktoro database object). The species name
converter (``convert_to_rkt_species``) supports PHREEQC databases; for others provide
``species_to_rkt_species_dict``.

PHREEQC databases
-----------------

PHREEQC databases distributed with Reaktoro, and the redox species they contain
(Reaktoro 2.13). A database supports redox calculations when elements exist in more than
one oxidation state; ``dissolved_O2`` and ``O2_fugacity`` also need O2(aq).

.. list-table::
   :header-rows: 1
   :widths: 25 12 13 50

   * - Database file
     - Aqueous species
     - O2(aq), H2(aq)
     - Redox couples
   * - ``phreeqc.dat``
     - 250
     - yes
     - Fe(II)/Fe(III), Mn, S(-II)/S(VI), N(-III)/N(0)/N(V), C(-IV)/C(IV)
   * - ``wateq4f.dat``
     - 378
     - yes
     - Fe, Mn (incl. MnO4-), S, N, C
   * - ``llnl.dat``
     - 1344
     - yes
     - Fe, Mn, S, N, C
   * - ``core10.dat``, ``Tipping_Hurley.dat``
     - 227, 378
     - yes
     - Fe, Mn, S, N, C
   * - ``minteq.dat``, ``minteq.v4.dat``
     - 1027, 1332
     - yes
     - Fe, Mn, S, N (no N2(aq))
   * - ``sit.dat``
     - 1424
     - yes
     - Fe, Mn, S, N (no N2(aq)), C
   * - ``iso.dat``
     - 2316
     - yes
     - Fe, S, N, C (no Mn)
   * - ``Amm.dat``
     - 250
     - yes
     - Fe, Mn, S, N (no NH4+), C
   * - ``pitzer.dat``
     - 48
     - **no**
     - **none**: Fe(II) and Mn(II) only, each element in one oxidation state
   * - ``frezchem.dat``, ``ColdChem.dat``
     - 21, 10
     - no
     - none

``pitzer.dat`` is the most accurate choice for brines (with ``ActivityModelPitzer``),
but has no redox chemistry, so pE has no effect and ``dissolved_O2`` and ``O2_fugacity``
can not be used. For redox sensitive waters use ``phreeqc.dat`` or ``wateq4f.dat`` with
``ActivityModelPhreeqc``.

Dissolved oxygen at saturation
------------------------------

Dissolved oxygen (mg/L) in equilibrium with air at 1 atm, computed with ``phreeqc.dat`` or
``wateq4f.dat`` and ``ActivityModelPhreeqc``, compared to standard tables for fresh water
in moist air (APHA / Benson and Krause, approximate values):

.. list-table::
   :header-rows: 1

   * - T (C)
     - f(O2) = 0.2123 bar (dry air)
     - f(O2) corrected for water vapor
     - Standard table
   * - 0
     - 14.91
     - 14.82
     - 14.6
   * - 10
     - 11.65
     - 11.51
     - 11.3
   * - 20
     - 9.50
     - 9.28
     - 9.09
   * - 30
     - 8.04
     - 7.70
     - 7.56
   * - 40
     - 7.04
     - 6.53
     - 6.41
   * - 60
     - 5.86
     - 4.71
     - ~4.7
   * - 20, 0.6 M NaCl (~seawater)
     - 8.27
     - 8.08
     - ~7.4

With the water vapor correction, Reaktoro is within ~2% of tabulated values from 0 to
60 C. A fixed f(O2) of 0.21 bar overestimates dissolved oxygen in warm water (by ~7% at
40 C, ~20% at 60 C), and the activity model underestimates salting out (~9% high at
seawater salinity, more for brines).
