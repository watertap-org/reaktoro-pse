Redox state of water
====================

Why a redox specification is needed
-----------------------------------

An aqueous state is fully defined by temperature, pressure, the total amounts of all
elements, and charge. Lab analyses report the elements other than H and O (as apparent
species, e.g. total Fe as ``Fe``), the amount of water, and pH. That leaves one degree of
freedom open: the balance between H and O beyond water, which is the **redox state**
(equivalently, the amount of electrons available to reduce species).

The redox state can be given as:

* **pE** (or Eh): the electron activity, :math:`pE = -\log a_{e^-}`,
  :math:`Eh = (\ln 10 \, RT / F) \, pE`
* **dissolved O2**: in oxygen controlled waters,
  :math:`O_2 + 4H^+ + 4e^- \rightleftharpoons 2H_2O`, so at fixed pH,
  :math:`pE = 20.78 - pH + \tfrac{1}{4} \log f_{O_2}` (25 C), and dissolved O2 is
  proportional to :math:`f_{O_2}` (Henry's law)
* **O2 fugacity**: equilibrium with a gas, e.g. air at 0.21 bar

Without one, the equilibrium problem has a family of solutions, and Reaktoro returns one
that depends on the solver path: often strongly reduced (nitrate converted to N2, O2 at
~1e-60), which makes results sensitive to ``epsilon`` and inconsistent between solves and
platforms.

Apparent species do not carry the redox state
---------------------------------------------

In flowsheets (e.g. WaterTAP MCAS), each apparent species represents the total amount of
an element: ``Fe_2+`` holds all Fe, whether it is Fe(II) or Fe(III) in the water. The
oxidation state split is not in the composition, it is set by the redox state of the
equilibrium. Assuming that apparent species charges encode the redox state is unsafe, as
they are an analysis convention, not a measurement of oxidation state. The redox state
has to be passed explicitly, as pE or dissolved oxygen, in the same way pH is.

Conservation between blocks
---------------------------

Within a ReaktoroBlock with a speciation block, all elements and charge are passed from
the speciation block to the property block, so the electron balance

.. math::

   Q = Z - b_H + 2 b_O

(charge :math:`Z`, and amounts of H and O, :math:`b_H` and :math:`b_O`) is conserved
(to ~1e-14 mol in practice). Chemical additions in the property block then change pE
and dissolved O2 consistently: e.g. acid at constant dissolved O2 raises pE.

Between blocks, a redox specification is passed, and its choice matters only where the
next block's inputs differ from the outlet state of the previous one:

* **mixing**: dissolved O2 adds linearly between streams, pE does not (and neither does
  pH, so mixed streams need their own equilibrium calculation)
* **units without equilibrium chemistry** (e.g. RO or ion exchange modeled with MCAS):
  dissolved O2 follows the mass balance, a passed through pE is not updated
* **pH set at the next speciation block** instead of by acid or base addition: holding pE
  while changing pH creates or destroys oxidant (:math:`pE + pH` is constant at fixed
  :math:`f_{O_2}`), holding dissolved O2 does not

So dissolved O2 is the better carrier for oxic waters. In anoxic waters (dissolved O2
~0), it is a poor specification: :math:`pE \propto \tfrac{1}{4} \log(DO)`, so tiny changes
in DO move pE by whole units, and it no longer describes the state, which is held by
couples such as Fe(II)/Fe(III), sulfide/sulfate, or ammonium/nitrate. pE is then the
better carrier.

Poorly poised waters
--------------------

When nothing holds the redox state (e.g. dissolved Mn(II) after all O2 and Mn(IV) solids
are gone), pE is set by trace species and is very sensitive to small changes. Equilibrium
calculations are numerically fragile in this regime regardless of which specification is
passed, and the result depends on the trace species considered. Avoid passing such
states between units where possible, e.g. by modeling the aeration step that fixes them,
or by including the redox couple that controls the water in practice.

Aeration
--------

Fixing O2 fugacity models equilibrium with an unlimited air reservoir. Dissolved O2 then
follows Henry's law, :math:`DO \approx K_H(T, P, \text{salinity}) \, f_{O_2}`:

* **temperature**: :math:`K_H` is computed by Reaktoro from the database
* **pH**: does not change O2 saturation (O2 does not hydrolyze), but changes the pE
  implied by air, and how strongly O2 oxidizes Fe and Mn
* **pressure**: :math:`f_{O_2}` is the O2 partial pressure of the gas in contact with water,
  :math:`f_{O_2} \approx 0.2095 \, (P_{air} - p_{H_2O}(T))`, not the liquid pressure
* **salinity**: salting out is only approximately captured by activity models

See :doc:`/technical_reference/databases` for predicted and tabulated saturation, and
:doc:`/how_to_guides/how_to_track_redox` for how to set this up.
