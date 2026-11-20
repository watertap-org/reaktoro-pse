How Reaktoro is embedded in Pyomo
=================================

.. image:: ../tutorials/notebooks/reaktoro_pse_api.png
   :width: 100%
   :alt: Reaktoro-PSE structure

Equilibrium problem
-------------------

Reaktoro finds the equilibrium amounts of all species :math:`n` by minimizing the Gibbs
energy of the system subject to mass and charge conservation:

.. math::

   \min_{n} \; G(n; T, P) \quad \text{s.t.} \quad A n = b, \quad n \ge 0

where :math:`A` is the formula matrix (atoms of each element and charge of each species)
and :math:`b` are the amounts of elements and charge. In flowsheet problems, not all of
:math:`b` is known: lab analyses report apparent species, and pH (rather than the total
amount of H) is measured. Reaktoro handles this with *equilibrium specifications*:
additional constraints (pH, pE, dissolved O2, temperature from enthalpy, ...) and
*open* elements or species (titrants) whose amounts are found to satisfy them.

Reaktoro-PSE builds these specifications from the block options:

* **apparent species** (``exact_speciation=False``): species amounts fix element amounts,
  except H and O (the amount of solvent H2O is fixed instead). pH and a redox
  specification close the two degrees of freedom H and O leave, and the charge balance
  ion (element) is opened to make the solution neutral.
* **exact speciation** (``exact_speciation=True``, or the property block after a
  speciation block): all element amounts, including H and O, are fixed by the inputs,
  and pH and pE are results.
* **chemical additions** add their elements to :math:`b` in the property block.
* **system state modifiers** replace temperature and pressure in the property block, or
  open the system to O2 (``O2_fugacity``).

Gray box model
--------------

The equilibrium problem is solved by Reaktoro, not by the Pyomo solver. Reaktoro-PSE wraps
it in a Pyomo `ExternalGreyBoxModel <https://pyomo.readthedocs.io/en/stable/explanation/modeling/network/greybox.html>`_:

.. math::

   y = f(u)

with inputs :math:`u` (composition, chemical additions, temperature, pressure, pH, ...)
and outputs :math:`y` (requested properties). For each evaluation by Ipopt, the gray box:

#. passes the current inputs to Reaktoro and solves the equilibrium problem,
#. evaluates the outputs from the chemical state,
#. returns the Jacobian :math:`\partial y / \partial u` from Reaktoro sensitivity
   derivatives (exact), chain rule (calculated), or numerical propagation (numeric),
   see :doc:`/technical_reference/outputs`,
#. returns an approximation of the Hessian of the Lagrangian of the gray box
   :math:`\sum_i \lambda_i \nabla^2 y_i`, since second derivatives are not available
   from Reaktoro (``hessian_options``, default SR1 per output).

Inputs and outputs are bound to user variables with equality constraints, so they can be
fixed, bounded, or linked to unit models as any other Pyomo variable. Inputs and outputs
are scaled with the variable scaling factors, and the Jacobian is scaled with
``jacobian_options``.

Speciation and property blocks
------------------------------

With ``build_speciation_block=True`` two equilibrium problems are solved: the speciation
block (apparent species at measured pH) and the property block (exact speciation with
chemical additions and state changes). With direct coupling (default), both are solved
in the same gray box evaluation, the true species amounts are passed from one to the
other in memory, and the Jacobian is combined:

.. math::

   \frac{\partial y}{\partial u} =
   \frac{\partial y}{\partial n_{\text{spec}}}
   \frac{\partial n_{\text{spec}}}{\partial u} +
   \frac{\partial y}{\partial u}\bigg|_{\text{property block}}

where :math:`n_{\text{spec}}` are the true species amounts from the speciation block.
Element amounts and charge (and therefore electrons) are conserved in the hand off, so
the redox state found in the speciation block is carried into the property block.

Reference
---------

Akkor, I., Amusat, O. O., Vecchiarelli, P., Gounaris, C. E., Knueven, B., & Dudchenko, A.
V. "Optimization of Desalination Systems with Detailed Water Chemistry through
Integration of Reaktoro in WaterTAP." ACS ES&T Engineering, vol. 6, no. 5, 2026,
pp. 1464-1477. `doi:10.1021/acsestengg.5c01008 <https://pubs.acs.org/doi/abs/10.1021/acsestengg.5c01008>`_
