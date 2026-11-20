Reaktoro-PSE
============

.. image:: images/reaktoro_pse.png
   :width: 350
   :alt: Reaktoro-PSE

Reaktoro-PSE configures `Reaktoro <https://reaktoro.org>`_ chemical equilibrium
calculations as gray box models in `Pyomo <https://pyomo.readthedocs.io/en/stable/>`_,
`IDAES <https://idaes-pse.readthedocs.io/en/stable/>`_, and
`WaterTAP <https://watertap.readthedocs.io/en/stable/>`_ flowsheets.

A ``ReaktoroBlock`` takes Pyomo variables for composition, temperature, pressure,
pH (and optionally pE or dissolved oxygen, and chemical additions), solves the
equilibrium problem in Reaktoro, and returns properties such as pH, scaling tendencies,
osmotic pressure, vapor pressure, enthalpy, or amounts of precipitated minerals,
together with their derivatives. This lets equation oriented solvers (e.g. Ipopt)
optimize flowsheets with detailed water chemistry.

Reaktoro-PSE is not a higher level API for Reaktoro: it sets up Reaktoro equilibrium
problems as blocks on Pyomo models and transfers Reaktoro data into Pyomo variables.

Cite this work
--------------

If you find this tool useful for your research, please cite the
`paper <https://pubs.acs.org/doi/abs/10.1021/acsestengg.5c01008>`_:

   Akkor, I., Amusat, O. O., Vecchiarelli, P., Gounaris, C. E., Knueven, B., &
   Dudchenko, A. V. "Optimization of Desalination Systems with Detailed Water Chemistry
   through Integration of Reaktoro in WaterTAP." ACS ES&T Engineering, vol. 6, no. 5,
   2026, pp. 1464-1477.

Contents
--------

.. toctree::
   :maxdepth: 2

   getting_started/index
   tutorials/index
   how_to_guides/index
   technical_reference/index
   background/index
   license
   copyright

Reporting issues and requesting features
----------------------------------------

Please open an issue on `GitHub <https://github.com/watertap-org/reaktoro-pse/issues>`_.
Include a minimal example that reproduces the issue (or demonstrates the feature) and,
for errors, the full traceback.

Indices
-------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
