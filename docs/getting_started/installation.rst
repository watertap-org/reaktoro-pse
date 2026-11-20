Installation
============

Requirements
------------

* A conda distribution that uses ``conda-forge``, e.g.
  `Miniforge <https://github.com/conda-forge/miniforge?tab=readme-ov-file#download>`_.
  Reaktoro and CyIpopt are distributed through conda-forge.
* Git, needed by `setuptools_scm <https://setuptools-scm.readthedocs.io/en/latest/>`_
  to set the version when installing from source.

Reaktoro-PSE is tested with:

* Python 3.9 through 3.12
* Reaktoro >= 2.12.3
* CyIpopt 1.4.1
* Pyomo >= 6.8.0
* idaes-pse >= 2.5.0

Installing
----------

.. code-block:: shell

   conda create --name reaktoro-pse --yes -c conda-forge python=3.12 cyipopt reaktoro
   conda activate reaktoro-pse
   pip install reaktoro-pse

To use Reaktoro-PSE with WaterTAP models (as in the tutorials), also install WaterTAP:

.. code-block:: shell

   pip install watertap

Solvers
-------

Reaktoro-PSE builds Pyomo gray box models, which are solved with
`CyIpopt <https://cyipopt.readthedocs.io/en/stable/>`_; other solvers have not been tested.
:func:`reaktoro_pse.core.util_classes.cyipopt_solver.get_cyipopt_watertap_solver`
returns CyIpopt with WaterTAP scaling and default options used throughout the examples.

For HSL linear solvers (MA27, MA57, etc., common in WaterTAP and IDAES) beyond the
default MUMPS, follow the `CyIpopt instructions <https://cyipopt.readthedocs.io/en/latest/install.html>`_.

Development installation
------------------------

.. code-block:: shell

   git clone https://github.com/watertap-org/reaktoro-pse.git
   cd reaktoro-pse
   conda create --name reaktoro-pse-dev --yes -c conda-forge python=3.12 cyipopt reaktoro pandoc
   conda activate reaktoro-pse-dev
   pip install -r requirements-dev.txt

Run the tests:

.. code-block:: shell

   pytest --pyargs reaktoro_pse --verbose

Format code before committing:

.. code-block:: shell

   black .

Build these docs (output in ``docs/_build/html``):

.. code-block:: shell

   cd docs
   make html

Set the ``SKIP_APIDOC`` environment variable to skip generating the module reference
for faster builds while writing pages. Code in the docs is tested with ``make doctest``.
