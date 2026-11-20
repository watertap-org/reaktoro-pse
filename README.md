
<img src="docs\images\reaktoro_pse.png" width="400"/>

## Overview
This is a package for configuring [Reaktoro](https://reaktoro.org/index.html) as a gray box model in [Pyomo](https://pyomo.readthedocs.io/en/stable/), [IDAES-PSE](https://idaes-pse.readthedocs.io/en/stable/), and [WaterTAP](https://watertap.readthedocs.io/en/stable/) modeling libraries. This package is not meant to replace or act as a higher level API for Reaktoro - it is only meant to enable setting up Reaktoro equilibrium problems as blocks on Pyomo models and automate transferring Reaktoro data into Pyomo variables.

**Documentation: [reaktoro-pse.readthedocs.io](https://reaktoro-pse.readthedocs.io/)** - installation, getting started, tutorials, how-to guides (inputs, outputs, chemical additions, redox tracking, scaling, parallel solves, troubleshooting), and the full options reference.

If you find this tool useful for your research, please consider citing the [paper](https://pubs.acs.org/doi/abs/10.1021/acsestengg.5c01008):
> Akkor, I., Amusat, O. O., Vecchiarelli, P., Gounaris, C. E., Knueven, B., & Dudchenko, A. V. “Optimization of Desalination Systems with Detailed Water Chemistry through Integration of Reaktoro in WaterTAP.” ACS ES&T Engineering, vol. 6, no. 5, 2026, pp. 1464–1477.

## Installation
Reaktoro and CyIpopt are distributed through conda-forge, use a conda distribution such as [Miniforge](https://github.com/conda-forge/miniforge?tab=readme-ov-file#download):

```sh
conda create --name reaktoro-pse --yes -c conda-forge python=3.12 cyipopt reaktoro
conda activate reaktoro-pse
pip install reaktoro-pse
```

**Familiarize yourself with Reaktoro options, especially databases, database files, and activity models. Reaktoro-PSE does not select these for you, and uses ideal activity models by default.** See [Choosing databases and activity models](https://reaktoro-pse.readthedocs.io/en/latest/getting_started/prerequisites.html).

## Tutorials and examples
* [Tutorials](https://github.com/watertap-org/reaktoro-pse/tree/main/src/reaktoro_pse/tutorials) (Jupyter notebooks): interacting with ReaktoroBlock, integration with WaterTAP RO and NF models, and a softening, acid addition, and desalination treatment train.
* [Examples](https://github.com/watertap-org/reaktoro-pse/tree/main/src/reaktoro_pse/examples): desalination, thermal precipitation, ion exchange, biogas combustion, redox and oxygen tracking, evaporative crystallization, and comparisons with PHREEQC.

## Compatibility
- Python 3.9 through 3.12
- Reaktoro>=2.12.3
- CyIpopt 1.4.1
- Pyomo>=6.8.0
- idaes-pse>=2.5.0

## Requesting new features or issues
Please include a minimal example using Reaktoro for your specific feature or issue request if possible. Please also include full traceback for errors that occur.

## For contributors
```sh
git clone https://github.com/watertap-org/reaktoro-pse.git
cd reaktoro-pse
conda create --name reaktoro-pse-dev --yes -c conda-forge python=3.12 cyipopt reaktoro pandoc
conda activate reaktoro-pse-dev
pip install -r requirements-dev.txt
```

Run tests, format code, and build docs:

```sh
pytest --pyargs reaktoro_pse --verbose
black .
make -C docs html   # output in docs/_build/html
```
