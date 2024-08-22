# reaktoro-pse introduction
### Overview
This is a package for configuring [reaktoro](https://reaktoro.org/index.html) as a gray box model in [pyomo](https://pyomo.readthedocs.io/en/stable/), [idaes-pse](https://idaes-pse.readthedocs.io/en/stable/), and [watertap](https://watertap.readthedocs.io/en/stable/) modeling libraries. This package is not meant to replace or act as a higher level API for reaktoro - it is only meant to enable setting up reaktoro equilibrium problems as blocks on pyomo models and automate transferring reaktoro data into pyomo variables. 

### Prerequisites
**The user must familiarize thyself with reaktoro options and usage, especially when selecting databases, database files, and activity models for aqueous, solid, and gas phases. This package does not automatically select any of these options, and will use ideal models by default.**

* Please refer here for information on [databases](https://reaktoro.org/tutorials/basics/loading-databases.html) (all are supported) 
* Please refer here for information on [activity models](https://reaktoro.org/tutorials/basics/specifying-activity-models.html) (all are supported, included chain operations, or passing in pre-configured activity models)

**By default the package will use:**

* **Database:** [PhreeqcDatabase](https://reaktoro.org/api/classReaktoro_1_1PhreeqcDatabase.html) 
* **Data file:** [pitzer.dat](https://reaktoro.org/api/classReaktoro_1_1PhreeqcDatabase.html) 
* **Aqueous activity mode:** [ActivityModelIdealAqueous](https://reaktoro.org/api/namespaceReaktoro.html#ae431d4c8a1f283910ae1cf35024091b8)
* **Gas activity mode:** [ActivityModelIdealGas](https://reaktoro.org/api/namespaceReaktoro.html#a7a0788a5a863d987a88b81303d80b427)
* **Solid activity mode:** [ActivityModelIdealSolution](https://reaktoro.org/api/namespaceReaktoro.html#a6581d5c0cde36cae6d9c46dbc32d56f8)
* **Ion exchange activity model:** [ActivityModelIonExchange](https://reaktoro.org/api/namespaceReaktoro.html#a6581d5c0cde36cae6d9c46dbc32d56f8)

### Inputs and outputs of the reaktoro blocks
The reaktoro blocks built by this package are designed to solve an equilibrium problem using user provided apparent species, temperature, pressure, and pH, which are broken down to base elements and equilibrated within reaktoro to provide exact speciation and equilibrium state. Using this state the block can return various information supported by reaktoro:

* [Chemical properties](https://reaktoro.org/api/classReaktoro_1_1ChemicalProps.html)
* [Aqueous properties](https://reaktoro.org/api/classReaktoro_1_1AqueousProps.html)
* Pyomo build properties, which are custom properties built in pyomo that use chemical properties or aqueous properties as inputs 

Note: Only properties that return a single values are supported, as they will be directly matched to a pyomo var, as such any functions that return arrays or strings are not supported. 

### Solvers 
To-date, this has only been tested with [cyipopt](https://cyipopt.readthedocs.io/en/stable/) - other solvers have not been tested. 
* For accesing MA27 solver (common to WaterTAP and IDAES-PSE) follow [these instructions](https://cyipopt.readthedocs.io/en/latest/install.html) for adding it to cyipopt

### Requesting new features or issues
Please include a minimal example using reaktoro for your specific feature or issue request if possible . 

## Compatibility

reaktoro-pse depends on the following packages and/or versions:

- Python 3.9 through 3.12
- Reaktoro 2.12
- CyIpopt 1.4.1
- Pyomo>=6.8.0
- idaes-pse>=2.5.0
- watertap>=1.0.0 - (required for watertap-cyipopt wrapper only)

## Getting started (for contributors)

### Prerequisites

- A Conda distribution compatible with `conda-forge`, e.g. [Miniforge](https://github.com/conda-forge/miniforge?tab=readme-ov-file#download)
- Git (needed by [setuptools_scm](https://setuptools-scm.readthedocs.io/en/latest/) to set the version dynamically during installation of the Python package distribution)

### Installation

```sh
git clone https://github.com/watertap-org/reaktoro-pse.git
cd reaktoro-pse
conda create --yes -c conda-forge --name reaktoro-pse-dev python=3.11 reaktoro=2.12.1 cyipopt=1.4.1
conda activate reaktoro-pse-dev
pip install -r requirements-dev.txt
```

### Running tests

```sh
conda activate reaktoro-pse-dev
pytest --pyargs reaktoro_pse --verbose
```

### Before committing

```sh
black .
```
