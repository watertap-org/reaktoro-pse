# reaktoro-pse introduction
## 1. Overview
This is a package for configuring [Reaktoro](https://reaktoro.org/index.html) as a gray box model in [Pyomo](https://pyomo.readthedocs.io/en/stable/), [IDAES-PSE](https://idaes-pse.readthedocs.io/en/stable/), and [WaterTAP](https://watertap.readthedocs.io/en/stable/) modeling libraries. This package is not meant to replace or act as a higher level API for Reaktoro - it is only meant to enable setting up Reaktoro equilibrium problems as blocks on Pyomo models and automate transferring Reaktoro data into Pyomo variables. 

## 2. Prerequisites
**The user must familiarize thyself with Reaktoro options and usage, especially when selecting Reaktoro provided databases, database files, and activity models for aqueous, solid, and gas phases. This package does not automatically select any of these options, and will use ideal models by default.**

* Please refer here for information on [Reaktoro supported databases](https://reaktoro.org/tutorials/basics/loading-databases.html) (all are supported) 
* Please refer here for information on [Reaktoro supported activity models](https://reaktoro.org/tutorials/basics/specifying-activity-models.html) (all are supported, included chain operations, or passing in pre-configured activity models)

**By default the package will use:**

* **Database:** [PhreeqcDatabase](https://reaktoro.org/api/classReaktoro_1_1PhreeqcDatabase.html) 
* **Data file:** [pitzer.dat](https://reaktoro.org/api/classReaktoro_1_1PhreeqcDatabase.html) 
* **Aqueous activity models:** [ActivityModelIdealAqueous](https://reaktoro.org/api/namespaceReaktoro.html#ae431d4c8a1f283910ae1cf35024091b8)
* **Gas activity models:** [ActivityModelIdealGas](https://reaktoro.org/api/namespaceReaktoro.html#a7a0788a5a863d987a88b81303d80b427)
* **Solid activity models:** [ActivityModelIdealSolution](https://reaktoro.org/api/namespaceReaktoro.html#a6581d5c0cde36cae6d9c46dbc32d56f8)
* **Ion exchange activity modells:** [ActivityModelIonExchange](https://reaktoro.org/api/namespaceReaktoro.html#a6581d5c0cde36cae6d9c46dbc32d56f8)

## 3. Inputs and outputs of the Reaktoro blocks
The Reaktoro blocks built by this package are designed to solve an equilibrium problem using user provided apparent species or true species, temperature, pressure, and pH, which are broken down to base elements and equilibrated within Reaktoro to provide exact speciation and equilibrium state. Using this state the block can return various information supported by Reaktoro:

* [Chemical properties](https://reaktoro.org/api/classReaktoro_1_1ChemicalProps.html)
* [Aqueous properties](https://reaktoro.org/api/classReaktoro_1_1AqueousProps.html)
* Pyomo build properties, which are custom properties built in Pyomo that use chemical properties or aqueous properties as inputs 

Note: Only Reaktoro properties that return a single floating point or real value are supported, as they will be directly matched to a [Pyomo Var](https://pyomo.readthedocs.io/en/stable/pyomo_modeling_components/Variables.html). Reaktoro functions that return arrays or strings are not supported.

## 4. Tutorials, Examples, and Comparisons
Currently, repo includes several tutorials and examples.

*Tutorials:*

1. [Demonstration of working with Reaktoro block that shows](https://github.com/avdudchenko/reaktoro-pse/blob/main/src/reaktoro_pse/tutorials/ReaktoroBlock_tutorial.ipynb)
    * How to balance feed charge with ReaktoroBlock
    * How to add indexed ReaktoroBlocks [WaterTAP RO1D model](https://watertap.readthedocs.io/en/stable/technical_reference/unit_models/reverse_osmosis_1D.html) for calculation of Osmotic pressure
 
*Examples:*

1. [Example of adding ReaktoroBlock to basic desalination problem](https://github.com/avdudchenko/reaktoro-pse/blob/main/src/reaktoro_pse/examples/simple_desalination.py) that demonstrates how to:
    * Setup up basic ReaktoroBlock
    * Calculate basic properties for Scaling Tendency, pH, and Osmotic pressure 
    * Optimize system pH for operation at target Scaling Tendency
2. [Example of thermal precipitation](https://github.com/avdudchenko/reaktoro-pse/blob/main/src/reaktoro_pse/examples/thermal_precipitation.py) that demonstrates how to:
    * Configure different database from default
    * Get enthalpy and vapor pressure from Reaktoro
    * Setup precipitation calculation 
    * Setup simulation for removal of Calcite over different temperatures and estimate required energy input
3. [Example of ion exchange calculations](https://github.com/avdudchenko/reaktoro-pse/blob/main/src/reaktoro_pse/examples/simple_ion_exchange.py) that demonstrates how to:
    * Set up ReaktoroBlock for charge neutralizing the feed composition 
    * Use outputs from speciation block as inputs into a second property block
    * Add Ion Exchange phase and species into ReaktoroBlock
    * Optimize addition of acid and bases for maximizing Calcium removal selectivity over Magnesium

*Comparisons of Reaktoro-pse to PhreeqC when simulating*

These comparisons further demonstrate how to setup Reaktoro-pse for each type of calculation.

1. [Water removal from solution (e.g evaporative processes)](https://github.com/avdudchenko/reaktoro-pse/blob/main/src/reaktoro_pse/examples/reaktoro_pse_to_phreeqc_comparison/water_removal_comparison.py)
2. [Vapor pressure calculation at different temperatures](https://github.com/avdudchenko/reaktoro-pse/blob/main/src/reaktoro_pse/examples/reaktoro_pse_to_phreeqc_comparison/vapor_pressure_comparison.py)
3. [Precipitation of mineral phases](https://github.com/avdudchenko/reaktoro-pse/blob/main/src/reaktoro_pse/examples/reaktoro_pse_to_phreeqc_comparison/precipitation_comparison.py)
4.  [Mixing of two solution](https://github.com/avdudchenko/reaktoro-pse/blob/main/src/reaktoro_pse/examples/reaktoro_pse_to_phreeqc_comparison/solution_mixing_comparison.py)

## 5. Solvers 
To-date, this has only been tested with [cyipopt](https://cyipopt.readthedocs.io/en/stable/) - other solvers have not been tested. 
* For accessing other [HSL linear solvers](http://www.hsl.rl.ac.uk/ipopt/) beyond MUMPS (such as MA27, MA57, etc. solver common to WaterTAP and IDAES-PSE) follow [these instructions](https://cyipopt.readthedocs.io/en/latest/install.html) for adding it to cyipopt.

## 6. Requesting new features or issues
Please include a minimal example using Reaktoro for your specific feature or issue request if possible. Please also include full traceback for errors the occur. 

## 7. Compatibility
Reaktoro-pse depends on the following packages and/or versions:

- Python 3.9 through 3.12
- Reaktoro 2.12
- CyIpopt 1.4.1
- Pyomo>=6.8.0
- idaes-pse>=2.5.0
- watertap>=1.0.0 - (required for watertap-cyipopt wrapper only)

## 8. Getting started (for contributors)

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
