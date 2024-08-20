# reaktoro-pse

## Compatibility

reaktoro-pse depends on the following packages and/or versions:

- Python 3.9 through 3.12
- Reaktoro 2.12
- CyIpopt 1.4.1
- Pyomo

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
