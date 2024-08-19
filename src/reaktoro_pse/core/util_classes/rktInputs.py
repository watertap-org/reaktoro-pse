from pyomo.environ import Var, units as pyunits
from pyomo.core.base.var import VarData

__author__ = "Alexander V. Dudchenko (SLAC)"


class rktInput:
    def __init__(self, var_name, pyomo_var=None):
        # TODO: Add more flexible check that user providd a pyomo variable or param
        self.varName = var_name
        self.tempValue = 0  # used during reaktoro solver or temporarty holding a value

        self.rktIndex = None  # tracking rkt input index
        self.jacobianIndex = None  # tracking rkt jacobian row index
        self.timeUnit = None
        self.mainUnit = None
        self.rktName = var_name
        if pyomo_var is not None:
            if isinstance(pyomo_var, (Var, VarData)) == False:
                raise TypeError(
                    "{varName} is not a pyomo variable, ensure its pyomo variable"
                )
            self.value = pyomo_var.value
            self.pyomoVar = pyomo_var
            self.check_unit()
        else:
            self.pyomoVar = None

    def update_values(self, update_temp=False):
        if self.pyomoVar is not None:
            self.value = self.pyomoVar.value
        else:
            self.value = None
        if update_temp:
            self.tempValue = self.value

    def get_value(self, update_temp=False):
        self.update_values(update_temp)
        return self.value

    def set_pyomoVar(self, pyomo_var):
        self.pyomoVar = pyomo_var

    def set_jacobian_index(self, idx):
        self.jacobianIndex = idx

    def get_jacobian_index(self):
        return self.jacobianIndex

    def set_rkt_index(self, idx):
        self.rktIndex = idx

    def set_rkt_input_name(self, name):
        self.rktName = name

    def get_rkt_input_name(self):
        return self.rktName

    def get_rkt_index(self):
        return self.rktIndex

    def set_pyomo_var_value(self, value):
        self.pyomoVar.value = value

    def get_pyomo_var_value(self):
        return self.pyomoVar.value

    def check_unit(self):
        """this checks if unit has a time component to it (e.g mol/s)
        reaktoro is "batch" and has no flow, so here we will isolate
        the primary mass unit from time unit and
        also convert them to string"""
        default_unit = str(pyunits.get_units(self.pyomoVar))
        if default_unit == "dimensionless" and self.varName not in [
            "pH",
            "temperature",
            "pressure",
        ]:
            self.mainUnit = "mol"
            self.timeUnit = None
        split_units = default_unit.split("/")
        if len(split_units) == 2:
            self.timeUnit = split_units[1]
            self.mainUnit = split_units[0]
        else:
            self.timeUnit = None
            self.mainUnit = default_unit


class rktInputs(dict):
    def __init__(self):
        self.speciesList = []
        self.rktInputList = []
        self.convertToRktSpecies = False
        self.compositionIsElements = False
        self.conversion_method = "default"

    def enable_rkt_species_conversion(self, convert=False, conversion_method="default"):
        self.convertToRktSpecies = convert
        self.conversion_method = conversion_method

    def set_composition_is_elements(self, comp_is_elements=False):
        self.compositionIsElements = comp_is_elements

    def __setitem__(self, var_name, var):
        if isinstance(var, rktInput):
            # reference rktInput being passed in
            rkt_input = var
        else:
            # create new one if its a pyomo var or vardata
            rkt_input = rktInput(var_name, var)

        self._set_species(var_name, rkt_input)
        return super().__setitem__(var_name, rkt_input)

    def _set_species(self, var_name, var):
        if var_name not in ["pH", "temperature", "pressure"]:
            if var_name not in self.speciesList:
                if self.convertToRktSpecies:
                    if self.conversion_method == "default":
                        var_name = specie_to_rkt_species(var_name)
                    elif isinstance(self.conversion_method, dict):
                        var_name = self.conversion_method[var_name]
                    else:
                        raise TypeError(
                            f"Conversion method of {type(self.conversion_method)} is not supported)"
                        )
                    super().__setitem__(var_name, var)
                self.speciesList.append(var_name)

    def __getitem__(self, var_name):
        var = super().__getitem__(var_name)
        var.update_values()
        return var


def specie_to_rkt_species(species):
    """basic function to convert speices to rkt names"""

    # TODO: needs to be better automated
    name_dict = {
        "-2": ["SO4", "CO3"],
        "-": [
            "Cl",
            "HCO3",
        ],
        "+": ["Na", "K"],
        "+2": ["Mg", "Ca", "Sr", "Ba"],
        "": ["H2O", "CO2"],
        "H4SiO4": ["Si", "SiO2"],
    }
    for charge, speciesList in name_dict.items():
        for spc in speciesList:
            if spc in species:
                return f"{spc}{charge}"

    raise TypeError(f"Species {species} not found, please add to conversion dict")
