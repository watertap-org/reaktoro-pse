from pyomo.environ import Var, units as pyunits
from pyomo.core.base.var import VarData

__author__ = "Alexander V. Dudchenko (SLAC)"


class RktInput:
    def __init__(self, var_name, pyomo_var=None):
        # TODO: Add more flexible check that user providd a pyomo variable or param
        self.var_name = var_name
        self.temp_value = 0  # used during reaktoro solver or temporarty holding a value

        self.rkt_index = None  # tracking rkt input index
        self.jacobian_index = None  # tracking rkt jacobian row index
        self.time_unit = None
        self.main_unit = None
        self.rkt_name = var_name
        if pyomo_var is not None:
            if isinstance(pyomo_var, (Var, VarData)) == False:
                raise TypeError(
                    "{var_name} is not a pyomo variable, ensure its pyomo variable"
                )
            self.value = pyomo_var.value
            self.pyomo_var = pyomo_var
            self.check_unit()
        else:
            self.pyomo_var = None

    def update_values(self, update_temp=False):
        if self.pyomo_var is not None:
            self.value = self.pyomo_var.value
        else:
            self.value = None
        if update_temp:
            self.temp_value = self.value

    def get_value(self, update_temp=False):
        self.update_values(update_temp)
        return self.value

    def set_pyomo_var(self, pyomo_var):
        self.pyomo_var = pyomo_var

    def set_jacobian_index(self, idx):
        self.jacobian_index = idx

    def get_jacobian_index(self):
        return self.jacobian_index

    def set_rkt_index(self, idx):
        self.rkt_index = idx

    def set_rkt_input_name(self, name):
        self.rkt_name = name

    def get_rkt_input_name(self):
        return self.rkt_name

    def get_rkt_index(self):
        return self.rkt_index

    def set_pyomo_var_value(self, value):
        self.pyomo_var.value = value

    def get_pyomo_var_value(self):
        return self.pyomo_var.value

    def check_unit(self):
        """this checks if unit has a time component to it (e.g mol/s)
        reaktoro is "batch" and has no flow, so here we will isolate
        the primary mass unit from time unit and
        also convert them to string"""
        default_unit = str(pyunits.get_units(self.pyomo_var))
        if default_unit == "dimensionless" and self.var_name not in [
            "pH",
            "temperature",
            "pressure",
        ]:
            self.main_unit = "mol"
            self.time_unit = None
        split_units = default_unit.split("/")
        if len(split_units) == 2:
            self.time_unit = split_units[1]
            self.main_unit = split_units[0]
        else:
            self.time_unit = None
            self.main_unit = default_unit


class RktInputs(dict):
    def __init__(self):
        self.species_list = []
        self.rkt_input_list = []
        self.convert_to_rkt_species = False
        self.composition_is_elements = False
        self.conversion_method = "default"

    def enable_rkt_species_conversion(self, convert=False, conversion_method="default"):
        self.convert_to_rkt_species = convert
        self.conversion_method = conversion_method

    def set_composition_is_elements(self, comp_is_elements=False):
        self.composition_is_elements = comp_is_elements

    def __setitem__(self, var_name, var):
        if isinstance(var, RktInput):
            # reference RktInput being passed in
            rkt_input = var
        else:
            # create new one if its a pyomo var or vardata
            rkt_input = RktInput(var_name, var)

        self._set_species(var_name, rkt_input)
        return super().__setitem__(var_name, rkt_input)

    def auto_convert_to_rkt_species(self, var_name):
        if self.convert_to_rkt_species:
            var_name = self.convert_rkt_species_fun(var_name)
        return var_name

    def convert_rkt_species_fun(self, var_name):
        if self.conversion_method == "default":
            var_name = specie_to_rkt_species(var_name)
        elif isinstance(self.conversion_method, dict):
            var_name = self.conversion_method[var_name]
        else:
            raise TypeError(
                f"Conversion method of {type(self.conversion_method)} is not supported)"
            )
        return var_name

    def _set_species(self, var_name, var):
        if var_name not in ["pH", "temperature", "pressure"]:
            if var_name not in self.species_list:
                if self.convert_to_rkt_species:
                    var_name = self.convert_rkt_species_fun(var_name)
                    super().__setitem__(var_name, var)
                self.species_list.append(var_name)

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
    for charge, species_list in name_dict.items():
        for spc in species_list:
            if spc in species:
                if charge == "H4SiO4":
                    return charge
                else:
                    return f"{spc}{charge}"

    raise TypeError(f"Species {species} not found, please add to conversion dict")