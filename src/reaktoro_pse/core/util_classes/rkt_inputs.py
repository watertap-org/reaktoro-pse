#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/reaktoro-pse/"
#################################################################################
from pyomo.environ import Var, value, units as pyunits
from pyomo.core.base.var import VarData


__author__ = "Alexander V. Dudchenko"


class RktInputTypes:
    """defines default types of inputs for reaktoro usage"""

    mol = "mol"
    K = "K"
    Pa = "Pa"
    pH = "pH"
    temperature = "temperature"
    enthalpy = "enthalpy"
    pressure = "pressure"
    dimensionless = "dimensionless"
    mass_units = ["kg", "mg"]
    system_state = "system_state"
    aqueous_phase = "aqueous_phase"
    gas_phase = "gas_phase"
    mineral_phase = "mineral_phase"
    ion_exchange_phase = "ion_exchange_phase"
    liquid_phase = "liquid_phase"
    solid_phase = "solid_phase"
    condensed_phase = "condensed_phase"
    supported_phases = [
        aqueous_phase,
        gas_phase,
        mineral_phase,
        solid_phase,
        ion_exchange_phase,
        condensed_phase,
        liquid_phase,
    ]
    non_species_types = [pH, enthalpy, pressure, temperature]


class RktInput:
    def __init__(self, var_name, pyomo_var=None):
        # TODO: Add more flexible check that user providd a pyomo variable or param
        self.var_name = var_name
        self.temp_value = 0  # used during reaktoro solver or temporarty holding a value

        self.rkt_index = None  # tracking rkt input index
        self.jacobian_index = None  # tracking rkt jacobian row index
        self.time_unit = None
        self.main_unit = None
        self.conversion_unit = None
        self.conversion_value = None
        self.required_unit = None
        self.rkt_name = var_name
        self.lower_bound = None
        self.input_type = None
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

    def get_input_type(self):
        return self.input_type

    def set_input_type(self, input_type):
        self.input_type = input_type

    def update_values(self, update_temp=False):
        if self.pyomo_var is not None:
            self.value = self.pyomo_var.value
        else:
            self.value = None
        if update_temp:
            self.set_temp_value(self.value)

    def set_temp_value(self, value):
        self.temp_value = value

    def get_temp_value(self):
        return self.temp_value

    def get_lower_bound(self):
        return self.lower_bound

    def set_lower_bound(self, value):
        self.lower_bound = value

    def get_value(self, update_temp=False, apply_conversion=False):
        self.update_values(update_temp)
        if apply_conversion:
            return value(self.get_pyomo_with_required_units())
        else:
            return self.value

    def get_pyomo_with_required_units(self):
        if self.conversion_value == None:
            return self.pyomo_var
        else:
            # print(
            #     pyunits.convert(
            #         self.pyomo_var / (self.conversion_value * self.conversion_unit),
            #         to_units=self.required_unit,
            #     )
            # )
            return pyunits.convert(
                self.pyomo_var / (self.conversion_value * self.conversion_unit),
                to_units=self.required_unit,
            )

    def set_unit_conversion(self, value, unit):
        self.conversion_unit = unit
        self.conversion_value = value

    def get_unit_conversion_value(self):
        return self.conversion_value

    def get_unit_conversion_units(self):
        return self.conversion_unit

    def get_required_unit(self):
        return self.required_unit

    def set_required_unit(self, main_unit):
        self.required_unit = pyunits.__getattr__(main_unit)
        if self.time_unit is not None:
            self.required_unit = self.required_unit / pyunits.__getattr__(
                self.time_unit
            )

    def set_pyomo_var(self, pyomo_var):
        self.pyomo_var = pyomo_var

    def get_pyomo_var(self):
        return self.pyomo_var

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
        if (
            default_unit == RktInputTypes.dimensionless
            and self.var_name not in RktInputTypes.non_species_types
        ):
            self.main_unit = RktInputTypes.mol
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

        self.rkt_input_list = []
        self.registered_phases = []
        self.species_list = {}
        self.all_species = []
        self.convert_to_rkt_species = {}
        self.composition_is_elements = {}
        self.conversion_method = {}
        for phase in RktInputTypes.supported_phases:
            self.species_list[phase] = []
            self.convert_to_rkt_species[phase] = False
            self.composition_is_elements[phase] = False
            self.conversion_method[phase] = "default"

    def enable_rkt_species_conversion(
        self, phase, convert=False, conversion_method="default"
    ):
        self.convert_to_rkt_species[phase] = convert
        self.conversion_method[phase] = conversion_method

    def set_composition_is_elements(self, phase, comp_is_elements=False):
        self.composition_is_elements[phase] = comp_is_elements

    def __setitem__(self, var_name, var):
        if isinstance(var, RktInput):
            # reference RktInput being passed in
            rkt_input = var
        else:
            # create new one if its a pyomo var or vardata
            rkt_input = RktInput(var_name, var)
        return super().__setitem__(var_name, rkt_input)

    def process_inputs(self):
        keys = list(self.keys())
        for key in keys:
            phase = self[key].get_input_type()
            self._set_species(key, self[key], phase)

    def auto_convert_to_rkt_species(self, var_name):
        if self.convert_to_rkt_species:
            var_name = self.convert_rkt_species_fun(var_name)
        return var_name

    def convert_rkt_species_fun(self, var_name, phase):
        if self.conversion_method[phase] == "default":
            var_name = specie_to_rkt_species(var_name)
        elif isinstance(self.conversion_method[phase], dict):
            var_name = self.conversion_method[phase][var_name]
        else:
            raise TypeError(
                f"Conversion method of {type(self.conversion_method[phase] )} is not supported)"
            )
        return var_name

    def _set_species(self, var_name, var, phase):
        if var_name not in RktInputTypes.non_species_types:
            if var_name not in self.species_list[phase]:
                if self.convert_to_rkt_species[phase]:
                    var_name = self.convert_rkt_species_fun(var_name, phase)
                    super().__setitem__(var_name, var)
                self.species_list[phase].append(var_name)
            if var_name not in self.all_species:
                self.all_species.append(var_name)

            if phase not in self.registered_phases:
                self.registered_phases.append(phase)

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
