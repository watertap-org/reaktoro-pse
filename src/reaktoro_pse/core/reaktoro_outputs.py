###############################################################################
# #################################################################################
# # WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# # through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# # National Renewable Energy Laboratory, and National Energy Technology
# # Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# # of Energy). All rights reserved.
# #
# # Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# # information, respectively. These files are also available online at the URL
# # "https://github.com/watertap-org/reaktoro-pse/"
# #################################################################################
###############################################################################
import reaktoro as rkt

import json
from reaktoro_pse.core.reaktoro_state import ReaktoroState
import reaktoro_pse.core.pyomo_property_writer.property_functions as propFuncs

# disabling warnings

__author__ = "Alexander Dudchenko"


""" class to setup output constraints, outputs, and jacobian reaktoro solver class"""


# NOTE: Consider refactoring to inherit RktInput class
class RktOutput:
    """general class to stor output metadata and vars"""

    def __init__(
        self,
        property_type,
        property_name,
        property_index=None,
        get_function=None,
        pyomo_var=None,
        value=None,
        stoichiometric_coeff=None,
        jacobian_type=None,
    ):
        self.property_type = property_type
        self.property_name = property_name  # properties from which to extract data
        self.property_index = property_index  # index if any
        if property_type != PropTypes.pyomo_built_prop:
            # function for getting reaktoro value
            self.set_get_function(get_function)
        else:
            self.set_poyomo_build_option(
                get_function
            )  # class that contains information for building pyomo constraints if any
        # pyomo var to reference if any - will be built if not user provided
        self.pyomo_var = pyomo_var
        self.value = value  # manually specified value
        self.jacobian_value = (
            None  # (might not be used) will be an array of jacobian values
        )
        self.stoichiometric_coeff = (
            stoichiometric_coeff  # for tracking stichometry if needed
        )
        self.jacobian_type = jacobian_type

    def get_value(self, prop_object, update_values=False):
        value = self.get_function(prop_object, self.property_name, self.property_index)
        # print(self.property_name, self.property_index, value)
        if update_values:
            self.value = value
        return value

    def set_poyomo_build_option(self, func):
        self.pyomo_build_options = func

    def set_get_function(self, func):
        self.get_function = func

    def set_property_type(self, prop):
        self.property_type = prop

    def set_pyomo_var(self, var):
        self.pyomo_var = var

    def get_pyomo_var(self):
        return self.pyomo_var

    def set_pyomo_var_value(self, value):
        self.pyomo_var.value = value

    def get_pyomo_var_value(self, value):
        return self.pyomo_var.value

    def set_jacobian_value(self, value):
        self.jacobian_value = value


class PyomoBuildOptions:
    def __init__(self):
        self.properties = {}  # creats dict of rkt Properties to get desired values
        self.options = (
            {}
        )  # creats dict for optijnal parmeters used during pyomo constraint - property specfic
        self.build_constraint_function = None

    def register_property(self, property_type, property_name, property_index=None):

        self.properties[(property_name, property_index)] = RktOutput(
            property_type=property_type,
            property_name=property_name,
            property_index=property_index,
        )

    def register_option(self, option, value):
        self.options[option] = value

    def register_build_function(self, function):
        self.build_constraint_function = function


class PyomoProperties:
    def __init__(self, reaktor_state, chem_props, aqueous_props):
        self.state = reaktor_state
        self.chem_props = chem_props
        self.aqueous_props = aqueous_props

    def scalingTendency(self, property_index):
        """build scaling tendencty - RKT has saturationIndex but no scalingIndex"""
        required_props = PyomoBuildOptions()
        required_props.register_property(
            PropTypes.aqueous_prop, "saturationIndex", property_index
        )
        required_props.register_build_function(
            propFuncs.build_scaling_tendency_constraint
        )
        return required_props

    def scalingTendencyDirect(self, property_index):
        """build pyomo constraint for scaling index calculations directly form chem props
        #TODO: Need to add check for database being used as only PhreeqC is really supported at the
        moment"""
        required_props = PyomoBuildOptions()
        ref_temp = 25  # degC
        ref_pressure = 1  # atm
        spec = self.aqueous_props.saturationSpecies().get(property_index)
        thermo_model = spec.standardThermoModel()
        pr = spec.props(ref_temp, "C", ref_pressure, "atm")
        specie_volume = float(pr.V0)  # returns auto diff/not usable with pyomo

        # get data from thermo prop
        jsp = thermo_model.params().dumpJson()
        jsp_dict = json.loads(jsp)
        if jsp_dict[0].get("PhreeqcLgK") is not None:
            required_props.register_option("logk_type", "Analytical")
            required_props.register_option("logk_paramters", jsp_dict[0]["PhreeqcLgK"])
        elif jsp_dict[0].get("VantHoff") is not None:
            required_props.register_option("logk_type", "VantHoff")
            required_props.register_option("logk_paramters", jsp_dict[0]["VantHoff"])
        else:
            raise NotImplemented(f"reaction type {jsp_dict} not supported")
        required_props.register_option("gas_constant", rkt.universalGasConstant)
        volume_reactants = 0
        for s, mol in spec.reaction().reactants():
            spec = self.state.system.species().get(s.name())
            thermo_model = spec.standardThermoModel()
            _pr = spec.props(ref_temp, "C", ref_pressure, "atm")
            volume_reactants += float(_pr.V0) * abs(mol)
            required_props.register_property(
                PropTypes.chem_prop, "speciesActivityLn", s.name()
            )
            required_props.properties[
                ("speciesActivityLn", s.name())
            ].stoichiometric_coeff = abs(
                mol
            )  # create on demand to track coefficients

        required_props.register_option(
            "delta_V", float(specie_volume - (volume_reactants))
        )
        required_props.register_property(PropTypes.chem_prop, "temperature")
        required_props.register_property(PropTypes.chem_prop, "pressure")
        required_props.register_build_function(
            propFuncs.build_direct_scaling_tendency_constraint
        )
        return required_props

    def osmoticPressure(self, property_index):
        """build osmoric pressure constraint, as its not available from reaktoro"""
        required_props = PyomoBuildOptions()
        required_props.register_property(
            PropTypes.chem_prop, "speciesStandardVolume", property_index
        )
        required_props.register_property(
            PropTypes.chem_prop, "speciesActivityLn", property_index
        )
        required_props.register_property(PropTypes.chem_prop, "temperature")
        required_props.register_build_function(propFuncs.build_osmotic_constraint)
        required_props.register_option("gas_constant", rkt.universalGasConstant)
        return required_props

    def pHDirect(self, property_index=None):
        """build direct pH caclautions from chem props"""
        required_props = PyomoBuildOptions()
        required_props.register_property(PropTypes.chem_prop, "speciesActivityLn", "H+")
        required_props.register_build_function(propFuncs.build_ph_constraint)
        return required_props

    def vaporPressure(self, property_index=None):
        """build direct pH caclautions from chem props"""
        required_props = PyomoBuildOptions()
        required_props.register_property(
            PropTypes.chem_prop, "speciesActivityLn", property_index
        )
        required_props.register_build_function(
            propFuncs.build_vapor_pressure_constraint
        )
        return required_props


class PropTypes:
    """define base property types"""

    chem_prop = "chemProp"
    aqueous_prop = "aqueousProp"
    pyomo_built_prop = "pyomoBuiltProperties"


class ReaktoroOutputSpec:
    def __init__(self, reaktor_state):
        self.state = reaktor_state
        if isinstance(self.state, ReaktoroState) == False:
            raise TypeError("Reator outputs require rektoroState class")
        self.supported_properties = {}
        self.supported_properties[PropTypes.chem_prop] = self.state.state.props()
        self.supported_properties[PropTypes.aqueous_prop] = rkt.AqueousProps(
            self.state.state.props()
        )
        self.supported_properties[PropTypes.pyomo_built_prop] = PyomoProperties(
            self.state,
            self.supported_properties[PropTypes.chem_prop],
            self.supported_properties[PropTypes.aqueous_prop],
        )
        self.supported_properties[PropTypes.pyomo_built_prop].osmoticPressure("Calcite")
        self.rkt_outputs = {}  # outputs that reaktoro needs to generate
        self.user_outputs = {}  # outputs user requests
        self.get_possible_indexes()

    def update_supported_props(self):
        self.state.state.props().update(self.state.state)
        self.supported_properties[PropTypes.aqueous_prop].update(
            self.state.state.props()
        )

    def evaluate_property(
        self, RktOutputObject, property_type=None, update_values_in_object=False
    ):
        """evaluating reaktoro output object, doing it here so we can
        provide custom property types -> this will be require for numerical derivatives

        Keywords:
        RktOutputObject -- output object that contains property info
        property_type -- either a propType, or supplied user property
        """

        if isinstance(RktOutputObject, RktOutput) == False:
            raise TypeError(
                "Provided object is not supported, pplease provide an rktOuput object"
            )
        if property_type is None:
            property_type = self.supported_properties[RktOutputObject.property_type]
        return RktOutputObject.get_value(property_type, update_values_in_object)

    def register_output(
        self,
        property_name,
        property_index=None,
        get_all_indexes=False,
        pyomo_var=None,
    ):
        """register a reaktoro output, couple it to property type.

        Keywords:
        property_name -- prop name (specieisActivityLn, pH etc)
        property_index -- prop index if any (H+, etc) (default: None)
        get_all_indexes -- if user want to get all possible indexs for specfied prop (default: False)
        pyomo_var -- pyomo var that should be used for the output of this property (optional: will be auto built) (default: None)
        """
        if get_all_indexes:
            self.get_all_indexes(property_name)
        else:
            property_type, get_function = self.get_prop_type(
                property_name, property_index
            )
            self.process_output(
                property_type=property_type,
                property_name=property_name,
                property_index=property_index,
                get_function=get_function,
                pyomo_var=pyomo_var,
            )

    def process_output(
        self,
        property_type,
        property_name,
        property_index=None,
        get_function=None,
        pyomo_var=None,
    ):
        # if property_index is None:
        #     index = property_name
        # else:
        index = (property_name, property_index)
        if property_type != PropTypes.pyomo_built_prop:
            self.user_outputs[index] = RktOutput(
                property_type=property_type,
                property_name=property_name,
                property_index=property_index,
                get_function=get_function,
                pyomo_var=pyomo_var,
            )
            if index not in self.rkt_outputs:
                self.rkt_outputs[index] = self.user_outputs[index]
        else:
            self.user_outputs[index] = RktOutput(
                property_type=property_type,
                property_name=property_name,
                property_index=property_index,
                get_function=get_function,
                pyomo_var=pyomo_var,
            )
            for index, prop in get_function.properties.items():
                """chcek if prop already exists if it does ont add it outputs
                otherwise overwrite it"""
                if index not in self.rkt_outputs:
                    self.rkt_outputs[index] = prop
                else:
                    get_function.properties[index] = self.rkt_outputs[index]

    def get_all_indexes(
        self,
        property_name,
    ):
        if "species" in property_name:
            for specie in self.species:
                property_type, get_function = self.get_prop_type(property_name, specie)
                self.process_output(
                    property_type=property_type,
                    property_name=property_name,
                    property_index=specie,
                    get_function=get_function,
                )
        elif "elements" in property_name:
            for element in self.elements:
                property_type, get_function = self.get_prop_type(property_name, element)
                self.process_output(
                    property_type=property_type,
                    property_name=property_name,
                    property_index=element,
                    get_function=get_function,
                )
        else:
            raise NotImplementedError(
                f"{property_name} is not supported for automatic indexing"
            )

    def get_prop_type(self, property_name, property_index=None):
        """this function will try differernt property types useing standard
        call functions to figure out property type (aquous, chem, etc) and how
        to get the actual value prop.value(), prop.value(index), prop.value(index).val()
        and so forth"""
        for supported_props, prop in self.supported_properties.items():
            for func_attempt in [
                self._get_prop_phase_name_val,
                self._get_prop_name,
                self._get_prop_name_val,
            ]:
                try:
                    if supported_props != PropTypes.pyomo_built_prop:
                        func_attempt(
                            prop,
                            property_name,
                            property_index,
                        )
                        return supported_props, func_attempt
                    else:
                        func_results = getattr(prop, property_name)(
                            property_index=property_index
                        )
                        for prop_key, obj in func_results.properties.items():
                            supported_prop, func_result = self.get_prop_type(
                                obj.property_name, obj.property_index
                            )
                            obj.set_get_function(func_result)
                            obj.set_property_type(supported_prop)

                        return supported_props, func_results
                except (TypeError, KeyError, AttributeError, RuntimeError):
                    pass

        raise NotImplementedError(
            f"The {property_name}, {property_index} is not supported at the moment"
        )

    def get_possible_indexes(self):
        """this gets possible indexes for when user wants to output all indexes for properties"""
        self.elements = [
            specie.symbol() for specie in self.state.state.system().elements()
        ]
        self.species = [specie.name() for specie in self.state.state.system().species()]
        self.saturation_species = [
            specie.name()
            for specie in self.supported_properties[
                PropTypes.aqueous_prop
            ].saturationSpecies()
        ]

    """ start of possible call function to extract values from reactoro properties"""

    def _get_prop_phase_name_val(self, prop_type, prop_name, prop_index):
        """get prop based on phase, used for chem_props.phaseProp"""
        value = getattr(prop_type.phaseProps(prop_index), prop_name)()

        return float(value)

    def _get_prop_name_val(self, prop_type, prop_name, prop_index=None):
        """get prop based on name/index and execute value call"""
        if prop_index is None:
            value = getattr(prop_type, prop_name)
        else:
            value = getattr(prop_type, prop_name)(prop_index)
        return float(value.val())

    def _get_prop_name(self, prop_type, prop_name, prop_index=None):
        """get prop based/index on name only"""
        if prop_index is None:
            value = getattr(prop_type, prop_name)()
        else:
            value = getattr(prop_type, prop_name)(prop_index)
        return float(value)
