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
from matplotlib.pylab import f
import reaktoro as rkt

import json
from reaktoro_pse.core.reaktoro_state import ReaktoroState
import reaktoro_pse.core.pyomo_property_writer.property_functions as propFuncs
from reaktoro_pse.core.util_classes.rkt_inputs import RktInputTypes
import copy
import idaes.logger as idaeslog
import math

_log = idaeslog.getLogger(__name__)
# disabling warnings

__author__ = "Alexander V. Dudchenko"


# class to setup output constraints, outputs, and jacobian reaktoro solver class


# NOTE: Consider refactoring to inherit RktInput class
class RktOutput:
    """general class to store output metadata and vars"""

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
        self.jacobian_index = (
            property_name,
            property_index,
        )  # index for jacobian if any
        self.get_function = None
        self.set_option_function(property_type, get_function)
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
        self.auto_scaled = False
        self.min_rkt_value = None
        self.max_rkt_value = None
        self.io_type = None
        self.conversion_function = None
        self.derivative_conversion_function = None
        self.derivative_conversion_value = 1
        self.converted_prop_type = None

    def get_value(self, prop_object, update_values=False, apply_der_conversion=True):
        value = self.get_function(prop_object, self.property_name, self.property_index)
        if self.min_rkt_value is not None and value < self.min_rkt_value:
            value = self.min_rkt_value
        if self.max_rkt_value is not None and value > self.max_rkt_value:
            value = self.max_rkt_value
        value, derivative_conversion_value = self.apply_conversions(value)
        if apply_der_conversion:
            self.derivative_conversion_value = derivative_conversion_value
        if update_values:
            self.value = value

        return value

    def get_derivative_conversion_factor(self):
        return self.derivative_conversion_value

    def apply_conversions(self, value):
        """apply conversion functions if any"""
        if self.conversion_function is not None:
            converted_value = self.conversion_function(value)
            derivative_conversion = self.derivative_conversion_function(value)
            return converted_value, derivative_conversion
        else:
            return value, 1

    def delete_pyomo_var(self):
        # self.update_values()
        del self.pyomo_var
        self.pyomo_var = None

    def remove_unpicklable_data(self):
        self.delete_pyomo_var()
        if self.property_type == PropTypes.pyomo_built_prop:
            del self.pyomo_build_options
            # self.get_function = None

    def set_option_function(self, property_type, get_function):
        if property_type != PropTypes.pyomo_built_prop:
            # function for getting reaktoro value
            self.set_get_function(get_function)
        else:
            self.set_poyomo_build_option(
                get_function
            )  # class that contains information for building pyomo constraints if any

    def set_poyomo_build_option(self, func):
        self.pyomo_build_options = func

    def set_pyo_value(self, value):
        self.pyomo_var.set_value(value)

    def set_get_function(self, func):
        self.get_function = func

    def set_property_type(self, prop):
        self.property_type = prop

    def set_pyomo_var(self, var):
        self.pyomo_var = var

    def get_pyomo_var(self):
        return self.pyomo_var

    def set_pyomo_var_value(self, value):
        self.value = value
        if self.pyomo_var is not None:
            self.pyomo_var.set_value(value)

    def get_pyomo_var_value(self):
        return self.pyomo_var.value

    def set_jacobian_value(self, value):
        self.jacobian_value = value

    def get_lb(self):
        return self.lb


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
        """build scaling tendency - RKT has saturationIndex but no scalingIndex"""
        required_props = PyomoBuildOptions()
        required_props.register_property(
            PropTypes.aqueous_prop, "saturationIndex", property_index
        )
        required_props.register_build_function(
            propFuncs.build_scaling_tendency_constraint
        )
        return required_props

    def chargeDirect(self, property_index):
        """build pyomo constraint for charge calculations directly form chem props"""
        required_props = PyomoBuildOptions()
        required_props.register_property(PropTypes.chem_prop, "charge", property_index)
        species = []
        for specie in self.state.state.system().species():
            name, charge = specie.name(), specie.charge()
            if charge != 0:
                species.append((name, charge))
                required_props.register_property(
                    PropTypes.chem_prop, "speciesAmount", name
                )
        required_props.register_option("species", species)
        required_props.register_build_function(propFuncs.build_direct_charge)
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
        # reference  https://help.syscad.net/PHREEQC_Reverse_Osmosis
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

    def alkalinityAsCaCO3(self, property_index=None):
        """build alkalinity and convert it to CaCO3 basis"""
        required_props = PyomoBuildOptions()
        required_props.register_property(PropTypes.aqueous_prop, "alkalinity")
        required_props.register_build_function(
            propFuncs.build_alkalinity_as_caco3_constraint
        )
        return required_props


class ConvertedPropTypes:
    def __init__(self, reaktor_state, chem_props, aqueous_props):
        self.state = reaktor_state
        self.chem_props = chem_props
        self.aqueous_props = aqueous_props

    def scalingTendency(self, property_index):
        """build scaling tendency - RKT has saturationIndex but no scalingIndex"""
        output = RktOutput(
            property_type=PropTypes.aqueous_prop,
            property_name="saturationIndex",
            property_index=property_index,
        )
        output.converted_prop_type = "scalingTendency"
        output.conversion_function = lambda x: 10 ** (x)
        output.derivative_conversion_function = lambda x: 10**x * math.log(10)
        return output

    def logSpeciesAmount(self, property_index):
        """build log species amount"""
        output = RktOutput(
            property_type=PropTypes.chem_prop,
            property_name="speciesAmount",
            property_index=property_index,
        )
        output.converted_prop_type = "logSpeciesAmount"
        output.conversion_function = lambda x: math.log10(x)
        output.derivative_conversion_function = lambda x: 1 / (x * math.log(10))
        return output


class PropTypes:
    """define base property types"""

    chem_prop = "chemProp"
    aqueous_prop = "aqueousProp"
    pyomo_built_prop = "pyomoBuiltProperties"
    converted_prop = (
        "convertedProp"  # used for converted properties, like pH, alkalinity, etc
    )


class ReaktoroOutputExport:
    def __init__(self):
        self.rkt_outputs = None
        self.user_outputs = None

    def copy_rkt_outputs(self, outputs):
        self.rkt_outputs = {}
        for key, obj in outputs.items():
            self.rkt_outputs[key] = RktOutput(
                obj.property_type,
                obj.property_name,
                obj.property_index,
                # get_function=obj.get_function,
                value=obj.value,
                jacobian_type=obj.jacobian_type,
            )
            self.rkt_outputs[key].min_rkt_value = obj.min_rkt_value
            self.rkt_outputs[key].max_rkt_value = obj.max_rkt_value
            self.rkt_outputs[key].converted_prop_type = obj.converted_prop_type
            self.rkt_outputs[key].remove_unpicklable_data()

    def copy_user_outputs(self, outputs):
        self.user_outputs = {}  # copy.deepcopy(outputs)
        for key, obj in outputs.items():
            self.user_outputs[key] = RktOutput(
                obj.property_type,
                obj.property_name,
                obj.property_index,
                # get_function=obj.get_function,
                value=obj.value,
                jacobian_type=obj.jacobian_type,
            )
            self.user_outputs[key].min_rkt_value = obj.min_rkt_value
            self.user_outputs[key].max_rkt_value = obj.max_rkt_value
            self.user_outputs[key].converted_prop_type = obj.converted_prop_type
            self.user_outputs[key].remove_unpicklable_data()


class ReaktoroOutputSpec:
    def __init__(self, reaktor_state):
        self.state = reaktor_state
        if isinstance(self.state, ReaktoroState) == False:
            raise TypeError("Reator outputs require rektoroState class")

        self.supported_properties = {}
        self.supported_properties[PropTypes.chem_prop] = self.state.state.props()

        if RktInputTypes.aqueous_phase in self.state.inputs.registered_phases:
            self.supported_properties[PropTypes.aqueous_prop] = rkt.AqueousProps(
                self.state.state.props()
            )
            aq_props = self.supported_properties[PropTypes.aqueous_prop]
        else:
            aq_props = None
        self.supported_properties[PropTypes.pyomo_built_prop] = PyomoProperties(
            self.state, self.supported_properties[PropTypes.chem_prop], aq_props
        )
        self.supported_properties[PropTypes.converted_prop] = ConvertedPropTypes(
            self.state, self.supported_properties[PropTypes.chem_prop], aq_props
        )

        self.rkt_outputs = {}  # outputs that reaktoro needs to generate
        self.user_outputs = {}  # outputs user requests
        self.output_limits = {}
        self.register_output_limits("speciesAmount", min_value=1e-16, max_value=None)
        self.get_possible_indexes()

    def update_supported_props(self):
        self.state.state.props().update(self.state.state)
        if RktInputTypes.aqueous_phase in self.state.inputs.registered_phases:
            self.supported_properties[PropTypes.aqueous_prop].update(
                self.state.state.props()
            )

    def evaluate_property(
        self,
        RktOutputObject,
        property_type=None,
        update_values_in_object=False,
        apply_der_conversion=True,
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
        return RktOutputObject.get_value(
            property_type, update_values_in_object, apply_der_conversion
        )

    def register_output(
        self,
        property_name,
        property_index=None,
        get_all_indexes=False,
        pyomo_var=None,
        ignore_indexes=None,
    ):
        """register a reaktoro output, couple it to property type.

        Keywords:
        property_name -- prop name (specieisActivityLn, pH etc)
        property_index -- prop index if any (H+, etc) (default: None)
        get_all_indexes -- if user want to get all possible indexs for specfied prop (default: False)
        pyomo_var -- pyomo var that should be used for the output of this property (optional: will be auto built) (default: None)
        """
        if get_all_indexes:
            self.get_all_indexes(property_name, ignore_indexes)
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

    def register_output_limits(self, property_name, min_value=None, max_value=None):
        """register output limits for the property, this will be used to limit the output values from reaktoro solve before passing them
        to pyomo variables, this is useful for scaling and limiting outputs to avoid numerical issues
        """

        self.output_limits[property_name] = {
            "min_rkt_value": min_value,
            "max_rkt_value": max_value,
        }

    def apply_output_limits(self, property_name, rkt_output):

        if property_name in self.output_limits:
            limits = self.output_limits[property_name]
            rkt_output.min_rkt_value = limits["min_rkt_value"]
            rkt_output.max_rkt_value = limits["max_rkt_value"]

    def process_output(
        self,
        property_type,
        property_name,
        property_index=None,
        get_function=None,
        pyomo_var=None,
    ):
        index = (property_name, property_index)
        if index not in self.user_outputs:
            prop_type = None
            if "specie" in property_name:
                prop_type = "specie"
            elif "element" in property_name:
                prop_type = "element"
            if property_type == PropTypes.pyomo_built_prop:
                self.user_outputs[index] = RktOutput(
                    property_type=property_type,
                    property_name=property_name,
                    property_index=property_index,
                    get_function=get_function,
                    pyomo_var=pyomo_var,
                )
                self.apply_output_limits(property_name, self.user_outputs[index])
                self.user_outputs[index].io_type = prop_type
                for index, prop in get_function.properties.items():
                    # check if prop already exists if it does nor add it outputs
                    # otherwise overwrite it
                    if index not in self.rkt_outputs:
                        self.rkt_outputs[index] = prop
                    else:
                        get_function.properties[index] = self.rkt_outputs[index]
                    self.apply_output_limits(
                        prop.property_name, self.rkt_outputs[index]
                    )

            elif property_type == PropTypes.converted_prop:
                # if converted prop, we need to get the converted prop type
                # and then set the get function
                self.user_outputs[index] = get_function
                self.user_outputs[index].set_pyomo_var(pyomo_var)
                self.apply_output_limits(
                    get_function.property_name, self.user_outputs[index]
                )
                if index not in self.rkt_outputs:
                    self.rkt_outputs[index] = self.user_outputs[index]
                    self.rkt_outputs[index].jacobian_index = (
                        get_function.property_name,
                        get_function.property_index,
                    )
            else:
                self.user_outputs[index] = RktOutput(
                    property_type=property_type,
                    property_name=property_name,
                    property_index=property_index,
                    get_function=get_function,
                    pyomo_var=pyomo_var,
                )

                self.apply_output_limits(property_name, self.user_outputs[index])
                self.user_outputs[index].io_type = prop_type
                if index not in self.rkt_outputs:
                    self.rkt_outputs[index] = self.user_outputs[index]
        else:
            _log.warning("Output {index}, already added!")

    def get_all_indexes(
        self,
        property_name,
        ignore_indexes,
    ):
        if "specie" in property_name.lower():
            for specie in self.species:
                if ignore_indexes is None or specie not in str(ignore_indexes):
                    property_type, get_function = self.get_prop_type(
                        property_name, specie
                    )
                    print(property_type, get_function)
                    self.process_output(
                        property_type=property_type,
                        property_name=property_name,
                        property_index=specie,
                        get_function=get_function,
                    )
        elif "element" in property_name.lower():
            for element in self.elements:
                if ignore_indexes is None or element not in str(ignore_indexes):
                    property_type, get_function = self.get_prop_type(
                        property_name, element
                    )
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
                    if supported_props == PropTypes.pyomo_built_prop:

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
                    elif supported_props == PropTypes.converted_prop:
                        func_results = getattr(prop, property_name)(
                            property_index=property_index
                        )
                        _, func_result = self.get_prop_type(
                            func_results.property_name,
                            func_results.property_index,
                        )
                        func_results.set_get_function(func_result)
                        return supported_props, func_results
                    else:
                        self._func_tester(
                            func_attempt,
                            prop,
                            property_name,
                            property_index,
                        )
                        return supported_props, func_attempt
                except (TypeError, KeyError, AttributeError, RuntimeError):
                    pass
        raise NotImplementedError(
            f"""The {property_name}, {property_index} was not found,
                its either not supported, or requested index is not in present.
            """
        )

    def get_possible_indexes(self):
        """this gets possible indexes for when user wants to output all indexes for properties"""
        self.elements = [
            specie.symbol() for specie in self.state.state.system().elements()
        ]
        self.species = [specie.name() for specie in self.state.state.system().species()]

        # self.saturation_species = [
        #     specie.name()
        #     for specie in self.supported_properties[
        #         PropTypes.aqueous_prop
        #     ].saturationSpecies()
        # ]

    def export_config(self):
        export_object = ReaktoroOutputExport()
        export_object.copy_rkt_outputs(self.rkt_outputs)
        export_object.copy_user_outputs(self.user_outputs)
        export_object.output_limits = copy.deepcopy(self.output_limits)
        return export_object

    def load_from_export_object(self, export_object):
        self.rkt_outputs = export_object.rkt_outputs
        self.user_outputs = export_object.user_outputs
        self.output_limits = export_object.output_limits

        def get_converted_function(obj):
            """get converted function for converted properties"""
            prop_info = getattr(
                self.supported_properties[PropTypes.converted_prop],
                obj.converted_prop_type,
            )(obj.property_index)
            obj.conversion_function = prop_info.conversion_function
            obj.derivative_conversion_function = (
                prop_info.derivative_conversion_function
            )
            print(obj, obj.conversion_function, obj.derivative_conversion_function)
            return obj.conversion_function

        for key, obj in self.rkt_outputs.items():
            property_type, get_function = self.get_prop_type(
                obj.property_name,
                obj.property_index,
            )
            assert property_type == obj.property_type
            print(key, property_type, obj.converted_prop_type)
            if obj.converted_prop_type is not None:
                get_converted_function(obj)
            obj.set_option_function(property_type, get_function)
        for key, obj in self.user_outputs.items():
            property_type, get_function = self.get_prop_type(
                obj.property_name,
                obj.property_index,
            )
            if obj.converted_prop_type is not None:
                get_converted_function(obj)

            assert property_type == obj.property_type
            obj.set_option_function(property_type, get_function)

    #### start of possible call function to extract values from reactoro properties #####

    def _func_tester(self, func, prop_type, prop_name, prop_index):
        """test function for reaktoro properties,
        The props can return errors due to not fully configured state, as such
        we want to accept those props as real, only if they specify explicit runtime error
        """
        try:
            value = func(prop_type, prop_name, prop_index)

        except RuntimeError as error:

            if "Unable to interpolate" not in str(error):
                raise error

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
