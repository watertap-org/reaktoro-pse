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
        self.value_clipped = False
        self.io_type = None
        self.calculate_value = None
        self.calculate_derivative_conversion = None
        self.calculation_options = None
        self.derivative = None

    def set_derivative(self, value):
        """set derivative value if available, used for jacobian calculations"""
        self.derivative = value

    def get_calculated_jacobian_value(self):
        """returns calculated jacobian value if available"""
        return self.calculation_options.calculate_derivative_conversion(
            self.calculation_options.properties
        )

    def compute_values(self, prop_object, supported_props):
        if self.calculation_options is None:
            return self.get_function(
                prop_object, self.property_name, self.property_index
            )
        else:
            for idx, prop in self.calculation_options.properties.items():
                prop.get_value(supported_props[prop.property_type], update_values=True)
            value = self.calculation_options.calculate_value(
                self.calculation_options.properties
            )
            return value

    def get_value(
        self,
        prop_object,
        update_values=False,
        supported_props=None,
    ):
        value = self.compute_values(prop_object, supported_props)
        self.value_clipped = False
        if update_values:
            self.value = value

        return value

    def delete_pyomo_var(self):
        # self.update_values()
        del self.pyomo_var
        self.pyomo_var = None

    def remove_unpicklable_data(self):
        self.delete_pyomo_var()
        if self.property_type == PropTypes.pyomo_built_prop:
            del self.pyomo_build_options
        if self.property_type == PropTypes.converted_prop:
            del self.calculation_options
            self.calculation_options = None
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


class PropOptions:
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

    def scalingTendencyPyomo(self, property_index):
        """build pyomo constraint for scaling index calculations directly form chem props
        #TODO: Need to add check for database being used as only PhreeqC is really supported at the
        moment"""
        required_props = PropOptions()
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
        system_species = [s.name() for s in self.state.state.system().species()]
        for s, mol in spec.reaction().reactants():
            if s.name() in system_species:
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

    def osmoticPressurePyomo(self, property_index):
        """build osmotic pressure constraint, as its not available from reaktoro"""
        # reference  https://help.syscad.net/PHREEQC_Reverse_Osmosis
        required_props = PropOptions()
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

    def vaporPressurePyomo(self, property_index=None):
        """build direct pH caclautions from chem props"""
        required_props = PropOptions()
        required_props.register_property(
            PropTypes.chem_prop, "speciesActivityLn", property_index
        )
        required_props.register_build_function(
            propFuncs.build_vapor_pressure_constraint
        )
        return required_props


class ConvertedPropTypes:
    def __init__(self, reaktor_state, chem_props, aqueous_props):
        self.state = reaktor_state
        self.chem_props = chem_props
        self.aqueous_props = aqueous_props

    def vaporPressure(self, property_index=None):
        """build vapor pressure"""
        output = PropOptions()
        output.register_property(
            PropTypes.chem_prop, "speciesActivityLn", property_index
        )
        output.calculate_value = (
            lambda x: math.exp(x["speciesActivityLn", property_index].value) * 101325
        )

        output.calculate_derivative_conversion = (
            lambda x: math.exp(x["speciesActivityLn", property_index].value)
            * 101325
            * x["speciesActivityLn", property_index].derivative
        )

        return output

    def osmoticPressure(self, property_index=None):
        """build osmotic pressure"""
        output = PropOptions()
        output.register_property(
            PropTypes.chem_prop, "speciesStandardVolume", property_index
        )
        output.register_property(
            PropTypes.chem_prop, "speciesActivityLn", property_index
        )
        output.register_property(PropTypes.chem_prop, "temperature")
        output.register_option("gas_constant", rkt.universalGasConstant)

        def calc_pressure(x):
            build_options = output.options
            return (
                -x[("speciesActivityLn", property_index)].value
                * build_options["gas_constant"]
                * x[("temperature", None)].value
                / x[("speciesStandardVolume", property_index)].value
            )

        def der_calc_pressure(x):
            build_options = output.options
            return sum(
                [
                    -x[("speciesActivityLn", property_index)].derivative
                    * build_options["gas_constant"]
                    * x[("temperature", None)].value
                    / x[("speciesStandardVolume", property_index)].value,
                    -x[("speciesActivityLn", property_index)].value
                    * build_options["gas_constant"]
                    * x[("temperature", None)].derivative
                    / x[("speciesStandardVolume", property_index)].value,
                    x[("speciesActivityLn", property_index)].value
                    * build_options["gas_constant"]
                    * x[("temperature", None)].value
                    / (x[("speciesStandardVolume", property_index)].value ** 2)
                    * x[("speciesStandardVolume", property_index)].derivative,
                ]
            )

        output.calculate_value = calc_pressure
        output.calculate_derivative_conversion = der_calc_pressure
        return output

    def elementAmount(self, property_index):
        """build element amount"""
        output = PropOptions()
        for mol, spc in self.state.element_to_species[property_index]:
            output.register_property(
                property_type=PropTypes.chem_prop,
                property_name="speciesAmount",
                property_index=spc,
            )
        output.calculate_value = lambda x: sum(
            mol * x["speciesAmount", spc].value
            for mol, spc in self.state.element_to_species[property_index]
        )
        output.calculate_derivative_conversion = lambda x: sum(
            mol * x["speciesAmount", spc].derivative
            for mol, spc in self.state.element_to_species[property_index]
        )
        return output

    def charge(self, property_index):
        """build element amount"""
        output = PropOptions()
        species = []
        for specie in self.state.state.system().species():
            name, charge = specie.name(), specie.charge()
            if charge != 0:
                species.append((charge, name))
                output.register_property(
                    property_type=PropTypes.chem_prop,
                    property_name="speciesAmount",
                    property_index=name,
                )
        output.calculate_value = lambda x: sum(
            charge * x["speciesAmount", spc].value for charge, spc in species
        )

        output.calculate_derivative_conversion = lambda x: sum(
            charge * x["speciesAmount", spc].derivative for charge, spc in species
        )
        return output

    def alkalinityAsCaCO3(self, property_index=None):
        """build alkalinity and convert it to CaCO3 basis"""
        output = PropOptions()
        output.register_property(
            property_type=PropTypes.aqueous_prop,
            property_name="alkalinity",
            property_index=None,
        )
        output.calculate_value = lambda x: (x["alkalinity", None].value * 100.09 * 1000)
        output.calculate_derivative_conversion = lambda x: (
            x["alkalinity", None].derivative * 100.09 * 1000
        )
        return output

    def scalingTendencySaturationIndex(self, property_index):
        """build scaling tendency - RKT has saturationIndex but no scalingIndex"""
        output = PropOptions()
        output.register_property(
            property_type=PropTypes.aqueous_prop,
            property_name="saturationIndex",
            property_index=property_index,
        )

        def calc_sat_value(x):
            try:
                return 10 ** (x["saturationIndex", property_index].value)
            except OverflowError:
                print("overflow error in scalingTendency calc_sat_value")
                return 1e100

        def calc_sat_dir(x):
            try:
                return (
                    x["saturationIndex", property_index].derivative
                    * (10 ** x["saturationIndex", property_index].value)
                    * math.log(10)
                )
            except OverflowError:
                print("overflow error in scalingTendency calc_sat_dir")
                return 1e100

        output.calculate_value = calc_sat_value
        output.calculate_derivative_conversion = calc_sat_dir
        return output

    def pH(self, property_index):
        """build log species amount"""
        output = PropOptions()
        output.register_property(
            property_type=PropTypes.chem_prop,
            property_name="speciesActivityLn",
            property_index="H+",
        )
        output.calculate_value = (
            lambda x: -1 * x["speciesActivityLn", "H+"].value / math.log(10)
        )
        output.calculate_derivative_conversion = (
            lambda x: x["speciesActivityLn", "H+"].derivative * -1 / math.log(10)
        )

        return output

    def scalingTendency(self, property_index):
        """build pyomo constraint for scaling index calculations directly form chem props
        #TODO: Need to add check for database being used as only PhreeqC is really supported at the
        moment"""
        # print("building sc direct")
        output = PropOptions()
        ref_temp = 25  # degC
        ref_pressure = 1  # atm
        spec = self.aqueous_props.saturationSpecies().get(property_index)
        reactant_species = spec.reaction().reactants()
        thermo_model = spec.standardThermoModel()
        pr = spec.props(ref_temp, "C", ref_pressure, "atm")
        specie_volume = float(pr.V0)  # returns auto diff/not usable with pyomo

        # get data from thermo prop
        jsp = thermo_model.params().dumpJson()
        jsp_dict = json.loads(jsp)
        not_implemented = False
        # TODO: need to add pE calculation to be able to calc scaling tendcies for these props.
        system_species = [s.name() for s in self.state.state.system().species()]
        all_species_exists = True
        for s, _ in reactant_species:
            if s.name() not in system_species:
                all_species_exists = False
                _log.warning(
                    f"Species {s.name()} not found in system species for index {property_index}, returning numerical scalingTendencySaturationIndex instead"
                )
        if isinstance(jsp_dict, list) and all_species_exists:
            if jsp_dict[0].get("PhreeqcLgK", None) is not None:
                output.register_option("logk_type", "Analytical")
                output.register_option("logk_paramters", jsp_dict[0]["PhreeqcLgK"])
            elif jsp_dict[0].get("VantHoff", None) is not None:
                output.register_option("logk_type", "VantHoff")
                output.register_option("logk_paramters", jsp_dict[0]["VantHoff"])
            else:
                not_implemented = True
        else:
            not_implemented = True
        if not_implemented:
            _log.warning(
                f"Exact derivatives for scaling tendency of {property_index} not implemented, returning numerical scalingTendencySaturationIndex instead"
            )
            return self.scalingTendencySaturationIndex(property_index)

        output.register_option("gas_constant", rkt.universalGasConstant)
        volume_reactants = 0
        for s, mol in reactant_species:
            if s.name() in system_species:
                spec = self.state.system.species().get(s.name())
                thermo_model = spec.standardThermoModel()
                _pr = spec.props(ref_temp, "C", ref_pressure, "atm")
                volume_reactants += float(_pr.V0) * (mol)
                output.register_property(
                    PropTypes.chem_prop, "speciesActivityLn", s.name()
                )
                output.properties[
                    ("speciesActivityLn", s.name())
                ].stoichiometric_coeff = mol  # create on demand to track coefficients

        output.register_option("delta_V", float(specie_volume - (volume_reactants)))
        output.register_property(PropTypes.chem_prop, "temperature")
        output.register_property(PropTypes.chem_prop, "pressure")

        def calc_scaling_tendency(x):
            build_options = output.options
            temperature_value = x["temperature", None].value
            if build_options["logk_type"] == "Analytical":
                A_params = build_options["logk_paramters"]
                log_k = [A_params["A1"]]
                # temp dependence for phreeqc
                if A_params["A2"] != 0:
                    log_k.append(A_params["A2"] * temperature_value)
                if A_params["A3"] != 0:
                    log_k.append(A_params["A3"] * temperature_value**-1)
                if A_params["A4"] != 0:
                    log_k.append(A_params["A4"] * math.log10(temperature_value))
                if A_params["A5"] != 0:
                    log_k.append(A_params["A5"] * temperature_value**-2)
                if A_params["A6"] != 0:
                    log_k.append(A_params["A6"] * temperature_value**2)

            if build_options["logk_type"] == "VantHoff":
                vfparams = build_options["logk_paramters"]
                log_k = [
                    vfparams["lgKr"]
                    - vfparams["dHr"]
                    / build_options["gas_constant"]
                    * (1 / temperature_value - 1 / vfparams["Tr"])
                    / math.log(10)
                ]
            # pressure dependence
            log_k.append(
                -(
                    build_options["delta_V"]
                    * (x[("pressure", None)].value - 101325)
                    / (math.log(10) * build_options["gas_constant"] * temperature_value)
                )
            )

            activities = []
            for key, obj in x.items():
                if "speciesActivityLn" in key:
                    activities.append(
                        obj.value * obj.stoichiometric_coeff / math.log(10)
                    )
            try:
                si = 10 ** (
                    sum(activities)
                    + sum(
                        log_k
                    )  # this is positive here and in log10 fom, so we add instead of subtract
                )
            except OverflowError:
                print("overflow error in scalingTendency calc_sat_dir")
                si = 1e100
            return si

        def calc_scaling_tendency_derivative(x):
            build_options = output.options
            temperature_value = x["temperature", None].value
            temperature_value_derivative = x["temperature", None].derivative
            if build_options["logk_type"] == "Analytical":
                # print("Analytical logk")
                A_params = build_options["logk_paramters"]
                log_k = [0]
                # temp dependence for phreeqc
                if A_params["A2"] != 0:
                    log_k.append(A_params["A2"] * temperature_value_derivative)
                if A_params["A3"] != 0:
                    log_k.append(
                        -A_params["A3"]
                        * temperature_value**-2
                        * temperature_value_derivative
                    )
                if A_params["A4"] != 0:
                    log_k.append(
                        A_params["A4"]
                        / (temperature_value * math.log(10))
                        * temperature_value_derivative
                    )
                if A_params["A5"] != 0:
                    log_k.append(
                        -2
                        * A_params["A5"]
                        * temperature_value**-3
                        * temperature_value_derivative
                    )
                if A_params["A6"] != 0:
                    log_k.append(
                        2
                        * A_params["A6"]
                        * temperature_value
                        * temperature_value_derivative
                    )

            if build_options["logk_type"] == "VantHoff":
                # print("VantHoff logk")
                vfparams = build_options["logk_paramters"]
                log_k = [
                    vfparams["dHr"]
                    / build_options["gas_constant"]
                    * (1 / temperature_value**2)
                    * temperature_value_derivative
                    / math.log(10)
                ]
            # pressure dependenance
            pressure_derivative = x[("pressure", None)].derivative
            log_k.append(
                -(
                    build_options["delta_V"]
                    * pressure_derivative
                    / (math.log(10) * build_options["gas_constant"] * temperature_value)
                )
            )
            log_k.append(
                (
                    build_options["delta_V"]
                    * (x[("pressure", None)].value - 101325)
                    / (
                        math.log(10)
                        * build_options["gas_constant"]
                        * temperature_value**2
                    )
                    * temperature_value_derivative
                )
            )
            activities = []
            for key, obj in x.items():
                if "speciesActivityLn" in key:
                    activities.append(
                        x[key].derivative * obj.stoichiometric_coeff / math.log(10)
                    )
            try:
                st = calc_scaling_tendency(x)
                sum_derivatives = sum(activities) + sum(
                    log_k
                )  # this is postive here and in log10 fom, so we add instead of subtract
                dir_st = sum_derivatives * st * math.log(10)
            except OverflowError:
                print("overflow error in scalingTendency calc_sat_dir")
                dir_st = 1e100
            return dir_st

        output.calculate_value = calc_scaling_tendency
        output.calculate_derivative_conversion = calc_scaling_tendency_derivative
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
            self.user_outputs[key].remove_unpicklable_data()


class ReaktoroOutputSpec:
    def __init__(self, reaktor_state=None):
        if reaktor_state is not None:
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
            self.prop_check_order = [
                PropTypes.converted_prop,
                PropTypes.pyomo_built_prop,
                PropTypes.chem_prop,
                PropTypes.aqueous_prop,
            ]
            self.rkt_outputs = {}  # outputs that reaktoro needs to generate
            self.user_outputs = {}  # outputs user requests

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
            property_type,
            update_values_in_object,
            self.supported_properties,
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

                self.user_outputs[index].io_type = prop_type
                for index, prop in get_function.properties.items():
                    # check if prop already exists if it does nor add it outputs
                    # otherwise overwrite it
                    if index not in self.rkt_outputs:
                        self.rkt_outputs[index] = prop
                    else:
                        get_function.properties[index] = self.rkt_outputs[index]

            elif property_type == PropTypes.converted_prop:
                # if converted prop, we need to get the converted prop type
                # and then set the get function

                self.user_outputs[index] = RktOutput(
                    property_type=property_type,
                    property_name=property_name,
                    property_index=property_index,
                    pyomo_var=pyomo_var,
                )

                self.user_outputs[index].calculation_options = get_function
                self.user_outputs[index].set_pyomo_var(pyomo_var)
                self.user_outputs[index].io_type = prop_type
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
        for supported_props in self.prop_check_order:
            prop = self.supported_properties.get(supported_props)
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
                        for prop_key, obj in func_results.properties.items():

                            supported_prop, func_result = self.get_prop_type(
                                obj.property_name, obj.property_index
                            )
                            obj.set_get_function(func_result)
                            obj.set_property_type(supported_prop)
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

    def export_config(self):
        export_object = ReaktoroOutputExport()
        export_object.copy_rkt_outputs(self.rkt_outputs)
        export_object.copy_user_outputs(self.user_outputs)
        return export_object

    def load_from_export_object(self, export_object):
        self.rkt_outputs = export_object.rkt_outputs
        self.user_outputs = export_object.user_outputs

        def get_converted_function(obj):
            """get converted function for converted properties"""
            prop_info = getattr(
                self.supported_properties[PropTypes.converted_prop],
                obj.property_name,
            )(obj.property_index)
            for prop_key, sub_obj in prop_info.properties.items():
                supported_prop, func_result = self.get_prop_type(
                    sub_obj.property_name, sub_obj.property_index
                )
                sub_obj.set_get_function(func_result)
                sub_obj.set_property_type(supported_prop)

            obj.calculation_options = prop_info

        for rkt_type in [self.rkt_outputs, self.user_outputs]:
            for key, obj in rkt_type.items():
                property_type, get_function = self.get_prop_type(
                    obj.property_name,
                    obj.property_index,
                )
                assert property_type == obj.property_type
                if obj.property_type == PropTypes.converted_prop:
                    get_converted_function(obj)
                else:
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
