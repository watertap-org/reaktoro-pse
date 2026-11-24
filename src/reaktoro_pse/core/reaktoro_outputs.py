#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
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
        property_sub_index=None,
        get_function=None,
        pyomo_var=None,
        value=None,
        stoichiometric_coeff=None,
        jacobian_type=None,
    ):
        self.property_type = property_type
        self.property_name = property_name  # properties from which to extract data
        self.property_index = property_index  # index if any
        self.property_sub_index = property_sub_index
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
                prop_object,
                self.property_name,
                self.property_index,
                self.property_sub_index,
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
        output.calculate_value = lambda x: (
            x["alkalinity", None].value * 100.09 * 1000 / 2
        )
        output.calculate_derivative_conversion = lambda x: (
            x["alkalinity", None].derivative * 100.09 * 1000 / 2
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

    def pE(self, property_index=None):
        """build pE from activities of a redox couple species (O2 or H2) with exact
        derivatives. Reaktoro does not expose an electron activity we can
        differentiate, but at equilibrium the couple formation reaction, e.g.
        O2 = 2 H2O - 4 H+ - 4 e-, gives

            log a(e-) = (log a(O2) - sum(nu_i log a_i) - log K(T, P)) / nu_e

        Only a couple above the epsilon floor is at equilibrium (O2 in oxidizing and
        H2 in reducing water), so the more abundant couple is used on each evaluation
        (both give the same value where both are present). Falls back to Reaktoro pE
        (numerical derivatives) when database has no couple species.
        property_index is not used (pE has no index)"""
        ref_temp = 25  # degC
        ref_pressure = 1  # atm
        system_species = [s.name() for s in self.state.state.system().species()]
        primaries = {
            s
            for s in system_species
            if not self.state.database.species().get(s).reaction().reactants()
        }
        couples = []
        for couple in ("O2", "H2"):
            if couple not in system_species:
                continue
            expanded = self._expand_reaction(couple, primaries)
            if expanded is None or "e-" not in expanded[0]:
                continue
            reactants, logk_terms = expanded
            nu_e = reactants.pop("e-")
            if not all(s in system_species for s in reactants):
                continue
            volume = float(
                self.state.system.species()
                .get(couple)
                .props(ref_temp, "C", ref_pressure, "atm")
                .V0
            )
            for s, mol in reactants.items():
                volume -= mol * float(
                    self.state.system.species()
                    .get(s)
                    .props(ref_temp, "C", ref_pressure, "atm")
                    .V0
                )
            couples.append(
                {
                    "couple": couple,
                    "reactants": reactants,
                    "nu_e": nu_e,
                    "logk_terms": logk_terms,
                    "delta_V": volume,
                }
            )
        if not couples:
            _log.warning(
                "No redox couple (O2 or H2) available for exact pE derivatives, using numerical pE"
            )
            # get_prop_type moves on to Reaktoro aqueous props pE
            raise AttributeError("no redox couple for exact pE")
        output = PropOptions()
        for c in couples:
            output.register_property(
                PropTypes.chem_prop, "speciesActivityLn", c["couple"]
            )
            for s in c["reactants"]:
                output.register_property(PropTypes.chem_prop, "speciesActivityLn", s)
        output.register_property(PropTypes.chem_prop, "temperature")
        output.register_property(PropTypes.chem_prop, "pressure")
        output.register_option("couples", couples)
        output.register_option("gas_constant", rkt.universalGasConstant)
        _log.info(
            f"pE is calculated from {[c['couple'] for c in couples]} redox couple(s)"
        )

        def get_couple(x):
            return max(
                output.options["couples"],
                key=lambda c: x[("speciesActivityLn", c["couple"])].value,
            )

        def calc_pE(x):
            c = get_couple(x)
            gas_constant = output.options["gas_constant"]
            temperature = x["temperature", None].value
            log_k = sum(
                mult * self._log_k_value(t, prm, temperature, gas_constant)
                for t, prm, mult in c["logk_terms"]
            ) - c["delta_V"] * (x[("pressure", None)].value - 101325) / (
                math.log(10) * gas_constant * temperature
            )
            log_a = (
                x[("speciesActivityLn", c["couple"])].value
                - sum(
                    mol * x[("speciesActivityLn", s)].value
                    for s, mol in c["reactants"].items()
                )
            ) / math.log(10)
            return -(log_a - log_k) / c["nu_e"]

        def calc_pE_derivative(x):
            c = get_couple(x)
            gas_constant = output.options["gas_constant"]
            temperature = x["temperature", None].value
            d_temperature = x["temperature", None].derivative
            d_log_k = sum(
                mult
                * self._log_k_derivative(
                    t, prm, temperature, d_temperature, gas_constant
                )
                for t, prm, mult in c["logk_terms"]
            )
            d_log_k -= (
                c["delta_V"]
                * x[("pressure", None)].derivative
                / (math.log(10) * gas_constant * temperature)
            )
            d_log_k += (
                c["delta_V"]
                * (x[("pressure", None)].value - 101325)
                / (math.log(10) * gas_constant * temperature**2)
                * d_temperature
            )
            d_log_a = (
                x[("speciesActivityLn", c["couple"])].derivative
                - sum(
                    mol * x[("speciesActivityLn", s)].derivative
                    for s, mol in c["reactants"].items()
                )
            ) / math.log(10)
            return -(d_log_a - d_log_k) / c["nu_e"]

        output.calculate_value = calc_pE
        output.calculate_derivative_conversion = calc_pE_derivative
        return output

    @staticmethod
    def _get_logk_term(thermo_model):
        """returns (logk_type, params) for a species standard thermo model, or None
        when the model type is not supported"""
        jsp_dict = json.loads(thermo_model.params().dumpJson())
        if not isinstance(jsp_dict, list) or jsp_dict[0] is None:
            return None
        if jsp_dict[0].get("PhreeqcLgK", None) is not None:
            return ("Analytical", jsp_dict[0]["PhreeqcLgK"])
        if jsp_dict[0].get("VantHoff", None) is not None:
            return ("VantHoff", jsp_dict[0]["VantHoff"])
        return None

    def _expand_reaction(self, species_name, keep, depth=0):
        """Formation reaction of species_name written in terms of species in keep
        (primary master species of the system) and e-, expanding secondary species
        recursively through their own formation reactions. Species that can not be
        expanded are kept as is.

        returns (reactants dict, list of (logk_type, params, multiplier)) or None"""
        spec = self.state.database.species().get(species_name)
        term = self._get_logk_term(spec.standardThermoModel())
        if term is None or depth > 10:
            return None
        reactants = {}
        logk_terms = [term + (1.0,)]
        for s, mol in spec.reaction().reactants():
            name, mol = s.name(), float(mol)
            sub = None
            if name not in keep and name != "e-" and s.reaction().reactants():
                sub = self._expand_reaction(name, keep, depth + 1)
            if sub is None:
                reactants[name] = reactants.get(name, 0.0) + mol
                continue
            sub_reactants, sub_terms = sub
            for sub_name, sub_mol in sub_reactants.items():
                reactants[sub_name] = reactants.get(sub_name, 0.0) + mol * sub_mol
            logk_terms += [(t, prm, mult * mol) for t, prm, mult in sub_terms]
        reactants = {k: v for k, v in reactants.items() if abs(v) > 1e-12}
        return reactants, logk_terms

    def _get_reaction_representations(self, property_index, system_species):
        """Reaction representations of a mineral that only use activities of system
        species, so scaling tendency has exact derivatives.

        The reaction is expanded to primary master species, as secondary species
        (e.g. H2S in oxic water) can sit at the Reaktoro epsilon floor, where their
        activity is not at equilibrium. Reaktoro does not expose an electron activity
        we can differentiate, so e- is replaced with a redox couple species (O2 or H2)
        through its formation reaction (e.g. O2 = 2 H2O - 4 H+ - 4 e-):

            reaction' = reaction - (nu_e / nu_e,couple) * couple reaction

        At equilibrium both couples give the same saturation index, but only one of
        them is above the epsilon floor (oxidizing water has O2, reducing water H2),
        so a representation is built for each couple and the one with the larger
        couple activity is used during evaluation.

        returns list of (couple or None, reactants dict, logk terms), empty if exact
        derivatives are not possible"""
        primaries = {
            s
            for s in system_species
            if not self.state.database.species().get(s).reaction().reactants()
        }
        expanded = self._expand_reaction(property_index, primaries)
        if expanded is None:
            return []
        reactants, logk_terms = expanded
        if "e-" not in reactants:
            if all(s in system_species for s in reactants):
                return [(None, reactants, logk_terms)]
            return []
        representations = []
        for couple in ("O2", "H2"):
            if couple not in system_species:
                continue
            couple_expanded = self._expand_reaction(couple, primaries)
            if couple_expanded is None or "e-" not in couple_expanded[0]:
                continue
            couple_reactants, couple_terms = couple_expanded
            factor = reactants["e-"] / couple_reactants["e-"]
            new_reactants = dict(reactants)
            for s, mol in couple_reactants.items():
                new_reactants[s] = new_reactants.get(s, 0.0) - factor * mol
            # couple species moves from product to reactant side of reaction
            new_reactants[couple] = new_reactants.get(couple, 0.0) + factor
            new_reactants.pop("e-")
            new_reactants = {s: v for s, v in new_reactants.items() if abs(v) > 1e-12}
            if not all(s in system_species for s in new_reactants):
                continue
            representations.append(
                (
                    couple,
                    new_reactants,
                    logk_terms
                    + [(t, prm, -factor * mult) for t, prm, mult in couple_terms],
                )
            )
        if representations:
            _log.info(
                f"Scaling tendency for {property_index} uses {[r[0] for r in representations]} redox couple(s) in place of e-"
            )
        return representations

    @staticmethod
    def _log_k_value(logk_type, params, temperature_value, gas_constant):
        """log10 K at temperature for PhreeqcLgK (Analytical) or VantHoff params"""
        if logk_type == "Analytical":
            log_k = params["A1"]
            if params["A2"] != 0:
                log_k += params["A2"] * temperature_value
            if params["A3"] != 0:
                log_k += params["A3"] * temperature_value**-1
            if params["A4"] != 0:
                log_k += params["A4"] * math.log10(temperature_value)
            if params["A5"] != 0:
                log_k += params["A5"] * temperature_value**-2
            if params["A6"] != 0:
                log_k += params["A6"] * temperature_value**2
            return log_k
        return params["lgKr"] - params["dHr"] / gas_constant * (
            1 / temperature_value - 1 / params["Tr"]
        ) / math.log(10)

    @staticmethod
    def _log_k_derivative(logk_type, params, temperature_value, dT, gas_constant):
        """derivative of log10 K given temperature derivative dT"""
        if logk_type == "Analytical":
            d_log_k = 0
            if params["A2"] != 0:
                d_log_k += params["A2"] * dT
            if params["A3"] != 0:
                d_log_k += -params["A3"] * temperature_value**-2 * dT
            if params["A4"] != 0:
                d_log_k += params["A4"] / (temperature_value * math.log(10)) * dT
            if params["A5"] != 0:
                d_log_k += -2 * params["A5"] * temperature_value**-3 * dT
            if params["A6"] != 0:
                d_log_k += 2 * params["A6"] * temperature_value * dT
            return d_log_k
        return (
            params["dHr"]
            / gas_constant
            * (1 / temperature_value**2)
            * dT
            / math.log(10)
        )

    def scalingTendency(self, property_index):
        """build scaling tendency (10**saturation index) directly from species
        activities, which have exact derivatives. Redox minerals (e- in reaction) are
        supported through a redox couple (see _get_reaction_representations), falls back
        to numerical saturation index when no representation is available
        #TODO: Need to add check for database being used as only PhreeqC is really supported at the
        moment"""
        output = PropOptions()
        ref_temp = 25  # degC
        ref_pressure = 1  # atm
        spec = self.aqueous_props.saturationSpecies().get(property_index)
        pr = spec.props(ref_temp, "C", ref_pressure, "atm")
        specie_volume = float(pr.V0)  # returns auto diff/not usable with pyomo

        system_species = [s.name() for s in self.state.state.system().species()]
        representations = self._get_reaction_representations(
            property_index, system_species
        )
        if not representations:
            _log.warning(
                f"Exact derivatives for scaling tendency of {property_index} not implemented, returning numerical scalingTendencySaturationIndex instead"
            )
            return self.scalingTendencySaturationIndex(property_index)

        # each representation: couple, reactant coefficients, log K terms, delta V
        reaction_options = []
        for couple, reactants, logk_terms in representations:
            volume_reactants = 0
            for s, mol in reactants.items():
                _pr = (
                    self.state.system.species()
                    .get(s)
                    .props(ref_temp, "C", ref_pressure, "atm")
                )
                volume_reactants += float(_pr.V0) * mol
                output.register_property(PropTypes.chem_prop, "speciesActivityLn", s)
            reaction_options.append(
                {
                    "couple": couple,
                    "reactants": reactants,
                    "logk_terms": logk_terms,
                    "delta_V": float(specie_volume - volume_reactants),
                }
            )
        output.register_option("reactions", reaction_options)
        output.register_option("gas_constant", rkt.universalGasConstant)
        output.register_property(PropTypes.chem_prop, "temperature")
        output.register_property(PropTypes.chem_prop, "pressure")

        def get_reaction(x):
            """reaction with the most abundant redox couple (smooth switch, as both
            give same value and derivatives at equilibrium)"""
            reactions = output.options["reactions"]
            if len(reactions) == 1:
                return reactions[0]
            return max(
                reactions,
                key=lambda r: x[("speciesActivityLn", r["couple"])].value,
            )

        def calc_log_si(x, reaction):
            gas_constant = output.options["gas_constant"]
            temperature_value = x["temperature", None].value
            log_k = sum(
                multiplier
                * self._log_k_value(logk_type, params, temperature_value, gas_constant)
                for logk_type, params, multiplier in reaction["logk_terms"]
            )
            # pressure dependence
            log_k -= (
                reaction["delta_V"]
                * (x[("pressure", None)].value - 101325)
                / (math.log(10) * gas_constant * temperature_value)
            )
            activities = sum(
                mol * x[("speciesActivityLn", s)].value / math.log(10)
                for s, mol in reaction["reactants"].items()
            )
            # log K is in formation direction here, so we add instead of subtract
            return activities + log_k

        def calc_scaling_tendency(x):
            try:
                return 10 ** calc_log_si(x, get_reaction(x))
            except OverflowError:
                print("overflow error in scalingTendency calc_sat_dir")
                return 1e100

        def calc_scaling_tendency_derivative(x):
            reaction = get_reaction(x)
            gas_constant = output.options["gas_constant"]
            temperature_value = x["temperature", None].value
            temperature_value_derivative = x["temperature", None].derivative
            d_log_k = sum(
                multiplier
                * self._log_k_derivative(
                    logk_type,
                    params,
                    temperature_value,
                    temperature_value_derivative,
                    gas_constant,
                )
                for logk_type, params, multiplier in reaction["logk_terms"]
            )
            # pressure dependence
            d_log_k -= (
                reaction["delta_V"]
                * x[("pressure", None)].derivative
                / (math.log(10) * gas_constant * temperature_value)
            )
            d_log_k += (
                reaction["delta_V"]
                * (x[("pressure", None)].value - 101325)
                / (math.log(10) * gas_constant * temperature_value**2)
                * temperature_value_derivative
            )
            d_activities = sum(
                mol * x[("speciesActivityLn", s)].derivative / math.log(10)
                for s, mol in reaction["reactants"].items()
            )
            try:
                st = 10 ** calc_log_si(x, reaction)
                return (d_activities + d_log_k) * st * math.log(10)
            except OverflowError:
                print("overflow error in scalingTendency calc_sat_dir")
                return 1e100

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
        property_sub_index=None,
        get_all_indexes=False,
        pyomo_var=None,
        ignore_indexes=None,
        multi_sub_index=False,
    ):
        """register a reaktoro output, couple it to property type.

        Keywords:
        property_name -- prop name (specieisActivityLn, pH etc)
        property_index -- prop index if any (H+, etc) (default: None)
        property_sub_index -- second index for dual-argument props
        get_all_indexes -- if user want to get all possible indexs for specfied prop (default: False)
        pyomo_var -- pyomo var that should be used for the output of this property (optional: will be auto built) (default: None)
        """
        if get_all_indexes:
            self.get_all_indexes(property_name, ignore_indexes)
        else:
            property_type, get_function = self.get_prop_type(
                property_name, property_index, property_sub_index
            )
            self.process_output(
                property_type=property_type,
                property_name=property_name,
                property_index=property_index,
                property_sub_index=property_sub_index,
                get_function=get_function,
                pyomo_var=pyomo_var,
                multi_sub_index=multi_sub_index,
            )

    def process_output(
        self,
        property_type,
        property_name,
        property_index=None,
        property_sub_index=None,
        get_function=None,
        pyomo_var=None,
        multi_sub_index=False,
    ):
        if multi_sub_index or property_sub_index is not None:
            index = (property_name, property_index, property_sub_index)
            multi_sub_index = True
        else:
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
                    property_sub_index=property_sub_index,
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
                    property_sub_index=property_sub_index,
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
                    property_sub_index=property_sub_index,
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

    def get_prop_type(
        self, property_name, property_index=None, property_sub_index=None
    ):
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
                            property_sub_index,
                        )
                        return supported_props, func_attempt
                except (TypeError, KeyError, AttributeError, RuntimeError):
                    pass
        raise NotImplementedError(
            f"""The {property_name}, {property_index}, {property_sub_index} was not found,
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

    def _func_tester(self, func, prop_type, prop_name, prop_index, property_sub_index):
        """test function for reaktoro properties,
        The props can return errors due to not fully configured state, as such
        we want to accept those props as real, only if they specify explicit runtime error
        """
        try:
            value = func(prop_type, prop_name, prop_index, property_sub_index)

        except RuntimeError as error:
            if "Unable to interpolate" not in str(error):
                raise error

    def _get_prop_phase_name_val(
        self, prop_type, prop_name, prop_index, property_sub_index=None
    ):
        """get prop based on phase, used for chem_props.phaseProp"""
        value = getattr(prop_type.phaseProps(prop_index), prop_name)()
        return float(value)

    def _get_reaktoro_function(
        self, prop_type, prop_name, prop_index=None, property_sub_index=None
    ):
        """get reaktoro function result, handling None, tuple, and single index"""
        if prop_index is None:
            return getattr(prop_type, prop_name)()
        elif property_sub_index is not None:
            return getattr(prop_type, prop_name)(prop_index, property_sub_index)
        else:
            return getattr(prop_type, prop_name)(prop_index)

    def _get_prop_name_val(
        self, prop_type, prop_name, prop_index=None, property_sub_index=None
    ):
        """get prop based on name/index and execute value call"""
        return float(
            self._get_reaktoro_function(
                prop_type, prop_name, prop_index, property_sub_index
            ).val()
        )

    def _get_prop_name(
        self, prop_type, prop_name, prop_index=None, property_sub_index=None
    ):
        """get prop based/index on name only"""
        return float(
            self._get_reaktoro_function(
                prop_type, prop_name, prop_index, property_sub_index
            )
        )
