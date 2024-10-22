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

import numpy as np
from reaktoro_pse.core.reaktoro_outputs import PropTypes
from reaktoro_pse.core.reaktoro_state import ReaktoroState
from reaktoro_pse.core.util_classes.rkt_inputs import RktInputTypes
from reaktoro_pse.core.reaktoro_outputs import (
    ReaktoroOutputSpec,
)

import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)

__author__ = "Alexander V. Dudchenko, Paul Vecchiarelli, Ben Knueven"


# class to setup jacobian for reaktoro


class JacType:
    average = "average"
    center_difference = "center_difference"
    exact = "exact"
    numeric = "numeric"


class JacboianRows:
    def __init__(self, rkt_state):
        # useing lists instead of dict to preserve order"""
        self.property = []
        self.property_index = []
        self.keys = []
        self.available_in_props = []
        self.get_functions = []
        self.absolute_values = []
        self.standard_keys = []
        self.get_rows(rkt_state)

    def get_rows(self, state):
        # we would need to manually to add gas and solid phase
        # this hsould not include things like Calcite
        _jac_species = [species.name() for species in state.system().species()]
        _jac_phases = [phase.name() for phase in state.system().phases()]
        # state temperature and pressure
        self.keys = [("temperature",), ("pressure",)]
        # _jac_species amounts * # _jac_species
        self.keys.extend([("speciesAmount", s) for s in _jac_species])
        # surface area of reacting surfaces (probably will be here if kinetics are used)
        # temperature of each phase * # _jac_phases
        self.keys.extend([("temperature", p) for p in _jac_phases])
        # pressure of each phase * # _jac_phases
        self.keys.extend([("pressure", p) for p in _jac_phases])
        # sum speciesAmounts * # _jac_phases
        self.keys.extend([("amount", p) for p in _jac_phases])
        # sum speciesMasses * # _jac_phases
        self.keys.extend([("mass", p) for p in _jac_phases])
        # mole fractions * # species
        self.keys.extend([("speciesMoleFraction", s) for s in _jac_species])
        # standard molar gibbs free energy * # species
        self.keys.extend([("speciesStandardGibbsEnergy", s) for s in _jac_species])
        # standard molar enthalpy * # species
        self.keys.extend([("speciesStandardEnthalpy", s) for s in _jac_species])
        # standard molar volumes * # species
        self.keys.extend([("speciesStandardVolume", s) for s in _jac_species])
        # temp derivative of standard molar volumes * # species
        self.keys.extend([("speciesStandardVolumeT", s) for s in _jac_species])
        # pres derivative of standard molar volumes * # species
        self.keys.extend([("speciesStandardVolumeP", s) for s in _jac_species])
        # standard isobaric molar heat capacities * # species
        self.keys.extend(
            [("speciesStandardHeatCapacityConstP", s) for s in _jac_species]
        )
        # corrective molar volume * # _jac_phases
        self.keys.extend([("correctiveMolarVolume", p) for p in _jac_phases])
        # temp derivative of corrective molar volume * # _jac_phases
        self.keys.extend([("correctiveMolarVolumeT", p) for p in _jac_phases])
        # pres derivative of corrective molar volume * # _jac_phases
        self.keys.extend([("correctiveMolarVolumeP", p) for p in _jac_phases])
        # mole frac derivative of corrective molar volume * # _jac_phases
        self.keys.extend(
            [("speciesCorrectiveMolarVolume", s) for s in _jac_species]
        )  # not sure what this is
        # corrective gibbs free energy * # _jac_phases
        self.keys.extend([("correctiveMolarGibbsEnergy", p) for p in _jac_phases])
        # corrective molar enthalpy * # _jac_phases
        self.keys.extend([("correctiveMolarEnthalpy", p) for p in _jac_phases])
        # corrective isobaric molar heat capacity * # _jac_phases
        self.keys.extend(
            [("correctiveMolarHeatCapacityConstP", p) for p in _jac_phases]
        )
        # activity coefficients (ln) * # species
        self.keys.extend([("speciesActivityCoefficientLn", s) for s in _jac_species])
        # activities (ln) * # species
        self.keys.extend([("speciesActivityLn", s) for s in _jac_species])
        # chemical potentials * # species
        self.keys.extend([("speciesChemicalPotential", s) for s in _jac_species])

        for idx, row in [(idx, row) for idx, row in enumerate(self.keys)]:
            self.available_in_props.append(False)
            if len(row) == 2:
                index = row[1]
                prop = row[0]
            else:
                prop = row[0]
                index = None
            self.absolute_values.append(0)
            self.property.append(prop)
            self.property_index.append(index)
            self.get_functions.append(None)
            self.standard_keys.append((prop, index))

    def get_index(self, key, index=None):
        if isinstance(key, tuple) == False:
            key = (key, index)
        idx = self.standard_keys.index(key)
        return idx

    def get_key(self, index):
        return self.keys[index]

    def get_property(self, key):
        idx = self.get_index(key)
        return self.property[idx], self.property_index[idx]

    def check_key(self, key, index=None):
        if (key, index) in self.standard_keys:
            return self.get_index((key, index))
        else:
            return False

    def set_availableInProps(self, key, available):
        idx = self.get_index(key)
        self.available_in_props[idx] = available

    def set_get_function(self, key, function):
        """set eval function"""
        idx = self.get_index(key)
        self.get_functions[idx] = function

    def compute_value(self, prop, key):
        idx = self.get_index(key)
        if self.get_functions[idx] is None:
            val = 0
        else:
            val = self.get_functions[idx](
                prop, self.property[idx], self.property_index[idx]
            )
        self.absolute_values[idx] = val

    def get_value(self, key):
        """eval jacobian using specified function other wise return 0"""
        idx = self.get_index(key)
        return self.absolute_values[idx]

    def display_available(self):
        available_dict = {}
        _log.info("-----displaying available jacobian values-----")
        for i, key in enumerate(self.standard_keys):
            _log.info(f"{key}: {self.available_in_props[i]}")
            available_dict[key] = self.available_in_props[i]
        _log.info("-----done-----")

        return available_dict


class ReaktoroJacobianExport:
    def __init__(self):
        self.der_step_size = None
        self.jacobian_type = None
        self.numerical_order = None


class ReaktoroJacobianSpec:
    def __init__(self, reaktor_state, reaktor_outputs):
        self.state = reaktor_state
        if isinstance(self.state, ReaktoroState) == False:
            raise TypeError("Reator jacobian require rektoroState class")
        self.output_specs = reaktor_outputs
        if isinstance(self.output_specs, ReaktoroOutputSpec) == False:
            raise TypeError("Reator outputs require ReaktoroOutputSpec class")
        self.jac_rows = JacboianRows(self.state.state)
        self.set_jacobian_type()
        self.configure_numerical_jacobian()
        self.check_existing_jacobian_props()

    def export_config(self):
        export_object = ReaktoroJacobianExport()
        export_object.der_step_size = self.der_step_size
        export_object.jacobian_type = self.jacobian_type
        export_object.numerical_order = self.numerical_order
        return export_object

    def load_from_export_object(self, export_object):
        self.configure_numerical_jacobian(
            jacobian_type=export_object.jacobian_type,
            order=export_object.numerical_order,
            step_size=export_object.der_step_size,
        )

    def configure_numerical_jacobian(
        self, jacobian_type="average", order=4, step_size=1e-8
    ):
        """Configure numerical derivate options

        Keywords:
        derivative_type -- type of derivative method (average, or center_difference)
        order -- order of derivative
        step_size -- numerical step size for approximating derivative
        """
        self.der_step_size = step_size
        self.jacobian_type = JacType.average
        self.numerical_order = order
        if jacobian_type == JacType.average:
            self.jacobian_type = JacType.average
            assert order % 2 == 0
            self.numerical_steps = np.arange(-order / 2, order / 2 + 1) * step_size
        if jacobian_type == JacType.center_difference:
            self.jacobian_type = JacType.center_difference
            self.center_diff_order(order)
        self.set_up_chem_and_aq_states()

    def set_up_chem_and_aq_states(self):
        self.chem_prop_states = []
        self.aqueous_prop_states = []
        for step in self.numerical_steps:
            self.chem_prop_states.append(rkt.ChemicalProps(self.state.state))
            if RktInputTypes.aqueous_phase in self.state.inputs.registered_phases:
                self.aqueous_prop_states.append(
                    rkt.AqueousProps(self.chem_prop_states[-1])
                )

    def update_states(self, new_states):
        for i in range(len(self.numerical_steps)):
            self.chem_prop_states[i].update(new_states[:, i])
            if RktInputTypes.aqueous_phase in self.state.inputs.registered_phases:
                self.aqueous_prop_states[i].update(self.chem_prop_states[i])

    def get_state_values(self, output_object):
        output_vals = []
        for i, step in enumerate(self.numerical_steps):
            if output_object.property_type == PropTypes.chem_prop:
                prop = self.chem_prop_states[i]
            elif output_object.property_type == PropTypes.aqueous_prop:
                prop = self.aqueous_prop_states[i]
            else:
                raise TypeError(
                    f"{output_object.property_type} not supported by numerical derivative method, please update"
                )

            val = self.output_specs.evaluate_property(output_object, prop)
            output_vals.append(val)
        return output_vals

    def set_jacobian_type(self):
        """function to check all inputs and identify if exact or numeric jacobian should be used"""
        _jac_types = []
        for key, obj in self.output_specs.rkt_outputs.items():
            jac_available = self.jac_rows.check_key(
                obj.property_name, obj.property_index
            )
            if jac_available is not False:
                obj.jacobian_type = JacType.exact
            else:
                obj.jacobian_type = JacType.numeric
            _jac_types.append(f"{key}: {obj.jacobian_type}")

    def display_jacobian_output_types(self):
        """used for displaying jac output types"""
        out_types = {}
        _log.info("-----displaying jacobian outputs and types-----")
        for key, obj in self.output_specs.rkt_outputs.items():
            _log.info(f"{key}: Jac type: {obj.jacobian_type}")
            out_types[key] = obj.jacobian_type
        _log.info("-----done-----")
        return out_types

    def check_existing_jacobian_props(self):
        """check if jacobian is available in current properties and identify function to use
        for getting values"""
        for idx, key in enumerate(self.jac_rows.standard_keys):
            prop, prop_index = self.jac_rows.get_property(key)
            try:
                prop_type, get_function = self.output_specs.get_prop_type(
                    prop, prop_index
                )
                if prop_type == PropTypes.chem_prop:
                    self.jac_rows.set_availableInProps(key, True)
                    self.jac_rows.set_get_function(key, get_function)
                else:
                    self.jac_rows.set_availableInProps(key, False)
            except NotImplementedError:
                self.jac_rows.set_availableInProps(key, False)

    def update_jacobian_absolute_values(self):
        """' used to update absolute values of jacobian - needed for numeric derivatives"""
        prop = self.output_specs.supported_properties[PropTypes.chem_prop]
        for jacIdx, key in enumerate(self.jac_rows.standard_keys):
            self.jac_rows.compute_value(prop, key)

    def process_jacobian_matrix(self, jacobian_matrix, input_index, input_value):
        """this function is used to pull out a specific column from the jacobian and also
        generate matrix for manually propagating derivatives
        Here we need to retain row order, as its same as input into chem properties"""
        jacobian_abs_matrix = []
        jacobian_dict = {}
        for jacIdx, key in enumerate(self.jac_rows.standard_keys):
            jacobian_abs_matrix.append([])
            jac_abs_value = self.jac_rows.get_value(key)
            jac_value = jacobian_matrix[jacIdx][input_index]
            jacobian_dict[key] = jac_value

            for step in self.numerical_steps:
                der_step = jac_value * input_value * step
                jacobian_abs_matrix[-1].append(jac_abs_value + der_step)
        return jacobian_dict, np.array(jacobian_abs_matrix)

    def get_jacobian(self, jacobian_matrix, input_object):
        input_index = input_object.get_jacobian_index()
        input_value = input_object.get_temp_value()
        jacobian_dict, jacobian_abs_matrix = self.process_jacobian_matrix(
            jacobian_matrix, input_index, input_value
        )
        self.update_states(jacobian_abs_matrix)

        output_jacobian = []
        for output_key, output_obj in self.output_specs.rkt_outputs.items():
            if output_obj.jacobian_type == JacType.exact:
                output_jacobian.append(jacobian_dict[output_key])
            else:
                values = self.get_state_values(output_obj)
                if JacType.average:
                    diff = np.diff(values)
                    step = np.diff(self.numerical_steps * input_value)
                    jac_val = np.average(diff / step)
                elif JacType.center_difference:
                    jac_val = np.array(values) * self.cdf_multipliers
                    jac_val = np.sum(jac_val) / (
                        self.rkt_aqueous_props_der_step * input_value
                    )
                output_jacobian.append(jac_val)
        return output_jacobian

    def center_diff_order(self, order):
        """paramterizes center different taylor series
        refer to https://en.wikipedia.org/wiki/Finite_difference_coefficient"""

        if order == 2:
            self.numerical_steps = (
                np.array(
                    [
                        -1,
                        1,
                    ]
                )
                * self.der_step_size
            )
            self.cdf_multipliers = np.array([-1 / 2, 1 / 2])
        if order == 4:
            self.numerical_steps = np.array([-2, -1, 1, 2]) * self.der_step_size
            self.cdf_multipliers = np.array([-1 / 12, 2 / 3, -2 / 3, 1 / 12])
        if order == 6:
            self.numerical_steps = np.array([-3, -2, -1, 1, 2, 3]) * self.der_step_size
            self.cdf_multipliers = np.array(
                [-1 / 60, 3 / 20, -3 / 4, 3 / 4, -3 / 20, 1 / 60]
            )
        if order == 8:
            self.numerical_steps = (
                np.array([-4, -3, -2, -1, 1, 2, 3, 4]) * self.der_step_size
            )
            self.cdf_multipliers = np.array(
                [1 / 280, -4 / 105, 1 / 5, -4 / 5, 4 / 5, -1 / 5, 4 / 105, -1 / 280]
            )
        if order == 10:
            self.numerical_steps = (
                np.array([-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]) * self.der_step_size
            )
            self.cdf_multipliers = (
                np.array([-2, 25, -150, 600, -2100, 2100, -600, 150, -25, 2]) / 2520
            )
