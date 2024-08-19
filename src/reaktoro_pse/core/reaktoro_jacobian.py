import reaktoro as rkt

import numpy as np
from reaktoro_pse.core.reaktoro_outputs import propTypes
from reaktoro_pse.core.reaktoro_state import reaktoroState
from reaktoro_pse.core.reaktoro_outputs import (
    reaktoroOutputSpec,
)

import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)

__author__ = "Alexander Dudchenko, Paul Vecchiarelli, Ben Knueven"


"""class to setup jacobian for reaktoro"""


class jacType:
    average = "average"
    center_difference = "center_difference"
    exact = "exact"
    numeric = "numeric"


class jacboianRows:
    def __init__(self, rkt_state):
        """useing lists instead of dict to preserve order"""
        self.property = []
        self.propertyIndex = []
        self.keys = []
        self.availableInProps = []
        self.getFunctions = []
        self.absoluteValues = []
        self.standardKeys = []
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
        self.keys.extend([("molarVolume", p) for p in _jac_phases])
        # temp derivative of corrective molar volume * # _jac_phases
        self.keys.extend([("molarVolumeT", p) for p in _jac_phases])
        # pres derivative of corrective molar volume * # _jac_phases
        self.keys.extend([("molarVolumeP", p) for p in _jac_phases])
        # mole frac derivative of corrective molar volume * # _jac_phases
        self.keys.extend(
            [("molarVolumeI", s) for s in _jac_species]
        )  # not sure what this is
        # corrective gibbs free energy * # _jac_phases
        self.keys.extend([("molarGibbsEnergy", p) for p in _jac_phases])
        # corrective molar enthalpy * # _jac_phases
        self.keys.extend([("molarEnthalpy", p) for p in _jac_phases])
        # corrective isobaric molar heat capacity * # _jac_phases
        self.keys.extend([("molarHeatCapacityConstP", p) for p in _jac_phases])
        # activity coefficients (ln) * # species
        self.keys.extend([("speciesActivityCoefficientLn", s) for s in _jac_species])
        # activities (ln) * # species
        self.keys.extend([("speciesActivityLn", s) for s in _jac_species])
        # chemical potentials * # species
        self.keys.extend([("speciesChemicalPotential", s) for s in _jac_species])

        for idx, row in [(idx, row) for idx, row in enumerate(self.keys)]:
            self.availableInProps.append(False)
            if len(row) == 2:
                index = row[1]
                prop = row[0]
            else:
                prop = row[0]
                index = None
            self.absoluteValues.append(0)
            self.property.append(prop)
            self.propertyIndex.append(index)
            self.getFunctions.append(None)
            self.standardKeys.append((prop, index))

    def get_index(self, key, index=None):
        if isinstance(key, tuple) == False:
            key = (key, index)
        idx = self.standardKeys.index(key)
        return idx

    def get_key(self, index):
        return self.keys[index]

    def get_property(self, key):
        idx = self.get_index(key)
        return self.property[idx], self.propertyIndex[idx]

    def check_key(self, key, index=None):
        if (key, index) in self.keys:
            return self.get_index((key, index))
        else:
            return False

    def set_availableInProps(self, key, available):
        idx = self.get_index(key)
        self.availableInProps[idx] = available

    def set_get_function(self, key, function):
        """set eval function"""
        idx = self.get_index(key)
        self.getFunctions[idx] = function

    def compute_value(self, prop, key):
        idx = self.get_index(key)
        if self.getFunctions[idx] is None:
            val = 0
        else:
            val = self.getFunctions[idx](
                prop, self.property[idx], self.propertyIndex[idx]
            )
        self.absoluteValues[idx] = val

    def get_value(self, key):
        """eval jacobian using specified function other wise return 0"""
        idx = self.get_index(key)
        return self.absoluteValues[idx]

    def display_available(self):
        available_dict = {}
        for i, key in enumerate(self.standardKeys):
            print(key, self.availableInProps[i])
            available_dict[key] = self.availableInProps[i]
        return available_dict


class reaktoroJacobianSpec:
    def __init__(self, reaktoroBase, rktOutputSpec):
        self.rktBase = reaktoroBase
        if isinstance(self.rktBase, reaktoroState) == False:
            raise TypeError("Reator jacobian require rektoroState class")
        self.rktOutputSpec = rktOutputSpec
        if isinstance(self.rktOutputSpec, reaktoroOutputSpec) == False:
            raise TypeError("Reator outputs require reaktoroOutputSpec class")
        self.jacRows = jacboianRows(self.rktBase.rktState)
        self.set_jacobian_type()
        self.configure_numerical_jacobian()
        self.check_existing_jacobian_props()

    def configure_numerical_jacobian(
        self, jacobian_type="average", order=4, step_size=1e-8
    ):
        """Configure numerical derivate options

        Keywords:
        derivative_type -- type of derivative method (average, or center_difference)
        order -- order of derivative
        step_size -- numerical step size for approximating derivative
        """
        self.derStepSize = step_size
        self.jacobianType = jacType.average
        if jacobian_type == jacType.average:
            self.jacobianType = jacType.average
            assert order % 2 == 0
            self.numericalSteps = np.arange(-order / 2, order / 2 + 1) * step_size
        if jacobian_type == jacType.center_difference:
            self.jacobianType = jacType.center_difference
            self.center_diff_order(order)
        self.set_up_chem_and_aq_states()

    def set_up_chem_and_aq_states(self):
        self.chemPropStates = []
        self.aqueousPropStates = []
        for step in self.numericalSteps:
            self.chemPropStates.append(rkt.ChemicalProps(self.rktBase.rktSystem))
            self.aqueousPropStates.append(rkt.AqueousProps(self.chemPropStates[-1]))

    def update_states(self, new_states):
        for i in range(len(self.numericalSteps)):
            self.chemPropStates[i].update(new_states[:, i])
            self.aqueousPropStates[i].update(self.chemPropStates[i])

    def get_state_values(self, outputObject):
        output_vals = []
        for i, step in enumerate(self.numericalSteps):
            if outputObject.propertyType == propTypes.chemProp:
                prop = self.chemPropStates[i]
            elif outputObject.propertyType == propTypes.aqueousProp:
                prop = self.aqueousPropStates[i]
            else:
                raise TypeError(
                    f"{outputObject.propertyType} not supported by numerical derivative method, please update"
                )

            val = self.rktOutputSpec.evaluate_property(outputObject, prop)
            output_vals.append(val)
        return output_vals

    def set_jacobian_type(self):
        """function to check all inputs and identify if exact or numeric jacobian should be used"""
        _jac_types = []
        for key, obj in self.rktOutputSpec.rktOutputs.items():
            jac_available = self.jacRows.check_key(obj.propertyName, obj.propertyIndex)
            if jac_available is not False:
                obj.jacobianType = jacType.exact
            else:
                obj.jacobianType = jacType.numeric
            _jac_types.append(f"{key}: {obj.jacobianType}")
        _log.info(f"Jacobian outputs are: {_jac_types}")

    def check_existing_jacobian_props(self):
        """check if jacobian is available in current properties and identfy function to use
        for getting values"""
        for idx, key in enumerate(self.jacRows.standardKeys):
            prop, prop_index = self.jacRows.get_property(key)
            try:
                prop_type, get_function = self.rktOutputSpec.get_prop_type(
                    prop, prop_index
                )
                if prop_type == propTypes.chemProp:
                    self.jacRows.set_availableInProps(key, True)
                    self.jacRows.set_get_function(key, get_function)
                else:
                    self.jacRows.set_availableInProps(key, False)
            except NotImplementedError:
                self.jacRows.set_availableInProps(key, False)

    def update_jacobian_absolute_values(self):
        """' used to update absolute values of jacobian - needed for numeric derivatives"""
        prop = self.rktOutputSpec.supported_properties[propTypes.chemProp]
        for jacIdx, key in enumerate(self.jacRows.standardKeys):
            self.jacRows.compute_value(prop, key)

    def process_jacobian_matrix(self, jacobian_matrix, input_index, input_value):
        """this function is used to pull out a specific column from the jacobian and also
        generate matrix for manually propagating derivatives
        Here we need to retain row order, as its same as input into chem properties"""
        jacobian_abs_matrix = []
        jacobian_dict = {}
        for jacIdx, key in enumerate(self.jacRows.standardKeys):
            jacobian_abs_matrix.append([])
            jac_abs_value = self.jacRows.get_value(key)
            jac_value = jacobian_matrix[jacIdx][input_index]
            jacobian_dict[key] = jac_value

            for step in self.numericalSteps:
                der_step = jac_value * input_value * step
                jacobian_abs_matrix[-1].append(jac_abs_value + der_step)
        return jacobian_dict, np.array(jacobian_abs_matrix)

    def get_jacobian(self, jacobian_matrix, inputObject):
        input_index = inputObject.jacobianIndex
        input_value = inputObject.tempValue
        jacobian_dict, jacobian_abs_matrix = self.process_jacobian_matrix(
            jacobian_matrix, input_index, input_value
        )
        self.update_states(jacobian_abs_matrix)

        output_jacobian = []
        for output_key, output_obj in self.rktOutputSpec.rktOutputs.items():
            if output_obj.jacobianType == jacType.exact:
                output_jacobian.append(jacobian_dict[output_key])
            else:
                values = self.get_state_values(output_obj)
                if jacType.average:
                    jac_val = np.average(
                        np.diff(values) / np.diff(self.numericalSteps * input_value)
                    )
                elif jacType.center_difference:
                    jac_val = np.array(values) * self.cdfMultipliers
                    jac_val = np.sum(jac_val) / (
                        self.rkt_aqueous_props_der_step * input_value
                    )
                output_jacobian.append(jac_val)
        return output_jacobian

    def center_diff_order(self, order):
        """paramterizes center different taylor series
        refer to https://en.wikipedia.org/wiki/Finite_difference_coefficient"""

        if order == 2:
            self.numericalSteps = (
                np.array(
                    [
                        -1,
                        1,
                    ]
                )
                * self.derStepSize
            )
            self.cdfMultipliers = np.array([-1 / 2, 1 / 2])
        if order == 4:
            self.numericalSteps = np.array([-2, -1, 1, 2]) * self.derStepSize
            self.cdfMultipliers = np.array([-1 / 12, 2 / 3, -2 / 3, 1 / 12])
        if order == 6:
            self.numericalSteps = np.array([-3, -2, -1, 1, 2, 3]) * self.derStepSize
            self.cdfMultipliers = np.array(
                [-1 / 60, 3 / 20, -3 / 4, 3 / 4, -3 / 20, 1 / 60]
            )
        if order == 8:
            self.numericalSteps = (
                np.array([-4, -3, -2, -1, 1, 2, 3, 4]) * self.derStepSize
            )
            self.cdfMultipliers = np.array(
                [1 / 280, -4 / 105, 1 / 5, -4 / 5, 4 / 5, -1 / 5, 4 / 105, -1 / 280]
            )
        if order == 10:
            self.numericalSteps = (
                np.array([-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]) * self.derStepSize
            )
            self.cdfMultipliers = (
                np.array([-2, 25, -150, 600, -2100, 2100, -600, 150, -25, 2]) / 2520
            )
