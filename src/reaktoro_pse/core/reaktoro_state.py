import reaktoro as rkt
from pyomo.environ import Var, units as pyunits

import reaktoro_pse.core.util_classes.rkt_inputs as RktInputs
from pyomo.core.base.var import IndexedVar, Var
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)

__author__ = "Alexander Dudchenko"


""" base class for configuring reaktoro states and solver"""


class ReaktoroState:
    def __init__(self):
        """initialize for all parameters need to build reaktor solver"""
        self.inputs = RktInputs.RktInputs()
        self.mineral_phases = []
        self.gas_phases = []
        self.ion_exchange_phases = []

    def register_inputs(
        self,
        composition,
        temperature=None,
        pressure=None,
        pH=None,
        convert_to_rkt_species=True,
        species_to_rkt_species_dict="default",
        composition_is_elements=False,
        composition_index=None,
        temperature_index=None,
        pressure_index=None,
        pH_index=None,
    ):
        """registers inputs

        Keyword arguments:
        composition -- dictionary or pyomo indexed block that contains apparent or elemental specie composition
        temperature -- pyomo var that contains temperature data (293.15K)
        pressure -- pyomo var that contains pressure data (default 1 atm)
        pH -- pyomo var that contains solution pH (default None)
        convert_to_rkt_species -- if set to True, specie input names are converted to RKT species names (set this to False if providing exact elemental input composition)
        species_to_rkt_species_dict -- dictionary that defines how to convert species to rkt species, "default" converter supports limited number of databases and files
        composition_index -- defines index for supplied input to use in configuring rkt inputs (e.g. input[(composition_index,specie)])
        temperature_index -- defines index for supplied input to use in configuring rkt inputs (e.g. input[temperature_index])
        pressure_index -- defines index for supplied input to use in configuring rkt inputs (e.g. input[pressure_index])
        pH_index -- defines index for supplied input to use in configuring rkt inputs (e.g. input[pH_index])
        """
        # unfold input for composition
        self.inputs.enable_rkt_species_conversion(
            convert_to_rkt_species, species_to_rkt_species_dict
        )
        self.inputs.set_composition_is_elements(composition_is_elements)
        for props, pyo_obj in composition.items():
            if composition_index is None or composition_index in props:
                if isinstance(props, str):
                    specie = props
                else:
                    specie = props[-1]
                self.inputs[specie] = pyo_obj

        self.inputs["temperature"] = self.process_input(
            temperature, temperature_index, default=293.15, units=pyunits.K
        )
        self.inputs["pressure"] = self.process_input(
            pressure, pressure_index, default=10000, units=pyunits.Pa
        )

        if pH is not None:
            self.inputs["pH"] = self.process_input(
                pH, pH_index, default=7, units=pyunits.dimensionless
            )
        """ generates rkt speciation"""
        self.aqueous_phase = rkt.AqueousPhase(rkt.speciate(self.inputs.species_list))

    def process_input(self, var, var_index, default=0, units=None):
        if var is None:
            var = Var(initialize=293.15, units=units)
            var.construct()
            var.fix()
            return var
        elif isinstance(var, (dict, IndexedVar)):
            for index, v in var.items():
                if var_index == index:
                    return v
        else:
            return var

    def _process_phase(self, phase, default_phase):
        if isinstance(phase, str):
            return default_phase(phase)
        else:
            return phase

    def register_mineral_phases(self, mineral_phases=[]):
        """register possible mineral phases"""
        if isinstance(mineral_phases, list) == False:
            mineral_phases = [mineral_phases]

        for mineral_phase in mineral_phases:
            self.mineral_phases.append(
                self._process_phase(mineral_phase, rkt.MineralPhase)
            )

    def register_gas_phases(self, gas_phases=[]):
        """register possible gas phases"""
        if isinstance(gas_phases, str):
            gas_phases = [gas_phases]
        for gas_phase in gas_phases:
            self.gas_phases.append(self._process_phase(gas_phase, rkt.GaseousPhase))

    def register_ion_exchange_phase(self, ion_phase=[]):
        """register possible ion exchange phases"""
        if isinstance(ion_phase, str):
            ion_phase = [ion_phase]
        self.ion_exchange_phases.append(
            self._process_phase(ion_phase, rkt.IonExchangePhase)
        )

    def set_database(self, dbtype="PhreeqcDatabase", database="pitzer.dat"):
        """set data base of reaktoro"""
        self.database = getattr(rkt, dbtype)(database)

    def _process_activity(self, activity_model, state_of_matter=None):
        def get_func(activity_model, state_of_matter):
            if isinstance(activity_model, str):
                if state_of_matter is None:
                    return getattr(rkt, activity_model)()
                else:
                    return getattr(rkt, activity_model)(state_of_matter)
            else:
                return activity_model

        if isinstance(activity_model, list):
            if state_of_matter is None:
                activity_model = rkt.chain(
                    *[get_func(am, state_of_matter) for am in activity_model]
                )
            else:
                activity_model = rkt.chain(
                    *[get_func(am, state_of_matter) for am in activity_model]
                )
        else:
            activity_model = get_func(activity_model, state_of_matter)
        return activity_model

    def set_aqueous_phase_activity_model(
        self, activity_model="ActivityModelIdealAqueous"
    ):
        """set activity model of aquous phases in reaktoro"""

        activity_model = self._process_activity(activity_model)
        self.aqueous_phase.set(activity_model)

    def set_gas_phase_activity_model(self, activity_model="ActivityModelIdealGas"):
        """set activity model of gas phases in reaktoro"""
        # TODO: add support for massking rkt initialized activty models directly
        activity_model = self._process_activity(activity_model)
        for phase in self.gas_phases:
            phase.set(activity_model)

    def set_mineral_phase_activity_model(
        self, activity_model="ActivityModelIdealSolution"
    ):
        """set activity model of mineral phases in reaktoro"""
        # TODO: add support for massking rkt initialized activty models directly
        activity_model = self._process_activity(
            activity_model, state_of_matter=rkt.StateOfMatter.Solid
        )
        for phase in self.mineral_phases:
            phase.set(activity_model)

    def set_ion_exchange_phase_activity_model(
        self, activity_model="ActivityModelIonExchange"
    ):
        """set activity model of mineral phases in reaktoro"""
        activity_model = self._process_activity(
            activity_model,
        )
        for phase in self.ion_exchange_phases:
            phase.set(activity_model)

    def build_state(self):
        """this will build reaktor states"""
        self.system = rkt.ChemicalSystem(
            self.database, self.aqueous_phase, *self.mineral_phases, *self.gas_phases
        )
        self.state = rkt.ChemicalState(self.system)
        self.set_rkt_state()

    def set_rkt_state(self):
        """sets initial rkt state using user provided inputs"""
        self.state.temperature(
            self.inputs["temperature"].get_value(),
            self.inputs["temperature"].main_unit,
        )
        self.state.pressure(
            self.inputs["pressure"].get_value(),
            self.inputs["pressure"].main_unit,
        )

        """ set apparant species if used """
        if self.inputs.composition_is_elements == False:
            for species in self.inputs.species_list:
                if species in self.inputs:  # user might not provide all
                    if self.inputs[species].get_value() != 0:
                        unit = self.inputs[species].main_unit
                        # print(species, self.inputs[species].get_value())
                        if unit == "dimensionless":
                            self.state.set(
                                species,
                                self.inputs[species].get_value(),
                                "mol",
                            )
                        else:
                            self.state.set(
                                species,
                                self.inputs[species].get_value(),
                                self.inputs[species].main_unit,
                            )

    def equilibrate_state(self):
        self.set_rkt_state()
        rkt.equilibrate(self.state)
        _log.info("Equilibrated successfully")
