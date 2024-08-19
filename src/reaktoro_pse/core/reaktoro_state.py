import reaktoro as rkt
from pyomo.environ import Var, units as pyunits

import reaktoro_pse.core.util_classes.rktInputs as rktInputs
from pyomo.core.base.var import IndexedVar, Var

__author__ = "Alexander Dudchenko"


""" base class for configuring reaktoro states and solver"""


class reaktoroState:
    def __init__(self):
        """initialize for all parameters need to build reaktor solver"""
        self.rktInputs = rktInputs.rktInputs()
        self.phases = []
        self.rktGasPhases = []

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
        self.rktInputs.enable_rkt_species_conversion(
            convert_to_rkt_species, species_to_rkt_species_dict
        )
        self.rktInputs.set_composition_is_elements(composition_is_elements)
        for props, pyo_obj in composition.items():
            if composition_index is None or composition_index in props:
                if isinstance(props, str):
                    specie = props
                else:
                    specie = props[-1]
                self.rktInputs[specie] = pyo_obj

        self.rktInputs["temperature"] = self.process_input(
            temperature, temperature_index, default=293.15, units=pyunits.K
        )
        self.rktInputs["pressure"] = self.process_input(
            pressure, pressure_index, default=10000, units=pyunits.Pa
        )

        if pH is not None:
            self.rktInputs["pH"] = self.process_input(
                pH, pH_index, default=7, units=pyunits.dimensionless
            )
        """ generates rkt speciation"""
        self.rktAqueousPhase = rkt.AqueousPhase(
            rkt.speciate(self.rktInputs.speciesList)
        )

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

    def register_mineral_phases(self, mineral_phases=[]):
        """register possible mineral phases"""
        if isinstance(mineral_phases, str):
            mineral_phases = [mineral_phases]

        for mineral_phase in mineral_phases:
            self.phases.append(rkt.MineralPhase(mineral_phase))

    def register_gas_phases(self, gas_phases=[]):
        """register possible gas phases"""
        if isinstance(gas_phases, str):
            gas_phases = [gas_phases]
        self.rktGasPhases.append(rkt.GaseousPhase(rkt.speciate(gas_phases)))

    def set_database(self, dbtype="PhreeqcDatabase", database="pitzer.dat"):
        """set data base of reaktoro"""
        self.rktDatabase = getattr(rkt, dbtype)(database)

    def set_activity_model(self, activity_model="ActivityModelPitzer"):
        """set activity model of reaktoro"""
        self.rktAqueousPhase.set(getattr(rkt, activity_model)())

    def build_state(self):
        """this will build reaktor states"""
        self.rktSystem = rkt.ChemicalSystem(
            self.rktDatabase, self.rktAqueousPhase, *self.phases, *self.rktGasPhases
        )
        self.rktState = rkt.ChemicalState(self.rktSystem)
        self.set_rkt_state()

    def set_rkt_state(self):
        """sets initial rkt state using user provided inputs"""
        self.rktState.temperature(
            self.rktInputs["temperature"].value, self.rktInputs["temperature"].mainUnit
        )
        self.rktState.pressure(
            self.rktInputs["pressure"].value, self.rktInputs["pressure"].mainUnit
        )

        """ set apparant species if used """
        if self.rktInputs.convertToRktSpecies:
            for species in self.rktInputs.speciesList:
                if species in self.rktInputs:  # user might not provide all
                    if self.rktInputs[species].value != 0:
                        self.rktState.set(
                            species,
                            self.rktInputs[species].value,
                            self.rktInputs[species].mainUnit,
                        )

    def equilibrate_state(self):
        self.set_rkt_state()
        rkt.equilibrate(self.rktState)
