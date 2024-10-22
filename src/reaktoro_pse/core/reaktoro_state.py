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
import sys
from numpy import isin
import reaktoro as rkt
from pyomo.environ import units as pyunits

import reaktoro_pse.core.util_classes.rkt_inputs as RktInputs
from reaktoro_pse.core.util_classes.rkt_inputs import RktInputTypes
from pyomo.core.base.var import IndexedVar
import idaes.logger as idaeslog
import copy

_log = idaeslog.getLogger(__name__)

__author__ = "Alexander V. Dudchenko"


class PhaseData:
    def __init__(self):
        """creates state for storing phase object information"""
        self.phase_list = None
        self.non_speciate_phase_list = None
        self.activity_model = None
        self.state_of_matter = None
        self.phase_function = None
        self.phase_list_mode = False
        self.speciate = True

    def update_species_list(
        self, species, phase_function, list_mode=False, speciate=True
    ):
        """updates list of speciesies for given phase, ensuring no duplication

        args:
            species - specie or list of specieis to register
            phase_function - reaktoro function to create phase
            list_mode - defines if we will need to create a list of phases or single phase
        """

        def _update_list(current_list, input_val):
            if isinstance(input_val, str):
                current_list = [input_val]

            elif isinstance(input_val, list):
                for i in input_val:
                    if i not in current_list:
                        current_list.append(i)
            else:
                current_list = input_val
            return current_list

        if speciate:
            if species != None and self.phase_list == None:
                self.phase_list = []
            self.phase_list = _update_list(self.phase_list, species)
        else:
            if species != None and self.non_speciate_phase_list == None:
                self.non_speciate_phase_list = []
            self.non_speciate_phase_list = _update_list(
                self.non_speciate_phase_list, species
            )
        self.phase_function = phase_function
        self.phase_list_mode = list_mode

    def set_activity_model(self, activity_model, state_of_matter=None):
        """sets activity model and it state"""
        self.activity_model = activity_model
        self.state_of_matter = state_of_matter


class PhaseManager:
    def __init__(self):
        self.registered_phases = {}
        self.exclude_species_list = []

    def register_phases_species(
        self, phase, species, phase_function, list_mode=False, speciate=True
    ):
        """this is used to add up all phases provided as inputs or registed phases by user,
        for example, user might provide CO2(g) as a gas input wants to track also N2(g) via register_gas_phase,
        this ensure we updates all the phases with all requested values

        args:
            phase - phase type supported by reaktoro-pse
            species - specie or list of specieis to register
            phase_function - reaktoro function to create phase
            list_mode - defines if we will need to create a list of phases or single phase
        """
        if phase not in self.registered_phases:
            self.registered_phases[phase] = PhaseData()

        self.registered_phases[phase].update_species_list(
            species, phase_function, list_mode, speciate
        )

    def set_activity_model(
        self, phase, activity_model, default_activity_model, state_of_matter=None
    ):
        """
        set activity model for reaktoro

        args:
            phase - phase type supported by reaktoro-pse
            activity_model - activity model (should be string or touple that contains options
            to be passed into function)
            default_activity_model - default activity model if one is not provided
            state_of_matter - defines state of matter for solids, otherwise none
        """
        if phase not in self.registered_phases:
            self.registered_phases[phase] = PhaseData()
        if activity_model is None:
            activity_model = default_activity_model
        self.registered_phases[phase].set_activity_model(
            activity_model, state_of_matter
        )

    def get_registered_phases(self, activate_database):
        """
        creates listof phases with applied activity models for generation of reaktoro state
        args:
            activate_database - database to use during creation of phases

        """

        activate_phase_list = []
        for phase, phase_object in self.registered_phases.items():
            if (
                phase_object.phase_list is not None
                or phase_object.non_speciate_phase_list is not None
            ):
                rkt_phase_object = self.create_rkt_phase(
                    activate_database,
                    phase_object.phase_function,
                    phase_object.phase_list,
                    phase_object.non_speciate_phase_list,
                    phase_object.phase_list_mode,
                )
                if isinstance(rkt_phase_object, list):
                    for rpo in rkt_phase_object:
                        self.apply_activity_model(
                            rpo,
                            phase_object.activity_model,
                            phase_object.state_of_matter,
                        )
                        activate_phase_list.append(rpo)
                else:
                    self.apply_activity_model(
                        rkt_phase_object,
                        phase_object.activity_model,
                        phase_object.state_of_matter,
                    )
                    activate_phase_list.append(rkt_phase_object)
        return activate_phase_list

    def apply_activity_model(
        self, rkt_phase_object, activity_model, state_of_matter=None
    ):
        """sets activity mode"""
        if activity_model is None:
            raise TypeError(f"Activity model for {rkt_phase_object} is not set.")
        rkt_activity_model_object = self._process_activity(
            activity_model, activity_model, state_of_matter
        )
        rkt_phase_object.set(rkt_activity_model_object)

    def create_rkt_phase(
        self,
        active_database,
        phase_function,
        phases_list,
        non_speciate_phase_list,
        phase_list_mode,
    ):
        """Function to remove species from speciation command"""
        if phase_list_mode:
            rkt_phase_list = []
            for phase in phases_list:
                if isinstance(phase, str):
                    rkt_phase_list.append(phase_function(phase))
                else:
                    rkt_phase_list.append(phase)
            if non_speciate_phase_list is not None:
                for phase in non_speciate_phase_list:
                    if isinstance(phase, str):
                        rkt_phase_list.append(phase_function(phase))
                    else:
                        rkt_phase_list.append(phase)
            return rkt_phase_list
        elif isinstance(phases_list, list):
            if isinstance(phases_list, (list, tuple)):
                phases_list = " ".join(phases_list)
            try:
                # try to spetiatiate
                phases = phase_function(rkt.speciate(phases_list))
            except TypeError:
                # if does nto spectiate, try passing string input directly
                phases = phase_function(phases_list)
            system = rkt.ChemicalSystem(active_database, phases)
            spc_list = []
            for specie in system.species():
                if specie.name() not in str(self.exclude_species_list):
                    spc_list.append(specie.name())
            if non_speciate_phase_list is not None:
                for phase in non_speciate_phase_list:
                    if phase not in spc_list:
                        spc_list.append(phase)
            return phase_function(" ".join(spc_list))
        elif isinstance(non_speciate_phase_list, list):
            return phase_function(" ".join(non_speciate_phase_list))
        else:
            return phases_list

    def _process_activity(
        self,
        active_database,
        activity_model,
        state_of_matter=None,
    ):
        """this will process give activity model if user provides a sting it will
        find it on reaktoro and initialize it either by it self or passing in state of matter argument
        if user provides initialized reaktoro activity model we use it directly."""

        def get_func(activity_model, state_of_matter):
            if isinstance(activity_model, str):
                if state_of_matter is None:
                    try:
                        return getattr(rkt, activity_model)()
                    except RuntimeError:
                        # might require database as input (e.g ActivityModelPhreeqc)
                        return getattr(rkt, activity_model)(active_database)

                else:
                    return getattr(rkt, activity_model)(state_of_matter)
            elif isinstance(activity_model, tuple):
                if isinstance(activity_model[1], (tuple, list)):
                    return getattr(rkt, activity_model[0])(*activity_model[1])
                else:
                    return getattr(rkt, activity_model[0])(activity_model[1])
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


class ReaktoroStateExport:
    def __init__(self):
        # defines required exports
        self.exclude_species_list = None
        # we can direcly pass original phase manager
        self.phase_manager = None
        # this will not have any of the pyomo vars
        self.inputs = RktInputs.RktInputs()
        self.exclude_species_list = None
        self.database_type = None
        self.database_file = None

    def copy_rkt_inputs(self, inputs):
        self.inputs = copy.deepcopy(inputs)
        for key, obj in self.inputs.items():
            # self.inputs[key] = copy.deepcopy(obj)
            self.inputs[key].delete_pyomo_var()


# base class for configuring reaktoro states and solver
class ReaktoroState:
    def __init__(self):
        """initialize all parameters need to build reaktor solver"""
        self.inputs = RktInputs.RktInputs()
        self.phase_manager = PhaseManager()
        self.exclude_species_list = []
        self.export_state = ReaktoroStateExport()
        self._inputs_not_processed = True

    def register_system_inputs(
        self,
        temperature=None,
        pressure=None,
        enthalpy=None,
        pH=None,
        temperature_index=None,
        pressure_index=None,
        pH_index=None,
        enthalpy_index=None,
    ):
        """registers inputs

        Keyword arguments:

        temperature -- pyomo var that contains temperature data (293.15K)
        pressure -- pyomo var that contains pressure data (default 1 atm)
        temperature_index -- defines index for supplied input to use in configuring rkt inputs (e.g. input[temperature_index])
        pressure_index -- defines index for supplied input to use in configuring rkt inputs (e.g. input[pressure_index])
        """
        if temperature is not None:
            self.inputs[RktInputTypes.temperature] = self.process_input(
                temperature, temperature_index
            )
            self.inputs[RktInputTypes.temperature].set_input_type(
                RktInputTypes.system_state
            )
        if pressure is not None:
            self.inputs[RktInputTypes.pressure] = self.process_input(
                pressure, pressure_index
            )
            self.inputs[RktInputTypes.pressure].set_input_type(
                RktInputTypes.system_state
            )
        if enthalpy is not None:
            self.inputs[RktInputTypes.enthalpy] = self.process_input(
                enthalpy, enthalpy_index
            )
            self.inputs[RktInputTypes.enthalpy].set_input_type(
                RktInputTypes.system_state
            )
        if pH is not None:
            self.inputs[RktInputTypes.pH] = self.process_input(pH, pH_index)
            self.inputs[RktInputTypes.pH].set_input_type(RktInputTypes.aqueous_phase)

    def set_input_options(
        self,
        phase_type,
        convert_to_rkt_species=True,
        species_to_rkt_species_dict="default",
        composition_is_elements=False,
    ):
        """registers inputs

        Keyword arguments:
        convert_to_rkt_species -- if set to True, specie input names are converted to RKT species names (set this to False if providing exact elemental input composition)
        species_to_rkt_species_dict -- dictionary that defines how to convert species to rkt species, "default" converter supports limited number of databases and files
        """

        self.inputs.enable_rkt_species_conversion(
            phase_type, convert_to_rkt_species, species_to_rkt_species_dict
        )
        self.inputs.set_composition_is_elements(phase_type, composition_is_elements)

    def register_aqueous_inputs(
        self,
        composition,
        composition_index=None,
    ):
        """registers inputs

        Keyword arguments:
        composition -- dictionary or pyomo indexed block that contains apparent or elemental specie composition
        pH -- pyomo var that contains solution pH (default None)
        composition_index -- defines index for supplied input to use in configuring rkt inputs (e.g. input[(composition_index,specie)])
        temperature_index -- defines index for supplied input to use in configuring rkt inputs (e.g. input[temperature_index])
        pressure_index -- defines index for supplied input to use in configuring rkt inputs (e.g. input[pressure_index])
        pH_index -- defines index for supplied input to use in configuring rkt inputs (e.g. input[pH_index])
        """
        self.register_inputs(
            composition, composition_index, RktInputTypes.aqueous_phase
        )

    def register_inputs(self, composition, composition_index, phase):
        """generic input registration method,
        unpacks composition (assumes its a dict or indexed var)
        if user provides index then extract values only for that index and
        phase to specify which phase input belongs to"""
        for props, pyo_obj in composition.items():
            if composition_index is None or composition_index in props:
                if isinstance(props, str):
                    specie = props
                else:
                    specie = props[-1]
                self.inputs[specie] = pyo_obj
                self.inputs[specie].set_input_type(phase)
        self._inputs_not_processed = True  # flag that inputs ver modified

    def register_species_to_exclude(self, species):
        """updates list of species to exclude from when creating phases"""
        if species != None:
            if isinstance(species, str):
                self.phase_manager.exclude_species_list.append(species)
            elif isinstance(species, list):
                for spc in species:
                    self.phase_manager.exclude_species_list.append(spc)
            else:
                raise TypeError(f"{species} is not supported, must be str or list")

    def register_gas_inputs(
        self,
        composition,
        composition_index=None,
    ):
        """registers inputs

        Keyword arguments:
        composition -- dictionary or pyomo indexed block that contains apparent or elemental specie composition
        composition_index -- defines index for supplied input to use in configuring rkt inputs (e.g. input[(composition_index,specie)])
        """
        # unfold input for composition
        self.register_inputs(composition, composition_index, RktInputTypes.gas_phase)

    def register_mineral_inputs(
        self,
        composition,
        composition_index=None,
    ):
        """registers inputs

        Keyword arguments:
        composition -- dictionary or pyomo indexed block that contains apparent or elemental specie composition
        composition_index -- defines index for supplied input to use in configuring rkt inputs (e.g. input[(composition_index,specie)])
        """
        # unfold input for composition
        self.register_inputs(
            composition, composition_index, RktInputTypes.mineral_phase
        )

    def register_liquid_inputs(
        self,
        composition,
        composition_index=None,
    ):
        """registers inputs

        Keyword arguments:
        composition -- dictionary or pyomo indexed block that contains apparent or elemental specie composition
        composition_index -- defines index for supplied input to use in configuring rkt inputs (e.g. input[(composition_index,specie)])
        """
        # unfold input for composition
        self.register_inputs(composition, composition_index, RktInputTypes.liquid_phase)

    def register_condensed_inputs(
        self,
        composition,
        composition_index=None,
    ):
        """registers inputs

        Keyword arguments:
        composition -- dictionary or pyomo indexed block that contains apparent or elemental specie composition
        composition_index -- defines index for supplied input to use in configuring rkt inputs (e.g. input[(composition_index,specie)])
        """
        # unfold input for composition
        self.register_inputs(
            composition, composition_index, RktInputTypes.condensed_phase
        )

    def register_solid_inputs(
        self,
        composition,
        composition_index=None,
    ):
        """registers inputs

        Keyword arguments:
        composition -- dictionary or pyomo indexed block that contains apparent or elemental specie composition
        composition_index -- defines index for supplied input to use in configuring rkt inputs (e.g. input[(composition_index,specie)])
        """
        # unfold input for composition
        self.register_inputs(composition, composition_index, RktInputTypes.solid_phase)

    def register_ion_exchange_inputs(
        self,
        composition,
        composition_index=None,
    ):
        """registers inputs

        Keyword arguments:
        composition -- dictionary or pyomo indexed block that contains apparent or elemental specie composition
        composition_index -- defines index for supplied input to use in configuring rkt inputs (e.g. input[(composition_index,specie)])
        """
        # unfold input for composition
        self.register_inputs(
            composition, composition_index, RktInputTypes.ion_exchange_phase
        )

    def verify_specie_units(self):
        for phase_type in self.inputs.registered_phases:
            for specie in self.inputs.species_list[phase_type]:
                mw, unit = self.get_molar_mass(specie)
                self.verify_unit(self.inputs[specie], mw, unit)

    def verify_unit(self, rkt_input_object, mw, mw_unit):
        """verify user provides mass or mol basis for species
        other wise do nothing"""
        if (
            rkt_input_object.main_unit == RktInputTypes.mol
            or rkt_input_object.main_unit == RktInputTypes.dimensionless
        ):
            rkt_input_object.set_unit_conversion(None, None)
        elif rkt_input_object.main_unit in RktInputTypes.mass_units:
            rkt_input_object.set_unit_conversion(mw, mw_unit)
            rkt_input_object.set_required_unit(RktInputTypes.mol)
        else:
            raise ValueError(
                f"Unit {rkt_input_object.main_unit} is not supported, only mol, kg, and mg is supported as main input unit (e.g. main_unit/time_unit)"
            )

    def get_molar_mass(self, comp):
        """get mw from database from specie or element"""
        if comp in self.database_species:
            return self.get_molar_mass_specie(comp)
        elif comp in self.database_elements:
            return self.get_molar_mass_element(comp)
        else:
            raise KeyError(f"Provided {comp} not found in species or elements in db")

    def get_molar_mass_specie(self, species_name):
        mw = self.database.species(species_name).molarMass()
        mw_unit = pyunits.kg / pyunits.mol
        return mw, mw_unit

    def get_molar_mass_element(self, element_name):
        mw = self.database.element(element_name).molarMass()
        mw_unit = pyunits.kg / pyunits.mol
        return mw, mw_unit

    def process_input(self, var, var_index):
        if isinstance(var, (dict, IndexedVar)):
            for index, v in var.items():
                if var_index == index:
                    return v
        else:
            return var

    def register_aqueous_phase(self, aqueous_phases=None):
        """register possible mineral phases"""
        if aqueous_phases != [] and aqueous_phases is not None:
            self.phase_manager.register_phases_species(
                RktInputTypes.aqueous_phase, aqueous_phases, rkt.AqueousPhase
            )

    def register_condensed_phase(self, condensed_phase=None):
        """register possible condensed phases"""
        if condensed_phase != [] and condensed_phase is not None:
            self.phase_manager.register_phases_species(
                RktInputTypes.condensed_phase, condensed_phase, rkt.CondensedPhases
            )

    def register_liquid_phase(self, liquid_phases=None):
        """register possible mineral phases"""
        if liquid_phases != [] and liquid_phases is not None:
            self.phase_manager.register_phases_species(
                RktInputTypes.liquid_phase, liquid_phases, rkt.LiquidPhase
            )

    def register_mineral_phases(self, mineral_phases=None):
        """register possible mineral phases"""
        self.phase_manager.register_phases_species(
            RktInputTypes.mineral_phase,
            mineral_phases,
            rkt.MineralPhase,
            list_mode=True,
        )

    def register_solid_phases(self, solid_phases=None):
        """register possible solid phases"""
        self.phase_manager.register_phases_species(
            RktInputTypes.solid_phase, solid_phases, rkt.SolidPhase, list_mode=True
        )

    def register_gas_phase(self, gas_phase=None, speciate=False):
        """register possible gas phases"""
        self.phase_manager.register_phases_species(
            RktInputTypes.gas_phase, gas_phase, rkt.GaseousPhase, speciate=speciate
        )

    def register_ion_exchange_phase(self, ion_exchange_phase=None, speciate=False):
        """register possible ion exchange phases"""
        self.phase_manager.register_phases_species(
            RktInputTypes.ion_exchange_phase,
            ion_exchange_phase,
            rkt.IonExchangePhase,
            speciate=speciate,
        )

    def set_database(self, dbtype="PhreeqcDatabase", database="pitzer.dat"):
        """set data base of reaktoro"""
        """ assume that if database is string we need to find and init in reaktoro
        other wise assume we received activated reaktoro db"""
        self.database_type = dbtype
        self.database_file = database

    def load_database(self):
        if isinstance(self.database_type, str):
            self.database = getattr(rkt, self.database_type)(self.database_file)
        else:
            self.database = self.database_type
        self.database_species = [specie.name() for specie in self.database.species()]
        self.database_elements = [
            element.symbol() for element in self.database.elements()
        ]

    def set_aqueous_phase_activity_model(
        self,
        activity_model=None,
    ):
        """set activity model of aqueous phases in reaktoro"""
        self.phase_manager.set_activity_model(
            RktInputTypes.aqueous_phase,
            activity_model,
            default_activity_model="ActivityModelIdealAqueous",
        )

    def set_liquid_phase_activity_model(
        self,
        activity_model=None,
    ):
        """set activity model of liquid phases in reaktoro"""

        self.phase_manager.set_activity_model(
            RktInputTypes.liquid_phase,
            activity_model,
            default_activity_model="ActivityModelIdealSolution",
        )

    def set_gas_phase_activity_model(self, activity_model=None):
        """set activity model of gas phases in reaktoro"""
        self.phase_manager.set_activity_model(
            RktInputTypes.gas_phase,
            activity_model,
            default_activity_model="ActivityModelIdealGas",
        )

    def set_condensed_phase_activity_model(self, activity_model=None):
        """set activity model of ccondensed phases in reaktoro"""
        self.phase_manager.set_activity_model(
            RktInputTypes.condensed_phase,
            activity_model,
            default_activity_model="ActivityModelIdealSolution",
            state_of_matter=rkt.StateOfMatter.Liquid,
        )

    def set_mineral_phase_activity_model(self, activity_model=None):
        """set activity model of mineral phases in reaktoro"""
        self.phase_manager.set_activity_model(
            RktInputTypes.mineral_phase,
            activity_model,
            state_of_matter=rkt.StateOfMatter.Solid,
            default_activity_model="ActivityModelIdealSolution",
        )

    def set_solid_phase_activity_model(self, activity_model=None):
        """set activity model of solid phases in reaktoro"""
        self.phase_manager.set_activity_model(
            RktInputTypes.solid_phase,
            activity_model,
            state_of_matter=rkt.StateOfMatter.Solid,
            default_activity_model="ActivityModelIdealSolution",
        )

    def set_ion_exchange_phase_activity_model(self, activity_model=None):
        """set activity model of mineral phases in reaktoro"""
        self.phase_manager.set_activity_model(
            RktInputTypes.ion_exchange_phase,
            activity_model,
            default_activity_model="ActivityModelIonExchange",
        )

    def process_registered_inputs(self):
        """this will process all inputs, this function must be called before setting activity models!"""
        if self._inputs_not_processed:
            """post-process inputs-this converts input specs to rkt and registers phases"""
            self.inputs.process_inputs()
            """ verify units"""
            self.verify_specie_units()
            """ register any added phases """
            self.register_all_input_phases()
            self._inputs_not_processed = False

    def register_all_input_phases(self):
        '''this will add phases to reaktoro state that are provided as inputs,
        this removes the need for user to manually specifies phase inputs if they are already
        part of explicit inputs.
        E.g. if user provided CO2 in gas phase, this will add CO2 as gas phase to reaktoro
        with out the need to call "register_gas_phase"'''

        # register liquid type phases
        if len(self.inputs.species_list[RktInputTypes.aqueous_phase]) > 0:
            self.register_aqueous_phase(
                self.inputs.species_list[RktInputTypes.aqueous_phase]
            )
        if len(self.inputs.species_list[RktInputTypes.liquid_phase]) > 0:
            self.register_liquid_phase(
                self.inputs.species_list[RktInputTypes.liquid_phase]
            )
        if len(self.inputs.species_list[RktInputTypes.condensed_phase]) > 0:
            self.register_condensed_phase(
                self.inputs.species_list[RktInputTypes.condensed_phase]
            )
        # register gas phases
        if len(self.inputs.species_list[RktInputTypes.mineral_phase]) > 0:
            self.register_mineral_phases(
                self.inputs.species_list[RktInputTypes.mineral_phase]
            )
        # register solid phase
        if len(self.inputs.species_list[RktInputTypes.solid_phase]) > 0:
            self.register_solid_phases(
                self.inputs.species_list[RktInputTypes.solid_phase]
            )
        if len(self.inputs.species_list[RktInputTypes.gas_phase]) > 0:
            self.register_gas_phase(self.inputs.species_list[RktInputTypes.gas_phase])
        # register ix phases
        if len(self.inputs.species_list[RktInputTypes.ion_exchange_phase]) > 0:
            self.register_ion_exchange_phase(
                self.inputs.species_list[RktInputTypes.ion_exchange_phase]
            )

    def build_state(self):
        # this will build reaktor states
        self.load_database()
        self.process_registered_inputs()
        phases = self.phase_manager.get_registered_phases(self.database)
        self.system = rkt.ChemicalSystem(
            self.database,
            *phases,
        )
        self.state = rkt.ChemicalState(self.system)
        self.set_rkt_state()

    def set_rkt_state(self):
        # sets initial rkt state using user provided inputs
        if self.inputs.get("temperature") is not None:
            self.state.temperature(
                self.inputs["temperature"].get_value(),
                self.inputs["temperature"].main_unit,
            )
        if self.inputs.get("pressure") is not None:
            self.state.pressure(
                self.inputs["pressure"].get_value(),
                self.inputs["pressure"].main_unit,
            )
        # set apparent species if used
        for phase in self.inputs.registered_phases:
            if self.inputs.composition_is_elements[phase] == False:
                for species in self.inputs.species_list[phase]:
                    if species in self.inputs:  # user might not provide all
                        if self.inputs[species].get_value() != 0:
                            unit = self.inputs[species].main_unit
                            if unit == "dimensionless":
                                # assume correct units are provided
                                self.state.set(
                                    species,
                                    self.inputs[species].get_value(
                                        apply_conversion=True
                                    ),
                                    "mol",
                                )
                            else:
                                self.state.set(
                                    species,
                                    self.inputs[species].get_value(
                                        apply_conversion=False
                                    ),
                                    self.inputs[species].main_unit,
                                )
        # _jac_phases = [phase.name() for phase in self.state.system().phases()]

    def equilibrate_state(self):
        self.set_rkt_state()
        rkt.equilibrate(self.state)
        _log.info("Equilibrated successfully")

    def export_config(self):
        export_object = ReaktoroStateExport()
        export_object.copy_rkt_inputs(self.inputs)
        export_object.exclude_species_list = self.exclude_species_list
        export_object.phase_manager = self.phase_manager
        export_object.database_file = self.database_file
        export_object.database_type = self.database_type
        return export_object

    def load_from_export_object(self, export_object):
        self.inputs = export_object.inputs
        self.phase_manager = export_object.phase_manager
        self.database_file = export_object.database_file
        self.database_type = export_object.database_type
