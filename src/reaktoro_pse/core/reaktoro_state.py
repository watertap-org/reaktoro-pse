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
from cmath import phase
from numpy import isin
import reaktoro as rkt
from pyomo.environ import units as pyunits
from sympy import comp

import reaktoro_pse.core.util_classes.rkt_inputs as RktInputs
from reaktoro_pse.core.util_classes.rkt_inputs import RktInputTypes
from pyomo.core.base.var import IndexedVar, Var
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)

__author__ = "Alexander Dudchenko"


""" base class for configuring reaktoro states and solver"""


class ReaktoroState:
    def __init__(self):
        """initialize all parameters need to build reaktor solver"""
        self.inputs = RktInputs.RktInputs()

        self.mineral_phase = []
        self.solid_phase = []
        self.liquid_phase = None
        self.condensed_phase = None
        self.aqueous_phase = None
        self.gas_phase = None
        self.ion_exchange_phase = None
        self.registered_phases = {}

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

    def _process_phase(self, phase, default_phase):
        """
        if phsae is list of speices convert it to single string
        other wise pass it in directly into phase object
        if phase object is not list, tuple or str assume its user initialized phase
        and use directly"""
        if isinstance(phase, (list, tuple)):
            return default_phase(" ".join(phase))
        elif isinstance(phase, str):
            return default_phase(phase)
        else:
            return phase

    def _process_phases(self, phase, default_phase, list_mode=False):
        # List mode is used for when we expect to build a list of phases, such as Mineral and Solid Phase objects
        if phase != [] and phase is not None:
            if isinstance(phase, (str)):
                phase = [phase]
            if list_mode:
                phases = []
                for p in phase:
                    phases.append(self._process_phase(p, default_phase))
                return phases
            else:
                return self._process_phase(phase, default_phase)
        else:
            return None

    def track_registered_phases(self, phase, inputs):
        """this is used to add up all phases provided as inputs or registed phases by user,
        for example, user might porbice CO2(g) as a gas input wants to track also N2(g) via register_gas_phase,
        this ensure we updates all the phases with all requested values"""
        if inputs is None:
            return inputs
        if phase not in self.registered_phases:
            self.registered_phases[phase] = []
        if isinstance(self.registered_phases[phase], (str, list)):
            if isinstance(inputs, str):
                inputs = [inputs]
            if isinstance(inputs, list):
                for i in inputs:
                    if i not in self.registered_phases[phase]:
                        self.registered_phases[phase].append(i)
            else:
                self.registered_phases[phase] = inputs
        return self.registered_phases[phase]

    def register_aqueous_phase(self, aqueous_phases=None):
        """register possible mineral phases"""
        if aqueous_phases != [] and aqueous_phases is not None:
            aqueous_phases = self.track_registered_phases(
                RktInputTypes.aqueous_phase, aqueous_phases
            )
            if isinstance(aqueous_phases, (list, str)):
                self.aqueous_phase = rkt.AqueousPhase(rkt.speciate(aqueous_phases))
            else:
                self.aqueous_phase = aqueous_phases

    def register_condensed_phase(self, condensed_phase=None):
        """register possible condensed phases"""
        if condensed_phase != [] and condensed_phase is not None:
            condensed_phase = self.track_registered_phases(
                RktInputTypes.condensed_phase, condensed_phase
            )
            if isinstance(condensed_phase, (list, str)):
                self.condensed_phase = rkt.CondensedPhases(
                    rkt.speciate(condensed_phase)
                )
            else:
                self.condensed_phase = condensed_phase

    def register_liquid_phase(self, liquid_phases=None):
        """register possible mineral phases"""
        if liquid_phases != [] and liquid_phases is not None:
            liquid_phases = self.track_registered_phases(
                RktInputTypes.liquid_phases, liquid_phases
            )
            if isinstance(liquid_phases, (list, str)):
                self.liquid_phase = rkt.LiquidPhase(rkt.speciate(liquid_phases))
            else:
                self.liquid_phase = liquid_phases

    def register_mineral_phases(self, mineral_phases=None):
        """register possible mineral phases"""
        mineral_phases = self.track_registered_phases(
            RktInputTypes.mineral_phase, mineral_phases
        )
        self.mineral_phase = self._process_phases(
            mineral_phases, rkt.MineralPhase, list_mode=True
        )

    def register_solid_phases(self, solid_phases=None):
        """register possible solid phases"""
        solid_phases = self.track_registered_phases(
            RktInputTypes.solid_phase, solid_phases
        )
        self.solid_phase = self._process_phases(
            solid_phases, rkt.SolidPhase, list_mode=True
        )

    def register_gas_phase(self, gas_phase=None):
        """register possible gas phases"""
        gas_phase = self.track_registered_phases(RktInputTypes.gas_phase, gas_phase)
        self.gas_phase = self._process_phases(
            gas_phase, rkt.GaseousPhase, list_mode=False
        )

    def register_ion_exchange_phase(self, ion_exchange_phase=None):
        """register possible ion exchange phases"""
        ion_exchange_phase = self.track_registered_phases(
            RktInputTypes.ion_exchange_phase, ion_exchange_phase
        )
        self.ion_exchange_phase = self._process_phases(
            ion_exchange_phase, rkt.IonExchangePhase, list_mode=False
        )

    def set_database(self, dbtype="PhreeqcDatabase", database="pitzer.dat"):
        """set data base of reaktoro"""
        """ assume that if database is string we need to find and init in reaktoro
        other wise assume we received activated reaktoro db"""
        if isinstance(dbtype, str):
            self.database = getattr(rkt, dbtype)(database)
        else:
            self.database = dbtype
        self.database_species = [specie.name() for specie in self.database.species()]
        self.database_elements = [
            element.symbol() for element in self.database.elements()
        ]

    def _process_activity(
        self, activity_model, state_of_matter=None, default_activity_model=None
    ):
        """this will process give activity model if user provides a strng it will
        find it on reaktoro and initialize it either by it self or passing in state of matter argument
        if user provides intialized reaktoro activity model we use it directly."""

        def get_func(activity_model, state_of_matter):
            if isinstance(activity_model, str):
                if state_of_matter is None:
                    return getattr(rkt, activity_model)()
                else:
                    return getattr(rkt, activity_model)(state_of_matter)
            else:
                return activity_model

        if activity_model is None:
            activity_model = default_activity_model
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
        self,
        activity_model=None,
    ):
        """set activity model of aqueous phases in reaktoro"""
        self.process_registered_inputs()
        if self.aqueous_phase is not None:
            activity_model = self._process_activity(
                activity_model, default_activity_model="ActivityModelIdealAqueous"
            )
            self.aqueous_phase.set(activity_model)

    def set_liquid_phase_activity_model(
        self,
        activity_model=None,
    ):
        """set activity model of liquid phases in reaktoro"""
        self.process_registered_inputs()
        if self.liquid_phase is not None:
            activity_model = self._process_activity(
                activity_model, default_activity_model="ActivityModelIdealSolution"
            )
            self.liquid_phase.set(activity_model)

    def set_gas_phase_activity_model(self, activity_model=None):
        """set activity model of gas phases in reaktoro"""
        self.process_registered_inputs()
        activity_model = self._process_activity(
            activity_model, default_activity_model="ActivityModelIdealGas"
        )
        if self.gas_phase is not None:
            self.gas_phase.set(activity_model)

    def set_condensed_phase_activity_model(self, activity_model=None):
        """set activity model of ccondensed phases in reaktoro"""
        self.process_registered_inputs()
        activity_model = self._process_activity(
            activity_model,
            default_activity_model="ActivityModelIdealSolution",
            state_of_matter=rkt.StateOfMatter.Liquid,
        )
        if self.condensed_phase is not None:
            self.condensed_phase.set(activity_model)

    def set_mineral_phase_activity_model(self, activity_model=None):
        """set activity model of mineral phases in reaktoro"""
        self.process_registered_inputs()
        activity_model = self._process_activity(
            activity_model,
            state_of_matter=rkt.StateOfMatter.Solid,
            default_activity_model="ActivityModelIdealSolution",
        )
        if self.mineral_phase is not None:
            for phase in self.mineral_phase:
                phase.set(activity_model)

    def set_solid_phase_activity_model(self, activity_model=None):
        """set activity model of solid phases in reaktoro"""
        self.process_registered_inputs()
        activity_model = self._process_activity(
            activity_model,
            state_of_matter=rkt.StateOfMatter.Solid,
            default_activity_model="ActivityModelIdealSolution",
        )
        if self.solid_phase is not None:
            for phase in self.solid_phase:
                phase.set(activity_model)

    def set_ion_exchange_phase_activity_model(self, activity_model=None):
        """set activity model of mineral phases in reaktoro"""
        self.process_registered_inputs()
        activity_model = self._process_activity(
            activity_model,
            default_activity_model="ActivityModelIonExchange",
        )
        if self.ion_exchange_phase is not None:
            self.ion_exchange_phase.set(activity_model)

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

        """ assemble all input species into single list 
        and pass into aqueous phase"""
        all_species = self.inputs.all_species
        """register liquid type phases"""
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

        """ register gas phases"""
        if len(self.inputs.species_list[RktInputTypes.mineral_phase]) > 0:
            self.register_mineral_phases(
                self.inputs.species_list[RktInputTypes.mineral_phase]
            )
        """ register solid phase"""
        if len(self.inputs.species_list[RktInputTypes.solid_phase]) > 0:
            self.register_solid_phases(
                self.inputs.species_list[RktInputTypes.solid_phase]
            )
        if len(self.inputs.species_list[RktInputTypes.gas_phase]) > 0:
            self.register_gas_phase(self.inputs.species_list[RktInputTypes.gas_phase])
        """ register ix phases"""
        if len(self.inputs.species_list[RktInputTypes.ion_exchange_phase]) > 0:
            self.register_ion_exchange_phase(
                self.inputs.species_list[RktInputTypes.ion_exchange_phase]
            )

    def build_state(self):
        """this will build reaktor states"""
        phases = []
        if self.aqueous_phase is not None:
            phases.append(self.aqueous_phase)
        if self.liquid_phase is not None:
            phases.append(self.liquid_phase)
        if self.condensed_phase is not None:
            phases.append(self.condensed_phase)
        if self.mineral_phase is not None:
            for m in self.mineral_phase:
                phases.append(m)
        if self.gas_phase is not None:
            phases.append(self.gas_phase)
        if self.ion_exchange_phase is not None:
            phases.append(self.ion_exchange_phase)
        self.system = rkt.ChemicalSystem(
            self.database,
            *phases,
        )
        self.state = rkt.ChemicalState(self.system)
        self.set_rkt_state()

    def set_rkt_state(self):
        """sets initial rkt state using user provided inputs"""
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
        """ set apparent species if used """
        for phase in self.inputs.registered_phases:
            if self.inputs.composition_is_elements[phase] == False:
                for species in self.inputs.species_list[phase]:
                    if species in self.inputs:  # user might not provide all
                        if self.inputs[species].get_value() != 0:
                            unit = self.inputs[species].main_unit
                            if unit == "dimensionless":
                                """assume correct units are provided"""
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
