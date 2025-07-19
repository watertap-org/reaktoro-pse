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

from idaes.core.base.process_base import declare_process_block_class, ProcessBlockData
import idaes.logger as idaeslog
import math

from pyomo.common.config import ConfigValue, IsInstance
from pyomo.core.base.var import IndexedVar
from pyomo.environ import Block, Constraint, Var, units as pyunits

from reaktoro_pse.core.reaktoro_state import ReaktoroState
from reaktoro_pse.core.reaktoro_inputs import (
    ReaktoroInputSpec,
)

from reaktoro_pse.core.util_classes.rkt_inputs import RktInputTypes
from reaktoro_pse.core.reaktoro_outputs import (
    ReaktoroOutputSpec,
)
from types import SimpleNamespace
from reaktoro_pse.core.reaktoro_jacobian import ReaktoroJacobianSpec
from reaktoro_pse.core.reaktoro_solver import ReaktoroSolver
from reaktoro_pse.core.reaktoro_block_builder import ReaktoroBlockBuilder
from reaktoro_pse.core.reaktoro_coupled_solver import ReaktoroCoupledSolver
from reaktoro_pse.parallel_tools.reaktoro_block_manager import ReaktoroBlockManager

from reaktoro_pse.reaktoro_block_config.jacobian_options import JacobianOptions
from reaktoro_pse.reaktoro_block_config.hessian_options import HessianOptions
from reaktoro_pse.reaktoro_block_config.reaktoro_solver_options import (
    ReaktoroSolverOptions,
)
from pyomo.contrib.pynumero.interfaces.external_grey_box import (
    ExternalGreyBoxModel,
)
from reaktoro_pse.reaktoro_block_config.input_options import (
    PhaseInput,
    SystemInput,
)

_log = idaeslog.getLogger(__name__)

__author__ = "Alexander V. Dudchenko"


@declare_process_block_class("ReaktoroBlock")
class ReaktoroBlockData(ProcessBlockData):
    CONFIG = ProcessBlockData.CONFIG()
    CONFIG.declare(
        RktInputTypes.aqueous_phase,
        PhaseInput().get_dict(
            include_pH=True, aqueous_phase=True, include_solvent_species=True
        ),
        # description="Defines aqueous phase properties",
        # doc="""Defines aqueous phase properties""",
    )
    CONFIG.declare(
        RktInputTypes.liquid_phase,
        PhaseInput().get_dict(
            include_pH=True, aqueous_phase=False, include_solvent_species=True
        ),
        # description="Defines aqueous liquid properties",
        # doc="""Defines aqueous liquid properties""",
    )
    CONFIG.declare(
        RktInputTypes.gas_phase,
        PhaseInput().get_dict(include_speciate_phase_component=True),
        # description="Defines aqueous gas properties",
        # doc="""Defines aqueous gas properties""",
    )
    CONFIG.declare(
        RktInputTypes.condensed_phase,
        PhaseInput().get_dict(),
        # description="Defines aqueous condensed phase properties",
        # doc="""Defines aqueous condensed phase properties""",
    )
    CONFIG.declare(
        RktInputTypes.mineral_phase,
        PhaseInput().get_dict(),
        # description="Defines mineral phase properties",
        # doc="""Defines aqueous mineral phase properties""",
    )
    CONFIG.declare(
        RktInputTypes.solid_phase,
        PhaseInput().get_dict(),
        # description="Defines aqueous liquid properties",
        # doc="""Defines aqueous liquid properties""",
    )
    CONFIG.declare(
        RktInputTypes.ion_exchange_phase,
        PhaseInput().get_dict(),
        # description="Defines ion exchange properties",
        # doc="""Defines ion exchange properties""",
    )
    CONFIG.declare(
        RktInputTypes.system_state,
        SystemInput().get_dict(),
        # description="Defines system state",
        # doc="""Defines system state""",
    )
    CONFIG.declare(
        RktInputTypes.system_state_modifier,
        SystemInput().get_dict(enable_pH=False),
        # description="Defines system state modifier (used if speciation block is built on property block)",
        # doc="""Defines system state modifier (used if speciation block is built on property block)
        #         e.g. define changes in system pressure, temperature or enthalpy
        #         Note: Changes in pH are not supported, please use acid or base addition to modulate pH
        #        """,
    )
    CONFIG.declare(
        "build_speciation_block",
        ConfigValue(
            default=False,
            domain=bool,
            description="Construct speciation block (use when exact composition is unknown and modified state is desired)",
            doc="""If enabled, will construct a speciation block to calculate initial equilibrium state 
             at specified pH. Any chemicals will be added to the property block, enable this when exact composition is unknown and the 
             property state is modified through addition of chemicals or formation of phases that modify final state""",
        ),
    )
    CONFIG.declare(
        "direct_speciation_to_property_block_coupling",
        ConfigValue(
            default=True,
            domain=bool,
            description="This mode will directly couple speciation block to property block",
            doc="""
            - If True this will build a single graybox model that contains both the speciation and property blocks, the true species from speciation block will
            be directly pass the property block, and chain rule will be used to compute the jacobian.
            - If False, this will build two separate graybox models, one for speciation and one for property block, the true species from speciation block
            will be passed to property block through pyomo constraints. This mode is less stable as property block can fail to due to solver taking steps that result
            in infeasible true species composition.
            """,
        ),
    )
    CONFIG.declare(
        "build_graybox_model",
        ConfigValue(
            default=True,
            domain=bool,
            description="This will build a graybox model",
            doc="""
            Use this option to disable building graybox models and generate a number of configured speciation blocks that can be passed
            together into a single property block in direct coupling mode. These configured Reaktoro model blocks can be passed into another ReaktoroBlock through
            'external_speciation_reaktoro_blocks' option.
            Use this to perform mixing operations, (e.g. speciate multiple streams and then mix them together in a single property block).
            """,
        ),
    )

    CONFIG.declare(
        "exclude_species_list",
        ConfigValue(
            default=None,
            domain=list,
            description="List of species to exclude from speciation in reaktoro",
            doc="""The species that should not be included when speciating aqueous, liquid, or condensed phases, this 
            will completely exclude the species from thermodynamic speciation. This is not the same as
            rkt.exclude. """,
        ),
    )
    CONFIG.declare(
        "speciation_block_exclude_species",
        ConfigValue(
            default=None,
            domain=list,
            description="List of species to not include in output of speciation block and input into property block",
            doc="""If provided then selected species will not be passed from speciation block and into property block. 
            This is useful for removing species that do not go above 0 and cause issues during solve. No warning or checks
            are performed by Reaktoro-PSE for scenario where the given species is above 0 at any time during solve. The user
            should manually verify that the specie did not go above 0 or expected value. """,
        ),
    )
    CONFIG.declare(
        "exact_speciation",
        ConfigValue(
            default=False,
            domain=bool,
            description="Defines if exact speciation is provided ",
            doc="If true, then the problem will constrain all elements to equal to provided inputs after dissolution",
        ),
    )

    CONFIG.declare(
        "build_speciation_block_with_phases",
        ConfigValue(
            default=False,
            domain=bool,
            description="Defines if gas or mineral phases should be included in speciation block when built",
            doc="""
            Defines if mineral and gas phases are added to speciation block if constructed. Generally, 
             it is not need to include gas and mineral phases when getting initial equilibrium state before 
             performing chemical addition. 
             Phases are always built when getting properties block. 
             Ion exchange phases are always included on all blocks
             """,
        ),
    )
    CONFIG.declare(
        "database_file",
        ConfigValue(
            default="pitzer.dat",
            domain=str,
            description="Data base file",
            doc="Defines which data base file to use from reaktoro",
        ),
    )
    CONFIG.declare(
        "database",
        ConfigValue(
            default="PhreeqcDatabase",
            description="Database to use with reaktoro (string, or initialized reaktoro database)",
            doc="""Defines which database to use in reaktoro.
            If providing database as string, then must also include database_file name as string, 
            Otherwise, can be supplied as initialized reaktoro database object e.g. database = reaktoro.PhreeqcDatabase('pitzer.dat')""",
        ),
    )

    CONFIG.declare(
        "dissolve_species_in_reaktoro",
        ConfigValue(
            default=True,
            domain=bool,
            description="Defines if species should be converted to elements using reaktoro or Pyomo",
            doc="""
            The equilibrium calculation requires element amounts as an input,
            since normally species are provided, they need to be converted to elements, this can be done 
            using pyomo constraints (Set to False) or reaktoro (Set to True).""",
        ),
    )

    CONFIG.declare(
        "assert_charge_neutrality",
        ConfigValue(
            default=True,
            domain=bool,
            description="Defines if charge neutrality should be maintained",
            doc="""Will enforce charge neutrality if set to true using specified charge neutrally ion""",
        ),
    )
    CONFIG.declare(
        "charge_neutrality_ion",
        ConfigValue(
            default="Cl",
            domain=str,
            description="How Reaktoro maintains charge neutrality (pH, ion, or element)",
            doc="""This will unfix specified ion during equilibrium calculations while enforcing charge==0 constraint
              in reaktoro, if exact speciation is provided, pH will be used by default""",
        ),
    )
    CONFIG.declare(
        "assert_charge_neutrality_on_property_block",
        ConfigValue(
            default=True,
            domain=bool,
            description="Defines if charge neutrality should be applied to property_block when build_speciation_block=True or when exact speciation is set to True",
            doc="""
            When user provides chemical inputs, its assumed that user wants to modify an equilibrated state, 
            as such a speciation and property block will be built. The speciation block will likely require charge balancing on ion basis (e.g Cl),
            while property block that received exact speciation do not necessarily need charge balance, 
            If this is set to True the property block will charge balance on pH
            """,
        ),
    )

    CONFIG.declare(
        "chemistry_modifier",
        ConfigValue(
            default=None,
            domain=IsInstance((dict, IndexedVar)),
            description="Chemicals being added to reaktoro",
            doc="""This should be a dictionary, or indexed var that defines type of chemical (key) and pyomo var that defines 
            amount being added""",
        ),
    )
    CONFIG.declare(
        "chemistry_modifier_indexed",
        ConfigValue(
            default=True,
            domain=bool,
            description="Chemical addition is indexed",
            doc="""
            Option that defines how to treat input variable when building indexed reaktoroBlock":
                - If true, the input has same indexing as block, and each indexed input will be passed into respective indexed reaktoroBlock
                - If false, all indexed blocks will get same input""",
        ),
    )
    CONFIG.declare(
        "register_new_chemistry_modifiers",
        ConfigValue(
            default=None,
            domain=dict,
            description="Defines speciation of chemistry modifiers if they are not available by default",
            doc="""
            Adds a new chemistry modifier with specific name and element breakdown 
            {'H2O_removal:{'H':2,'O':1},
            'CaO':'Ca':1,'O':1}.""",
        ),
    )
    CONFIG.declare(
        "outputs",
        ConfigValue(
            default=None,
            domain=IsInstance((dict, IndexedVar)),
            description="Defines outputs that should be returned by reaktoro",
            doc="""
            This will configure all outputs that should be produced by the block, the block will either use 
            provided pyomo var as output, or create a var on the block. Provide as a dictionary or indexed var where keys are 
            desired properties. For example:
            {
            (scalingTendency, Calcite):m.calcite_var,- will use provided var for output
            pH: None, - will create a new var as m.reaktoroBlock.outputs[(pH,None)]
            "speciesAmount":True - this will force reaktor to return all species
            }
            """,
        ),
    )
    CONFIG.declare(
        "reaktoro_block_manager",
        ConfigValue(
            default=None,
            domain=IsInstance(ReaktoroBlockManager),
            description="Reaktoro block manager for parallelizing reaktoro solves",
            doc="""
            Option to provide a reaktoro block manager which would manage all blocks built on a model and
            allow their parallel execution. When using ReaktorBlockManager, make sure to run
            ReaktoroBlockManager.build_reaktoro_blocks() after constructing all of ReaktoroBlocks and 
            before doing any initialization calls, or interacting with ReaktoroBlocks themself .""",
        ),
    )
    CONFIG.declare(
        "external_speciation_reaktoro_blocks",
        ConfigValue(
            default=None,
            domain=IsInstance(list),
            description="Provides a list of external reaktoro blocks for speciation",
            doc="""
            Used to for mixing different speciation blocks when useing "species_direct" options
            Build external reaktoro block and then pass them into this option. .""",
        ),
    )
    CONFIG.declare("jacobian_options", JacobianOptions().get_dict())
    CONFIG.declare("hessian_options", HessianOptions().get_dict())
    CONFIG.declare(
        "reaktoro_solve_options",
        ReaktoroSolverOptions().get_dict(advanced_options=True),
    )
    CONFIG.declare(
        "reaktoro_presolve_options",
        ReaktoroSolverOptions().get_dict(presolve_options=True),
    )

    def build(self):
        super().build()
        # configure state
        self.direct_coupling_mode = False
        self.external_block_mode = False
        if self.config.build_speciation_block:
            # create speciation block and then property block
            if self.config.direct_speciation_to_property_block_coupling:
                self.direct_coupling_mode = True
                self.speciation_block = SimpleNamespace()
            else:
                self.speciation_block = Block()
            self.build_rkt_state(self.speciation_block, speciation_block=True)
            self.build_rkt_inputs(self.speciation_block, speciation_block=True)
            self.build_rkt_outputs(self.speciation_block, speciation_block=True)
            self.build_rkt_jacobian(self.speciation_block)
            self.build_rkt_solver(self.speciation_block, speciation_block=True)

            self.build_gray_box(self.speciation_block, speciation_block=True)

            self.build_rkt_state(
                self,
                speciation_block=False,
                speciation_block_built=True,
            )
            self.build_rkt_inputs(
                self, speciation_block=False, speciation_block_built=True
            )
        elif (
            self.config.build_speciation_block == False
            and self.config.external_speciation_reaktoro_blocks is not None
        ):
            self.external_block_mode = True
            self.direct_coupling_mode = True
            self.speciation_block = SimpleNamespace()
            self.build_gray_box(self.speciation_block, speciation_block=True)
            self.build_rkt_state(
                self,
                speciation_block=False,
                speciation_block_built=True,
            )
            self.build_rkt_inputs(
                self, speciation_block=False, speciation_block_built=True
            )
        else:
            # create property block only
            self.build_rkt_state(self)
            self.build_rkt_inputs(self)
        if (
            self.config.external_speciation_reaktoro_blocks is not None
            and self.direct_coupling_mode == False
        ):
            raise TypeError(
                "Provide external_speciation_reaktoro_blocks only when using 'species_direct' option"
            )
        self.build_rkt_outputs(self)
        self.build_rkt_jacobian(self)
        self.build_rkt_solver(self)
        if self.config.build_graybox_model:
            self.build_gray_box(self)

    def build_rkt_state(
        self,
        block,
        speciation_block=False,
        speciation_block_built=False,
    ):
        """This will build our reaktoro state specified block.
        The keyword arguments are for automatic configuration of speciation and property blocks

        Keywords:
        block -- pyomo block to build the model on
        speciation_block -- sets to build speciation block which by default does not
            add mineral or gas phases (unless configured to do so), and only outputs speciesAmount from the rktModel
        speciation_block_built -- this should only be True if speciation block was built, in this case the model
            will build user requested outputs, not provide user supplied pH and disable charge neutrality

        Defining which indices should be indexed
         if block is not indexed, index is None
         if specific input is set to be not Indexed the index will be set to None

        Function to return index only when requested by user and
        when we specify if speciation_block_built was built
        """

        def building_prop_block_after_speciation():
            if speciation_block == False and speciation_block_built:
                return True
            else:
                return False

        def get_indexing(index_enabled, speciation_block_built=False):
            if index_enabled == False or speciation_block_built:
                return None
            else:
                return self.index()

        # return False when building property block after building speciation block
        #  (e.g. When speciation_block=True and speciation_block_built=False) """

        def return_false_option(input_option):
            if building_prop_block_after_speciation():
                return False
            else:
                return input_option

        # return None when building property block after building speciation block
        #  (e.g. When speciation_block=True and speciation_block_built=False) """

        def return_none_option(input_option):
            if building_prop_block_after_speciation():
                return None
            else:
                return input_option

        def return_empty_dict_option(input_option):
            if (
                building_prop_block_after_speciation()
                and self.external_block_mode == False
            ):
                return {}
            else:
                return input_option

        def get_prop_state(initial_state, modified_state, variable):
            if getattr(modified_state, variable) is None:
                return initial_state
            else:
                return modified_state

        def get_phases(phase_type):
            if (
                building_prop_block_after_speciation()
                and getattr(self.config, phase_type).phase_components is None
            ):
                return self.speciation_block.rkt_state.phase_manager.registered_phases[
                    phase_type
                ].phase_list
            else:
                return getattr(self.config, phase_type).phase_components

        block.rkt_state = ReaktoroState()

        # setup database
        block.rkt_state.set_database(
            dbtype=self.config.database, database=self.config.database_file
        )

        # setup input for different phase options
        for phase in RktInputTypes.supported_phases:
            options = getattr(self.config, phase)
            if building_prop_block_after_speciation():
                composition_is_elements = False
            else:
                composition_is_elements = options.composition_is_elements
            block.rkt_state.set_input_options(
                phase,
                convert_to_rkt_species=return_false_option(
                    options.convert_to_rkt_species
                ),
                species_to_rkt_species_dict=options.species_to_rkt_species_dict,
                composition_is_elements=composition_is_elements,
            )

        # setup system inputs
        if building_prop_block_after_speciation():
            temperature_state = get_prop_state(
                self.config.system_state,
                self.config.system_state_modifier,
                "temperature",
            )
            pressure_state = get_prop_state(
                self.config.system_state, self.config.system_state_modifier, "pressure"
            )
            enthalpy_state = get_prop_state(
                self.config.system_state, self.config.system_state_modifier, "enthalpy"
            )
            pH_state = self.config.system_state
        else:
            # use initial state
            temperature_state = self.config.system_state
            pressure_state = self.config.system_state
            enthalpy_state = self.config.system_state
            pH_state = self.config.system_state
        config_ph = return_none_option(pH_state.pH)
        config_ph_index = get_indexing(pH_state.pH_indexed)

        # setup aqueous inputs
        aqueous_input_composition = self.config.aqueous_phase.composition
        liquid_input_composition = self.config.liquid_phase.composition
        condensed_input_composition = self.config.condensed_phase.composition

        block.rkt_state.register_system_inputs(
            temperature=temperature_state.temperature,
            pressure=pressure_state.pressure,
            enthalpy=enthalpy_state.enthalpy,
            enthalpy_index=get_indexing(enthalpy_state.enthalpy_indexed),
            temperature_index=get_indexing(temperature_state.temperature_indexed),
            pressure_index=get_indexing(pressure_state.pressure_indexed),
            pH=config_ph,
            pH_index=config_ph_index,
        )

        if building_prop_block_after_speciation():
            # we need to ensure when we provide initial input compo into
            # speciation block we don't have extremely high ion concentration
            # these value swill be overwritten during initialization anyway

            for ion, obj in self.speciation_block.outputs.items():
                if self.config.aqueous_phase.fixed_solvent_specie in ion:
                    obj.set_value(obj.value * 10)
                else:
                    obj.set_value(obj.value / 1000)

            if aqueous_input_composition is not {}:
                aqueous_input_composition = self.speciation_block.outputs
                liquid_input_composition = {}
                condensed_input_composition = {}
            elif liquid_input_composition is not {}:
                aqueous_input_composition = {}
                liquid_input_composition = self.speciation_block.outputs
                condensed_input_composition = {}
            elif condensed_input_composition is not {}:
                aqueous_input_composition = {}
                liquid_input_composition = {}
                condensed_input_composition = self.speciation_block.outputs
            else:
                raise ValueError(
                    "Speciation block requires that either liquid or aqueous phase is provided"
                )
        else:
            self.speciation_input_scale = False
            aqueous_input_composition = self.config.aqueous_phase.composition
        block.rkt_state.register_species_to_exclude(self.config.exclude_species_list)
        block.rkt_state.register_aqueous_inputs(
            composition=aqueous_input_composition,
            composition_index=get_indexing(
                self.config.aqueous_phase.composition_indexed, speciation_block_built
            ),
        )
        block.rkt_state.register_liquid_inputs(
            composition=liquid_input_composition,
            composition_index=get_indexing(
                self.config.liquid_phase.composition_indexed, speciation_block_built
            ),
        )
        block.rkt_state.register_condensed_inputs(
            composition=condensed_input_composition,
            composition_index=get_indexing(
                self.config.condensed_phase.composition_indexed, speciation_block_built
            ),
        )
        block.rkt_state.register_gas_inputs(
            composition=return_empty_dict_option(self.config.gas_phase.composition),
            composition_index=get_indexing(self.config.gas_phase.composition_indexed),
        )
        block.rkt_state.register_mineral_inputs(
            composition=return_empty_dict_option(self.config.mineral_phase.composition),
            composition_index=get_indexing(
                self.config.mineral_phase.composition_indexed
            ),
        )
        block.rkt_state.register_ion_exchange_inputs(
            composition=return_empty_dict_option(
                self.config.ion_exchange_phase.composition
            ),
            composition_index=get_indexing(
                self.config.ion_exchange_phase.composition_indexed
            ),
        )
        # register phases
        if speciation_block == False or self.config.build_speciation_block_with_phases:
            # dont add phases if we are speciating
            block.rkt_state.register_aqueous_phase(get_phases("aqueous_phase"))
            block.rkt_state.register_liquid_phase(get_phases("liquid_phase"))
            block.rkt_state.register_solid_phases(get_phases("solid_phase"))
            block.rkt_state.register_condensed_phase(
                get_phases("condensed_phase"),
            )
            block.rkt_state.register_gas_phase(
                get_phases("gas_phase"), self.config.gas_phase.speciate_phase_component
            )
            block.rkt_state.register_mineral_phases(get_phases("mineral_phase"))
            block.rkt_state.register_ion_exchange_phase(
                get_phases("ion_exchange_phase")
            )

        # setup activity models - if no phases present they will do nothing
        block.rkt_state.set_aqueous_phase_activity_model(
            self.config.aqueous_phase.activity_model
        )
        block.rkt_state.set_liquid_phase_activity_model(
            self.config.liquid_phase.activity_model
        )
        block.rkt_state.set_gas_phase_activity_model(
            self.config.gas_phase.activity_model
        )
        block.rkt_state.set_condensed_phase_activity_model(
            self.config.condensed_phase.activity_model
        )
        block.rkt_state.set_solid_phase_activity_model(
            self.config.solid_phase.activity_model
        )
        block.rkt_state.set_mineral_phase_activity_model(
            self.config.mineral_phase.activity_model
        )
        block.rkt_state.set_ion_exchange_phase_activity_model(
            self.config.ion_exchange_phase.activity_model
        )
        # build state
        block.rkt_state.build_state()

    def build_rkt_inputs(
        self,
        block,
        speciation_block=False,
        speciation_block_built=False,
    ):
        """this will build our rkt inputs specified block.
        The keyword arguments are for automatic configuration of speciation and property blocks

        Keywords:
        block -- pyomo block to build the model on
        speciation_block -- sets to build speciation block which by default does not
            add mineral or gas phases (unless configured to do so), and only outputs speciesAmount from the rktModel
        speciation_block_built -- this should only be True if speciation block was built, in this case the model
            will build user requested outputs, not provide user supplied pH and disable charge neutrality
        """

        # get index for chemicals - refer to build_rkt_state on indexing notes
        chemistry_modifier_indexed = self.index()
        if self.config.chemistry_modifier_indexed == False:
            chemistry_modifier_indexed = None

        block.rkt_inputs = ReaktoroInputSpec(block.rkt_state)

        # add chemical only if its not a speciation block (normal mode)
        if speciation_block == False:
            block.rkt_inputs.register_modifier(
                self.config.register_new_chemistry_modifiers
            )
            if self.config.chemistry_modifier is not None:
                block.rkt_inputs.register_chemistry_modifiers(
                    self.config.chemistry_modifier,
                    index=chemistry_modifier_indexed,
                )
            block.rkt_inputs.register_open_species(
                self.config.reaktoro_solve_options.open_species_on_property_block
            )
        else:
            block.rkt_inputs.register_open_species(
                self.config.reaktoro_solve_options.open_species_on_speciation_block
            )

        # register solvents and elements to open for aqueous and liquid phase
        block.rkt_inputs.register_fixed_solvent_specie(
            RktInputTypes.aqueous_phase, self.config.aqueous_phase.fixed_solvent_specie
        )
        block.rkt_inputs.register_fixed_solvent_specie(
            RktInputTypes.liquid_phase, self.config.liquid_phase.fixed_solvent_specie
        )
        block.rkt_inputs.register_free_elements(self.config.aqueous_phase.free_element)
        block.rkt_inputs.register_free_elements(self.config.liquid_phase.free_element)
        # register charge neutrality
        # default charge balance is for non exact speciation on property block
        # if this is not speciation block and we are not given exact speciation
        # the charge balance on user defined ion
        exact_speciation = self.config.exact_speciation
        assert_charge_neutrality = self.config.assert_charge_neutrality
        ion_for_balancing = self.config.charge_neutrality_ion

        # if this is not speciation block and we are given exact speciation
        # we charge balance on pH
        if speciation_block == False and (exact_speciation or speciation_block_built):
            # only do so if we have 'H+' in species
            ion_for_balancing = "pH"
            if (
                "H+" not in block.rkt_state.database_species
                or self.config.assert_charge_neutrality_on_property_block == False
            ):
                assert_charge_neutrality = False
            exact_speciation = True
        block.rkt_inputs.register_charge_neutrality(
            assert_neutrality=assert_charge_neutrality,
            ion=ion_for_balancing,
        )
        block.rkt_inputs.configure_specs(
            dissolve_species_in_rkt=self.config.dissolve_species_in_reaktoro,
            exact_speciation=exact_speciation,
        )
        block.rkt_inputs.build_input_specs()

    def build_rkt_outputs(self, block, speciation_block=False):
        """this will build rkt outputs specified block.
        The keyword arguments are for automatic configuration of speciation and property blocks

        Keywords:
        block -- pyomo block to build the model on
        speciation_block -- sets to build speciation block which by default does not
            add mineral or gas phases (unless configured to do so), and only outputs speciesAmount from the rktModel
        """

        # configure outputs
        index = self.index()

        block.rkt_outputs = ReaktoroOutputSpec(block.rkt_state)

        if self.config.outputs is None:
            raise ValueError("Outputs must be provided!")
        if speciation_block:
            # when speciating we only want species amounts as output
            block.rkt_outputs.register_output(
                "speciesAmount",
                get_all_indexes=True,
                ignore_indexes=self.config.speciation_block_exclude_species,
            )

        else:

            # build user requested outputs
            for output_key, output_var in self.config.outputs.items():
                if index is None or index in output_key:
                    if isinstance(output_key, tuple):
                        if len(output_key) == 2:
                            output_key, output_prop = output_key
                        if len(output_key) == 3:
                            _, output_key, output_prop = output_key
                    else:
                        output_prop = None
                    if isinstance(output_var, bool):
                        if (
                            self.config.exact_speciation == False
                            and output_key == "speciesAmount"
                        ):
                            ignore_species = (
                                self.config.speciation_block_exclude_species
                            )
                        else:
                            ignore_species = None
                        block.rkt_outputs.register_output(
                            output_key,
                            get_all_indexes=True,
                            ignore_indexes=ignore_species,
                        )
                    elif isinstance(output_var, list):
                        block.rkt_outputs.register_output(
                            output_key,
                            get_all_indexes=True,
                            ignore_indexes=output_var,
                        )
                    else:
                        block.rkt_outputs.register_output(
                            output_key, output_prop, pyomo_var=output_var
                        )

    def convert_outputs_to_dict(self):
        """This will convert the outputs to a dictionary"""
        if isinstance(self.config.outputs, dict) == False:
            self.config.outputs = {
                key: self.config.outputs[key] for key in self.config.outputs
            }

    def build_rkt_jacobian(self, block):
        """this will build rkt jacboian on specified block.
        The keyword arguments are for automatic configuration of speciation and property blocks

        Keywords:
        block -- pyomo block to build the model on
        """
        # config outputs
        block.rkt_jacobian = ReaktoroJacobianSpec(block.rkt_state, block.rkt_outputs)
        block.rkt_jacobian.configure_numerical_jacobian(
            jacobian_type=self.config.jacobian_options.numerical_type,
            order=self.config.jacobian_options.numerical_order,
            step_size=self.config.jacobian_options.numerical_step,
        )

    def build_rkt_solver(self, block, speciation_block=False):
        """this will build rkt outputs specified block.
        The keyword arguments are for automatic configuration of speciation and property blocks

        Keywords:
        block -- pyomo block to build the model on
        speciation_block -- sets to build speciation block which by default does not
            add mineral or gas phases (unless configured to do so), and only outputs speciesAmount from the rktModel
        """

        # config solver
        name = str(self)
        if speciation_block:
            name = f"{name}_speciation_block"
        block.rkt_solver = ReaktoroSolver(
            block.rkt_state,
            block.rkt_inputs,
            block.rkt_outputs,
            block.rkt_jacobian,
            block_name=name,
            maximum_failed_solves=self.config.reaktoro_solve_options.max_reaktoro_failed_solves,
        )
        block.rkt_solver.set_system_bounds(
            self.config.system_state.temperature_bounds,
            self.config.system_state.pressure_bounds,
        )

        if speciation_block:
            presolve = self.config.reaktoro_presolve_options.presolve_speciation_block
        else:
            presolve = self.config.reaktoro_presolve_options.presolve_property_block

        block.rkt_solver.set_solver_options(
            tolerance=self.config.reaktoro_solve_options.solver_tolerance,
            epsilon=self.config.reaktoro_solve_options.epsilon,
            presolve=presolve,
            presolve_tolerance=self.config.reaktoro_presolve_options.solver_tolerance,
            presolve_epsilon=self.config.reaktoro_presolve_options.epsilon,
            max_iters=self.config.reaktoro_solve_options.max_iterations,
            presolve_max_iters=self.config.reaktoro_presolve_options.max_iterations,
            hessian_type=self.config.hessian_options.hessian_type,
            bfgs_initialization_type=self.config.hessian_options.bfgs_initialization_type,
            bfgs_init_min_hessian_value=self.config.hessian_options.bfgs_init_min_hessian_value,
            bfgs_init_max_hessian_value=self.config.hessian_options.bfgs_init_max_hessian_value,
            bfgs_init_const_hessian_value=self.config.hessian_options.bfgs_init_const_hessian_value,
            bfgs_hessian_memory=self.config.hessian_options.bfgs_hessian_memory,
            bfgs_epsilon=self.config.hessian_options.bfgs_epsilon,
        )

    def build_gray_box(self, block, speciation_block=False):
        """this will build rkt outputs specified block.
        The keyword arguments are for automatic configuration of speciation and property blocks

        Keywords:
        block -- pyomo block to build the model on
        """
        # build block

        if (
            self.config.reaktoro_block_manager is not None
            and hasattr(self, "speciation_block_data") == False
        ):
            self.speciation_block_data = []
        if self.direct_coupling_mode:
            if speciation_block:
                if (
                    hasattr(block, "rkt_solver") == False
                    and self.config.external_speciation_reaktoro_blocks is not None
                ):
                    working_block = self.config.external_speciation_reaktoro_blocks[0]
                    block.rkt_state = working_block.rkt_state
                    block.rkt_inputs = working_block.rkt_inputs
                    block.rkt_outputs = working_block.rkt_outputs
                    block.rkt_jacobian = working_block.rkt_jacobian
                    block.rkt_solver = working_block.rkt_solver
                    self.config.external_speciation_reaktoro_blocks.pop(0)
                elif (
                    hasattr(block, "rkt_solver") == False
                    and self.config.external_speciation_reaktoro_blocks is None
                ):
                    raise ValueError(
                        "Direct coupling mode is enabled, but no speciation block is provided."
                    )

                solver_blocks = [block.rkt_solver]

                if self.config.external_speciation_reaktoro_blocks is not None:
                    for spc in self.config.external_speciation_reaktoro_blocks:
                        if spc != block:
                            solver_blocks.append(spc.rkt_solver)
                self.coupled_solver = ReaktoroCoupledSolver(solver_blocks)
                block.outputs = self.coupled_solver.outputs

                solver = None
                if self.config.reaktoro_block_manager is not None:
                    self.speciation_block_data.append(
                        {
                            "state": block.rkt_state,
                            "inputs": block.rkt_inputs,
                            "outputs": block.rkt_outputs,
                            "jacobian": block.rkt_jacobian,
                            "solver": block.rkt_solver,
                        }
                    )
                    if self.config.external_speciation_reaktoro_blocks is not None:
                        for spc in self.config.external_speciation_reaktoro_blocks:
                            self.speciation_block_data.append(
                                {
                                    "state": spc.rkt_state,
                                    "inputs": spc.rkt_inputs,
                                    "outputs": spc.rkt_outputs,
                                    "jacobian": spc.rkt_jacobian,
                                    "solver": spc.rkt_solver,
                                }
                            )
            else:
                self.coupled_solver.register_property_solver(block.rkt_solver)
                solver = self.coupled_solver
                prop_block_data = {
                    "state": block.rkt_state,
                    "inputs": block.rkt_inputs,
                    "outputs": block.rkt_outputs,
                    "jacobian": block.rkt_jacobian,
                    "solver": block.rkt_solver,
                }
        else:
            solver = block.rkt_solver
            prop_block_data = None
        if solver is not None:
            block.rkt_block_builder = ReaktoroBlockBuilder(
                block, solver, build_on_init=False
            )
            block.rkt_block_builder.configure_jacobian_scaling(
                jacobian_scaling_type=self.config.jacobian_options.scaling_type,
                user_scaling=self.config.jacobian_options.user_scaling,
                jacobian_scaling_bounds=self.config.jacobian_options.jacobian_scale_bounds,
                update_jacobian_scale_every_solve=self.config.jacobian_options.update_jacobian_scale_every_solve,
                jacobian_scaling_bounds_output_based=self.config.jacobian_options.jacobian_scaling_bounds_output_based,
            )

            if self.config.reaktoro_block_manager is not None:
                if self.direct_coupling_mode:
                    managed_block = self.config.reaktoro_block_manager.register_block(
                        inputs=solver.input_specs,
                        outputs=solver.output_specs,
                        builder=block.rkt_block_builder,
                        speciation_block=self.speciation_block_data,
                        property_block=prop_block_data,
                    )
                else:
                    managed_block = self.config.reaktoro_block_manager.register_block(
                        state=block.rkt_state,
                        inputs=block.rkt_inputs,
                        outputs=block.rkt_outputs,
                        jacobian=block.rkt_jacobian,
                        solver=solver,
                        builder=block.rkt_block_builder,
                    )
                block.managed_block = managed_block
            else:
                block.rkt_block_builder.build_reaktoro_block()

    # TODO: Update to provide output location (e.g. StringIO)
    def display_jacobian_outputs(self):
        """Displays jacobian output types"""
        if self.config.build_speciation_block and self.direct_coupling_mode == False:
            _log.info("-----Displaying information for speciation block ------")
            self.speciation_block.rkt_jacobian.display_jacobian_output_types()
        _log.info("-----Displaying information for property block ------")
        self.rkt_jacobian.display_jacobian_output_types()

    # TODO:# Update to provide output location (e.g. StringIO)
    def display_jacobian_scaling(self):
        """Displays jacobian scaling"""
        jacobian_scaling = {}
        if self.config.build_speciation_block and self.direct_coupling_mode == False:
            _log.info("-----Displaying information for speciation block ------")
            jac_scale = (
                self.speciation_block.rkt_block_builder.display_jacobian_scaling()
            )
            jacobian_scaling["speciation_block"] = jac_scale
        _log.info("-----Displaying information for property block ------")
        jac_scale = self.rkt_block_builder.display_jacobian_scaling()
        jacobian_scaling["property_block"] = jac_scale
        return jacobian_scaling

    # TODO:# Update to provide output location (e.g. StringIO)
    def display_reaktoro_state(self):
        """Displays reaktoro state"""
        if self.config.build_speciation_block and self.direct_coupling_mode == False:
            _log.info("-----Displaying information for speciation block ------")
            self.speciation_block.rkt_block_builder.display_state()
        _log.info("-----Displaying information for property block ------")
        self.rkt_block_builder.display_state()

    def update_block_scaling(self, update_on_speciation_block=True):
        """This will update the scaling of input and output constraints on all blocks, can be
        helpful when chemistry changes dramatically between solves.
        Keywords:
        update_on_speciation_block -- if scaling should be also updated on speciation block if built.
        """
        if (
            self.config.build_speciation_block
            and update_on_speciation_block
            and self.direct_coupling_mode == False
        ):
            self.speciation_block.rkt_block_builder.initialize_input_variables_and_constraints(
                use_default_scaling=False
            )
            self.speciation_block.rkt_block_builder.set_output_vars_and_scale(
                use_default_scaling=False
            )
        self.rkt_block_builder.initialize_input_variables_and_constraints(
            use_default_scaling=False
        )
        self.rkt_block_builder.set_output_vars_and_scale(use_default_scaling=False)

    def update_jacobian_scaling(
        self, user_scaling_dict=None, set_on_speciation_block=True
    ):
        """This will recalculate jacobian scaling and update with user provided values
        Keywords:
        user_scaling_dict -- Dictionary that contains jacobian keys and scaling block
        set_on_speciation_block -- if scaling should be also set on speciation block if built.
        """
        if (
            self.config.build_speciation_block
            and set_on_speciation_block
            and self.direct_coupling_mode == False
        ):
            self.speciation_block.rkt_block_builder.set_jacobian_scaling()
            self.speciation_block.rkt_block_builder.set_user_jacobian_scaling(
                user_scaling_dict
            )
        self.rkt_block_builder.set_jacobian_scaling()
        self.rkt_block_builder.set_user_jacobian_scaling(user_scaling_dict)

    # TODO: Update to use new initialization method https://idaes-pse.readthedocs.io/en/stable/reference_guides/initialization/developing_initializers.html?highlight=Initializer
    def initialize(self):
        """Initialize reaktoro blocks"""
        if (
            self.config.reaktoro_presolve_options.presolve_during_initialization
            or self.config.reaktoro_presolve_options.presolve_speciation_block
            or self.config.reaktoro_presolve_options.presolve_property_block
        ):
            presolve = True
        else:
            presolve = False
        if self.config.build_speciation_block and self.direct_coupling_mode == False:
            _log.info(f"---initializing speciation block {str(self)}----")
            self.speciation_block.rkt_block_builder.initialize(presolve)
        _log.info(f"---initializing property block {str(self)}----")
        self.rkt_block_builder.initialize(presolve)

    def deactivate(self):
        """Deactivates all constraints and grayboxes"""
        super().deactivate()
        for v in self.component_data_objects(Constraint):
            v.deactivate()
        for v in self.component_data_objects(ExternalGreyBoxModel):
            v.deactivate()

    def activate(self):
        """Activates all constraints and grayboxes"""
        super().activate()
        for v in self.component_data_objects(Constraint):
            v.activate()
        for v in self.component_data_objects(ExternalGreyBoxModel):
            v.activate()
