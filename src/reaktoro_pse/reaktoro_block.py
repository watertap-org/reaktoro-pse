from idaes.core.base.process_base import declare_process_block_class, ProcessBlockData
from pyomo.common.config import ConfigValue, IsInstance
from pyomo.core.base.var import IndexedVar, Var, VarData

from reaktoro_pse.core.reaktoro_state import ReaktoroState
from reaktoro_pse.core.reaktoro_inputs import (
    ReaktoroInputSpec,
)

from reaktoro_pse.core.reaktoro_outputs import (
    ReaktoroOutputSpec,
)
from reaktoro_pse.core.reaktoro_jacobian import ReaktoroJacobianSpec, JacType
from reaktoro_pse.core.reaktoro_solver import (
    ReaktoroSolver,
)

from reaktoro_pse.core.reaktoro_block_builder import (
    ReaktoroBlockBuilder,
    JacScalingTypes,
)
import reaktoro as rkt
from reaktoro_pse.core.reaktoro_gray_box import HessTypes
from pyomo.environ import Block
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)

__author__ = "Alexander Dudchenko"


@declare_process_block_class("ReaktoroBlock")
class ReaktoroBlockData(ProcessBlockData):
    CONFIG = ProcessBlockData.CONFIG()
    CONFIG.declare(
        "build_speciation_block",
        ConfigValue(
            default=False,
            domain=bool,
            description="Construct speciation block (use when exact composition is not known and modified state is desired)",
            doc="""If enabled, will construct a speciation block to calculate initial equilibrium state 
             at specified pH. Any chemicals will be added to the property block, enable this when exact composition is know known and the 
             property state is modified through addition of chemicals or formation of phases that modify final state""",
        ),
    )
    CONFIG.declare(
        "composition",
        ConfigValue(
            default=None,
            domain=IsInstance((dict, IndexedVar)),
            description="Input composition to reaktoro block",
            doc="An input dictionary, or IndexedVar that contains species their amounts",
        ),
    )
    CONFIG.declare(
        "composition_indexed",
        ConfigValue(
            default=True,
            domain=bool,
            description="composition is indexed",
            doc="""Option that defines how to treat input variable when building indexed reaktoroBlock":
                - If true, the input has same indexing as block, and each indexed input will be passed into respective indexed reaktoroBlock
                - If false, all indexed blocks will get same input""",
        ),
    )
    CONFIG.declare(
        "temperature",
        ConfigValue(
            default=None,
            domain=IsInstance((VarData, Var, dict, IndexedVar)),
            description="Input temperature for reaktoro block",
            doc="Var or IndexedVar that references system temperature",
        ),
    )
    CONFIG.declare(
        "temperature_indexed",
        ConfigValue(
            default=True,
            domain=bool,
            description="Temperature is indexed",
            doc="""Option that defines how to treat input variable when building indexed reaktoroBlock":
                - If true, the input has same indexing as block, and each indexed input willbe passed into respective indexed reaktoroBlock
                - If false, all indexed blocks will get same input""",
        ),
    )
    CONFIG.declare(
        "pressure",
        ConfigValue(
            default=None,
            domain=IsInstance((VarData, Var, dict, IndexedVar)),
            description="Input pressure for reaktoro block",
            doc="Var or IndexedVar that references system pressure",
        ),
    )
    CONFIG.declare(
        "pressure_indexed",
        ConfigValue(
            default=True,
            domain=bool,
            description="Pressure is indexed",
            doc="""Option that defines how to treat input variable when building indexed reaktoroBlock":
                - If true, the input has same indexing as block, and each indexed input will be passed into respective indexed reaktoroBlock
                - If false, all indexed blocks will get same input""",
        ),
    )
    CONFIG.declare(
        "pH",
        ConfigValue(
            default=None,
            domain=IsInstance((VarData, Var, dict, IndexedVar)),
            description="Input pH for reaktoro block",
            doc="Var or IndexedVar that references system pH",
        ),
    )
    CONFIG.declare(
        "pH_indexed",
        ConfigValue(
            default=True,
            domain=bool,
            description="pH is indexed",
            doc="""Option that defines how to treat input variable when building indexed reaktoroBlock":
                - If true, the input has same indexing as block, and each indexed input will be passed into respective indexed reaktoroBlock
                - If false, all indexed blocks will get same input""",
        ),
    )
    CONFIG.declare(
        "aqueous_solvent_specie",
        ConfigValue(
            default="H2O",
            domain=str,
            description="Defines aqueous specie to use when speciating system",
            doc="""When speciating, the H2O is fixed, while H and O is unfixed to allow system to equilibrate,
            this should be same specie as being passed in composition, it will be automatically translated to reaktoro notation
            if enabled""",
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
        "mineral_phases",
        ConfigValue(
            default=None,
            description="List or str of mineral phases",
            doc="""
            Mineral phases supported by selected data base
            Accepts:
                string name for phase
                reaktoro.MineralPhase initialized object (e.g reaktoro.MineralPhase('Calcite'))
                List of strings or reaktor.MineralPhases
            """,
        ),
    )
    CONFIG.declare(
        "gas_phases",
        ConfigValue(
            default=None,
            description="List or str for gas phases",
            doc="""
            Gas phases supported by selected data base
            Accepts:
                string name for phase
                reaktoro.GaseousPhase initialized object (e.g reaktoro.GaseousPhase('H2O(g)'))
                List of strings or reaktor.GaseousPhase
                """,
        ),
    )
    CONFIG.declare(
        "ion_exchange_phases",
        ConfigValue(
            default=None,
            description="List or str for Ion exchange phases",
            doc="""Ion exchange phases supported by selected data base
            Accepts:
                string name for phase
                reaktoro.IonExchangePhase initialized object (e.g reaktoro.IonExchangePhase('X'))
                List of strings or reaktor.IonExchangePhase
                """,
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
        "convert_to_rkt_species",
        ConfigValue(
            default=False,
            domain=bool,
            description="Defines if provided species should be converted to RKT notation",
            doc="Enable conversion provided species names to reaktoro names - (currently supports PhreeqC database)",
        ),
    )
    CONFIG.declare(
        "species_to_rkt_species_dict",
        ConfigValue(
            default="default",
            domain=IsInstance((dict, str)),
            description="Dictionary for translating user supplied species to RKT species specific to selected database",
            doc="""
                Dictionary that connects user species to reaktoro data base species (please reference chosen data base in question):
                Dictionary should have following structure {user_species:reaktoro_database_specie}""",
        ),
    )
    CONFIG.declare(
        "composition_is_elements",
        ConfigValue(
            default=False,
            domain=bool,
            description="Defines if provided composition is elements and not species",
            doc="Defines if provided composition is elements and not species",
        ),
    )
    CONFIG.declare(
        "aqueous_phase_activity_model",
        ConfigValue(
            default="ActivityModelIdealAqueous",
            description="Activity model for aqueous phase",
            doc="""
            Defines which activity model to use for aqueous phase in reaktoro
            Accepts:
            String name of activity model
            Initialized reaktoro.ActivityModel (e.g. reaktoro.ActivityModelIdealAqueous())
            Chain of reaktoro activity models (reaktoro.chain(reaktoro.ActivityModelA, reaktoro.ActivityModelB))
            List of activity models strings or reaktoro.ActivityModel initialized objects.
            """,
        ),
    )
    CONFIG.declare(
        "gas_phase_activity_model",
        ConfigValue(
            default="ActivityModelIdealGas",
            description="Activity model for gas phase",
            doc="""Defines which activity model to use in reaktoro
            Accepts:
            String name of activity model
            Initialized reaktoro.ActivityModel (e.g. reaktoro.ActivityModelIdealAqueous())
            Chain of reaktoro activity models (reaktoro.chain(reaktoro.ActivityModelA, reaktoro.ActivityModelB))
            List of activity models strings or reaktoro.ActivityModel initialized objects.""",
        ),
    )
    CONFIG.declare(
        "mineral_phase_activity_model",
        ConfigValue(
            default="ActivityModelIdealSolution",
            description="Activity model for mineral phase",
            doc="""Defines which activity model to use in reaktoro
            Accepts:
            String name of activity model
            Initialized reaktoro.ActivityModel (e.g. reaktoro.ActivityModelIdealAqueous())
            Chain of reaktoro activity models (reaktoro.chain(reaktoro.ActivityModelA, reaktoro.ActivityModelB))
            List of activity models strings or reaktoro.ActivityModel initialized objects.""",
        ),
    )
    CONFIG.declare(
        "ion_exchange_phase_activity_model",
        ConfigValue(
            default="ActivityModelIonExchange",
            description="Activity model for mineral phase",
            doc="""Defines which activity model to use in reaktoro
            Accepts:
            String name of activity model
            Initialized reaktoro.ActivityModel (e.g. reaktoro.ActivityModelIdealAqueous())
            Chain of reaktoro activity models (reaktoro.chain(reaktoro.ActivityModelA, reaktoro.ActivityModelB))
            List of activity models strings or reaktoro.ActivityModel initialized objects.""",
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
            domain=str,
            description="Database to use with reaktoro",
            doc="Defines which database to use in reaktoro",
        ),
    )
    CONFIG.declare(
        "dissolve_species_in_reaktoro",
        ConfigValue(
            default=True,
            domain=bool,
            description="Defines if species should be converted to elements using reaktoro or pyomo",
            doc="""The equilibrium calculation requires element amounts as an input,
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
            description="Ion to use for maintaining charge neutrality",
            doc="""This will unfix specified ion during equilibrium calculations while enforcing charge==0 constraint in reaktoro""",
        ),
    )
    CONFIG.declare(
        "assert_charge_neutrality_on_all_blocks",
        ConfigValue(
            default=False,
            domain=bool,
            description="Defines if charge neutrality should be applied to both speciation and property block",
            doc="""When user provides chemical inputs, its assumed that user wants to modify an equilibrated state, 
            as such a speciation and property block will be built. Charge neutrality would only be applied on speciation block if enabled
            if this option is set to True then charge neutrality will also be applied on property block """,
        ),
    )

    CONFIG.declare(
        "chemical_addition",
        ConfigValue(
            default=None,
            domain=IsInstance((dict, IndexedVar)),
            description="Chemicals being added to reaktoro",
            doc="""This should be a dictionary, or indexed var that defines type of chemical (key) and pyomo var that defines 
            amount being added""",
        ),
    )
    CONFIG.declare(
        "chemical_addition_indexed",
        ConfigValue(
            default=True,
            domain=bool,
            description="Chemical addition is indexed",
            doc="""Option that defines how to treat input variable when building indexed reaktoroBlock":
                - If true, the input has same indexing as block, and each indexed input willbe passed into respective indexed reaktoroBlock
                - If false, all indexed blocks will get same input""",
        ),
    )
    CONFIG.declare(
        "chemical_speciation",
        ConfigValue(
            default=None,
            domain=dict,
            description="Defines if species should be converted to elements using reaktoro or pyomo",
            doc="""The equilibrium calculation requires element amounts as an input,
            since normally species are provided, they need to be converted to elements, this can be done 
            useing pyomo constraints (Set to False) or reaktoro (Set to True).""",
        ),
    )
    CONFIG.declare(
        "outputs",
        ConfigValue(
            default=None,
            domain=IsInstance((dict, IndexedVar)),
            description="Defines outputs that should be returned by reaktoro",
            doc="""This will configure all outputs that should be produced by the block, the block will either use 
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
        "numerical_jac_type",
        ConfigValue(
            default=JacType.average,
            domain=IsInstance((str, JacType)),
            description="Defines method for numerical jacobian approximations",
            doc="""Derivatives for many of the properties in reaktro are not directly available, 
            thus we numerically propagate derivatives from chemical state to methods for estimation of these properties. 
            Two methods are available, average and center_difference
                - average methods takes defined number of derivatives by numerical_jac_order from center points and gets the average of them
                - center_difference methods applies classical taylor difference approximation methods 
            In theory the two should yield same result- but due to round off errors the average method might provide better error dampening. 

            """,
        ),
    )
    CONFIG.declare(
        "numerical_jac_order",
        ConfigValue(
            default=2,
            domain=int,
            description="Defines order of numerical jacobian",
            doc="""This will define how many points to discretize the derivate over 
             - for numerical_jac_type==average - the number of points = order/2 + 1 
             - for numerical_jac_type==center_difference - the number of points equals the order
            """,
        ),
    )
    CONFIG.declare(
        "numerical_jac_step",
        ConfigValue(
            default=1e-6,
            domain=float,
            description="Defines the step to use for numerical descritiazaiton",
            doc="""This will define how small of a step to use for numerical derivative propagation which takes
             the absolute chemical property and multiplies it by chemical property derivative multiplied by step 
                chemical_property_step=chemical_property_absolute_value*chemical_property_derivative*step
            """,
        ),
    )
    CONFIG.declare(
        "jacobian_scaling_type",
        ConfigValue(
            default=JacScalingTypes.variable_scaling,
            domain=IsInstance((str, JacScalingTypes)),
            description="Defines how to scale Jacobian matrix",
            doc="""
            Defines methods for jacobian scaling:
            - if option is no_scaling, jacobian scale will == 1 for all outputs
            - if option is 'variable_scaling' will use output variable scaling factors
            - if option is jacobian_matrix will use actual jac matrix to calculate scaling factors
            - if user_scaling is not None then uses user provided scaling
            """,
        ),
    )
    CONFIG.declare(
        "jacobian_user_scaling",
        ConfigValue(
            default=None,
            domain=IsInstance((float, list, dict)),
            description="Manual scaling factors for jacobian",
            doc="""
            Applies user provided jacobian scaling values:
             - either single value that will be applied to all outputs in jacobian
             - array applied across jacobian
             - dict that specifics output and scaling factor to which apply scaling, (variable_scaling will be applied to non specified outputs)
                e.g. {output_name:scaling_factor} applies to specific jac output 
            """,
        ),
    )
    CONFIG.declare(
        "tolerance",
        ConfigValue(
            default=1e-8,
            domain=float,
            description="Tolerance for reaktoro solver",
            doc="""Tolerance for primary reaktoro solver""",
        ),
    )
    CONFIG.declare(
        "epsilon",
        ConfigValue(
            default=1e-32,
            domain=float,
            description="Epsilon for reaktoro solver",
            doc="""Defines what is considered to be 0 for ion composition""",
        ),
    )
    CONFIG.declare(
        "max_iters",
        ConfigValue(
            default=400,
            domain=int,
            description="Maximum number of iterations for reaktoro solver",
            doc="""The maximum number of iterations for reaktoro solver""",
        ),
    )
    CONFIG.declare(
        "presolve_during_initialization",
        ConfigValue(
            default=False,
            domain=bool,
            description="Option to pre-solve to low tolerance first, before primary solve but only during initialization",
            doc="""In some cases reaktoro might fail to solve to high tolerance first,
            a presolve at low tolerance can enable the reaktoro solve to high tolerance, this will only presolve during initialization""",
        ),
    )
    CONFIG.declare(
        "presolve_property_block",
        ConfigValue(
            default=False,
            domain=bool,
            description="Option to pre-solve to low tolerance first on main property block, before primary solve",
            doc="""In some cases reaktoro might fail to solve to high tolerance first,
            a presolve at low tolerance can enable the reaktoro solve to high tolerance""",
        ),
    )
    CONFIG.declare(
        "presolve_speciation_block",
        ConfigValue(
            default=False,
            domain=bool,
            description="Option to pre-solve to low tolerance for specification block, before primary solve",
            doc="""In some cases reaktoro might fail to solve to high tolerance first,
            a presolve at low tolerance can enable the reaktoro solve to high tolerance""",
        ),
    )
    CONFIG.declare(
        "presolve_max_iters",
        ConfigValue(
            default=400,
            domain=int,
            description="Presolve maximum number of iterations for reaktoro solver",
            doc="""The maximum number of iterations for reaktoro pre-solver""",
        ),
    )
    CONFIG.declare(
        "presolve_tolerance",
        ConfigValue(
            default=1e-4,
            domain=float,
            description="Tolerance for reaktoro pre-solve",
            doc="""Tolerance for reaktoro pre-solver""",
        ),
    )
    CONFIG.declare(
        "presolve_epsilon",
        ConfigValue(
            default=1e-32,
            domain=float,
            description="Tolerance for reaktoro pre-solve",
            doc="""Defines what is considered to be 0 for ion composition""",
        ),
    )
    CONFIG.declare(
        "hessian_type",
        ConfigValue(
            default="Jt.J",
            domain=IsInstance((str, HessTypes)),
            description="Hessian type to use for reaktor gray box",
            doc="""Hessian type to use, some might provide better stability
                options (Jt.J, BFGS, BFGS-mod,BFGS-damp, BFGS-ipopt""",
        ),
    )
    CONFIG.declare(
        "open_species_on_property_block",
        ConfigValue(
            default=None,
            domain=IsInstance((str, list)),
            description="Registers species to open to optimization, this can help with solvability of some problems",
            doc="""Registers species (or list of species) to open to optimization and write empty constraint for,
            this can help with solvability of some problems, but can
            lead to unexpected results depending on database, activity coefficients, and inputs chosen.
            This generally should not be left as None""",
        ),
    )
    CONFIG.declare(
        "open_species_on_speciation_block",
        ConfigValue(
            default=None,
            domain=IsInstance((str, list)),
            description="Registers species to open to optimization, this can help with solvability of some problems",
            doc="""Registers species (or list of species) to open to optimization and write empty constraint for,
            this can help with solvability of some problems, but can
            lead to unexpected results depending on database, activity coefficients, and inputs chosen.
            This generally should not be needed""",
        ),
    )

    def build(self):
        super().build()
        """ configure state"""
        if self.config.build_speciation_block:
            """create spectitation block and then property block"""
            self.speciation_block = Block()
            self.build_rkt_state(self.speciation_block, speciation_block=True)
            self.build_rkt_inputs(self.speciation_block, speciation_block=True)
            self.build_rkt_outputs(self.speciation_block, speciation_block=True)
            self.build_rkt_jacobian(self.speciation_block)
            self.build_rkt_solver(self.speciation_block, speciation_block=True)
            self.build_gray_box(self.speciation_block)

            self.build_rkt_state(
                self,
                speciation_block=False,
                speciation_block_built=True,
                input_composition=self.speciation_block.outputs,
            )
            self.build_rkt_inputs(
                self, speciation_block=False, speciation_block_built=True
            )
            self.build_rkt_outputs(self, speciation_block=False)
            self.build_rkt_jacobian(self)
            self.build_rkt_solver(self, speciation_block=False)
            self.build_gray_box(self)
        else:
            """create property block only"""
            self.build_rkt_state(self)
            self.build_rkt_inputs(self)
            self.build_rkt_outputs(self)
            self.build_rkt_jacobian(self)
            self.build_rkt_solver(self)
            self.build_gray_box(self)

    def build_rkt_state(
        self,
        block,
        speciation_block=False,
        speciation_block_built=False,
        input_composition=None,
    ):
        """this will buld our rkt state specified block.
        The keyword arguments are for automatic configuration of speciation and property blocks

        Keywords:
        block -- pyomo block to build the model on
        speciation_block -- sets to build speciation block which by default does not
            add mineral or gas phases (unless configured to do so), and only outputs speciesAmount from the rktModel
        speciation_block_built -- this should only be True if speciation block was built, in this case the model
            will build user requested outputs, not provide user supplied pH and disable charge neutrality
        input_composition -- if alternate input composition should be used rather one provided by config (used to pass in
            speciesAmount output from speciation_block)

        """

        """ defining which indexes shold be indexed 
         if block is not indexed, index is None 
         if specific input is set to be not Indexed the index will be set to None"""

        index = self.index()
        composition_indexed = index
        temperature_indexed = index
        pressure_indexed = index
        pH_indexed = index

        if self.config.composition_indexed == False:
            composition_indexed = None
        if self.config.temperature_indexed == False:
            temperature_indexed = None
        if self.config.pressure_indexed == False:
            pressure_indexed = None
        if self.config.pH_indexed == False:
            pH_indexed = None
        if input_composition is None:
            input_composition = self.config.composition

        convert_to_rkt_species = self.config.convert_to_rkt_species
        pH = self.config.pH
        if speciation_block == False:
            if speciation_block_built:
                """specitition block already used pH to calcualte
                state and we are getitng back exact species as inputs"""
                pH = None
                convert_to_rkt_species = False
                composition_indexed = None  # passing speciation output directly
        block.rkt_state = ReaktoroState()
        block.rkt_state.register_inputs(
            composition=input_composition,
            temperature=self.config.temperature,
            pressure=self.config.pressure,
            pH=pH,
            convert_to_rkt_species=convert_to_rkt_species,
            species_to_rkt_species_dict=self.config.species_to_rkt_species_dict,
            composition_is_elements=self.config.composition_is_elements,
            composition_index=composition_indexed,
            temperature_index=temperature_indexed,
            pressure_index=pressure_indexed,
            pH_index=pH_indexed,
        )

        """ register phases"""
        if speciation_block == False or self.config.build_speciation_block_with_phases:
            """dont add phases if we are speciating"""
            if self.config.gas_phases is not None:
                block.rkt_state.register_gas_phases(self.config.gas_phases)
            if self.config.mineral_phases is not None:
                block.rkt_state.register_mineral_phases(self.config.mineral_phases)
        if self.config.ion_exchange_phases is not None:
            block.rkt_state.register_ion_exchange_phase(self.config.ion_exchange_phases)
            block.rkt_state.register_ion_exchange_phase(
                self.config.ion_exchange_phase_activity_model
            )
        """ setup activity models - if no phases present they will do nothing """
        block.rkt_state.set_aqueous_phase_activity_model(
            self.config.aqueous_phase_activity_model
        )
        block.rkt_state.set_gas_phase_activity_model(
            self.config.gas_phase_activity_model
        )
        block.rkt_state.set_mineral_phase_activity_model(
            self.config.mineral_phase_activity_model
        )

        """ setup database """
        block.rkt_state.set_database(
            dbtype=self.config.database, database=self.config.database_file
        )
        """ build state """
        block.rkt_state.build_state()

    def build_rkt_inputs(
        self,
        block,
        speciation_block=False,
        speciation_block_built=False,
    ):
        """this will buld our rkt inputs specified block.
        The keyword arguments are for automatic configuration of speciation and property blocks

        Keywords:
        block -- pyomo block to build the model on
        speciation_block -- sets to build speciation block which by default does not
            add mineral or gas phases (unless configured to do so), and only outputs speciesAmount from the rktModel
        speciation_block_built -- this should only be True if speciation block was built, in this case the model
            will build user requested outputs, not provide user supplied pH and disable charge neutrality
        """
        """ get index for chemicals - refer to build_rkt_state on indexing notes"""
        chemical_addition_indexed = self.index()
        if self.config.chemical_addition_indexed == False:
            chemical_addition_indexed = None

        block.rkt_inputs = ReaktoroInputSpec(block.rkt_state)

        """ add chemical only if its not a specitation block (normal mode)"""
        if speciation_block == False:
            if self.config.chemical_addition is not None:
                block.rkt_inputs.register_chemical_additions(
                    self.config.chemical_addition, index=chemical_addition_indexed
                )
            if self.config.chemical_speciation is not None:
                for chemical, speciation in self.config.chemical_speciation.items():
                    if isinstance(chemical, tuple):
                        chemical = chemical[-1]
                    block.rkt_inputs.register_chemical(chemical, speciation)
            block.rkt_inputs.register_open_species(
                self.config.open_species_on_property_block
            )
        else:
            block.rkt_inputs.register_open_species(
                self.config.open_species_on_speciation_block
            )

        """ register aqueous solvent phase"""
        block.rkt_inputs.register_aqueous_solvent(self.config.aqueous_solvent_specie)

        """ register char neutrality"""
        if (
            speciation_block_built == False
            or self.config.assert_charge_neutrality_on_all_blocks
        ):
            """only ensure charge neutrality when doing first calculation"""
            block.rkt_inputs.register_charge_neutrality(
                assert_neutrality=self.config.assert_charge_neutrality,
                ion=self.config.charge_neutrality_ion,
            )
            block.rkt_inputs.configure_specs(
                dissolve_species_in_rkt=self.config.dissolve_species_in_reaktoro,
                exact_speciation=self.config.exact_speciation,
            )
        else:
            """if we have built a speciation block, the feed should be charge neutral and
            exact speciation is provided"""
            block.rkt_inputs.register_charge_neutrality(
                assert_neutrality=False, ion=self.config.charge_neutrality_ion
            )

            block.rkt_inputs.configure_specs(
                dissolve_species_in_rkt=self.config.dissolve_species_in_reaktoro,
                exact_speciation=True,
            )

    def build_rkt_outputs(self, block, speciation_block=False):
        """this will build rkt outputs specified block.
        The keyword arguments are for automatic configuration of speciation and property blocks

        Keywords:
        block -- pyomo block to build the model on
        speciation_block -- sets to build speciation block which by default does not
            add mineral or gas phases (unless configured to do so), and only outputs speciesAmount from the rktModel
        """

        """ configure outputs """
        index = self.index()

        block.rkt_outputs = ReaktoroOutputSpec(block.rkt_state)
        if self.config.outputs is None:
            raise ValueError("Outputs must be provided!")
        if speciation_block:
            """when speciating we only want species amounts as output"""
            block.rkt_outputs.register_output("speciesAmount", get_all_indexes=True)
        else:
            """build user requested outputs"""
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
                        block.rkt_outputs.register_output(
                            output_key, get_all_indexes=output_var
                        )
                    else:
                        block.rkt_outputs.register_output(
                            output_key, output_prop, pyomo_var=output_var
                        )

    def build_rkt_jacobian(self, block):
        """this will build rkt jacboian on specified block.
        The keyword arguments are for automatic configuration of speciation and property blocks

        Keywords:
        block -- pyomo block to build the model on
        """
        """ config outputs """
        block.rkt_jacobian = ReaktoroJacobianSpec(block.rkt_state, block.rkt_outputs)
        block.rkt_jacobian.configure_numerical_jacobian(
            jacobian_type=self.config.numerical_jac_type,
            order=self.config.numerical_jac_order,
            step_size=self.config.numerical_jac_step,
        )

    def build_rkt_solver(self, block, speciation_block=False):
        """this will build rkt outputs specified block.
        The keyword arguments are for automatic configuration of speciation and property blocks

        Keywords:
        block -- pyomo block to build the model on
        speciation_block -- sets to build speciation block which by default does not
            add mineral or gas phases (unless configured to do so), and only outputs speciesAmount from the rktModel
        """

        """config solver """
        name = str(self)
        if speciation_block:
            name = f"{name}_speciation_block"
        block.rktSolver = ReaktoroSolver(
            block.rkt_state,
            block.rkt_inputs,
            block.rkt_outputs,
            block.rkt_jacobian,
            block_name=name,
        )
        if speciation_block:
            presolve = self.config.presolve_speciation_block
        else:
            presolve = self.config.presolve_property_block
        block.rktSolver.set_solver_options(
            tolerance=self.config.tolerance,
            epsilon=self.config.epsilon,
            presolve=presolve,
            presolve_tolerance=self.config.presolve_tolerance,
            presolve_epsilon=self.config.presolve_epsilon,
            max_iters=self.config.max_iters,
            presolve_max_iters=self.config.presolve_max_iters,
            hessian_type=self.config.hessian_type,
        )

    def build_gray_box(self, block):
        """this will build rkt outputs specified block.
        The keyword arguments are for automatic configuration of speciation and property blocks

        Keywords:
        block -- pyomo block to build the model on
        """
        """ build block"""
        scaling = self.config.jacobian_user_scaling
        scaling_type = self.config.jacobian_scaling_type
        block.rkt_block_builder = ReaktoroBlockBuilder(
            block, block.rktSolver, build_on_init=False
        )
        block.rkt_block_builder.configure_jacobian_scaling(
            jacobian_scaling_type=scaling_type, user_scaling=scaling
        )
        block.rkt_block_builder.build_reaktoro_block()

    def display_jacobian_outputs(self):
        if self.config.build_speciation_block:
            _log.info("-----Displaying information for speciation block ------")
            self.speciation_block.rkt_jacobian.display_jacobian_output_types()
        _log.info("-----Displaying information for property block ------")
        self.rkt_jacobian.display_jacobian_output_types()

    def display_jacobian_scaling(self):
        jacobian_scaling = {}
        if self.config.build_speciation_block:
            _log.info("-----Displaying information for speciation block ------")
            jac_scale = (
                self.speciation_block.rkt_block_builder.display_jacobian_scaling()
            )
            jacobian_scaling["speciation_block"] = jac_scale
        _log.info("-----Displaying information for property block ------")
        jac_scale = self.rkt_block_builder.display_jacobian_scaling()
        jacobian_scaling["property_block"] = jac_scale
        return jacobian_scaling

    def set_jacobian_scaling(self, user_scaling_dict, speciation_block=False):
        if speciation_block:
            self.speciation_block.rkt_block_builder.set_user_jacobian_scaling(
                user_scaling_dict
            )
        else:
            self.rkt_block_builder.set_user_jacobian_scaling(user_scaling_dict)

    def initialize(self):
        if (
            self.config.presolve_during_initialization
            or self.config.presolve_speciation_block
            or self.config.presolve_property_block
        ):
            presolve = True
        else:
            presolve = False
        if self.config.build_speciation_block:
            self.speciation_block.rkt_block_builder.initialize(presolve)
        self.rkt_block_builder.initialize(presolve)
