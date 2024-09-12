from ast import Is
from pyomo.common.config import ConfigValue, IsInstance, ConfigDict

from pyomo.core.base.var import IndexedVar, Var, VarData


class PhaseInput:
    def __init__(self):
        pass

    def get_dict(
        self, include_pH=False, aqueous_phase=False, include_solvent_species=False
    ):
        phase_input = ConfigDict()
        phase_input.declare(
            "composition",
            ConfigValue(
                default=dict,
                domain=IsInstance((dict, IndexedVar)),
                description="Input composition to Reaktoro block",
                doc="An input dictionary, or IndexedVar that contains species and their amounts",
            ),
        )
        phase_input.declare(
            "composition_indexed",
            ConfigValue(
                default=True,
                domain=bool,
                description="composition is indexed",
                doc="""Option that defines how to treat input variable when building indexed ReaktoroBlock":
                    - If true, the input has same indexing as block, and each indexed input will be passed into respective indexed ReaktoroBlock
                    - If false, all indexed blocks will get same input""",
            ),
        )
        phase_input.declare(
            "phase_components",
            ConfigValue(
                default=None,
                description="List or str of phase components",
                doc="""
            Phases supported by selected data base
            Accepts:
                component name for phase
                list of phases or elements making up the phases
                single or list of initialized Reaktoro phase object (e.g. Reaktoro.MineralPhase('Calcite'))
            """,
            ),
        )
        phase_input.declare(
            "convert_to_rkt_species",
            ConfigValue(
                default=False,
                domain=bool,
                description="Defines if provided species should be converted to RKT notation",
                doc="Enable conversion provided species names to Reaktoro names - (currently supports PhreeqC database)",
            ),
        )
        phase_input.declare(
            "species_to_rkt_species_dict",
            ConfigValue(
                default="default",
                domain=IsInstance((dict, str)),
                description="Dictionary for translating user supplied species to RKT species specific to selected database",
                doc="""
                    Dictionary that connects user species to Reaktoro data base species (please reference chosen data base in question):
                    Dictionary should have following structure {user_species: reaktoro_database_specie}. 
                    e.g.
                    {'Na':'Na+',
                    'Ca':'Ca+2'}
                    """,
            ),
        )
        phase_input.declare(
            "composition_is_elements",
            ConfigValue(
                default=False,
                domain=bool,
                description="Defines if provided composition is elements and not species",
                doc="Defines if provided composition is elements and not species",
            ),
        )
        phase_input.declare(
            "activity_model",
            ConfigValue(
                default=None,
                description="Activity model for phase",
                doc="""
            Defines which activity model to use for aqueous phase in Reaktoro
            Accepts:
            String name of activity model
            Initialized Reaktoro.ActivityModel (e.g. reaktoro.ActivityModelIdealAqueous())
            Chain of Reaktoro activity models (reaktoro.chain(reaktoro.ActivityModelA, reaktoro.ActivityModelB))
            List of activity models strings or reaktoro.ActivityModel initialized objects.
            """,
            ),
        )
        if include_solvent_species:
            if aqueous_phase:
                default = "H2O"
            else:
                default = None
            phase_input.declare(
                "fixed_solvent_specie",
                ConfigValue(
                    default=default,
                    domain=str,
                    description="Defines solvent specie to fix when speciating system",
                    doc="""
                    When speciating, the exact amount of all elements is rarely known, as such fixing exact amount of solvent 
                    provides a simpler and more stable alterative.
                    For aqueous systems its amount of H2O, for organic system it would be primary solvent specie. If omitted the equilibrium will
                    be found assuming all species provide appropriate mass balance.
                    Providing this option will open all elements in solvent to optimization (e.g. if H2O is provided amount of O nad H will 
                    be open in equilibrium calculations)
                    if enabled
                    """,
                ),
            ),
            phase_input.declare(
                "free_element",
                ConfigValue(
                    default=None,
                    domain=IsInstance((str, list, tuple)),
                    description="Defines free element to unfix during speciation",
                    doc="""
                    When speciating exact amount all elements might not be known, as such freeing one to be found 
                    might be required - for example in Aqueous systems the total amount of oxygen might be known due to measurement of 
                    all oxygen containing species, but amount of H might not be due to required pH balance. Although in general even amount of 
                    Oxygen is unknown due to various non-measured species in solution contain oxygen (for example NaOH) as such its 
                    if enabled
                    """,
                ),
            )
        return phase_input


class SystemInput:
    def __init__(self):
        pass

    def get_dict(self):
        system_input = ConfigDict()
        system_input.declare(
            "temperature",
            ConfigValue(
                default=None,
                domain=IsInstance((VarData, Var, dict, IndexedVar)),
                description="Input temperature for reaktoro block",
                doc="Var or IndexedVar that references system temperature",
            ),
        )
        system_input.declare(
            "temperature_indexed",
            ConfigValue(
                default=True,
                domain=bool,
                description="Temperature is indexed",
                doc="""Option that defines how to treat input variable when building indexed reaktoroBlock":
                    - If true, the input has same indexing as block, and each indexed input will be passed into respective indexed reaktoroBlock
                    - If false, all indexed blocks will get same input""",
            ),
        )
        system_input.declare(
            "pressure",
            ConfigValue(
                default=None,
                domain=IsInstance((VarData, Var, dict, IndexedVar)),
                description="Input pressure for reaktoro block",
                doc="Var or IndexedVar that references system pressure",
            ),
        )
        system_input.declare(
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
        system_input.declare(
            "enthalpy",
            ConfigValue(
                default=None,
                domain=IsInstance((VarData, Var, dict, IndexedVar)),
                description="Input enthalpy for reaktoro block",
                doc="Var or IndexedVar that references system enthalpy",
            ),
        )
        system_input.declare(
            "enthalpy_indexed",
            ConfigValue(
                default=True,
                domain=bool,
                description="enthalpy is indexed",
                doc="""Option that defines how to treat input variable when building indexed reaktoroBlock":
                    - If true, the input has same indexing as block, and each indexed input will be passed into respective indexed reaktoroBlock
                    - If false, all indexed blocks will get same input""",
            ),
        )
        system_input.declare(
            "pH",
            ConfigValue(
                default=None,
                domain=IsInstance((VarData, Var, dict, IndexedVar)),
                description="Input pH for reaktoro block",
                doc="Var or IndexedVar that references system pH",
            ),
        )
        system_input.declare(
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
        system_input.declare(
            "temperature_bounds",
            ConfigValue(
                default=(10, 10000),
                domain=tuple,
                description="Bounds for temperature, only needed if we are temperature is being solved for ",
                doc="""Defines bounds used by reaktoro solver for temperature""",
            ),
        )
        system_input.declare(
            "pressure_bounds",
            ConfigValue(
                default=(1, 10000),
                domain=tuple,
                description="Bounds for pressure, only needed if we are setting pressure to be unknown",
                doc="""Defines bounds used by reaktoro solver for pressure""",
            ),
        )
        return system_input
