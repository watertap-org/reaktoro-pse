from pyomo.common.config import ConfigValue, IsInstance, ConfigDict

from pyomo.core.base.var import IndexedVar, Var, VarData


class PhaseInput:
    def __init__(self):
        pass

    def get_dict(self, include_pH=False, include_aqueous_solvent_species=False):
        phase_input = ConfigDict()
        phase_input.declare(
            "composition",
            ConfigValue(
                default=dict,
                domain=IsInstance((dict, IndexedVar)),
                description="Input composition to reaktoro block",
                doc="An input dictionary, or IndexedVar that contains species and their amounts",
            ),
        )
        phase_input.declare(
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
                single or list of initialized reaktoro phase object (e.g. reaktoro.MineralPhase('Calcite'))
            """,
            ),
        )
        phase_input.declare(
            "convert_to_rkt_species",
            ConfigValue(
                default=False,
                domain=bool,
                description="Defines if provided species should be converted to RKT notation",
                doc="Enable conversion provided species names to reaktoro names - (currently supports PhreeqC database)",
            ),
        )
        phase_input.declare(
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
            Defines which activity model to use for aqueous phase in reaktoro
            Accepts:
            String name of activity model
            Initialized reaktoro.ActivityModel (e.g. reaktoro.ActivityModelIdealAqueous())
            Chain of reaktoro activity models (reaktoro.chain(reaktoro.ActivityModelA, reaktoro.ActivityModelB))
            List of activity models strings or reaktoro.ActivityModel initialized objects.
            """,
            ),
        )
        if include_pH:
            phase_input.declare(
                "pH",
                ConfigValue(
                    default=None,
                    domain=IsInstance((VarData, Var, dict, IndexedVar)),
                    description="Input pH for reaktoro block",
                    doc="Var or IndexedVar that references system pH",
                ),
            )
            phase_input.declare(
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
        if include_aqueous_solvent_species:
            phase_input.declare(
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
                    - If true, the input has same indexing as block, and each indexed input willbe passed into respective indexed reaktoroBlock
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
        return system_input
