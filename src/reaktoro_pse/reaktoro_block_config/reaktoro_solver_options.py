from pyomo.common.config import ConfigValue, IsInstance, ConfigDict
from reaktoro_pse.core.reaktoro_jacobian import ReaktoroJacobianSpec, JacType

from reaktoro_pse.core.reaktoro_gray_box import HessTypes


class ReaktoroSolverOptions:
    def __init__(self):
        pass

    def get_dict(self, presolve_options=False, advanced_options=False):
        CONFIG = ConfigDict()

        CONFIG.declare(
            "solver_tolerance",
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
                description="epsilon for reaktoro solver",
                doc="""Defines what is considered to be 0 for ion composition""",
            ),
        )
        CONFIG.declare(
            "max_iterations",
            ConfigValue(
                default=400,
                domain=int,
                description="Maximum number of iterations for reaktoro solver",
                doc="""The maximum number of iterations for reaktoro solver""",
            ),
        )
        if presolve_options:
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
                    description="Option to pre-solve to low tolerance first on main property block, before primary solve",
                    doc="""In some cases reaktoro might fail to solve to high tolerance first,
                    a presolve at low tolerance can enable the reaktoro solve to high tolerance""",
                ),
            )
        if advanced_options:

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
        return CONFIG
