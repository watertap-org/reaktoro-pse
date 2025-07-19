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

from pyomo.common.config import ConfigValue, IsInstance, ConfigDict


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
                description="Tolerance for Reaktoro solver",
                doc="""Tolerance for primary Reaktoro solver""",
            ),
        )
        CONFIG.declare(
            "epsilon",
            ConfigValue(
                default=1e-32,
                domain=float,
                description="epsilon for Reaktoro solver",
                doc="""Defines what is considered to be 0 for ion composition""",
            ),
        )
        CONFIG.declare(
            "max_iterations",
            ConfigValue(
                default=400,
                domain=int,
                description="Maximum number of iterations for Reaktoro solver",
                doc="""The maximum number of iterations for Reaktoro solver""",
            ),
        )
        CONFIG.declare(
            "max_reaktoro_failed_solves",
            ConfigValue(
                default=2,
                domain=int,
                description="Number of attempts to re-solve Reaktoro block when running inside a solver",
                doc="""Defines number of tries Reaktoro block can fail to solve when running inside solver. When Reaktoro fails
                it will raise CyIpoptEvaluationError, this defines how many sequential raises are allowed before terminating.""",
            ),
        )
        if presolve_options:
            CONFIG.declare(
                "presolve_during_initialization",
                ConfigValue(
                    default=False,
                    domain=bool,
                    description="Option to pre-solve to low tolerance first, before primary solve but only during initialization",
                    doc="""In some cases Reaktoro might fail to solve to high tolerance first,
                a presolve at low tolerance can enable the Reaktoro solve to high tolerance, this will only presolve during initialization""",
                ),
            )
            CONFIG.declare(
                "presolve_property_block",
                ConfigValue(
                    default=False,
                    domain=bool,
                    description="Option to pre-solve to low tolerance first on main property block, before primary solve",
                    doc="""In some cases Reaktoro might fail to solve to high tolerance first,
                    a presolve at low tolerance can enable the Reaktoro solve to high tolerance""",
                ),
            )
            CONFIG.declare(
                "presolve_speciation_block",
                ConfigValue(
                    default=False,
                    domain=bool,
                    description="Option to pre-solve to low tolerance first on main property block, before primary solve",
                    doc="""In some cases Reaktoro might fail to solve to high tolerance first,
                    a presolve at low tolerance can enable the Reaktoro solve to high tolerance""",
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
