#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/reaktoro-pse/"
#################################################################################

from pyomo.common.config import ConfigValue, IsInstance, ConfigDict, In


class ReaktoroSolverOptions:
    def __init__(self):
        pass

    def get_dict(self, presolve_options=False, advanced_options=False):
        CONFIG = ConfigDict()

        CONFIG.declare(
            "solver_tolerance",
            ConfigValue(
                default=1e-9,
                domain=float,
                description="Tolerance for Reaktoro solver",
                doc="""Tolerance for primary Reaktoro solver""",
            ),
        )
        CONFIG.declare(
            "auto_scale_constraints",
            ConfigValue(
                default=True,
                domain=bool,
                description="Option to automatically scale Reaktoro constraints.",
                doc="""If True, Reaktoro constraints will be automatically scaled based on their values.
                These constraints are written in Reaktoro solver it self, and not pyomo, as such they only impact the
                Reaktoro solver and precision. This should ensure that all elements are solved to same tolerance.""",
            ),
        )
        CONFIG.declare(
            "epsilon",
            ConfigValue(
                default=1e-40,
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
            CONFIG.declare(
                "species_amount_zero_offset",
                ConfigValue(
                    default=1e-4,
                    domain=float,
                    description="Scaled offset added to speciesAmount outputs",
                    doc="""Offset added to user provided speciesAmount output variables, in scaled units
                    (output_var = reaktoro_amount + species_amount_zero_offset / scaling_factor of output_var).
                    When a species or phase is not present Reaktoro returns its amount as epsilon (e.g. 1e-40)
                    with zero derivatives. If the output variable (or a variable linked to it through
                    equality constraints) has a lower bound of 0, it then sits at epsilon from the bound, which
                    makes the problem degenerate and can prevent Ipopt from converging dual infeasibility.
                    The offset keeps these variables away from the bound, at the cost of a mass balance error of
                    offset/scaling_factor. Relies on user provided scaling factors on output variables, without
                    them the scaling factor is computed from the output value and the offset has no effect for
                    species that are not present. Set to 0 to disable. Species amounts passed from a
                    speciation block to the property block are never offset.""",
                ),
            )
            CONFIG.declare(
                "species_amount_offset_type",
                ConfigValue(
                    default="smooth_floor",
                    domain=In(["additive", "smooth_floor", None]),
                    description="How species_amount_zero_offset is applied",
                    doc="""With a = species_amount_zero_offset / scaling_factor of output_var:
                    - additive: output_var = reaktoro_amount + a, every species is shifted by a
                    - smooth_floor: output_var = sqrt(reaktoro_amount**2 + a**2), species that are not
                    present sit at a, present species (amount n >> a) are shifted by only ~a**2/(2n).
                    This allows a larger offset (e.g. 1e-4) that keeps absent phases further from
                    their bounds (bound multipliers scale as barrier parameter / a) with a small error
                    for species that are present.
                    - None: no offset is applied (same as species_amount_zero_offset = 0)
                    The offset is applied to the gray box outputs (values and jacobian), output
                    constraints stay linear (output_var == gray box output). a is computed from the
                    current scaling factor of output_var on every gray box evaluation, so scaling
                    changes after initialization take effect on the next solve.""",
                ),
            )
        return CONFIG
