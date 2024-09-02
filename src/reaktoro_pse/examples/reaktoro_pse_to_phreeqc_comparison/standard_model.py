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
from reaktoro_pse.reaktoro_block import ReaktoroBlock

from pyomo.environ import (
    ConcreteModel,
    Var,
    Constraint,
    units as pyunits,
)
from watertap.core.solvers import get_solver
from pyomo.util.calc_var_value import calculate_variable_from_constraint
import idaes.core.util.scaling as iscale

import reaktoro as rkt

__author__ = "Alexander Dudchenko"


def build_modification_example(water_comp):
    m = ConcreteModel()
    m.feed_composition = Var(
        water_comp.keys(),
        initialize=1,
        units=pyunits.mol / pyunits.s,
    )
    for key, value in water_comp.items():
        m.feed_composition[key].fix(value)
    m.feed_temperature = Var(initialize=273.15 + 20, units=pyunits.K)
    m.feed_temperature.fix()
    m.feed_pressure = Var(
        initialize=10e5, units=pyunits.Pa
    )  # 10 bar used in phreeqc feed pressure
    m.feed_pressure.fix()
    # pressure.construct()
    m.feed_pH = Var(initialize=7.8, bounds=(4, 12), units=pyunits.dimensionless)
    m.feed_pH.fix(7.8)  # feed pH used in phreeqc sim
    m.acid_addition = Var(initialize=0.0, units=pyunits.mol / pyunits.s)
    m.acid_addition.fix()
    m.base_addition = Var(initialize=0.0, units=pyunits.mol / pyunits.s)
    m.base_addition.fix()
    m.lime_addition = Var(initialize=0.0, units=pyunits.mol / pyunits.s)
    m.lime_addition.fix()
    m.modified_properties_water_removal = Var(
        initialize=0,
        units=pyunits.mol / pyunits.s,
    )
    m.water_recovery = Var(
        initialize=0.0,
        bounds=(0.0, 0.9),
        units=pyunits.dimensionless,
    )
    m.water_recovery.fix()
    m.eq_water_fow = Constraint(
        expr=m.water_recovery
        == -m.modified_properties_water_removal / m.feed_composition["H2O"]
    )
    return m


def add_standard_properties(m):
    m.modified_properties = Var(
        [
            ("scalingTendency", "Calcite"),
            ("scalingTendency", "Gypsum"),
            ("pH", None),
            ("osmoticPressure", "H2O"),
        ],
        initialize=1,
    )
    m.eq_modified_properties = ReaktoroBlock(
        aqueous_phase={
            "composition": m.feed_composition,
            "convert_to_rkt_species": True,
            "activity_model": rkt.ActivityModelPitzer(),
        },
        system_state={
            "temperature": m.feed_temperature,
            "pressure": m.feed_pressure,
            "pH": m.feed_pH,
        },
        outputs=m.modified_properties,
        chemistry_modifier={
            "HCl": m.acid_addition,
            "H2O_evaporation": m.modified_properties_water_removal,
            "NaOH": m.base_addition,
        },
        dissolve_species_in_reaktoro=False,
        # we can use default converter as its defined for default database (Phreeqc and pitzer)
        # we are modifying state and must speciate inputs before adding acid to find final prop state.
        build_speciation_block=True,
        reaktoro_solve_options={"open_species_on_property_block": ["H+", "OH-"]},
        jacobian_options={
            "user_scaling": {
                ("saturationIndex", "Calcite"): 1,
                ("saturationIndex", "Gypsum"): 1,
                ("pH", None): 1,
                ("speciesActivityLn", "H2O"): 1,
            },
        },
    )
    scale_model(m)


def scale_model(m):
    for key in m.feed_composition:
        iscale.set_scaling_factor(
            m.feed_composition[key], 1 / m.feed_composition[key].value
        )
    iscale.set_scaling_factor(m.water_recovery, 1)
    iscale.set_scaling_factor(m.acid_addition, 1 / 0.001)
    iscale.set_scaling_factor(m.base_addition, 1 / 0.001)


def initialize(m):
    calculate_variable_from_constraint(
        m.modified_properties_water_removal, m.eq_water_fow
    )
    m.eq_modified_properties.initialize()
    solve(m)


def solve(m):
    cy_solver = get_solver(solver="cyipopt-watertap")
    result = cy_solver.solve(m, tee=True)
    return result
