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
    Objective,
    Constraint,
    units as pyunits,
)

from watertap.core.solvers import get_solver
from pyomo.util.calc_var_value import calculate_variable_from_constraint

import idaes.core.util.scaling as iscale


import reaktoro as rkt


# This examples demonstrates how Reaktoro graybox can be used to estimates
# properties in desalination process.

# This example demonstrates how to:
# (1) Setup up basic ReaktoroBlock
# (2) Calculate basic properties for Scaling Tendency, pH, and Osmotic pressure
# (3) Optimize system pH for operation at target Scaling Tendency
# Key assumptions:
# Assumes that process concentrating the feed does not alter the pH.
# This might be a good assumptions for process such as RO, but might be a poor
# assumption for evaporative processes.


def main():
    m = build_simple_desal()
    m_open = build_simple_desal(True)
    initialize(m)
    initialize(m_open)
    setup_optimization(m)
    setup_optimization(m_open)
    print("---result with out extra open species---")
    solve(m)
    print("---result with open extra open species---")
    solve(m_open)
    return m, m_open


def build_simple_desal(open_species=False):
    m = ConcreteModel()
    m.feed_composition = Var(
        ["H2O", "Mg", "Na", "Cl", "SO4", "Ca", "HCO3"],
        initialize=1,
        units=pyunits.mol / pyunits.s,
    )
    m.feed_composition.construct()
    m.feed_composition["H2O"].fix(55)
    m.feed_composition["Mg"].fix(0.01)
    m.feed_composition["Na"].fix(0.025)
    m.feed_composition["Cl"].fix(0.025)
    m.feed_composition["Ca"].fix(0.001)
    m.feed_composition["HCO3"].fix(0.01)
    m.feed_composition["SO4"].fix(0.01)
    m.feed_temperature = Var(initialize=293.15, units=pyunits.K)
    m.feed_temperature.fix()
    m.feed_pressure = Var(initialize=1e5, units=pyunits.Pa)
    m.feed_pressure.fix()
    m.feed_pH = Var(initialize=7, bounds=(4, 12), units=pyunits.dimensionless)
    m.feed_pH.fix()
    m.acid_addition = Var(initialize=0.00001, units=pyunits.mol / pyunits.s)
    m.acid_addition.fix()
    m.desal_composition = Var(
        list(m.feed_composition.keys()),
        initialize=1,
        units=pyunits.mol / pyunits.s,
    )
    m.desal_product_flow = Var(
        initialize=1,
        units=pyunits.mol / pyunits.s,
    )
    m.water_recovery = Var(
        initialize=0.5,
        bounds=(0.1, 0.9),
        units=pyunits.dimensionless,
    )
    m.water_recovery.fix()
    m.desal_osmotic_pressure = Var(initialize=1e5, units=pyunits.Pa)

    # Note how we declare different output variables and assemble them into
    # single dict to pass into reaktoro block, this enables user to
    # mix outputs from different pyomo variables in their models and directly bind them as
    # outputs in reaktoro model
    m.desal_scaling = Var(
        [
            ("scalingTendency", "Calcite"),
            ("scalingTendency", "Gypsum"),
        ],
        initialize=1,
    )
    m.desal_pH = Var(initialize=1)

    m.desal_properties = {}
    for key, obj in m.desal_scaling.items():
        m.desal_properties[key] = obj
    m.desal_properties[("pH", None)] = m.desal_pH
    m.desal_properties[("osmoticPressure", "H2O")] = m.desal_osmotic_pressure

    m.desal_properties[("pH", None)].setlb(4)

    @m.Constraint(list(m.feed_composition.keys()))
    def eq_desal_composition(fs, key):
        if key == "H2O":
            return (
                m.desal_composition["H2O"]
                == m.feed_composition["H2O"] - m.desal_product_flow
            )
        else:
            return m.desal_composition[key] == m.feed_composition[key]

    m.eq_water_recovery = Constraint(
        expr=m.water_recovery == m.desal_product_flow / m.feed_composition["H2O"]
    )
    if open_species:

        # Example for opening additional species to improve
        # reaktoro solver stability,
        # Note how this does not alter reaktoro output results as
        # the H+ and OH- is already constrained by total H amount in property
        # block, as such the DOFs are still zero.

        # However, this can result in incorrect speciation for some databases, so please use with caution.

        species_to_open = ["H+", "OH-"]
    else:
        species_to_open = None
    m.eq_desal_properties = ReaktoroBlock(
        aqueous_phase={
            "composition": m.desal_composition,
            "convert_to_rkt_species": True,
            "activity_model": rkt.ActivityModelPitzer(),
        },
        system_state={
            "temperature": m.feed_temperature,
            "pressure": m.feed_pressure,
            "pH": m.feed_pH,
        },
        outputs=m.desal_properties,
        chemistry_modifier={"HCl": m.acid_addition},
        dissolve_species_in_reaktoro=True,
        # we can use default converter as its defined for default database (Phreeqc and pitzer)
        # we are modifying state and must speciate inputs before adding acid to find final prop state.
        build_speciation_block=True,
        reaktoro_solve_options={"open_species_on_property_block": species_to_open},
    )
    scale_model(m)
    return m


def scale_model(m):
    for key in m.feed_composition:
        iscale.set_scaling_factor(
            m.feed_composition[key], 1 / m.feed_composition[key].value
        )
        iscale.set_scaling_factor(
            m.desal_composition[key], 1 / m.feed_composition[key].value
        )
    iscale.set_scaling_factor(m.water_recovery, 1)
    iscale.set_scaling_factor(m.acid_addition, 1 / 0.001)


def initialize(m):
    calculate_variable_from_constraint(m.desal_product_flow, m.eq_water_recovery)
    for key in m.eq_desal_composition:
        calculate_variable_from_constraint(
            m.desal_composition[key], m.eq_desal_composition[key]
        )
    m.eq_desal_properties.initialize()
    solve(m)


def setup_optimization(m):
    m.objective = Objective(expr=(1 - m.water_recovery) * 10 + m.acid_addition)
    m.desal_properties[("scalingTendency", "Calcite")].setub(1)
    m.desal_properties[("scalingTendency", "Gypsum")].setub(1)
    m.water_recovery.unfix()
    m.acid_addition.unfix()


def display_results(m):
    print("result")
    print(f"""Osmotic pressure {m.desal_osmotic_pressure.value} (Pa)""")
    for key, obj in m.desal_scaling.items():
        print(f"{key}, {obj.value}")
    print(f"feed pH {m.feed_pH.value}, desal pH {m.desal_pH.value}")


def solve(m):
    cy_solver = get_solver(solver="cyipopt-watertap")
    # cy_solver.options["max_iter"] = 200
    # only enable if avaialbe !
    # cy_solver.options["linear_solver"] = "ma27"
    result = cy_solver.solve(m, tee=True)
    display_results(m)
    return result


if __name__ == "__main__":
    main()
