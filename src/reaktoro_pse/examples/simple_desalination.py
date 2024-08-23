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
from pyomo.common.modeling import unique_component_name

import pyomo.environ as pyo
import reaktoro as rkt


def main():
    m = build_simple_desal()
    initialize(m)
    setup_optimization(m)
    solve(m)
    m.display()
    return m


def build_simple_desal():
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
    # pressure.construct()
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

    """ Note how we declare different output variables and assemble them into
    single dict to pass into reaktoro block, this enables user to 
    mix outputs from different pyomo varabiles in thier models and directly bind them as 
    outputs in reaktoro model"""
    m.desal_scaling = Var(
        [
            ("scalingTendency", "Calcite"),
            ("scalingTendency", "Gypsum"),
        ],
        initialize=1,
    )
    m.desal_ph = Var(initialize=1)
    m.osmotic_pressure = Var(initialize=1)

    m.desal_properties = {}
    for key, obj in m.desal_scaling.items():
        m.desal_properties[key] = obj
    m.desal_properties[("pH", None)] = m.desal_ph
    m.desal_properties[("osmoticPressure", "H2O")] = m.osmotic_pressure

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

    m.eq_desal_properties = ReaktoroBlock(
        composition=m.desal_composition,
        temperature=m.feed_temperature,
        pressure=m.feed_pressure,
        pH=m.feed_pH,
        outputs=m.desal_properties,
        chemical_addition={"HCl": m.acid_addition},
        aqueous_phase_activity_model=rkt.ActivityModelPitzer(),
        dissolve_species_in_reaktoro=True,
        # we can use default converter as its defined for default database
        convert_to_rkt_species=True,
        # we are modifying state so lets speciate inputs before adding acid to find final prop state.
        build_speciation_block=True,
        open_species_on_property_block=["H+"],
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
    m.display()


def setup_optimization(m):
    m.objective = Objective(expr=(1 - m.water_recovery) * 10 + m.acid_addition)
    m.desal_properties[("scalingTendency", "Calcite")].setub(1)
    m.desal_properties[("scalingTendency", "Gypsum")].setub(1)
    m.water_recovery.unfix()
    m.acid_addition.unfix()


# def solve(m):
#     # add dummy objective, if needed
#     _dummy_objective = None
#     n_obj = 0
#     for c in m.component_data_objects(pyo.Objective, active=True):
#         n_obj += 1
#     # Add an objective if there isn't one
#     if n_obj == 0:
#         _dummy_objective = pyo.Objective(expr=0)
#         name = unique_component_name(m, "objective")
#         m.add_component(name, _dummy_objective)
#     cy_solver = get_solver(solver="cyipopt")
#     # cy_solver.options["max_iter"] = 200
#     # only enable if avaialbe !
#     # cy_solver.options["linear_solver"] = "ma27"
#     try:
#         result = cy_solver.solve(m, tee=True)
#     finally:
#         # delete the dummy objective
#         if _dummy_objective is not None:
#             m.del_component(_dummy_objective)
#     return result


def solve(m):
    cy_solver = get_solver(solver="cyipopt-watertap")
    # cy_solver.options["max_iter"] = 200
    # only enable if avaialbe !
    # cy_solver.options["linear_solver"] = "ma27"
    result = cy_solver.solve(m, tee=True)
    return result


if __name__ == "__main__":
    main()
