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
from pyomo.common.modeling import unique_component_name

import pyomo.environ as pyo
import reaktoro as rkt

from reaktoro_pse.reaktoro_block_config import jacobian_options

"""
This examples demonstrates how Reaktoro graybox can be used to estimates combustion perofrming an optimization 
on following reaktoro example: https://reaktoro.org/applications/biomass-gasification/biomass-gasification.html

This example demonstrates how to:
(1) Setup up basic ReaktoroBlock
(2) Calculate exhaust temperature from combustion
(3) Find air to fuel ratio needed to get to 2000 K

 
"""


def main():
    m = build_biogas()
    initialize(m)
    m.air_to_fuel_ratio.unfix()
    m.exhaust_temperature.fix(2000)
    solve(m)
    return m


def build_biogas(open_species=False):
    m = ConcreteModel()
    m.fuel = Var(
        ["Fuel"],
        initialize=1,
        units=pyunits.mol / pyunits.s,
    )
    m.fuel["Fuel"].fix(1)

    m.air = Var(
        ["O2", "N2"],
        initialize=1,
        units=pyunits.mol / pyunits.s,
    )
    m.air["O2"].value = 0.025
    m.air["N2"].value = 0.075

    m.feed_pressure = Var(initialize=1e5, units=pyunits.Pa)
    m.feed_pressure.fix()

    m.air_to_fuel_ratio = Var(initialize=3, units=pyunits.dimensionless)
    m.air_to_fuel_ratio.fix()

    m.exhaust_temperature = Var(initialize=293.15, units=pyunits.K)
    m.mass_gas = Var(initialize=0, units=pyunits.mol / pyunits.s)

    m.combusted_ratio_g_to_c = Var(initialize=1, units=pyunits.dimensionless)
    m.heat_duty = Var(initialize=-455, units=pyunits.kJ / pyunits.s)

    m.outputs = {
        ("temperature", None): m.exhaust_temperature,
        ("amount", "GaseousPhase"): m.mass_gas,
    }
    # Add fuel to nasa-cee database
    db = rkt.NasaDatabase("nasa-cea")

    massC = 49.30  # g/mol
    massH = 5.5  # g/mol
    massO = 45.2  # g/mol
    nC = massC / 12.011
    nH = massH / 1.00797
    nO = massO / 15.9994
    a = nH / nC
    b = nO / nC
    formula = rkt.ChemicalFormula(f"CH{a}O{b}")
    Mfuel = formula.molarMass() * 1000
    Mair = 28.850334

    factor = Mfuel / Mair
    HHV = 18.933  # kJ/g
    m.HHV = Var(initialize=HHV, units=pyunits.kJ / pyunits.g)
    m.HHV.fix()
    h0CO2 = -393.522  # kJ/g
    h0H2O = -285.83  # kJ/g

    # add fuel to db
    stmodelparams = rkt.StandardThermoModelParamsConstant()
    stmodelparams.G0 = 1.0e3
    stmodel = rkt.StandardThermoModelConstant(stmodelparams)
    species = rkt.Species()
    species = species.withName("Fuel")
    species = species.withElements(formula.elements())
    species = species.withAggregateState(rkt.AggregateState.CondensedPhase)
    species = species.withStandardThermoModel(stmodel)
    db.addSpecies(species)

    m.oxygen_flow_eq = Constraint(
        expr=m.air["O2"] == m.air_to_fuel_ratio * factor * 0.21
    )
    m.nitrogen_flow_eq = Constraint(
        expr=m.air["N2"] == m.air_to_fuel_ratio * factor * 0.79
    )

    m.heat_duty_eq = Constraint(
        expr=m.heat_duty == m.HHV * Mfuel + h0CO2 + 0.5 * h0H2O * a
    )

    m.eq_combustion = ReaktoroBlock(
        condensed_phase={
            "composition": m.fuel,
            "convert_to_rkt_species": False,
            "phase_components": ["Fuel", "H2O(l)", "C(gr)"],
        },
        gas_phase={
            "composition": m.air,
            "convert_to_rkt_species": False,
            "phase_components": rkt.GaseousPhase(
                rkt.speciate("C H O N")
            ),  # Expects all possible gasses
        },
        system_state={
            "pressure": m.feed_pressure,
            "enthalpy": m.heat_duty,
            "temperature_bounds": (750, 3000),
        },
        outputs=m.outputs,
        dissolve_species_in_reaktoro=True,
        assert_charge_neutrality=False,
        database=db,
        build_speciation_block=False,
        reaktoro_solve_options={"solver_tolerance": 1e-10},
        # exact_speciation=False,
        # jacobian_options={"user_scaling": {("temperature", None): 1000}},
    )

    m.eq_gtoc = Constraint(
        expr=m.combusted_ratio_g_to_c
        == m.mass_gas / (m.fuel["Fuel"] + m.air["O2"] + m.air["N2"])
    )
    """ input constraint for summing  to elements """
    for key, d in m.eq_combustion.rkt_inputs.constraint_dict.items():
        print(key, d)
    scale_model(m)
    return m


def scale_model(m):
    calculate_variable_from_constraint(m.air["O2"], m.oxygen_flow_eq)
    calculate_variable_from_constraint(m.air["N2"], m.nitrogen_flow_eq)
    calculate_variable_from_constraint(m.heat_duty, m.heat_duty_eq)
    for key in m.fuel:
        iscale.set_scaling_factor(m.fuel[key], 1 / m.fuel[key].value)
    for key in m.air:
        iscale.set_scaling_factor(m.air[key], 1 / 10)
    iscale.set_scaling_factor(m.air_to_fuel_ratio, 1 / 10)
    iscale.set_scaling_factor(m.mass_gas, 1)
    iscale.set_scaling_factor(m.exhaust_temperature, 1 / 1000)
    iscale.set_scaling_factor(m.feed_pressure, 1 / 10000)
    iscale.constraint_scaling_transform(m.eq_gtoc, 1)
    iscale.constraint_scaling_transform(
        m.oxygen_flow_eq, 1 / 10
    )  # / (m.air["O2"].value))
    iscale.constraint_scaling_transform(
        m.nitrogen_flow_eq, 1 / 10
    )  #  / (m.air["N2"].value))


def initialize(m):
    m.eq_combustion.initialize()
    m.eq_combustion.display_jacobian_scaling()
    m.eq_combustion.display_reaktoro_state()
    solve(m)


def display_results(m):
    print("result")
    print(f"""Gas temperature {m.exhaust_temperature.value} (K)""")
    print(f"ratio of gas to total mass {m.combusted_ratio_g_to_c.value}")
    print(f"gas mols {m.mass_gas.value}")
    print(f"air to fuel ratio {m.air_to_fuel_ratio.value}")
    print(f"flow of air O2 {m.air['O2'].value}, N2 {m.air['N2'].value}")
    print(f"flow of fuel {m.fuel['Fuel'].value}")


def solve(m):
    cy_solver = get_solver(solver="cyipopt-watertap")
    cy_solver.options["max_iter"] = 20
    # only enable if avaialbe !
    # cy_solver.options["linear_solver"] = "ma27"
    result = cy_solver.solve(m, tee=True)
    display_results(m)
    return result


if __name__ == "__main__":
    main()
