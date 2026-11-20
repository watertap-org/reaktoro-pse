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
from reaktoro_pse.reaktoro_block import ReaktoroBlock
from reaktoro_pse.core.util_classes.cyipopt_solver import (
    get_cyipopt_watertap_solver,
)
from pyomo.environ import (
    ConcreteModel,
    Var,
    assert_optimal_termination,
    value,
    units as pyunits,
)
from pyomo.util.calc_var_value import calculate_variable_from_constraint
import idaes.core.util.scaling as iscale

# This example shows how to track dissolved oxygen (DO) through a treatment train and
# how to equilibrate water with air in a process step.
#
# Feed composition is given as total amounts of each element (as is typical for lab data),
# with H and O left open, so one more specification is needed to define the redox state.
# Here it is the measured dissolved O2 (system_state "dissolved_O2"), an alternative to pE.
#
# The feed is nearly anoxic groundwater with dissolved Mn(II), softened with lime:
#   - Sealed softening: the only oxidant for Mn(II) is the small amount of incoming DO, so
#     only a small fraction of Mn precipitates as Birnessite (MnO2), and the DO is consumed.
#   - Softening open to air (system_state_modifier "O2_fugacity" = 0.21 bar): O2 is taken
#     up until the water is saturated, and Mn is fully oxidized and removed.
#
# After aerated softening, the water is neutralized with HCl. The redox state leaving
# softening is passed into the neutralization step the same way as pH, by connecting
# DO leaving softening (speciesAmount of O2) to the dissolved_O2 input of the next step.
#
# Note on choosing the redox state to pass between steps:
#   - DO is a good choice for oxic water: it is conserved, and changes with aeration.
#   - When oxidants are consumed (as in sealed softening, DO ~1e-32 mol/s), DO near 0 is a
#     poor specification (pE depends on log(DO)); pass pE instead.
#   - Water with nothing left to hold its redox state (e.g. Mn(II) after all oxidants and
#     Mn(IV) solids are gone) is poorly poised: its pE is set by trace species, and adding
#     chemicals to it can be numerically fragile regardless of which redox state is passed.
#     This is why the sealed case stops after softening.
#
# Note on O2_fugacity (bar) as a model of aeration:
#   Fixing O2 fugacity models equilibrium with an unlimited air reservoir (open basin),
#   DO then follows Henry's law, DO ~ K_H(T, P, salinity) * f_O2.
#   - Temperature: K_H is computed by Reaktoro at the block temperature from the database.
#   - pH: does not change O2 saturation (O2 does not hydrolyze), but it does change the pE
#     implied by air (pE ~ 20.75 - pH + 0.25 * log10(f_O2)) and how strongly O2 oxidizes
#     Mn and Fe.
#   - Pressure: f_O2 is the O2 partial pressure of the gas in contact with water, not the
#     liquid pressure, so it does not follow the block pressure. For air,
#     f_O2 ~ 0.2095 * (P_air - p_H2O(T)), where P_air drops with elevation
#     (~ 1 atm * exp(-elevation / 8400 m)) and p_H2O(T) is the water vapor pressure
#     (gas at the water surface is saturated with vapor).
#   - Partial aeration: DO is linear in f_O2, so an aerator reaching e.g. 80% of
#     saturation can be modeled as 0.8 * f_O2 of air.
#   - Salinity: salting out of O2 is captured only approximately by the activity model
#     (pitzer.dat does not contain O2).
#
#   DO at saturation (mg/L) with phreeqc.dat or wateq4f.dat and ActivityModelPhreeqc,
#   at 1 atm, compared to standard tables for fresh water in moist air
#   (APHA / Benson and Krause, approximate values):
#
#   T (C)  | f_O2 = 0.2123 bar | f_O2 corrected for | standard
#          | (dry air)         | water vapor        | table
#   -------+-------------------+--------------------+---------
#   0      | 14.91             | 14.82              | 14.6
#   10     | 11.65             | 11.51              | 11.3
#   20     | 9.50              | 9.28               | 9.09
#   30     | 8.04              | 7.70               | 7.56
#   40     | 7.04              | 6.53               | 6.41
#   60     | 5.86              | 4.71               | ~4.7
#   20, 0.6 M NaCl (~seawater) | 8.27   | 8.08      | ~7.4
#
#   With the water vapor correction, Reaktoro is within ~2% of tabulated values from 0 to
#   60 C. Fixed f_O2 of 0.21 overestimates DO in warm water (by ~7% at 40 C, ~20% at
#   60 C), and the activity model underestimates salting out (~9% high DO at seawater
#   salinity, more for brines).

feed = {
    "H2O": 55.5,  # mol/s, ~1 kg/s
    "Ca": 2e-3,
    "Mg": 1e-3,
    "Na": 3e-3,
    "Cl": 3e-3,
    "HCO3": 5e-3,
    "SO4": 1e-3,
    "Mn": 1e-5,  # ~0.55 mg/L
}
feed_DO = 3e-7  # mol/s, ~0.01 mg/L, water is nearly anoxic
air_O2_fugacity = 0.21  # bar
minerals = ["Calcite", "Birnessite"]
O2_mw = 32e3  # mg/mol, water flow is ~1 L/s so mol/s * mw gives ~mg/L


def main():
    results = {}
    for aerate in [False, True]:
        m = build_model(aerate=aerate)
        initialize(m)
        solve(m)
        results[aerate] = get_results(m)
        display_results(results[aerate], aerate)
    return results


def build_model(aerate=True):
    m = ConcreteModel()
    m.aerate = aerate
    m.feed_composition = Var(
        list(feed.keys()), initialize=feed, units=pyunits.mol / pyunits.s
    )
    m.feed_composition.fix()
    m.feed_temperature = Var(initialize=293.15, units=pyunits.K)
    m.feed_temperature.fix()
    m.feed_pressure = Var(initialize=101325, units=pyunits.Pa)
    m.feed_pressure.fix()
    m.feed_pH = Var(initialize=7.2, units=pyunits.dimensionless)
    m.feed_pH.fix()
    # measured dissolved oxygen in feed, defines its redox state
    m.feed_dissolved_O2 = Var(initialize=feed_DO, units=pyunits.mol / pyunits.s)
    m.feed_dissolved_O2.fix()
    # O2 in air the softening step is open to (bar)
    m.O2_fugacity = Var(initialize=air_O2_fugacity, units=pyunits.dimensionless)
    m.O2_fugacity.fix()

    """ lime softening """
    m.lime_dose = Var(initialize=3e-3, units=pyunits.mol / pyunits.s)
    m.lime_dose.fix()
    m.softening_outputs = Var(
        [("pH", None), ("pE", None), ("speciesAmount", "O2")]
        + [("speciesAmount", mineral) for mineral in minerals],
        initialize=1e-5,
    )
    m.softening_outputs[("pH", None)].set_value(10)
    m.softening_block = ReaktoroBlock(
        aqueous_phase={
            "composition": m.feed_composition,
            "convert_to_rkt_species": True,
            "activity_model": "ActivityModelPhreeqc",
        },
        mineral_phase={"phase_components": minerals},
        system_state={
            "temperature": m.feed_temperature,
            "pressure": m.feed_pressure,
            "pH": m.feed_pH,
            "dissolved_O2": m.feed_dissolved_O2,
        },
        # open softening step to air, O2 is taken up (or stripped) until
        # equilibrium with given O2 fugacity
        system_state_modifier=({"O2_fugacity": m.O2_fugacity} if aerate else {}),
        chemistry_modifier={"CaO": m.lime_dose},
        outputs=m.softening_outputs,
        database_file="wateq4f.dat",
        build_speciation_block=True,
    )
    if aerate:
        build_neutralization(m)
    scale_model(m)
    return m


def build_neutralization(m):
    """softened water composition, precipitated solids are removed"""
    m.softened_composition = Var(
        list(feed.keys()), initialize=feed, units=pyunits.mol / pyunits.s
    )
    precipitated = {
        "Ca": m.softening_outputs[("speciesAmount", "Calcite")],
        "HCO3": m.softening_outputs[("speciesAmount", "Calcite")],
        "Mn": m.softening_outputs[("speciesAmount", "Birnessite")],
    }
    added = {"Ca": m.lime_dose}

    @m.Constraint(list(feed.keys()))
    def eq_softened_composition(m, ion):
        return m.softened_composition[ion] == (
            m.feed_composition[ion]
            + added.get(ion, 0)
            - precipitated.get(ion, 0) * pyunits.mol / pyunits.s
        )

    """ neutralization, pH and DO leaving softening are passed in """
    m.acid_dose = Var(initialize=1e-3, units=pyunits.mol / pyunits.s)
    m.neutralized_outputs = Var(
        [("pH", None), ("pE", None), ("speciesAmount", "O2")], initialize=1
    )
    m.neutralized_outputs[("pH", None)].set_value(8.5)
    m.neutralization_block = ReaktoroBlock(
        aqueous_phase={
            "composition": m.softened_composition,
            "convert_to_rkt_species": True,
            "activity_model": "ActivityModelPhreeqc",
        },
        system_state={
            "temperature": m.feed_temperature,
            "pressure": m.feed_pressure,
            "pH": m.softening_outputs[("pH", None)],
            "dissolved_O2": m.softening_outputs[("speciesAmount", "O2")],
        },
        chemistry_modifier={"HCl": m.acid_dose},
        outputs=m.neutralized_outputs,
        database_file="wateq4f.dat",
        build_speciation_block=True,
    )


def scale_model(m):
    iscale.set_scaling_factor(m.feed_dissolved_O2, 1 / feed_DO)
    iscale.set_scaling_factor(m.lime_dose, 1 / m.lime_dose.value)
    for key in m.feed_composition:
        iscale.set_scaling_factor(
            m.feed_composition[key], 1 / m.feed_composition[key].value
        )
    if m.aerate:
        for key in m.feed_composition:
            sf = 1 / m.feed_composition[key].value
            iscale.set_scaling_factor(m.softened_composition[key], sf)
            iscale.constraint_scaling_transform(m.eq_softened_composition[key], sf)
        iscale.set_scaling_factor(m.acid_dose, 1 / m.acid_dose.value)


def initialize(m):
    m.softening_block.initialize()
    if m.aerate:
        for ion in m.eq_softened_composition:
            calculate_variable_from_constraint(
                m.softened_composition[ion], m.eq_softened_composition[ion]
            )
        m.neutralization_block.initialize()


def solve(m, final_pH=8.5):
    if m.aerate:
        # find acid dose that neutralizes softened water to final pH
        m.acid_dose.unfix()
        m.neutralized_outputs[("pH", None)].fix(final_pH)
    cy_solver = get_cyipopt_watertap_solver()
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    return result


def get_results(m):
    results = {
        "feed DO (mg/L)": value(m.feed_dissolved_O2) * O2_mw,
        "softened pH": value(m.softening_outputs[("pH", None)]),
        "softened pE": value(m.softening_outputs[("pE", None)]),
        "softened DO (mg/L)": value(m.softening_outputs[("speciesAmount", "O2")])
        * O2_mw,
        "Mn removed (%)": value(
            m.softening_outputs[("speciesAmount", "Birnessite")]
            / m.feed_composition["Mn"]
        )
        * 100,
        "Calcite (mol/s)": value(m.softening_outputs[("speciesAmount", "Calcite")]),
    }
    if m.aerate:
        results["acid dose (mol/s)"] = value(m.acid_dose)
        results["final pH"] = value(m.neutralized_outputs[("pH", None)])
        results["final pE"] = value(m.neutralized_outputs[("pE", None)])
        results["final DO (mg/L)"] = (
            value(m.neutralized_outputs[("speciesAmount", "O2")]) * O2_mw
        )
    return results


def display_results(results, aerate):
    print(f"---- softening {'open to air' if aerate else 'sealed'} ----")
    for key, val in results.items():
        print(f"{key}: {val:.4g}")


if __name__ == "__main__":
    main()
