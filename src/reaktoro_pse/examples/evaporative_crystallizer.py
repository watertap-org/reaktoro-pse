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
from pyomo.environ import (
    ConcreteModel,
    Var,
    Constraint,
    Objective,
    Param,
    value,
    check_optimal_termination,
    units as pyunits,
)
from pyomo.util.calc_var_value import calculate_variable_from_constraint
import idaes.core.util.scaling as iscale

from reaktoro_pse.reaktoro_block import ReaktoroBlock
from reaktoro_pse.core.util_classes.cyipopt_solver import (
    get_cyipopt_watertap_solver,
)

# This example demonstrates why species_amount_zero_offset is needed when Reaktoro
# outputs of mineral phases feed bounded variables in the flowsheet.
#
# An evaporative crystallizer concentrates seawater. A simple (naive) model computes
# the evaporation and the liquor composition from a mass balance, Reaktoro is only used to
# compute how much of each of 17 minerals precipitates. Most minerals are not present at
# any given point, Reaktoro returns their amount as epsilon (~1e-40) with a zero jacobian
# row, and they are linked to solids variables with a lower bound of 0.
#
# Without an offset, every absent mineral pins its solids variable at its bound, Ipopt
# keeps it a tiny distance from the bound, and the bound multipliers grow as
# barrier parameter / distance. The optimization (maximize evaporation, subject to a
# slurry solids limit, with lime and soda ash dosing) then fails to converge.
# With the default offset (smooth_floor, 1e-4 in scaled units) absent minerals sit at
# a = 1e-4 / scaling factor away from the bound, present minerals are shifted by only
# ~a**2/(2n), and the same optimization converges.
#
# Key assumptions:
# - pH is fixed (pH control), temperature 25 C, all water removed is pure vapor
# - dosed CaO and Na2CO3 react in the crystallizer

# seawater, mol/s for ~1 kg/s of water
FEED = {
    "H2O": 55.5,
    "Na": 0.469,
    "Cl": 0.546,
    "Mg": 0.0528,
    "SO4": 0.0282,
    "Ca": 0.0103,
    "K": 0.0102,
    "HCO3": 0.00206,
}
# mineral stoichiometry used in the naive mass balance
MINERALS = {
    "Calcite": {"Ca": 1, "HCO3": 1},
    "Dolomite": {"Ca": 1, "Mg": 1, "HCO3": 2},
    "Magnesite": {"Mg": 1, "HCO3": 1},
    "Gypsum": {"Ca": 1, "SO4": 1, "H2O": 2},
    "Anhydrite": {"Ca": 1, "SO4": 1},
    "Glauberite": {"Na": 2, "Ca": 1, "SO4": 2},
    "Polyhalite": {"K": 2, "Mg": 1, "Ca": 2, "SO4": 4, "H2O": 2},
    "Halite": {"Na": 1, "Cl": 1},
    "Mirabilite": {"Na": 2, "SO4": 1, "H2O": 10},
    "Thenardite": {"Na": 2, "SO4": 1},
    "Epsomite": {"Mg": 1, "SO4": 1, "H2O": 7},
    "Hexahydrite": {"Mg": 1, "SO4": 1, "H2O": 6},
    "Bloedite": {"Na": 2, "Mg": 1, "SO4": 2, "H2O": 4},
    "Kainite": {"K": 1, "Mg": 1, "Cl": 1, "SO4": 1, "H2O": 3},
    "Sylvite": {"K": 1, "Cl": 1},
    "Carnallite": {"K": 1, "Mg": 1, "Cl": 3, "H2O": 6},
    "Bischofite": {"Mg": 1, "Cl": 2, "H2O": 6},
}


def main(species_amount_zero_offset=1e-4, species_amount_offset_type="smooth_floor"):
    """build, initialize at 50% evaporation and optimize evaporation, returns
    (model, optimization converged)"""
    m = build_crystallizer(species_amount_zero_offset, species_amount_offset_type)
    initialize(m)
    setup_optimization(m)
    converged = solve(m)
    display_results(m)
    return m, converged


def build_crystallizer(
    species_amount_zero_offset=1e-4, species_amount_offset_type="smooth_floor"
):
    m = ConcreteModel()
    m.feed = Var(list(FEED), initialize=FEED, units=pyunits.mol / pyunits.s)
    m.feed.fix()
    m.temperature = Var(initialize=298.15, units=pyunits.K)
    m.temperature.fix()
    m.pressure = Var(initialize=1e5, units=pyunits.Pa)
    m.pressure.fix()
    m.pH = Var(initialize=8.0)
    m.pH.fix()
    m.evaporation = Var(initialize=0.5, bounds=(0, 0.999))
    m.evaporation.fix()
    # concentrate entering the crystallizer
    m.concentrate = Var(
        list(FEED), initialize=FEED, bounds=(0, None), units=pyunits.mol / pyunits.s
    )
    # precipitated solids from Reaktoro, lower bound of 0
    m.solids = Var(
        [("speciesAmount", k) for k in MINERALS], initialize=0, bounds=(0, None)
    )
    m.dose = Var(
        ["CaO", "Na2CO3"],
        initialize=1e-5,
        bounds=(1e-5, None),
        units=pyunits.mol / pyunits.s,
    )
    m.dose.fix()
    m.liquor = Var(
        list(FEED), initialize=FEED, bounds=(0, None), units=pyunits.mol / pyunits.s
    )

    @m.Constraint(list(FEED))
    def eq_concentrate(b, i):
        if i == "H2O":
            return b.concentrate[i] == b.feed[i] * (1 - b.evaporation)
        return b.concentrate[i] == b.feed[i]

    @m.Constraint(list(FEED))
    def eq_liquor(b, i):
        dosed = {
            "Ca": b.dose["CaO"],
            "Na": 2 * b.dose["Na2CO3"],
            "HCO3": b.dose["Na2CO3"],
        }.get(i, 0)
        return b.liquor[i] == b.concentrate[i] + dosed - sum(
            st.get(i, 0) * b.solids["speciesAmount", k] for k, st in MINERALS.items()
        )

    m.crystallizer = ReaktoroBlock(
        aqueous_phase={
            "composition": m.concentrate,
            "convert_to_rkt_species": True,
            "activity_model": "ActivityModelPitzer",
        },
        mineral_phase={"phase_components": list(MINERALS)},
        system_state={
            "temperature": m.temperature,
            "pressure": m.pressure,
            "pH": m.pH,
        },
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        outputs=m.solids,
        chemistry_modifier=m.dose,
        hessian_options={"hessian_type": "SR1"},
        reaktoro_solve_options={
            "species_amount_zero_offset": species_amount_zero_offset,
            "species_amount_offset_type": species_amount_offset_type,
        },
    )
    scale_model(m)
    return m


def scale_model(m):
    for i, v in FEED.items():
        sf = 1 / v
        for var in (m.feed[i], m.concentrate[i], m.liquor[i]):
            iscale.set_scaling_factor(var, sf)
        iscale.constraint_scaling_transform(m.eq_concentrate[i], sf)
        iscale.constraint_scaling_transform(m.eq_liquor[i], sf)
    # the offset is applied in scaled units, so the scaling of solids sets the offset size,
    # scale each mineral by the amount of its most limiting ion in the feed
    for k, st in MINERALS.items():
        limiting = min(FEED[i] / n for i, n in st.items() if i != "H2O")
        iscale.set_scaling_factor(m.solids["speciesAmount", k], 1 / limiting)


def initialize(m):
    for i in FEED:
        calculate_variable_from_constraint(m.concentrate[i], m.eq_concentrate[i])
    m.crystallizer.initialize()
    for i in FEED:
        calculate_variable_from_constraint(m.liquor[i], m.eq_liquor[i])
    if not solve(m):
        raise RuntimeError("Crystallizer failed to initialize")


def setup_optimization(m, solids_limit=0.005):
    """maximize evaporation, subject to slurry solids limit (mol solids per mol of water in
    liquor), with lime and soda ash dosing as decision variables"""
    m.solids_limit = Param(initialize=solids_limit, mutable=True)
    m.eq_slurry = Constraint(
        expr=sum(m.solids.values()) <= m.solids_limit * m.liquor["H2O"]
    )
    for v in m.dose.values():
        v.unfix()
        iscale.set_scaling_factor(v, 1e3)
    m.objective = Objective(expr=-10 * m.evaporation + m.dose["CaO"] + m.dose["Na2CO3"])
    m.evaporation.unfix()


def solve(m, max_iter=300, tee=True):
    solver = get_cyipopt_watertap_solver(max_iter=max_iter)
    result = solver.solve(m, tee=tee)
    return check_optimal_termination(result)


def display_results(m):
    print(f"evaporation: {value(m.evaporation):.4f}")
    print(
        f"CaO dose: {value(m.dose['CaO']):.3e}, Na2CO3 dose: {value(m.dose['Na2CO3']):.3e}"
    )
    # absent minerals sit at the offset, use raw reaktoro amounts to find present minerals
    raw = m.crystallizer.rkt_solver.output_specs.user_outputs
    for key, var in m.solids.items():
        if raw[key].value > 1e-12:
            print(f"{key[1]} precipitated: {var.value:.4e} mol/s")
    if m.crystallizer.rkt_block_builder.species_amount_offset_index:
        # error introduced by the offset, per species and per element (logged), warns if
        # a present mineral or an element is changed by more than 1%
        m.crystallizer.check_species_amount_offsets()


if __name__ == "__main__":
    print("--- without species amount offset ---")
    _, converged = main(species_amount_zero_offset=0, species_amount_offset_type=None)
    print(f"converged: {converged}")
    print("--- with default species amount offset ---")
    _, converged = main()
    print(f"converged: {converged}")
