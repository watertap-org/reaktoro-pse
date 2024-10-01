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
import pyomo.environ as pyo
import reaktoro as rkt

__author__ = "Alexander V. Dudchenko"

# This examples demonstrates how Reaktoro graybox can be used to estimates removal of specific ion through use of ion exchange material.

# This example shows how to:
# (1) Set up ReaktoroBlock for charge neutralizing the feed composition
# (2) Use outputs from speciation block as inputs into a second property block
# (3) Add Ion Exchange phase and species into ReaktoroBlock
# (4) Optimize addition of acid and bases for maximizing Calcium removal selectivity over Magnesium


def main():
    m = build_simple_desal()
    initialize(m)
    setup_optimization(m)
    solve(m)
    return m


def build_simple_desal():
    m = ConcreteModel()
    m.feed_composition = Var(
        ["H2O", "Mg", "Na", "Cl", "SO4", "Ca", "HCO3"],
        initialize=1,
        units=pyunits.mol / pyunits.s,
    )
    m.treated_composition = Var(
        ["H2O", "Mg", "Na", "Cl", "SO4", "Ca", "HCO3"],
        initialize=1,
        units=pyunits.mol / pyunits.s,
    )
    m.removal_percent = Var(
        ["H2O", "Mg", "Na", "Cl", "SO4", "Ca", "HCO3"],
        initialize=1,
        units=pyunits.dimensionless,
    )
    m.feed_composition.construct()
    m.feed_composition["H2O"].fix(55)
    m.feed_composition["Mg"].fix(0.1)
    m.feed_composition["Na"].fix(0.25)
    m.feed_composition["Cl"].fix(0.25)
    m.feed_composition["Ca"].fix(0.01)
    m.feed_composition["HCO3"].fix(0.01)
    m.feed_composition["SO4"].fix(0.01)
    m.feed_temperature = Var(initialize=293.15, units=pyunits.K)
    m.feed_temperature.fix()
    m.feed_pressure = Var(initialize=1e5, units=pyunits.Pa)
    m.feed_pressure.fix()
    m.feed_pH = Var(initialize=7, bounds=(4, 12), units=pyunits.dimensionless)
    m.feed_pH.fix()
    m.treated_pH = Var(initialize=7, bounds=(0, 14), units=pyunits.dimensionless)
    m.Ca_to_Mg_selectivity = Var(initialize=1, units=pyunits.dimensionless)

    # acid addition
    m.acid_addition = Var(initialize=0.00001, units=pyunits.mol / pyunits.s)
    m.acid_addition.fix()
    # base addition
    m.base_addition = Var(initialize=0.00001, units=pyunits.mol / pyunits.s)
    m.base_addition.fix()
    m.feed_charge = Var(initialize=0.00001, units=pyunits.mol / pyunits.s)
    m.treated_feed_charge = Var(initialize=0.00001, units=pyunits.mol / pyunits.s)
    # Clean ion exchange
    m.ion_exchange_material = Var(
        ["NaX", "CaX2", "MgX2"],
        initialize=0,
        units=pyunits.mol / pyunits.s,
    )
    m.used_ion_exchange_material = Var(
        ["NaX", "CaX2", "MgX2"],  # , "X-"],
        initialize=0,
        units=pyunits.mol / pyunits.s,
    )
    # We have a resin that is fresh and primarily contains Na balancing ions
    # note that each NaX/CaX2/MgX2 is really a site on a resin material rather then
    # total mass of resin
    m.ion_exchange_material["NaX"].fix(0.4)
    m.ion_exchange_material["CaX2"].fix(1e-5)
    m.ion_exchange_material["MgX2"].fix(1e-5)

    # We will build a block to charge neutralize the feed and adjust apparent species
    # to achieve this and then build a separate block to do ion exchange calculation
    m.eq_speciation_block = ReaktoroBlock(
        system_state={
            "temperature": m.feed_temperature,
            "pressure": m.feed_pressure,
            "pH": m.feed_pH,
        },
        aqueous_phase={
            "composition": m.feed_composition,
            "convert_to_rkt_species": True,
            "activity_model": rkt.ActivityModelPitzer(),
        },
        outputs={
            ("charge", None): m.feed_charge,
            "speciesAmount": True,
        },
        dissolve_species_in_reaktoro=True,
        assert_charge_neutrality=False,
        build_speciation_block=False,
    )

    # build the IX block

    # Note how we do not privde feed pH any more, as it will be estimated
    # from our specitated and charge balanced block aad adjusted due to addition of resin - which will
    # impact both charge balance and final pH of solution

    m.eq_speciation_block.display_reaktoro_state()
    m.eq_speciation_block.outputs.display()
    m.eq_ix_properties = ReaktoroBlock(
        aqueous_phase={
            "composition": m.eq_speciation_block.outputs,
            "convert_to_rkt_species": False,  # already exact
            "activity_model": rkt.ActivityModelPitzer(),
        },
        ion_exchange_phase={
            "composition": m.ion_exchange_material,
            "activity_model": rkt.ActivityModelIonExchangeGainesThomas(),
            "convert_to_rkt_species": False,  # already exact
        },
        system_state={"temperature": m.feed_temperature, "pressure": m.feed_pressure},
        outputs={
            "speciesAmount": True,
            ("pH", None): m.treated_pH,
            ("charge", None): m.treated_feed_charge,
        },
        chemistry_modifier={"NaOH": m.base_addition, "HCl": m.acid_addition},
        dissolve_species_in_reaktoro=True,
        assert_charge_neutrality=False,
        reaktoro_solve_options={
            "solver_tolerance": 1e-12,
            "open_species_on_property_block": ["OH-"],
        },
        # we do not need to re-speciate.
        exact_speciation=True,
        build_speciation_block=False,
    )
    m.eq_ix_properties.display_reaktoro_state()
    # assert False
    # Currently ReaktoroBlock does not support
    # automatic conversion of output True species to Apparent species, instead
    # we can use lower level core api to get input to output conversion dictionary
    # for our Apparent species to True/Exact species
    conversion_dict = m.eq_ix_properties.rkt_inputs.constraint_dict
    for key, speciation in conversion_dict.items():
        print(key, speciation)
    apparent_element_to_species = {
        "Na": "Na",
        "HCO3": "C",
        "Mg": "Mg",
        "Ca": "Ca",
        "SO4": "S",
        "Cl": "Cl",
    }

    # Here we write constraints to estimate treated composition

    @m.Constraint(list(m.treated_composition.keys()))
    def eq_treated_comp(fs, ion):
        if ion == "H2O":
            return (
                m.treated_composition[ion]
                == m.eq_ix_properties.outputs[("speciesAmount", ion)]
                * pyunits.mol
                / pyunits.s
            )
        if ion != "H2O":
            sum_ion = []
            element = apparent_element_to_species[ion]
            for mols, specie in conversion_dict[element]:
                if "X" not in specie and "NaOH" not in specie and "HCl" not in specie:
                    sum_ion.append(
                        mols * m.eq_ix_properties.outputs[("speciesAmount", specie)]
                    )
            return m.treated_composition[ion] == sum(sum_ion) * pyunits.mol / pyunits.s

    # Track changes to IX resin

    @m.Constraint(list(m.used_ion_exchange_material.keys()))
    def eq_used_ix(fs, ion):
        return (
            m.used_ion_exchange_material[ion]
            == m.eq_ix_properties.outputs[("speciesAmount", ion)]
        )

    # Track percent removal of each Apparent specie

    @m.Constraint(list(m.treated_composition.keys()))
    def eq_removal(fs, ion):
        return (
            m.removal_percent[ion]
            == (m.treated_composition[ion] - m.feed_composition[ion])
            / m.feed_composition[ion]
            * 100
        )

    # Calculate Ca to Mg selectivity
    m.eq_selectivity = Constraint(
        expr=m.Ca_to_Mg_selectivity == m.removal_percent["Ca"] / m.removal_percent["Mg"]
    )
    scale_model(m)
    return m


def scale_model(m):
    for key in m.feed_composition:
        iscale.set_scaling_factor(
            m.feed_composition[key], 1 / m.feed_composition[key].value
        )
        iscale.set_scaling_factor(
            m.treated_composition[key], 1 / m.feed_composition[key].value
        )
        iscale.constraint_scaling_transform(
            m.eq_treated_comp[key], 1 / m.feed_composition[key].value
        )
        iscale.constraint_scaling_transform(m.eq_removal[key], 1 / 100)
        iscale.set_scaling_factor(m.removal_percent[key], 1 / 100)
    for key in m.eq_used_ix:
        iscale.set_scaling_factor(m.used_ion_exchange_material[key], 10)
        iscale.constraint_scaling_transform(m.eq_used_ix[key], 10)
    iscale.constraint_scaling_transform(m.eq_selectivity, 10)
    iscale.set_scaling_factor(m.acid_addition, 1e2)
    iscale.set_scaling_factor(m.base_addition, 1e2)
    iscale.set_scaling_factor(m.Ca_to_Mg_selectivity, 10)
    iscale.set_scaling_factor(m.feed_charge, 1e8)
    iscale.set_scaling_factor(m.treated_feed_charge, 1e8)


def initialize(m):
    m.eq_speciation_block.initialize()
    m.eq_ix_properties.initialize()
    m.eq_ix_properties.display_jacobian_scaling()
    m.eq_ix_properties.display_reaktoro_state()
    for key in m.treated_composition:
        calculate_variable_from_constraint(
            m.treated_composition[key], m.eq_treated_comp[key]
        )
    # Unfix our feed Cl and fix charge to zero so we can charge neutralize the feed
    m.feed_composition["Cl"].unfix()
    m.feed_charge.fix(0)
    solve(m)


def setup_optimization(m):
    # Find resin amount and base/acid addition that maximize calcium selectivity
    m.objective = Objective(
        expr=(
            (1 / m.Ca_to_Mg_selectivity) * 10
            + m.base_addition * 10
            + m.acid_addition * 10
        )
    )
    m.base_addition.unfix()
    m.acid_addition.fix()
    m.removal_percent["Mg"].setub(-10)
    m.removal_percent["Ca"].setub(-10)


def display_results(m):
    print("result")
    for key, obj in m.treated_composition.items():
        print(f"{key} feed {m.feed_composition[key].value}, treated {obj.value}")
    print(f"Ca to Mg eq_selectivity {m.Ca_to_Mg_selectivity.value}")
    for key, obj in m.removal_percent.items():
        print(f"Change in {key} = {obj.value} %")
    for key, obj in m.used_ion_exchange_material.items():
        print(f"Mols of {key} = {obj.value} in final resin")
    print(f"feed pH {m.feed_pH.value}, treated pH {m.treated_pH.value}")
    print(f"HCl added {m.acid_addition.value}")
    print(f"NaOH added {m.base_addition.value}")
    print(f"Treated solution charge {m.treated_feed_charge.value}")


def solve(m):
    cy_solver = get_solver(solver="cyipopt-watertap")
    cy_solver.options["max_iter"] = 100
    # only enable if avaialbe !
    # cy_solver.options["linear_solver"] = "ma27"
    result = cy_solver.solve(m, tee=True)
    display_results(m)
    return result


if __name__ == "__main__":
    main()
