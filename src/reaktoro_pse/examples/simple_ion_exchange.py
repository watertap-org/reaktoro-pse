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

"""
This examples demonstrates how reaktoro graybox can be used to estimates remove of specific ion through use of ion exchange material
as well as how to find optimal amount of ion exchange material to maximize removal of target ion. 

Key assumptions:
TBD
"""


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
    m.treated_pH = Var(initialize=7, bounds=(0, 12), units=pyunits.dimensionless)
    m.Ca_to_Mg_selectivity = Var(initialize=1, units=pyunits.dimensionless)

    # acid addition
    m.acid_addition = Var(initialize=0.00001, units=pyunits.mol / pyunits.s)
    m.acid_addition.fix()
    # base addition
    m.base_addition = Var(initialize=0.00001, units=pyunits.mol / pyunits.s)
    m.base_addition.fix()
    m.feed_charge = Var(initialize=0.00001, units=pyunits.mol / pyunits.s)

    # Clean ion exchange
    m.ion_exchange_material = Var(
        ["NaX", "CaX2", "MgX2"],
        initialize=0,
        units=pyunits.mol / pyunits.s,
    )
    m.used_ion_exchange_material = Var(
        ["NaX", "CaX2", "MgX2", "X-"],
        initialize=0,
        units=pyunits.mol / pyunits.s,
    )
    # We have a resin that is fresh and primarily contains Na balancing ions
    # note that each NaX/CaX2/MgX2 is realyl a site on a resin material rather then
    # total mass of resin
    m.ion_exchange_material["NaX"].fix(1)
    m.ion_exchange_material["CaX2"].fix(1e-5)
    m.ion_exchange_material["MgX2"].fix(1e-5)

    """ We will build a block to charge neutralize the feed and adjust apparent species 
    to achieve this and separate block to do ion exchange equilibration"""
    m.eq_speciation_block = ReaktoroBlock(
        composition=m.feed_composition,
        temperature=m.feed_temperature,
        pressure=m.feed_pressure,
        pH=m.feed_pH,
        outputs={("charge", None): m.feed_charge},
        aqueous_phase_activity_model=rkt.ActivityModelPitzer(),
        dissolve_species_in_reaktoro=True,
        assert_charge_neutrality=False,
        # we can use default converter as its defined for default database (Phreeqc and pitzer)
        convert_to_rkt_species=True,
        # we are modifying state and must speciate inputs before adding acid to find final prop state.
        build_speciation_block=True,
        open_species_on_property_block=["H+", "OH-"],
    )

    """combine all inputs"""
    m.input_dict = {}
    for key, obj in m.feed_composition.items():
        m.input_dict[key] = obj
    for key, obj in m.ion_exchange_material.items():
        m.input_dict[key] = obj
    """ this will use charge neutralized feed and perform ion exchange calculations"""
    m.eq_ix_properties = ReaktoroBlock(
        composition=m.input_dict,
        temperature=m.feed_temperature,
        pressure=m.feed_pressure,
        pH=m.feed_pH,
        outputs={
            "speciesAmount": True,
            ("pH", None): m.treated_pH,
        },
        aqueous_phase_activity_model=rkt.ActivityModelPitzer(),
        ion_exchange_phase=["NaX", "KX", "CaX2"],
        ion_exchange_phase_activity_model=rkt.ActivityModelIonExchangeGainesThomas(),
        chemistry_modifier={"NaOH": m.base_addition, "HCl": m.acid_addition},
        dissolve_species_in_reaktoro=True,
        assert_charge_neutrality=False,
        # we can use default converter as its defined for default database (Phreeqc and pitzer)
        convert_to_rkt_species=True,
        species_to_rkt_species_dict={  # WE need to supply our own conversion dict as default does not support X
            "MgX2": "MgX2",
            "CaX2": "CaX2",
            "NaX": "NaX",
            "H2O": "H2O",
            "Mg": "Mg+2",
            "Na": "Na+",
            "Ca": "Ca+2",
            "Cl": "Cl-",
            "HCO3": "HCO3-",
            "SO4": "SO4-2",
        },
        # we are modifying state and must speciate inputs before adding acid to find final prop state.
        build_speciation_block=True,
        # jacobian_scaling_type="no_scaling",
    )

    """ currently ReaktoroBlock does not support 
    automatic conversion of True species to apparent species with out getting direct element amounts, instead
    we can use lower level core api to get input to output conversion dictionary for our apparant species to exact species"""
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

    @m.Constraint(list(m.used_ion_exchange_material.keys()))
    def eq_used_ix(fs, ion):
        return (
            m.used_ion_exchange_material[ion]
            == m.eq_ix_properties.outputs[("speciesAmount", ion)]
        )

    @m.Constraint(list(m.treated_composition.keys()))
    def eq_removal(fs, ion):
        return (
            m.removal_percent[ion]
            == (m.treated_composition[ion] - m.feed_composition[ion])
            / m.feed_composition[ion]
            * 100
        )

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
    iscale.constraint_scaling_transform(m.eq_selectivity, 1)
    iscale.set_scaling_factor(m.acid_addition, 1e5)
    iscale.set_scaling_factor(m.base_addition, 1e5)
    iscale.set_scaling_factor(m.Ca_to_Mg_selectivity, 1)
    iscale.set_scaling_factor(m.feed_charge, 1)


def initialize(m):
    m.eq_speciation_block.initialize()
    m.eq_ix_properties.initialize()
    for key in m.treated_composition:
        calculate_variable_from_constraint(
            m.treated_composition[key], m.eq_treated_comp[key]
        )
    m.feed_composition["Cl"].unfix()
    m.feed_charge.fix(0)
    solve(m)


def setup_optimization(m):
    """Find resin amount and base/acid addition that maximize calcium selectivity"""
    m.objective = Objective(
        expr=(1 / m.Ca_to_Mg_selectivity) + m.base_addition + m.acid_addition
    )
    m.base_addition.unfix()
    m.acid_addition.unfix()
    m.removal_percent["Mg"].setub(0)
    m.removal_percent["Ca"].setub(-20)


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
