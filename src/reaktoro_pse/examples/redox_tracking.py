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
from reaktoro_pse.core.util_classes.cyipopt_solver import (
    get_cyipopt_watertap_solver,
)
from pyomo.environ import (
    ConcreteModel,
    Var,
    Objective,
    Constraint,
    assert_optimal_termination,
    units as pyunits,
)

from pyomo.util.calc_var_value import calculate_variable_from_constraint

import idaes.core.util.scaling as iscale

from reaktoro_pse.parallel_tools.reaktoro_block_manager import (
    ReaktoroBlockManager,
)

# This example demonstrates importance of tracking pE in larger databases such as waterq4f or pitzer,
# redox species are typically Fe/Mn, but can include H2O, depending on how database is setup.
# In these cases its critical to pass pE from one block to the next to ensure true speciation is resolved.


def main():
    m = build_modifer_blocks()
    initialize(m)
    solve(m)
    results_array = log_results(m)
    for r in [0.01]:
        m.acid_modifier.fix(r)
        results_array = log_results(m, results_array)
        solve(m)
    print(results_array)
    return results_array


def build_modifer_blocks(parallel_mode=False):
    m = ConcreteModel()
    m.feed_composition = Var(
        ["H2O", "Mg", "Na", "Fe", "Mn", "Cl", "SO4", "Ca", "HCO3"],
        initialize=1,
        units=pyunits.mol / pyunits.s,
    )
    m.feed_composition.construct()
    m.feed_composition["H2O"].fix(55)
    m.feed_composition["Mg"].fix(0.01)
    m.feed_composition["Mn"].fix(0.01)
    m.feed_composition["Fe"].fix(0.01)
    m.feed_composition["Na"].fix(0.025)
    m.feed_composition["Cl"].fix(0.025)
    m.feed_composition["Ca"].fix(0.001)
    m.feed_composition["HCO3"].fix(0.01)
    m.feed_composition["SO4"].fix(0.01)
    m.feed_temperature = Var(initialize=293.15, units=pyunits.K)
    m.feed_temperature.fix()
    m.feed_pressure = Var(initialize=1e5, units=pyunits.Pa)
    m.feed_pressure.fix()

    m.pE = Var(initialize=4, units=pyunits.dimensionless)
    m.pE.fix()
    m.pH = Var(initialize=7, bounds=(4, 12), units=pyunits.dimensionless)
    m.pE.fix()

    m.acid_modifier = Var(
        initialize=1e-5, bounds=(1e-10, 1e-2), units=pyunits.mol / pyunits.s
    )
    m.acid_modifier.fix(0.00001)
    m.modified_properties_water_removal = Var(
        initialize=1e-8, units=pyunits.mol / pyunits.s
    )
    m.modified_properties_water_removal.fix(10)
    m.acid_outputs = Var(
        [("pH", None), ("pE", None)], initialize=7, units=pyunits.dimensionless
    )
    m.scaling_no_pE_outputs = Var(
        [
            ("pH", None),
            ("pE", None),
            ("scalingTendency", "Calcite"),
            ("scalingTendency", "Gypsum"),
        ],
        initialize=1,
    )
    m.scaling_with_pE_outputs = Var(
        [
            ("pH", None),
            ("pE", None),
            ("scalingTendency", "Calcite"),
            ("scalingTendency", "Gypsum"),
        ],
        initialize=1,
    )
    if parallel_mode:
        m.parallel_block_manager = ReaktoroBlockManager()
    else:
        m.parallel_block_manager = None
    m.eq_acidifier_block = ReaktoroBlock(
        aqueous_phase={
            "composition": m.feed_composition,
            "convert_to_rkt_species": True,
            "activity_model": "ActivityModelPhreeqc",
        },
        database_file="wateq4f.dat",
        system_state={
            "temperature": m.feed_temperature,
            "pressure": m.feed_pressure,
            "pH": m.pH,
            "pE": m.pE,
        },
        outputs=m.acid_outputs,
        chemistry_modifier={
            "HCl": m.acid_modifier,
        },
        dissolve_species_in_reaktoro=True,
        # we can use default converter as its defined for default database (Phreeqc and pitzer)
        # we are modifying state and must speciate inputs before adding acid to find final prop state.
        build_speciation_block=True,
        reaktoro_block_manager=m.parallel_block_manager,
    )
    m.eq_scaling_no_pE = ReaktoroBlock(
        aqueous_phase={
            "composition": m.feed_composition,
            "convert_to_rkt_species": True,
            "activity_model": "ActivityModelPhreeqc",
        },
        database_file="wateq4f.dat",
        system_state={
            "temperature": m.feed_temperature,
            "pressure": m.feed_pressure,
            "pH": m.acid_outputs[("pH", None)],
        },
        outputs=m.scaling_no_pE_outputs,
        chemistry_modifier={
            "HCl": m.acid_modifier,
            "H2O_evaporation": m.modified_properties_water_removal,
        },
        dissolve_species_in_reaktoro=True,
        # we can use default converter as its defined for default database (Phreeqc and pitzer)
        # we are modifying state and must speciate inputs before adding acid to find final prop state.
        build_speciation_block=True,
        reaktoro_block_manager=m.parallel_block_manager,
    )
    m.eq_scaling_with_pE = ReaktoroBlock(
        aqueous_phase={
            "composition": m.feed_composition,
            "convert_to_rkt_species": True,
            "activity_model": "ActivityModelPhreeqc",
        },
        database_file="wateq4f.dat",
        system_state={
            "temperature": m.feed_temperature,
            "pressure": m.feed_pressure,
            "pH": m.acid_outputs[("pH", None)],
            "pE": m.acid_outputs[("pE", None)],
        },
        outputs=m.scaling_with_pE_outputs,
        chemistry_modifier={
            "HCl": m.acid_modifier,
            "H2O_evaporation": m.modified_properties_water_removal,
        },
        dissolve_species_in_reaktoro=True,
        # we can use default converter as its defined for default database (Phreeqc and pitzer)
        # we are modifying state and must speciate inputs before adding acid to find final prop state.
        build_speciation_block=True,
        reaktoro_block_manager=m.parallel_block_manager,
    )
    if parallel_mode:
        m.parallel_block_manager.build_reaktoro_blocks()
    scale_model(m)
    return m


def scale_model(m):
    for key in m.feed_composition:
        iscale.set_scaling_factor(
            m.feed_composition[key], 1 / m.feed_composition[key].value
        )
    iscale.set_scaling_factor(m.acid_modifier, 1 / 0.0001)
    iscale.set_scaling_factor(m.modified_properties_water_removal, 1 / 0.0001)


def initialize(m):
    m.eq_acidifier_block.initialize()
    m.eq_scaling_no_pE.initialize()
    m.eq_scaling_with_pE.initialize()


def display_results(m):
    for key, obj in m.scaling_with_pE_outputs.items():
        print(
            f"Output for {key} has value of {obj.value} with pE tracking and {m.scaling_no_pE_outputs[key].value} with out pE tracking"
        )


def log_results(m, result_array=None):
    if result_array is None:
        result_array = {"inputs": {}, "no_pe": {}, "with_pe": {}}
        result_array["inputs"]["pH"] = []
        result_array["inputs"]["acid_addition"] = []
        result_array["inputs"]["pE"] = []
        result_array["inputs"]["water_removal"] = []
        for key, obj in m.scaling_with_pE_outputs.items():
            result_array["no_pe"][key] = []
            result_array["with_pe"][key] = []
    if result_array is not None:
        result_array["inputs"]["pH"].append(float(m.pH.value))
        result_array["inputs"]["acid_addition"].append(float(m.acid_modifier.value))
        result_array["inputs"]["pE"].append(float(m.pE.value))
        result_array["inputs"]["water_removal"].append(
            float(m.modified_properties_water_removal.value)
        )
        for key, obj in m.scaling_with_pE_outputs.items():
            result_array["no_pe"][key].append(float(obj.value))
            result_array["with_pe"][key].append(
                float(m.scaling_no_pE_outputs[key].value)
            )
    return result_array


def solve(m):
    cy_solver = get_cyipopt_watertap_solver()
    result = cy_solver.solve(m, tee=True)
    display_results(m)
    assert_optimal_termination(result)
    return result


if __name__ == "__main__":
    main()
