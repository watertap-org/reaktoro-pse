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
    Constraint,
    assert_optimal_termination,
    units as pyunits,
)
from reaktoro_pse.parallel_tools.reaktoro_block_manager import (
    ReaktoroBlockManager,
)

from pyomo.util.calc_var_value import calculate_variable_from_constraint

import idaes.core.util.scaling as iscale

__author__ = "Alexander V. Dudchenko"

# This examples demonstrates how reaktoro graybox can be used to enthalpy and water vapor pressure.

# NOTE: For water vapor calculations, pay attention to speciation and assumptions. Please
# refere to these two discussions:

# https://github.com/reaktoro/reaktoro/discussions/398
# https://github.com/reaktoro/reaktoro/discussions/285

# This example demonstrates:
# (1) How to configure different database from default
# (2) How to get enthalpy and vapor pressure from Reaktoro
# (3) How to setup precipitation calculation
# (4) Setup simulation for removal of Calcite over different temperatures and estimate required energy input

# Key assumptions:
# Assumes that process concentrating the feed does not alter the pH.
# This might be a good assumptions for process such as RO, but might be a poor
# assumption for evaporative processes.


def main(hess_type="BFGS_mod"):
    m = build_simple_precipitation(hess_type=hess_type)
    initialize(m)

    return m


def build_simple_precipitation(hess_type=None, parallel_mode=False):
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
    m.feed_composition["Cl"].fix(0.03)
    m.feed_composition["Ca"].fix(0.002)
    m.feed_composition["HCO3"].fix(0.01)
    m.feed_composition["SO4"].fix(0.02)
    m.feed_temperature = Var(initialize=293.15, units=pyunits.K)
    m.feed_temperature.fix()
    m.feed_pressure = Var(initialize=101325, units=pyunits.Pa)
    m.feed_pressure.fix()
    m.feed_pH = Var(initialize=7, bounds=(4, 12), units=pyunits.dimensionless)
    m.feed_pH.fix()
    m.precipitator_composition = Var(
        list(m.feed_composition.keys()),
        initialize=1,
        units=pyunits.mol / pyunits.s,
    )
    m.sludge_water_content = Var(initialize=0.8)
    m.sludge_water_content.fix()
    m.treated_composition = Var(
        list(m.feed_composition.keys()),
        initialize=1,
        units=pyunits.mol / pyunits.s,
    )

    m.sludge_composition = Var(
        list(m.feed_composition.keys()),
        initialize=1,
        units=pyunits.mol / pyunits.s,
    )
    m.precipitator_temperature = Var(
        initialize=273.15 + 50, bounds=(273.15, 273.15 + 99), units=pyunits.K
    )
    m.precipitator_temperature.fix()
    m.Q_heating = Var(initialize=0, units=pyunits.J / pyunits.s)

    # we only need enthalpy - can also request output pH, and pass it to precipitator
    # but dont need to - assume the temperature is not impacting pH
    m.feed_properties = Var(
        [
            ("molarEnthalpy", None),
            ("vaporPressure", "H2O(g)"),
        ],
        initialize=1,
    )
    m.precipitation_properties = Var(
        [
            ("speciesAmount", "Calcite"),
            ("speciesAmount", "Anhydrite"),
            ("molarEnthalpy", None),
            ("pH", None),
            ("vaporPressure", "H2O(g)"),
        ],
        initialize=1e-5,
    )
    m.treated_properties = Var(
        [
            ("molarEnthalpy", None),
        ],
        initialize=1,
    )

    reactant_dict = {
        "Ca": (1, "Calcite"),
        "HCO3": (1, "Calcite"),
        "Ca": (1, "Anhydrite"),
        "SO4": (1, "Anhydrite"),
    }
    m.eq_Q_heating = Constraint(
        expr=m.Q_heating
        == (
            m.precipitation_properties[("molarEnthalpy", None)]
            * sum([obj for key, obj in m.precipitator_composition.items()])
            - m.feed_properties[("molarEnthalpy", None)]
            * sum([obj for key, obj in m.feed_composition.items()])
        )
    )

    # connecting feed to precipitator composition
    @m.Constraint(list(m.feed_composition.keys()))
    def eq_precipitator_composition(fs, key):
        if key == "H2O":
            return m.precipitator_composition["H2O"] == m.feed_composition["H2O"]
        else:
            return m.precipitator_composition[key] == m.feed_composition[key]

    @m.Constraint(list(m.feed_composition.keys()))
    def eq_treated_composition(fs, key):
        return (
            m.precipitator_composition[key] - m.sludge_composition[key]
            == m.treated_composition[key]
        )

    # calculate sludge composition - we are not tracking solids that form, but rather
    # apparat species that would make them in addition those present in aqueous phase of the
    # sludge.

    @m.Constraint(list(m.feed_composition.keys()))
    def eq_sludge_composition(fs, key):
        if key == "H2O":
            # water flow is percent of total solids
            solid_mass = [
                m.precipitation_properties[("speciesAmount", "Calcite")],
                m.precipitation_properties[("speciesAmount", "Anhydrite")],
            ]
            return (
                m.precipitator_composition["H2O"]
                * m.sludge_water_content
                * sum(solid_mass)
                == m.sludge_composition["H2O"]
            )
        elif key in reactant_dict:
            return (
                m.precipitation_properties[("speciesAmount", reactant_dict[key][1])]
                * reactant_dict[key][0]
                + m.sludge_composition["H2O"] * m.treated_composition[key]
                == m.sludge_composition[key] * m.treated_composition["H2O"]
            )
        else:
            return (
                m.sludge_composition["H2O"] * m.treated_composition[key]
                == m.sludge_composition[key] * m.treated_composition["H2O"]
            )

    # We have to use super critical database to get enthalpy data as
    # our default PhreeqC data base with pitzer data file does not contain
    # enthalpy information - please refer to reaktoro documentation on supported data bases

    # we also need to define an ion translation dicionary for this data base
    # as default translator only support PhreeqCdata base with pitzer data file notation
    #  - once again refer to reaktoro documentation on specific data base and species to define translation dictionary
    #  - this dict should connect the name of species you are supplying to name of species in the databases
    #  in the data base file

    translation_dict = {
        "H2O": "H2O(aq)",
        "Mg": "Mg+2",
        "Na": "Na+",
        "Cl": "Cl-",
        "SO4": "SO4-2",
        "Ca": "Ca+2",
        "HCO3": "HCO3-",
    }
    # We will exclude species here that are present in zero concentration and we do not expect to participate
    exclude_species_list = [
        "S2O5-2",
        "S2O6-2",
        "S2O8-2",
        "S5O6-2",
        "HO2-",
        "HClO2(aq)",
        "HClO(aq)",
        "H2S2O4(aq)",
        "H2O2(aq)",
        "ClO-",
        "ClO2-",
        "ClO3-",
        "ClO4-",
    ]
    # note how we included nitrogen as one of gas species, this will prevent
    # PengRobinson EOS from forcing all of the water into vapor phase (refer to NOTE above)"""

    if parallel_mode:
        m.parallel_block_manager = ReaktoroBlockManager()
    else:
        m.parallel_block_manager = None
    if hess_type is None:
        hess_options = {}
    else:
        hess_options = {"hessian_type": hess_type}
    m.eq_feed_properties = ReaktoroBlock(
        system_state={
            "temperature": m.feed_temperature,
            "pressure": m.feed_pressure,
            "pH": m.feed_pH,
        },
        aqueous_phase={
            "composition": m.feed_composition,
            "convert_to_rkt_species": True,
            "species_to_rkt_species_dict": translation_dict,
            "activity_model": "ActivityModelPitzer",
        },
        outputs=m.feed_properties,
        mineral_phase={"phase_components": ["Calcite", "Anhydrite"]},
        gas_phase={
            "phase_components": ["H2O(g)", "N2(g)"],
            "activity_model": "ActivityModelRedlichKwong",
        },
        database="SupcrtDatabase",  # need to specify new data base to use
        database_file="supcrtbl",
        dissolve_species_in_reaktoro=True,
        exclude_species_list=exclude_species_list,
        reaktoro_block_manager=m.parallel_block_manager,
        build_speciation_block=False,
        assert_charge_neutrality_on_property_block=True,
        hessian_options=hess_options,
    )

    # """ need to get precipitator enthalpy to find required power input """
    m.eq_precipitation_properties = ReaktoroBlock(
        system_state={
            "temperature": m.precipitator_temperature,
            "pressure": m.feed_pressure,
            "pH": m.feed_pH,
        },
        aqueous_phase={
            "composition": m.precipitator_composition,
            "convert_to_rkt_species": True,
            "species_to_rkt_species_dict": translation_dict,
            "activity_model": "ActivityModelPitzer",
        },
        outputs=m.precipitation_properties,
        mineral_phase={"phase_components": ["Calcite", "Anhydrite"]},
        gas_phase={
            "phase_components": ["H2O(g)", "N2(g)"],
            "activity_model": "ActivityModelRedlichKwong",
        },
        database="SupcrtDatabase",  # need to specify new data base to use
        database_file="supcrtbl",
        dissolve_species_in_reaktoro=True,
        build_speciation_block=True,
        exclude_species_list=exclude_species_list,
        reaktoro_block_manager=m.parallel_block_manager,
        assert_charge_neutrality_on_property_block=True,
        hessian_options=hess_options,
    )

    m.eq_treated_properties = ReaktoroBlock(
        system_state={
            "temperature": m.precipitator_temperature,
            "pressure": m.feed_pressure,
            "pH": m.precipitation_properties[("pH", None)],
        },
        aqueous_phase={
            "composition": m.treated_composition,
            "convert_to_rkt_species": True,
            "species_to_rkt_species_dict": translation_dict,
            "activity_model": "ActivityModelPitzer",
        },
        outputs=m.treated_properties,
        mineral_phase={"phase_components": ["Calcite", "Anhydrite"]},
        gas_phase={
            "phase_components": ["H2O(g)", "N2(g)"],
            "activity_model": "ActivityModelRedlichKwong",
        },
        database="SupcrtDatabase",  # need to specify new data base to use
        database_file="supcrtbl",
        dissolve_species_in_reaktoro=True,
        exclude_species_list=exclude_species_list,
        reaktoro_block_manager=m.parallel_block_manager,
        assert_charge_neutrality_on_property_block=True,
        build_speciation_block=False,
        hessian_options=hess_options,
    )
    # assert False
    if parallel_mode:
        m.parallel_block_manager.build_reaktoro_blocks()
    scale_model(m)
    return m


def scale_model(m):
    for key in m.feed_composition:
        iscale.set_scaling_factor(
            m.feed_composition[key], 1 / m.feed_composition[key].value
        )
        iscale.set_scaling_factor(
            m.precipitator_composition[key], 1 / m.feed_composition[key].value
        )
        iscale.set_scaling_factor(
            m.sludge_composition[key], 1 / m.feed_composition[key].value
        )
        iscale.set_scaling_factor(
            m.treated_composition[key], 1 / m.feed_composition[key].value
        )
        iscale.constraint_scaling_transform(
            m.eq_sludge_composition[key], 1 / m.feed_composition[key].value
        )
        iscale.constraint_scaling_transform(
            m.eq_treated_composition[key], 1 / m.feed_composition[key].value
        )
        iscale.constraint_scaling_transform(
            m.eq_precipitator_composition[key], 1 / m.feed_composition[key].value
        )
    iscale.set_scaling_factor(m.feed_pH, 1)
    iscale.set_scaling_factor(
        m.precipitation_properties[("speciesAmount", "Calcite")], 1e3
    )
    iscale.set_scaling_factor(
        m.precipitation_properties[("speciesAmount", "Anhydrite")], 1e3
    )

    me_scale = 1e-2

    iscale.set_scaling_factor(m.feed_properties[("molarEnthalpy", None)], me_scale)
    iscale.set_scaling_factor(
        m.precipitation_properties[("molarEnthalpy", None)], me_scale
    )
    iscale.set_scaling_factor(m.treated_properties[("molarEnthalpy", None)], me_scale)
    iscale.set_scaling_factor(m.feed_temperature, 1 / 273)
    iscale.set_scaling_factor(m.precipitator_temperature, 1 / 273)
    iscale.set_scaling_factor(m.Q_heating, 1e-3)
    iscale.constraint_scaling_transform(m.eq_Q_heating, 1e-3)


def initialize(m):
    # propagate feed to precipitation comp
    for key in m.eq_precipitator_composition:
        calculate_variable_from_constraint(
            m.precipitator_composition[key], m.eq_precipitator_composition[key]
        )
    # initialize feed and precipitation properties
    # This will also get us initial precipitation amounts"""
    m.eq_feed_properties.initialize()
    m.eq_feed_properties.display_jacobian_scaling()
    m.eq_precipitation_properties.initialize()
    m.eq_precipitation_properties.display_jacobian_scaling()
    # get sludge flow volume first
    calculate_variable_from_constraint(
        m.sludge_composition["H2O"], m.eq_sludge_composition["H2O"]
    )

    # we wrote lazy formulation and cant explicitly calculate
    # treated or sludge ion composition, lets for initialization assume
    # that treated comp is same as feed composition and use that to estimate
    # initial sludge comp
    for key, obj in m.feed_composition.items():
        if key == "H2O":
            calculate_variable_from_constraint(
                m.treated_composition[key], m.eq_treated_composition[key]
            )
        else:
            m.treated_composition[key].value = obj.value

    m.eq_treated_properties.initialize()
    solve(m)


def solve(
    m,
):
    cy_solver = get_cyipopt_watertap_solver(max_iter=50)
    result = cy_solver.solve(m, tee=True)
    display_results(m)
    assert_optimal_termination(result)
    return result


def display_results(m):
    print("result")
    m.eq_precipitation_properties.display_reaktoro_state()
    print(
        f"Feed temp {m.feed_temperature.value-273.15}, precipitator temp {m.precipitator_temperature.value-273.15}"
    )
    print(
        f"""Feed vapor pressure {m.feed_properties[("vaporPressure", "H2O(g)")].value} (Pa)"""
    )
    print(
        f"""Precipitator vapor pressure {m.precipitation_properties[("vaporPressure", "H2O(g)")].value} (Pa),"""
    )
    print(f"Q heating {m.Q_heating.value/1000} kJ/s")
    print(
        f"Specific heat input {m.Q_heating.value/1000/(m.precipitator_temperature.value-m.feed_temperature.value)/(55*18.015/1000)} kJ/K/kg"
    )

    print(
        f'Calcite precipitation {m.precipitation_properties[("speciesAmount", "Calcite")].value} mol/s'
    )
    print(f'precipitator pH {m.precipitation_properties[("pH", None)].value}')


if __name__ == "__main__":
    main()
