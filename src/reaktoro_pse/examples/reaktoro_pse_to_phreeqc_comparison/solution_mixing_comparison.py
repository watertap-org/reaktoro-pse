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
import reaktoro_pse.examples.reaktoro_pse_to_phreeqc_comparison.standard_model as standardModel
import reaktoro_pse.examples.reaktoro_pse_to_phreeqc_comparison.comparison_utils as compUtils

from reaktoro_pse.reaktoro_block import ReaktoroBlock
from pyomo.environ import (
    ConcreteModel,
    Var,
    units as pyunits,
)
import idaes.core.util.scaling as iscale
import reaktoro as rkt

__author__ = "Alexander V. Dudchenko"


# This examples compares reaktoro_pse implementation to phreeqcinwt calculation mixing two
# solutions together based on mixing ratio - this is a common phenomena for recycling process.
# Unlike other examples, we will manually do speciation tracking for each feed, and pass exact speciation into
# reaktoro- the ReaktoroBlock does not automatically support mixing so lower level interaction is required.

# Key assumptions:
# Removing water impacts pH (note in simple_desalination example its assumed removing water does not alter pH)
# No solids form


def main(save_fig=False, show_fig=True):
    phreeqc_config = compUtils.get_phreeqc_data()
    m = build_modification_example(
        phreeqc_config["feed_comp"],
        phreeqc_config["feed_comp_2"],
        phreeqc_config["ph_2"],
    )
    add_standard_properties(m)
    initialize(m)
    m.display()
    # assert False
    reaktoro_output_dict = {}
    reaktoro_output_dict["mix_ratio_sweep"] = {}
    for wr in phreeqc_config["mix_ratio"]:
        m.ratio.fix(wr)
        standardModel.solve(m)
        compUtils.get_reaktoro_solved_outputs(
            m, reaktoro_output_dict["mix_ratio_sweep"]
        )
    errors = compUtils.plot_data_sets(
        phreeqc_config["mix_ratio"],
        phreeqc_config,
        reaktoro_output_dict,
        "mix_ratio_sweep",
        "Mixing ratio (-)",
        show_fig=show_fig,
        save_fig=save_fig,
    )
    print(errors)
    return errors


def main_wateq4f(save_fig=False, show_fig=True):
    phreeqc_config = compUtils.get_phreeqc_data(data_type="phreeqc_data_waterq4f.json")
    m = build_modification_example(
        phreeqc_config["feed_comp"],
        phreeqc_config["feed_comp_2"],
        phreeqc_config["ph_2"],
        phreeqc_config["pE_2"],
    )
    add_standard_properties(
        m, database="wateq4f.dat", activity_model="ActivityModelPhreeqc"
    )
    initialize(m)
    m.display()
    # assert False
    reaktoro_output_dict = {}
    reaktoro_output_dict["mix_ratio_sweep"] = {}
    for wr in phreeqc_config["mix_ratio"]:
        m.ratio.fix(wr)
        standardModel.solve(m)
        compUtils.get_reaktoro_solved_outputs(
            m, reaktoro_output_dict["mix_ratio_sweep"]
        )
    errors = compUtils.plot_data_sets(
        phreeqc_config["mix_ratio"],
        phreeqc_config,
        reaktoro_output_dict,
        "mix_ratio_sweep",
        "Mixing ratio (-)",
        show_fig=show_fig,
        save_fig=save_fig,
    )
    print(errors)
    return errors


def build_modification_example(water_comp, water_comp_2, pH2, pE2=None):
    m = ConcreteModel()
    m.feed_composition = Var(
        water_comp.keys(),
        initialize=1,
        units=pyunits.mol / pyunits.s,
    )
    for key, value in water_comp.items():
        m.feed_composition[key].fix(value)
    water_comp_2["H2O"] = water_comp["H2O"]
    print(water_comp_2)
    m.feed_temperature = Var(initialize=273.15 + 20, units=pyunits.K)
    m.feed_temperature.fix()
    m.feed_pressure = Var(
        initialize=10e5, units=pyunits.Pa
    )  # 10 bar used in phreeqc feed pressure
    m.feed_pressure.fix()
    # pressure.construct()
    m.feed_pH = Var(initialize=7.8, bounds=(4, 12), units=pyunits.dimensionless)
    m.feed_pH.fix(7.8)  # feed pH used in phreeqc sim
    if pE2 is not None:
        m.pE = Var(initialize=4, units=pyunits.dimensionless)
        m.pE.fix(4)
        m.pE_2 = Var(initialize=pE2, units=pyunits.dimensionless)
        m.pE_2.fix(pE2)
    m.feed_composition_2 = Var(
        water_comp_2.keys(),
        initialize=1,
        units=pyunits.mol / pyunits.s,
    )
    m.ratio = Var(initialize=1, units=pyunits.dimensionless)
    m.ratio.fix(1)
    for key, value in water_comp_2.items():
        m.feed_composition_2[key].fix(value)
    m.ratio_feed_composition_2 = Var(
        water_comp_2.keys(),
        initialize=1,
        units=pyunits.mol / pyunits.s,
    )
    for key, value in water_comp_2.items():
        m.ratio_feed_composition_2[key].value = value

    @m.Constraint(list(m.feed_composition.keys()))
    def eq_ratio_feed_composition_2(fs, key):
        return m.ratio_feed_composition_2[key] == m.feed_composition_2[key] * m.ratio

    m.feed_2_pH = Var(initialize=pH2, bounds=(4, 12), units=pyunits.dimensionless)
    m.feed_2_pH.fix(pH2)
    return m


def add_standard_properties(
    m,
    database="pitzer.dat",
    activity_model="ActivityModelPitzer",
):
    m.modified_properties = Var(
        [
            ("scalingTendency", "Calcite"),
            ("scalingTendency", "Gypsum"),
            ("pH", None),
            ("pE", None),
            ("osmoticPressure", "H2O"),
        ],
        initialize=1,
    )

    # building reaktoro blocks to speciate the two
    # different feeds, we can use {"speciesAmount", True} to get all possible
    # exact species as an output from reaktoro with out knowing what they
    # are apriori.
    if m.find_component("pE") is None:
        pE = None
        pE_2 = None
    else:
        pE_2 = m.pE_2
        pE = m.pE
    m.eq_feed_properties = ReaktoroBlock(
        aqueous_phase={
            "composition": m.feed_composition,
            "convert_to_rkt_species": True,
            "activity_model": activity_model,
        },
        database_file=database,
        system_state={
            "temperature": m.feed_temperature,
            "pressure": m.feed_pressure,
            "pH": m.feed_pH,
            "pE": pE,
        },
        outputs={"speciesAmount": True},  # get exact speciation for the feed
        dissolve_species_in_reaktoro=True,
        build_speciation_block=False,
        build_graybox_model=False,
    )
    m.eq_mixed_properties = ReaktoroBlock(
        aqueous_phase={
            "composition": m.ratio_feed_composition_2,
            "convert_to_rkt_species": True,
            "activity_model": activity_model,
        },
        database_file=database,
        system_state={
            "temperature": m.feed_temperature,
            "pressure": m.feed_pressure,
            "pH": m.feed_2_pH,
            "pE": pE_2,
        },
        outputs=m.modified_properties,  # get exact speciation for the feed
        build_speciation_block=True,
        external_speciation_reaktoro_blocks=[m.eq_feed_properties],
    )

    scale_model(m)


def scale_model(m):
    for key in m.feed_composition:
        iscale.set_scaling_factor(
            m.feed_composition[key], 1 / m.feed_composition[key].value
        )
        iscale.set_scaling_factor(
            m.feed_composition_2[key], 1 / m.feed_composition[key].value
        )
        iscale.set_scaling_factor(
            m.ratio_feed_composition_2[key], 1 / m.feed_composition[key].value
        )
    iscale.set_scaling_factor(m.ratio, 1)


def initialize(m):

    m.eq_mixed_properties.initialize()
    standardModel.solve(m)


if __name__ == "__main__":
    main_wateq4f()
