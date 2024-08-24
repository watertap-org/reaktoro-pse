import reaktoro_pse.examples.reaktoro_pse_to_phreeqc_comparison.standard_model as standardModel
import reaktoro_pse.examples.reaktoro_pse_to_phreeqc_comparison.comparison_utils as compUtils
import reaktoro as rkt
from reaktoro_pse.reaktoro_block import ReaktoroBlock
from pyomo.environ import (
    Var,
)


__author__ = "Alexander Dudchenko"


"""
This examples compares reaktoro_pse implementation to phreeqcinwt calculation of precipitation when adding lime.

Key assumptions:
Lime is added as pure CaO (100wt%)
"""


def main():
    phreeqc_config = compUtils.get_phreeqc_data()
    m = standardModel.build_modification_example(phreeqc_config["feed_comp"])
    add_mineral_properties(m)
    standardModel.initialize(m)
    reaktoro_output_dict = {}
    reaktoro_output_dict["CaO_sweep"] = {}
    for lime in phreeqc_config["CaO"]:
        m.lime_addition.fix(lime)
        standardModel.solve(m)
        compUtils.get_reaktoro_solved_outputs(m, reaktoro_output_dict["CaO_sweep"])
        m.display()
    compUtils.plot_data_sets(
        phreeqc_config["CaO"],
        phreeqc_config,
        reaktoro_output_dict,
        "CaO_sweep",
        "Lime addition (mol)",
    )


def add_mineral_properties(m):
    """getting vapor pressure, we
    the system presure is unkown,(as we don't know vapor pressure),
    and we are nod modifying state"""
    m.modified_properties = Var(
        [
            ("speciesAmount", "Calcite"),
            ("speciesAmount", "Gypsum"),
            ("scalingTendency", "Calcite"),
            ("scalingTendency", "Gypsum"),
            ("pH", None),
        ],
        initialize=1,
    )
    m.feed_pressure.fix(1e5)  # fixing at 1 bar
    m.eq_modified_properties = ReaktoroBlock(
        composition=m.feed_composition,
        temperature=m.feed_temperature,
        pressure=m.feed_pressure,
        pH=m.feed_pH,
        outputs=m.modified_properties,
        aqueous_phase_activity_model=rkt.ActivityModelPitzer(),
        mineral_phases=["Calcite", "Gypsum"],
        chemical_addition={"CaO": m.lime_addition},
        dissolve_species_in_reaktoro=False,
        # we can use default converter as its defined for default database (Phreeqc and pitzer)
        convert_to_rkt_species=True,
        # we are modifying state and must speciate inputs before adding acid to find final prop state.
        build_speciation_block=True,
        jacobian_user_scaling={
            ("speciesAmount", "Calcite"): 1e5,
            ("speciesAmount", "Gypsum"): 1e5,
            ("pH", None): 1,
            ("saturationIndex", "Calcite"): 1,
            ("saturationIndex", "Gypsum"): 1,
        },
    )


if __name__ == "__main__":
    main()
