import reaktoro_pse.examples.reaktoro_pse_to_phreeqc_comparison.standard_model as standardModel
import reaktoro_pse.examples.reaktoro_pse_to_phreeqc_comparison.comparison_utils as compUtils
import reaktoro as rkt
from reaktoro_pse.reaktoro_block import ReaktoroBlock
import idaes.core.util.scaling as iscale
from pyomo.environ import (
    Var,
)


__author__ = "Alexander Dudchenko"


"""
This examples compares reaktoro_pse implementation to phreeqcinwt calculation of precipitation when adding lime.

Key assumptions:
Lime is added as pure CaO (100wt%)
"""


def main(save_fig=False, show_fig=True):
    phreeqc_config = compUtils.get_phreeqc_data()
    m = standardModel.build_modification_example(phreeqc_config["feed_comp"])
    add_mineral_properties(m)
    m.display()
    # assert False
    standardModel.initialize(m)
    reaktoro_output_dict = {}
    reaktoro_output_dict["CaO_sweep"] = {}
    for lime in phreeqc_config["CaO"]:
        m.lime_addition.fix(lime)
        print(m.lime_addition.value)
        standardModel.solve(m)
        compUtils.get_reaktoro_solved_outputs(m, reaktoro_output_dict["CaO_sweep"])

    errors = compUtils.plot_data_sets(
        phreeqc_config["CaO"],
        phreeqc_config,
        reaktoro_output_dict,
        "CaO_sweep",
        "Lime addition (mol)",
        show_fig=show_fig,
        save_fig=save_fig,
    )
    return errors


def add_mineral_properties(m):
    """getting vapor pressure, we
    the system presure is unkown,(as we don't know vapor pressure),
    and we are nod modifying state"""
    m.modified_properties = Var(
        [
            ("speciesAmount", "Calcite"),
            ("scalingTendency", "Calcite"),
            ("pH", None),
        ],
        initialize=1,
    )
    iscale.set_scaling_factor(m.modified_properties[("speciesAmount", "Calcite")], 1e5)
    iscale.set_scaling_factor(m.modified_properties[("pH", None)], 1 / 10)
    # iscale.set_scaling_factor(m.modefied_properties[("speciesAmount", "Calcite")], 1e-5)
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
        presolve_property_block=True,
        jacobian_user_scaling={
            ("speciesAmount", "Calcite"): 1,
            ("speciesAmount", "Gypsum"): 1,
            ("pH", None): 0.1,
            ("saturationIndex", "Calcite"): 1,
            ("saturationIndex", "Gypsum"): 1,
        },
    )


if __name__ == "__main__":
    main()
