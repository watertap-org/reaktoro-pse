import reaktoro_pse.examples.reaktoro_pse_to_phreeqc_comparison.standard_model as standardModel
import reaktoro_pse.examples.reaktoro_pse_to_phreeqc_comparison.comparison_utils as compUtils


__author__ = "Alexander Dudchenko"


"""
This examples compares reaktoro_pse implementation to phreeqcinwt calculation of solution
state after removing a specified amount of water (imitating desalination, or evaporation process).

Key assumptions:
Removing water impacts pH (note in simple_desalination example its assumed removing water does not alter pH)
No solids form
"""


def main(save_fig=False, show_fig=True):
    phreeqc_config = compUtils.get_phreeqc_data()
    m = standardModel.build_modification_example(phreeqc_config["feed_comp"])
    standardModel.add_standard_properties(m)
    standardModel.initialize(m)
    reaktoro_output_dict = {}
    reaktoro_output_dict["water_removal_percent_sweep"] = {}
    for wr in phreeqc_config["water_removal_percent"]:
        m.water_recovery.fix(wr / 100)
        standardModel.solve(m)
        compUtils.get_reaktoro_solved_outputs(
            m, reaktoro_output_dict["water_removal_percent_sweep"]
        )
    errors = compUtils.plot_data_sets(
        phreeqc_config["water_removal_percent"],
        phreeqc_config,
        reaktoro_output_dict,
        "water_removal_percent_sweep",
        "Water evaporation (%)",
        show_fig=show_fig,
        save_fig=save_fig,
    )
    return errors


if __name__ == "__main__":
    main()
