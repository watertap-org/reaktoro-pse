import reaktoro_pse.examples.reaktoro_pse_to_phreeqc_comparison.comparison_utils as compUtils
import reaktoro_pse.examples.reaktoro_pse_to_phreeqc_comparison.standard_model as standardModel

__author__ = "Alexander Dudchenko"


"""
This examples compares reaktoro_pse implementation to phreeqcinwt calculation of solution
state adding acid or base. 

Key assumptions:

Uses 100wt% HCl and 100wt% NaOH
No solids form

"""


def main(save_fig=False, show_fig=True):
    phreeqc_config = compUtils.get_phreeqc_data()
    m = standardModel.build_modification_example(phreeqc_config["feed_comp"])
    standardModel.add_standard_properties(m)
    standardModel.initialize(m)
    reaktoro_output_dict = {}
    reaktoro_output_dict["HCl_sweep"] = {}
    m.water_recovery.fix(0)
    for wr in phreeqc_config["HCl"]:

        m.acid_addition.fix(wr)
        standardModel.solve(m)
        compUtils.get_reaktoro_solved_outputs(m, reaktoro_output_dict["HCl_sweep"])
    errors_hcl = compUtils.plot_data_sets(
        phreeqc_config["HCl"],
        phreeqc_config,
        reaktoro_output_dict,
        "HCl_sweep",
        "HCl dose (mol)",
        show_fig=show_fig,
        save_fig=save_fig,
    )
    reaktoro_output_dict["NaOH_sweep"] = {}
    m.acid_addition.fix(0)
    m.water_recovery.fix(0)
    for wr in phreeqc_config["NaOH"]:
        m.base_addition.fix(wr)
        standardModel.solve(m)
        compUtils.get_reaktoro_solved_outputs(m, reaktoro_output_dict["NaOH_sweep"])
    errors_naoh = compUtils.plot_data_sets(
        phreeqc_config["NaOH"],
        phreeqc_config,
        reaktoro_output_dict,
        "NaOH_sweep",
        "NaOH dose (mol)",
        show_fig=show_fig,
        save_fig=save_fig,
    )
    return errors_hcl, errors_naoh


if __name__ == "__main__":
    main()
