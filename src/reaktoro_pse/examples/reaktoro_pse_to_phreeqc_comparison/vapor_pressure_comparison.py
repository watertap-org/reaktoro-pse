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
import reaktoro as rkt
from reaktoro_pse.reaktoro_block import ReaktoroBlock
from pyomo.environ import (
    Var,
)

__author__ = "Alexander Dudchenko"

"""
This examples compares reaktoro_pse implementation to phreeqcinwt for calculation of vapor pressure at different temperatures. 

NOTE: For water vapor calculations, pay attention to speciation and assumptions. Please
refer to these two discussions:
https://github.com/reaktoro/reaktoro/discussions/398
https://github.com/reaktoro/reaktoro/discussions/285

Key assumptions:
The calculation is down with out fixed volume, or pressure and as such an equilibrium will find a pressure and volume. 
No solids form due to heating. 

General note:
Its strongly recommended that user explores their equilibrium problem in reaktoro directly for any vapor calculations 
to understand the changes the system would under go before using reaktoro-pse.
"""


def main(save_fig=False, show_fig=True):
    phreeqc_config = compUtils.get_phreeqc_data()

    m = standardModel.build_modification_example(phreeqc_config["feed_comp"])
    add_vapor_pressure_properties(m)
    standardModel.initialize(m)
    reaktoro_output_dict = {}
    reaktoro_output_dict["temperature_sweep"] = {}
    for t in phreeqc_config["temperature"]:
        m.feed_temperature.fix(273.15 + t)
        standardModel.solve(m)
        compUtils.get_reaktoro_solved_outputs(
            m, reaktoro_output_dict["temperature_sweep"]
        )
    errors = compUtils.plot_data_sets(
        phreeqc_config["temperature"],
        phreeqc_config,
        reaktoro_output_dict,
        "temperature_sweep",
        "Temperature (C)",
        show_fig=show_fig,
        save_fig=save_fig,
    )
    print(errors)
    return errors


def add_vapor_pressure_properties(m):
    """getting vapor pressure, we
    the system presure is unkown,(as we don't know vapor pressure),
    and we are nod modifying state"""
    m.modified_properties = Var(
        [
            # ("scalingTendency", "Calcite"),
            # ("scalingTendency", "Gypsum"),
            ("vaporPressure", "H2O(g)"),
        ],
        initialize=1,
    )
    m.feed_pressure.fix(1e5)
    """ note how we included nitrogen as one of gas species, this will prevent 
    PengRobinson EOS from forcing all of the water into vapor phase (refer to NOTE above)"""
    m.eq_modified_properties = ReaktoroBlock(
        composition=m.feed_composition,
        temperature=m.feed_temperature,
        pressure=m.feed_pressure,
        pH=m.feed_pH,
        outputs=m.modified_properties,
        aqueous_phase_activity_model=rkt.ActivityModelPitzer(),
        gas_phase=["H2O(g)", "Ntg(g)"],
        gas_phase_activity_model=rkt.ActivityModelPengRobinsonPhreeqc(),
        dissolve_species_in_reaktoro=True,
        # we can use default converter as its defined for default database (Phreeqc and pitzer)
        convert_to_rkt_species=True,
        # we are modifying state and must speciate inputs before adding acid to find final prop state.
        build_speciation_block=True,
        jacobian_user_scaling={
            ("speciesActivityLn", "H2O(g)"): 1,
        },
    )


if __name__ == "__main__":
    main()
