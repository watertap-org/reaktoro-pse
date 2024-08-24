import reaktoro_pse.examples.reaktoro_pse_to_phreeqc_comparison.standard_model as standardModel
import reaktoro_pse.examples.reaktoro_pse_to_phreeqc_comparison.comparison_utils as compUtils
import reaktoro as rkt
from reaktoro_pse.reaktoro_block import ReaktoroBlock
from pyomo.environ import (
    Var,
)

__author__ = "Alexander Dudchenko"

"""
This examples compares reaktoro_pse implementation to phreeqcinwt calculation of solution
state for calculation of vapor pressure at different temperatures. 

A key thing to note when calculating vapor with reaktoro is that when useing ActivityModelPengRobinsonPhreeqc it does not differentiate 
correctly between liquid and aqueous state at low vapor pressures, resulting in large portioning of 
water in "vapor phase". this leads to nonsensical estimates of scaling potential and possibly other aqueous properties.
But properties like Vapor pressure are reasonably estimated.
This in in general means, thats vapor pressure should be estimated separately from aqueous/solid properties. 

Please refer to discussion on this here: https://github.com/reaktoro/reaktoro/discussions/285

Key assumptions:
The calculation is down with out fixed volume, or pressure and as such an equilibrium will find a pressure and volume. 
No solids form due to heating. 

General note:
Its strongly recommended that user explores their equilibrium problem in reaktoro directly for any vapor calculations 
to understand the changes the system would under go before using reaktoro-pse.
"""


def main():
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
        m.display()
    compUtils.plot_data_sets(
        phreeqc_config["temperature"],
        phreeqc_config,
        reaktoro_output_dict,
        "temperature_sweep",
        "Temperature (C)",
    )


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
    m.feed_pressure.fix(1e5)  # fixing at 1 bar
    m.eq_modified_properties = ReaktoroBlock(
        composition=m.feed_composition,
        temperature=m.feed_temperature,
        pressure=m.feed_pressure,
        pH=m.feed_pH,
        outputs=m.modified_properties,
        aqueous_phase_activity_model=rkt.ActivityModelPitzer(),
        gas_phases=["H2O(g)"],
        gas_phase_activity_model=rkt.ActivityModelPengRobinsonPhreeqc(),
        dissolve_species_in_reaktoro=True,
        # we can use default converter as its defined for default database (Phreeqc and pitzer)
        convert_to_rkt_species=True,
        # we are modifying state and must speciate inputs before adding acid to find final prop state.
        build_speciation_block=True,
        jacobian_user_scaling={
            ("speciesActivityLn", "H2O(g)"): 1,
            # ("saturationIndex", "Calcite"): 1,
            # ("saturationIndex", "Gypsum"): 1,
        },
        # presolve_property_block=True,
        # presolve_tolerance=1e-3,
        # presolve_epsilon=1e-16,
    )


if __name__ == "__main__":
    main()
