#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/reaktoro-pse/"
#################################################################################
import numpy as np
import pytest
import idaes.core.util.scaling as iscale
from reaktoro_pse.examples import (
    biogas_combustion,
    evaporative_crystallizer,
    simple_desalination,
    simple_ion_exchange,
    thermal_precipitation,
    redox_tracking,
    oxygen_tracking,
)


def test_oxygen_tracking():
    results = oxygen_tracking.main()
    sealed, aerated = results[False], results[True]
    # sealed softening, incoming DO (3e-7 mol/s) oxidizes 2x as much Mn (1e-5 mol/s feed)
    assert pytest.approx(sealed["Mn removed (%)"], 1e-2) == 6.0
    assert sealed["softened DO (mg/L)"] < 1e-20
    # aerated softening saturates water with O2 (~9.4 mg/L at 20 C) and removes all Mn
    assert pytest.approx(aerated["softened DO (mg/L)"], 1e-2) == 9.382
    assert pytest.approx(aerated["Mn removed (%)"], 1e-3) == 100
    # DO is carried into neutralization
    assert pytest.approx(aerated["final DO (mg/L)"], 1e-2) == 9.382
    assert pytest.approx(aerated["final pH"], 1e-3) == 8.5
    assert pytest.approx(aerated["acid dose (mol/s)"], 1e-2) == 5.21e-4
    assert pytest.approx(aerated["final pE"], 1e-2) == 12.56


def test_redox_couples():
    results_array = redox_tracking.main()
    expected_result = {
        "inputs": {
            "pH": [7.315729083913665, 7.315729083913665],
            "acid_addition": [1e-05, 0.01],
            "pE": [4.0, 4.0],
            "water_removal": [10.0, 10.0],
        },
        # we do not test no_pe as it will change besed on OS due to pE being unbounded!
        # "no_pe": {
        #     ("pH", None): [7.297370329017587, 7.297370329017587],
        #     ("pE", None): [4.082175889008317, 4.082175889008317],
        #     ("scalingTendency", "Calcite"): [0.474175277894885, 0.474175277894885],
        #     ("scalingTendency", "Gypsum"): [0.03851168496475971, 0.03851168496475971],
        # },
        "with_pe": {
            ("pH", None): [7.297370328995919, 7.297370328995919],
            ("pE", None): [4.082175888930677, 4.082175888930677],
            ("scalingTendency", "Calcite"): [0.4741752782085243, 0.4741752782085243],
            ("scalingTendency", "Gypsum"): [0.03851168498111663, 0.03851168498111663],
        },
    }

    print(results_array)
    for key, item in expected_result.items():
        for subkey, subitem in item.items():
            for i, v in enumerate(subitem):
                assert pytest.approx(results_array[key][subkey][i], 1e-2) == v


@pytest.mark.parametrize(
    "hess_type",
    [
        "ZeroHessian",
        "GaussNewton",
        "LBFGS",
        "BFGS",
        # "BFGS_mod", does not work on this example
        "BFGS_damp",
        "BFGS_ipopt",
    ],
)
def test_desal(hess_type):
    m = simple_desalination.main(hess_type=hess_type)

    assert (
        pytest.approx(m.desal_properties[("scalingTendency", "Gypsum")].value, 1e-3)
        == 0.604051223942643
    )
    assert (
        pytest.approx(m.desal_properties[("osmoticPressure", "H2O")].value, 1e-1)
        == 1548396.415543
    )

    assert pytest.approx(m.desal_properties[("pH", None)].value, 1e-3) == 6.284055
    assert pytest.approx(m.water_recovery.value, 1e-3) == 0.899999
    assert pytest.approx(m.acid_addition.value, 1e-3) == 0.003043


@pytest.mark.parametrize("hess_type", ["SR1", "BFGS"])
@pytest.mark.parametrize("ipopt_multiplier_init", ["default", "least_squares"])
def test_hessian_reinitialized_on_new_solve(
    monkeypatch, hess_type, ipopt_multiplier_init
):
    """quasi-Newton matrices hold curvature of the previous solve, which is stale once the
    problem changes (e.g. a new stage of a sweep), so they are re-initialized at the start
    of every Ipopt solve. The trigger is the Hessian evaluation PyNumero does with zero
    duals when it builds the NLP (to get the sparsity pattern), so it does not depend on
    how Ipopt initializes its multipliers (constr_mult_init_max)"""
    from pyomo.contrib.pynumero.interfaces.pyomo_grey_box_nlp import (
        PyomoNLPWithGreyBoxBlocks,
    )
    from reaktoro_pse.core.util_classes.cyipopt_solver import (
        get_cyipopt_watertap_solver,
    )
    from pyomo.environ import assert_optimal_termination

    m = simple_desalination.build_simple_desal(hess_type)
    simple_desalination.initialize(m)
    simple_desalination.setup_optimization(m)
    calc = m.eq_desal_properties.reaktoro_model.get_external_model().hessian_calculator

    solve_number = [0]
    original_nlp_init = PyomoNLPWithGreyBoxBlocks.__init__

    def counting_nlp_init(self, *args, **kwargs):
        solve_number[0] += 1
        original_nlp_init(self, *args, **kwargs)

    monkeypatch.setattr(PyomoNLPWithGreyBoxBlocks, "__init__", counting_nlp_init)

    reinitializations = [0]
    original_initial_matrix = calc.get_initial_matrix

    def counting_initial_matrix(*args, **kwargs):
        reinitializations[0] += 1
        return original_initial_matrix(*args, **kwargs)

    monkeypatch.setattr(calc, "get_initial_matrix", counting_initial_matrix)

    # (solve number, duals non zero, matrices re-initialized during this evaluation)
    evaluations = []
    original_get_hessian = calc.get_hessian

    def recording_get_hessian(input_values, jacobian, dual_multipliers):
        before = reinitializations[0]
        result = original_get_hessian(input_values, jacobian, dual_multipliers)
        evaluations.append(
            (
                solve_number[0],
                bool(np.any(np.asarray(dual_multipliers) != 0)),
                reinitializations[0] > before,
            )
        )
        return result

    monkeypatch.setattr(calc, "get_hessian", recording_get_hessian)

    solver_options = None
    if ipopt_multiplier_init == "least_squares":
        # always start from least squares multiplier estimate (never reset to zero)
        solver_options = {"constr_mult_init_max": 1e20}
    solver = get_cyipopt_watertap_solver(solver_options=solver_options)
    # optimization, then re-solve with tighter scaling limit, warm started from optimum
    for calcite_limit in [1.0, 0.8]:
        m.desal_properties[("scalingTendency", "Calcite")].setub(calcite_limit)
        result = solver.solve(m)
        assert_optimal_termination(result)
    assert solve_number[0] == 2
    for solve in [1, 2]:
        solve_evaluations = [e for e in evaluations if e[0] == solve]
        weighted = [e for e in solve_evaluations if e[1]]
        assert len(weighted) > 0
        # matrices are rebuilt before the first evaluation that uses them (non zero duals)
        assert weighted[0][
            2
        ], f"solve {solve}: first weighted hessian used stale matrices"
        # after that, the solve keeps updating the same matrices
        assert not any(e[2] for e in weighted[1:])


@pytest.mark.parametrize(
    "hess_type",
    [
        "ZeroHessian",
        "GaussNewton",
        "LBFGS",
        "BFGS",
        "BFGS_mod",
        "BFGS_damp",
        "BFGS_ipopt",
    ],
)
def test_thermal_precipt(hess_type):
    m = thermal_precipitation.main(hess_type=hess_type)
    assert (
        pytest.approx(
            m.precipitation_properties[("speciesAmount", "Calcite")].value, 1e-2
        )
        == 8.38746210798396e-05
    )
    assert (
        pytest.approx(
            m.precipitation_properties[("vaporPressure", "H2O(g)")].value, 1e-3
        )
        == 12166.985206776764
    )
    assert (
        pytest.approx(m.precipitation_properties[("pH", None)].value, 1e-3)
        == 6.965357182573612
    )
    assert pytest.approx(m.precipitation_properties[("pE", None)].value, 1e-1) == -1.37
    assert pytest.approx(m.Q_heating.value, abs=4e4) == 123.99859118830717 * 1000
    assert pytest.approx(m.precipitator_temperature.value, 1e-3) == 273.15 + 50


@pytest.mark.parametrize(
    "hess_type",
    [
        "ZeroHessian",
        "GaussNewton",
        "LBFGS",
        "BFGS",
        "BFGS_mod",
        "BFGS_damp",
        "BFGS_ipopt",
    ],
)
def test_ion_exchange(hess_type):
    m = simple_ion_exchange.main(hess_type=hess_type)

    assert pytest.approx(m.removal_percent["Mg"].value, 1e-3) == -35.54130924283
    assert pytest.approx(m.removal_percent["Ca"].value, 1e-3) == -79.15299911033
    assert pytest.approx(m.treated_pH.value, 1e-3) == 13.374349619456911
    assert pytest.approx(m.base_addition.value, abs=1e-3) == 0.31967185413270405


@pytest.mark.parametrize(
    "hess_type",
    [
        "ZeroHessian",
        "GaussNewton",
        "LBFGS",
        "BFGS",
        "BFGS_mod",
        # "BFGS_damp", # does not work on this example
        "BFGS_ipopt",
    ],
)
def test_biogas(hess_type):
    m = biogas_combustion.main(hess_type=hess_type)

    assert pytest.approx(m.air_to_fuel_ratio.value, 1e-3) == 3.8751662012681587
    assert pytest.approx(m.exhaust_temperature.value, 1e-3) == 2000


def test_evaporative_crystallizer_needs_species_amount_offset():
    """absent minerals pin bounded solids vars at their bound without an offset, with the
    default offset they sit at a = offset / scaling factor away from it and the
    optimization converges.
    Without an offset, whether the optimization fails or only needs many more iterations
    depends on platform (e.g. 41 iterations on linux vs not converged in 300 on windows),
    so it is not asserted, the pinning that causes it is."""

    def scaled_absent_solids(m):
        """scaled distance from bound of solids vars for minerals absent in reaktoro"""
        raw = m.crystallizer.rkt_solver.output_specs.user_outputs
        return [
            var.value * iscale.get_scaling_factor(var)
            for key, var in m.solids.items()
            if raw[key].value < 1e-30
        ]

    # without offset, absent minerals are pinned at their bound of 0 (Ipopt keeps them
    # only ~1e-11 away in scaled units)
    m = evaporative_crystallizer.build_crystallizer(
        species_amount_zero_offset=0, species_amount_offset_type=None
    )
    evaporative_crystallizer.initialize(m)
    absent = scaled_absent_solids(m)
    assert len(absent) > len(m.solids) / 2
    assert max(absent) < 1e-8

    # with default offset, absent minerals sit at offset (1e-4 scaled) from their bound,
    # and optimization converges
    m = evaporative_crystallizer.build_crystallizer()
    evaporative_crystallizer.initialize(m)
    absent = scaled_absent_solids(m)
    assert len(absent) > len(m.solids) / 2
    assert pytest.approx(min(absent), 1e-2) == 1e-4
    evaporative_crystallizer.setup_optimization(m)
    assert evaporative_crystallizer.solve(m)
    assert pytest.approx(m.evaporation.value, 1e-3) == 0.9122
    raw = m.crystallizer.rkt_solver.output_specs.user_outputs
    present = {key[1] for key in m.solids if raw[key].value > 1e-12}
    assert present == {"Magnesite", "Gypsum", "Halite"}
    a = 1e-4 / iscale.get_scaling_factor(m.solids["speciesAmount", "Calcite"])
    # absent minerals sit at the offset, present minerals are shifted by ~a**2/(2n)
    assert pytest.approx(m.solids["speciesAmount", "Calcite"].value, 1e-3) == a
    halite = ("speciesAmount", "Halite")
    a_halite = 1e-4 / iscale.get_scaling_factor(m.solids[halite])
    n = raw[halite].value
    # (1% tolerance covers ipopt constraint violation on a ~1e-7 shift)
    assert pytest.approx(m.solids[halite].value - n, 1e-2) == a_halite**2 / (2 * n)
