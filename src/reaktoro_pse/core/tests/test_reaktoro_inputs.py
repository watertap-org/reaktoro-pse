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
import reaktoro_pse.core.reaktoro_inputs as RktInputspec

from pyomo.environ import Var, units as pyunits
from reaktoro_pse.core.tests.test_reaktoro_state import (
    build_rkt_state_with_species,
    build_rkt_state_with_species_no_ph,
)

import pickle

__author__ = "Alexander V. Dudchenko (SLAC)"


def test_with_rkt_sum(build_rkt_state_with_species):
    """testing if we can construct all constraints correctly if summing with
    rkt, so inputs should be apprant species, with single empty constraint for CO3-2"""
    m, rkt_state = build_rkt_state_with_species
    rkt_state.build_state()
    rkt_state.equilibrate_state()
    rkt_input = RktInputspec.ReaktoroInputSpec(rkt_state)
    rkt_input.configure_specs(dissolve_species_in_rkt=True)
    rkt_input.build_input_specs()
    input_names = rkt_input.equilibrium_specs.namesControlVariables()
    print(input_names)
    input_constraints = rkt_input.equilibrium_specs.namesConstraints()
    print(input_constraints)
    expected_names = [
        "[Cl]",
        "[C]",
        "[Na]",
        "[Mg]",
        "[Ca]",
        "[H2O]",
        "[O]",
        "[H]",
        "[H+]",
    ]

    expected_constraints = [
        "charge",
        "C_constraint",
        "Na_constraint",
        "Mg_constraint",
        "Ca_constraint",
        "H2O_constraint",
        "O_dummy_constraint",
        "H_dummy_constraint",
        "pH",
    ]

    for en in expected_names:
        assert en in input_names
    for ec in expected_constraints:
        assert ec in input_constraints
    assert len(input_names) == len(expected_names)
    assert len(expected_constraints) == len(input_constraints)
    assert len(input_names) == len(input_constraints)
    expected_inputs = [
        "temperature",
        "pressure",
        "pH",
        "CO3-2",
        "CO2",
        "H2O",
        "Na+",
        "Mg+2",
        "Ca+2",
    ]
    for ei in expected_inputs:
        assert ei in rkt_input.rkt_inputs.keys()
        assert ei in rkt_input.user_inputs.keys()
    rkt_expected_inputs = [
        "inputCO3-2",
        "inputCO2",
        "H2O",
        "inputNa+",
        "inputMg+2",
        "inputCa+2",
    ]
    for ei in rkt_expected_inputs:
        key = ei.replace("input", "")
        assert ei == rkt_input.rkt_inputs[key].get_rkt_input_name()
    assert len(rkt_input.rkt_inputs.keys()) == len(expected_inputs)
    expected_con_dict = {
        "C": [(1.0, "CO3-2"), (1.0, "CO2")],
        "Na": [(1, "Na+")],
        "Ca": [(1, "Ca+2")],
        "Mg": [(1, "Mg+2")],
    }
    expected_active_species = ["CO3-2", "CO2", "Na+", "Mg+2", "Ca+2"]
    for ion, ecd in expected_con_dict.items():
        rkt_ecd = rkt_input.constraint_dict[ion]
        # order might change....
        expected_ecds = len(ecd)
        counted_ecds = 0
        for i, (mol, ion) in enumerate(ecd):
            if mol == rkt_ecd[i][0] and ion == rkt_ecd[i][1]:
                counted_ecds += 1
        assert counted_ecds == expected_ecds
    for eas in expected_active_species:
        assert eas in rkt_input.active_species
    assert len(rkt_input.active_species) == len(expected_active_species)
    assert len(rkt_input.constraint_dict) == len(expected_con_dict)

    """ check charge neutrality enabled"""
    assert rkt_input.assert_charge_neutrality == True
    assert rkt_input.neutrality_ion == "Cl"


def test_with_rkt_sum_no_ph(build_rkt_state_with_species_no_ph):
    """testing if we can construct all constraints correctly if summing with
    rkt, so inputs should be apprant species, with single empty constraint for CO3-2"""
    m, rkt_state = build_rkt_state_with_species_no_ph
    rkt_state.build_state()
    rkt_state.equilibrate_state()
    rkt_input = RktInputspec.ReaktoroInputSpec(rkt_state)
    # lime = Var(initialize=0.01, units=pyunits.mol / pyunits.s)
    # lime.construct()

    # rkt_input.register_chemistry_modifier("CaO", lime)
    rkt_input.register_charge_neutrality(False)
    rkt_input.configure_specs(dissolve_species_in_rkt=True)
    rkt_input.build_input_specs()
    input_names = rkt_input.equilibrium_specs.namesControlVariables()
    input_constraints = rkt_input.equilibrium_specs.namesConstraints()
    print("inputtn", input_names)
    print("input_c", input_constraints)
    expected_names = ["[C]", "[Na]", "[Mg]", "[Cl]", "[Ca]", "[H2O]", "[H]", "[O]"]
    expected_constraints = [
        "C_constraint",
        "Na_constraint",
        "Mg_constraint",
        "Cl_constraint",
        "Ca_constraint",
        "H2O_constraint",
        "H_dummy_constraint",
        "O_dummy_constraint",
    ]

    for en in expected_names:
        assert en in input_names
    for ec in expected_constraints:
        assert ec in input_constraints
    assert len(input_names) == len(expected_names)
    assert len(expected_constraints) == len(input_constraints)
    assert len(input_names) == len(input_constraints)
    expected_inputs = [
        "temperature",
        "pressure",
        "H2O",
        "CO3-2",
        "Na+",
        "Mg+2",
        "Ca+2",
        "Cl-",
    ]
    for ei in expected_inputs:
        assert ei in rkt_input.rkt_inputs.keys()
        # assert ei in rkt_input.userInputs.keys()
    assert len(rkt_input.rkt_inputs.keys()) == len(expected_inputs)
    expected_con_dict = {
        "C": [(1.0, "CO3-2")],
        "Na": [(1, "Na+")],
        "Mg": [(1, "Mg+2")],
        "Ca": [(1, "Ca+2")],
        "Cl": [(1, "Cl-")],
    }
    expected_active_species = ["CO3-2", "Na+", "Mg+2", "Ca+2", "Cl-"]
    for ion, ecd in expected_con_dict.items():
        rkt_ecd = rkt_input.constraint_dict[ion]
        # order might change....
        expected_ecds = len(ecd)
        counted_ecds = 0
        for i, (mol, ion) in enumerate(ecd):
            if mol == rkt_ecd[i][0] and ion == rkt_ecd[i][1]:
                counted_ecds += 1
        assert counted_ecds == expected_ecds
    for eas in expected_active_species:
        assert eas in rkt_input.active_species
    assert len(rkt_input.active_species) == len(expected_active_species)
    assert len(rkt_input.constraint_dict) == len(expected_con_dict)

    """ check charge neutrality enabled"""
    assert rkt_input.assert_charge_neutrality == False
    assert rkt_input.neutrality_ion == "Cl"


def test_with_pyomo_sum(build_rkt_state_with_species):
    """testing if we can construct all constraints correctly but summing in pyomo
    so rkt inputs are elment species only"""
    m, rkt_state = build_rkt_state_with_species
    rkt_state.build_state()
    rkt_state.equilibrate_state()
    rkt_input = RktInputspec.ReaktoroInputSpec(rkt_state)
    rkt_input.configure_specs(dissolve_species_in_rkt=False)
    rkt_input.build_input_specs()
    input_names = rkt_input.equilibrium_specs.namesControlVariables()
    input_constraints = rkt_input.equilibrium_specs.namesConstraints()
    print(input_names)
    print(input_constraints)
    expected_names = [
        "[Cl]",
        "[C]",
        "[Na]",
        "[Mg]",
        "[Ca]",
        "[H2O]",
        "[H]",
        "[O]",
        "[H+]",
    ]
    expected_constraints = [
        "charge",
        "C_constraint",
        "Na_constraint",
        "Mg_constraint",
        "Ca_constraint",
        "H2O_constraint",
        "H_dummy_constraint",
        "O_dummy_constraint",
        "pH",
    ]

    for en in expected_names:
        assert en in input_names
    for ec in expected_constraints:
        assert ec in input_constraints
    assert len(input_names) == len(expected_names)
    assert len(expected_constraints) == len(input_constraints)
    assert len(input_names) == len(input_constraints)
    expected_inputs = [
        "H2O",
        "Mg+2",
        "Mg",
        "Na+",
        "Na",
        "Cl-",
        "Cl",
        "Ca+2",
        "Ca",
        "CO3-2",
        "HCO3",
        "CO2",
        "temperature",
        "pressure",
        "pH",
    ]
    for ei in expected_inputs:
        assert ei in rkt_input.user_inputs.keys()
    assert len(rkt_input.user_inputs.keys()) == len(expected_inputs)
    expected_inputs = ["temperature", "pressure", "pH", "C", "H2O", "Na", "Mg", "Ca"]
    for ei in expected_inputs:
        assert ei in rkt_input.rkt_inputs.keys()
    assert len(rkt_input.rkt_inputs.keys()) == len(expected_inputs)
    expected_con_dict = {
        "C": [(1.0, "CO3-2"), (1.0, "CO2")],
        "Na": [(1, "Na+")],
        "Ca": [(1, "Ca+2")],
        "Mg": [(1, "Mg+2")],
    }
    expected_active_species = ["CO3-2", "CO2", "Na+", "Mg+2", "Ca+2"]
    for ion, ecd in expected_con_dict.items():
        rkt_ecd = rkt_input.constraint_dict[ion]
        # order might change....
        expected_ecds = len(ecd)
        counted_ecds = 0
        for i, (mol, ion) in enumerate(ecd):
            if mol == rkt_ecd[i][0] and ion == rkt_ecd[i][1]:
                counted_ecds += 1
        assert counted_ecds == expected_ecds
    for eas in expected_active_species:
        assert eas in rkt_input.active_species
    assert len(rkt_input.active_species) == len(expected_active_species)
    assert len(rkt_input.constraint_dict) == len(expected_con_dict)


def test_with_chemical(build_rkt_state_with_species):
    """testing if we can add chemicals"""
    m, rkt_state = build_rkt_state_with_species
    rkt_state.build_state()
    rkt_state.equilibrate_state()
    rkt_input = RktInputspec.ReaktoroInputSpec(rkt_state)

    lime = Var(initialize=0.01, units=pyunits.mol / pyunits.s)
    lime.construct()

    rkt_input.register_chemistry_modifier("CaO", lime)
    rkt_input.configure_specs(dissolve_species_in_rkt=True)
    rkt_input.build_input_specs()
    expected_con_dict = {
        "C": [(1.0, "CO3-2"), (1.0, "CO2")],
        "Na": [(1, "Na+")],
        "Ca": [(1, "Ca+2"), (1, "CaO")],
        "Mg": [(1, "Mg+2")],
    }
    expected_active_species = ["CaO", "CO3-2", "CO2", "Na+", "Mg+2", "Ca+2"]
    print(rkt_input.constraint_dict)
    for ion, ecd in expected_con_dict.items():
        rkt_ecd = rkt_input.constraint_dict[ion]
        # order might change....
        expected_ecds = len(ecd)
        counted_ecds = 0
        for i, (mol, ion) in enumerate(ecd):
            if mol == rkt_ecd[i][0] and ion == rkt_ecd[i][1]:
                counted_ecds += 1
        assert counted_ecds == expected_ecds
    for eas in expected_active_species:
        assert eas in rkt_input.active_species
    assert len(rkt_input.active_species) == len(expected_active_species)
    assert len(rkt_input.constraint_dict) == len(expected_con_dict)


def test_input_with_pickle_copy(build_rkt_state_with_species):
    """testing if we can add chemicals"""
    m, rkt_state = build_rkt_state_with_species
    rkt_state.build_state()
    rkt_state.equilibrate_state()
    rkt_input = RktInputspec.ReaktoroInputSpec(rkt_state)

    lime = Var(initialize=0.01, units=pyunits.mol / pyunits.s)
    lime.construct()

    rkt_input.register_chemistry_modifier("CaO", lime)
    rkt_input.configure_specs(dissolve_species_in_rkt=True)
    rkt_input.build_input_specs()
    export_object = rkt_input.export_config()
    pickled_object = pickle.dumps(export_object)
    unpickled_object = pickle.loads(pickled_object)

    new_rkt_input = RktInputspec.ReaktoroInputSpec(rkt_state)
    new_rkt_input.load_from_export_object(unpickled_object)
    new_rkt_input.build_input_specs()
    expected_con_dict = {
        "C": [(1.0, "CO3-2"), (1.0, "CO2")],
        "Na": [(1, "Na+")],
        "Ca": [(1, "Ca+2"), (1, "CaO")],
        "Mg": [(1, "Mg+2")],
    }
    expected_active_species = ["CaO", "CO3-2", "CO2", "Na+", "Mg+2", "Ca+2"]
    print(new_rkt_input.constraint_dict)
    for ion, ecd in expected_con_dict.items():
        rkt_ecd = new_rkt_input.constraint_dict[ion]
        # order might change....
        expected_ecds = len(ecd)
        counted_ecds = 0
        for i, (mol, ion) in enumerate(ecd):
            if mol == rkt_ecd[i][0] and ion == rkt_ecd[i][1]:
                counted_ecds += 1
        assert counted_ecds == expected_ecds
    for eas in expected_active_species:
        assert eas in new_rkt_input.active_species
    assert len(new_rkt_input.active_species) == len(expected_active_species)
    assert len(new_rkt_input.constraint_dict) == len(expected_con_dict)
