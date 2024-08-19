import pytest
import reaktoro_pse.core.reaktoro_state as rktState
import reaktoro_pse.core.reaktoro_inputs as rktInputSpec

from pyomo.environ import Var, units as pyunits

from reaktoro_pse.core.tests.test_reaktoro_state import (
    build_rkt_state_with_species,
    build_rkt_state_with_species_no_ph,
)

__author__ = "Alexander V. Dudchenko (SLAC)"


def test_with_rkt_sum(build_rkt_state_with_species):
    """testing if we can construct all constraints correctly if summing with
    rkt, so inputs should be apprant species, with single empty constraint for CO3-2"""
    m, rkt_state = build_rkt_state_with_species
    rkt_state.build_state()
    rkt_state.equilibrate_state()
    rkt_input = rktInputSpec.reaktoroInputSpec(rkt_state)
    rkt_input.configure_specs(dissolve_species_in_rkt=True)
    input_names = rkt_input.rktEquilibriumSpecs.namesControlVariables()
    print(input_names)
    input_constraints = rkt_input.rktEquilibriumSpecs.namesConstraints()
    print(input_constraints)
    expected_names = [
        "[Cl-]",
        "[C]",
        "[O]",
        "[Na]",
        "[Mg]",
        "[Ca]",
        "[H2O]",
        "[CO3-2]",
        "[CO2]",
        "[Ca+2]",
        "[HCO3-]",
        "[Mg+2]",
        "[MgCO3]",
        "[MgOH+]",
        "[Na+]",
        "[OH-]",
        "[H+]",
    ]
    expected_constraints = [
        "charge",
        "C_constraint",
        "O_constraint",
        "Na_constraint",
        "Mg_constraint",
        "Ca_constraint",
        "H2O_constraint",
        "CO3-2_constraint",
        "CO2_constraint",
        "Ca+2_constraint",
        "HCO3-_constraint",
        "Mg+2_constraint",
        "MgCO3_constraint",
        "MgOH+_constraint",
        "Na+_constraint",
        "OH-_constraint",
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
        "Na",
        "Mg",
        "Ca",
    ]
    for ei in expected_inputs:
        assert ei in rkt_input.rktInputs.keys()
        assert ei in rkt_input.userInputs.keys()
    rkt_expected_inputs = [
        "inputCO3-2",
        "inputCO2",
        "inputH2O",
        "inputNa",
        "inputMg",
        "inputCa",
    ]
    for ei in rkt_expected_inputs:
        key = ei.replace("input", "")
        assert ei == rkt_input.rktInputs[key].get_rkt_input_name()
    assert len(rkt_input.rktInputs.keys()) == len(expected_inputs)
    expected_con_dict = {
        "C": [(1.0, "CO3-2"), (1.0, "CO2")],
        "O": [(1.0, "H2O"), (3.0, "CO3-2"), (2.0, "CO2")],
        "Na": [(1, "Na")],
        "Ca": [(1, "Ca")],
        "Mg": [(1, "Mg")],
    }
    expected_active_species = ["CO3-2", "CO2", "H2O", "Na", "Mg", "Ca"]
    for ion, ecd in expected_con_dict.items():
        rkt_ecd = rkt_input.constraintDict[ion]
        # order might change....
        expected_ecds = len(ecd)
        counted_ecds = 0
        for i, (mol, ion) in enumerate(ecd):
            if mol == rkt_ecd[i][0] and ion == rkt_ecd[i][1]:
                counted_ecds += 1
        assert counted_ecds == expected_ecds
    for eas in expected_active_species:
        assert eas in rkt_input.rktActiveSpecies
    assert len(rkt_input.rktActiveSpecies) == len(expected_active_species)
    assert len(rkt_input.constraintDict) == len(expected_con_dict)

    """ check charge neutrality enabled"""
    assert rkt_input.assertChargeNeutrality == True
    assert rkt_input.neutralityIon == "Cl"


def test_with_rkt_sum_no_ph(build_rkt_state_with_species_no_ph):
    """testing if we can construct all constraints correctly if summing with
    rkt, so inputs should be apprant species, with single empty constraint for CO3-2"""
    m, rkt_state = build_rkt_state_with_species_no_ph
    rkt_state.build_state()
    rkt_state.equilibrate_state()
    rkt_input = rktInputSpec.reaktoroInputSpec(rkt_state)
    # lime = Var(initialize=0.01, units=pyunits.mol / pyunits.s)
    # lime.construct()

    # rkt_input.register_chemical_addition("CaO", lime)
    rkt_input.assert_charge_neutrality(False)
    rkt_input.configure_specs(dissolve_species_in_rkt=True)
    input_names = rkt_input.rktEquilibriumSpecs.namesControlVariables()
    input_constraints = rkt_input.rktEquilibriumSpecs.namesConstraints()
    print("inputtn", input_names)
    print("input_c", input_constraints)
    expected_names = [
        "[H]",
        "[C]",
        "[O]",
        "[Na]",
        "[Mg]",
        "[Cl]",
        "[Ca]",
        "[H+]",
        "[H2O]",
        "[CO3-2]",
        "[CO2]",
        "[Ca+2]",
        "[Cl-]",
        "[HCO3-]",
        "[Mg+2]",
        "[MgCO3]",
        "[MgOH+]",
        "[Na+]",
        "[OH-]",
    ]
    expected_constraints = [
        "H_constraint",
        "C_constraint",
        "O_constraint",
        "Na_constraint",
        "Mg_constraint",
        "Cl_constraint",
        "Ca_constraint",
        "H+_constraint",
        "H2O_constraint",
        "CO3-2_constraint",
        "CO2_constraint",
        "Ca+2_constraint",
        "Cl-_constraint",
        "HCO3-_constraint",
        "Mg+2_constraint",
        "MgCO3_constraint",
        "MgOH+_constraint",
        "Na+_constraint",
        "OH-_constraint",
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
        "Na",
        "Mg",
        "Ca",
        "Cl",
    ]
    for ei in expected_inputs:
        assert ei in rkt_input.rktInputs.keys()
        # assert ei in rkt_input.userInputs.keys()
    assert len(rkt_input.rktInputs.keys()) == len(expected_inputs)
    expected_con_dict = {
        "H": [(2.0, "H2O")],
        "C": [(1.0, "CO3-2")],
        "O": [(1.0, "H2O"), (3.0, "CO3-2")],
        "Na": [(1, "Na")],
        "Mg": [(1, "Mg")],
        "Ca": [(1, "Ca")],
        "Cl": [(1, "Cl")],
    }
    expected_active_species = ["H2O", "CO3-2", "Na", "Mg", "Ca", "Cl"]
    for ion, ecd in expected_con_dict.items():
        rkt_ecd = rkt_input.constraintDict[ion]
        # order might change....
        expected_ecds = len(ecd)
        counted_ecds = 0
        for i, (mol, ion) in enumerate(ecd):
            if mol == rkt_ecd[i][0] and ion == rkt_ecd[i][1]:
                counted_ecds += 1
        assert counted_ecds == expected_ecds
    for eas in expected_active_species:
        assert eas in rkt_input.rktActiveSpecies
    assert len(rkt_input.rktActiveSpecies) == len(expected_active_species)
    assert len(rkt_input.constraintDict) == len(expected_con_dict)

    """ check charge neutrality enabled"""
    assert rkt_input.assertChargeNeutrality == False
    assert rkt_input.neutralityIon == "Cl"


def test_with_pyomo_sum(build_rkt_state_with_species):
    """testing if we can construct all constraints correctly but summing in pyomo
    so rkt inputs are elment species only"""
    m, rkt_state = build_rkt_state_with_species
    rkt_state.build_state()
    rkt_state.equilibrate_state()
    rkt_input = rktInputSpec.reaktoroInputSpec(rkt_state)
    rkt_input.configure_specs(dissolve_species_in_rkt=False)
    input_names = rkt_input.rktEquilibriumSpecs.namesControlVariables()
    input_constraints = rkt_input.rktEquilibriumSpecs.namesConstraints()
    print(input_names)
    print(input_constraints)
    expected_names = [
        "[Cl-]",
        "[C]",
        "[O]",
        "[Na]",
        "[Mg]",
        "[Ca]",
        "[H2O]",
        "[CO3-2]",
        "[CO2]",
        "[Ca+2]",
        "[HCO3-]",
        "[Mg+2]",
        "[MgCO3]",
        "[MgOH+]",
        "[Na+]",
        "[OH-]",
        "[H+]",
    ]
    expected_constraints = [
        "charge",
        "C_constraint",
        "O_constraint",
        "Na_constraint",
        "Mg_constraint",
        "Ca_constraint",
        "H2O_constraint",
        "CO3-2_constraint",
        "CO2_constraint",
        "Ca+2_constraint",
        "HCO3-_constraint",
        "Mg+2_constraint",
        "MgCO3_constraint",
        "MgOH+_constraint",
        "Na+_constraint",
        "OH-_constraint",
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
        assert ei in rkt_input.userInputs.keys()
    assert len(rkt_input.userInputs.keys()) == len(expected_inputs)
    expected_inputs = ["temperature", "pressure", "pH", "C", "O", "Na", "Mg", "Ca"]
    for ei in expected_inputs:
        assert ei in rkt_input.rktInputs.keys()
    assert len(rkt_input.rktInputs.keys()) == len(expected_inputs)
    expected_con_dict = {
        "C": [(1.0, "CO3-2"), (1.0, "CO2")],
        "O": [(1.0, "H2O"), (3.0, "CO3-2"), (2.0, "CO2")],
        "Na": [(1, "Na")],
        "Ca": [(1, "Ca")],
        "Mg": [(1, "Mg")],
    }
    expected_active_species = ["CO3-2", "CO2", "H2O", "Na", "Mg", "Ca"]
    for ion, ecd in expected_con_dict.items():
        rkt_ecd = rkt_input.constraintDict[ion]
        # order might change....
        expected_ecds = len(ecd)
        counted_ecds = 0
        for i, (mol, ion) in enumerate(ecd):
            if mol == rkt_ecd[i][0] and ion == rkt_ecd[i][1]:
                counted_ecds += 1
        assert counted_ecds == expected_ecds
    for eas in expected_active_species:
        assert eas in rkt_input.rktActiveSpecies
    assert len(rkt_input.rktActiveSpecies) == len(expected_active_species)
    assert len(rkt_input.constraintDict) == len(expected_con_dict)


def test_with_chemical(build_rkt_state_with_species):
    """testing if we can add chemicals"""
    m, rkt_state = build_rkt_state_with_species
    rkt_state.build_state()
    rkt_state.equilibrate_state()
    rkt_input = rktInputSpec.reaktoroInputSpec(rkt_state)

    lime = Var(initialize=0.01, units=pyunits.mol / pyunits.s)
    lime.construct()

    rkt_input.register_chemical_addition("CaO", lime)
    rkt_input.configure_specs(dissolve_species_in_rkt=True)
    expected_con_dict = {
        "C": [(1.0, "CO3-2"), (1.0, "CO2")],
        "O": [(1.0, "H2O"), (3.0, "CO3-2"), (2.0, "CO2"), (1, "CaO")],
        "Na": [(1, "Na")],
        "Ca": [(1, "Ca"), (1, "CaO")],
        "Mg": [(1, "Mg")],
    }
    expected_active_species = ["CaO", "CO3-2", "CO2", "H2O", "Na", "Mg", "Ca"]
    for ion, ecd in expected_con_dict.items():
        rkt_ecd = rkt_input.constraintDict[ion]
        # order might change....
        expected_ecds = len(ecd)
        counted_ecds = 0
        for i, (mol, ion) in enumerate(ecd):
            if mol == rkt_ecd[i][0] and ion == rkt_ecd[i][1]:
                counted_ecds += 1
        assert counted_ecds == expected_ecds
    for eas in expected_active_species:
        assert eas in rkt_input.rktActiveSpecies
    assert len(rkt_input.rktActiveSpecies) == len(expected_active_species)
    assert len(rkt_input.constraintDict) == len(expected_con_dict)
