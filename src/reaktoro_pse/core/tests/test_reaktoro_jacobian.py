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
import pytest
from reaktoro_pse.core.reaktoro_jacobian import (
    ReaktoroJacobianSpec,
    JacType,
)
from reaktoro_pse.core.reaktoro_outputs import (
    ReaktoroOutputSpec,
)
from reaktoro_pse.core.tests.test_reaktoro_state import (
    build_rkt_state_with_species,
)
import numpy as np

__author__ = "Alexander V. Dudchenko (NETL)"


@pytest.fixture
def build_standard_state(build_rkt_state_with_species):
    m, rkt_state = build_rkt_state_with_species
    rkt_state.register_mineral_phases("Calcite")
    rkt_state.set_mineral_phase_activity_model()
    rkt_state.build_state()
    rkt_state.equilibrate_state()
    rkt_outputs = ReaktoroOutputSpec(rkt_state)

    rkt_outputs.register_output("speciesAmount", get_all_indexes=True)
    rkt_outputs.register_output("pH")
    rkt_jacobian = ReaktoroJacobianSpec(rkt_state, rkt_outputs)
    return rkt_jacobian


def test_available_and_eval_jacobian(build_standard_state):
    """testing setting out outputs"""
    rkt_jac = build_standard_state
    rkt_jac.update_jacobian_absolute_values()
    av_dict = rkt_jac.jac_rows.display_available()
    print(av_dict)
    expected_dict = {
        ("temperature", None): True,
        ("pressure", None): True,
        ("speciesAmount", "H+"): True,
        ("speciesAmount", "H2O"): True,
        ("speciesAmount", "CO3-2"): True,
        ("speciesAmount", "CO2"): True,
        ("speciesAmount", "Ca+2"): True,
        ("speciesAmount", "Cl-"): True,
        ("speciesAmount", "HCO3-"): True,
        ("speciesAmount", "Mg+2"): True,
        ("speciesAmount", "MgCO3"): True,
        ("speciesAmount", "MgOH+"): True,
        ("speciesAmount", "Na+"): True,
        ("speciesAmount", "OH-"): True,
        ("speciesAmount", "Calcite"): True,
        ("temperature", "AqueousPhase"): True,
        ("temperature", "Calcite"): True,
        ("pressure", "AqueousPhase"): True,
        ("pressure", "Calcite"): True,
        ("amount", "AqueousPhase"): True,
        ("amount", "Calcite"): True,
        ("mass", "AqueousPhase"): True,
        ("mass", "Calcite"): True,
        ("speciesMoleFraction", "H+"): True,
        ("speciesMoleFraction", "H2O"): True,
        ("speciesMoleFraction", "CO3-2"): True,
        ("speciesMoleFraction", "CO2"): True,
        ("speciesMoleFraction", "Ca+2"): True,
        ("speciesMoleFraction", "Cl-"): True,
        ("speciesMoleFraction", "HCO3-"): True,
        ("speciesMoleFraction", "Mg+2"): True,
        ("speciesMoleFraction", "MgCO3"): True,
        ("speciesMoleFraction", "MgOH+"): True,
        ("speciesMoleFraction", "Na+"): True,
        ("speciesMoleFraction", "OH-"): True,
        ("speciesMoleFraction", "Calcite"): True,
        ("speciesStandardGibbsEnergy", "H+"): True,
        ("speciesStandardGibbsEnergy", "H2O"): True,
        ("speciesStandardGibbsEnergy", "CO3-2"): True,
        ("speciesStandardGibbsEnergy", "CO2"): True,
        ("speciesStandardGibbsEnergy", "Ca+2"): True,
        ("speciesStandardGibbsEnergy", "Cl-"): True,
        ("speciesStandardGibbsEnergy", "HCO3-"): True,
        ("speciesStandardGibbsEnergy", "Mg+2"): True,
        ("speciesStandardGibbsEnergy", "MgCO3"): True,
        ("speciesStandardGibbsEnergy", "MgOH+"): True,
        ("speciesStandardGibbsEnergy", "Na+"): True,
        ("speciesStandardGibbsEnergy", "OH-"): True,
        ("speciesStandardGibbsEnergy", "Calcite"): True,
        ("speciesStandardEnthalpy", "H+"): True,
        ("speciesStandardEnthalpy", "H2O"): True,
        ("speciesStandardEnthalpy", "CO3-2"): True,
        ("speciesStandardEnthalpy", "CO2"): True,
        ("speciesStandardEnthalpy", "Ca+2"): True,
        ("speciesStandardEnthalpy", "Cl-"): True,
        ("speciesStandardEnthalpy", "HCO3-"): True,
        ("speciesStandardEnthalpy", "Mg+2"): True,
        ("speciesStandardEnthalpy", "MgCO3"): True,
        ("speciesStandardEnthalpy", "MgOH+"): True,
        ("speciesStandardEnthalpy", "Na+"): True,
        ("speciesStandardEnthalpy", "OH-"): True,
        ("speciesStandardEnthalpy", "Calcite"): True,
        ("speciesStandardVolume", "H+"): True,
        ("speciesStandardVolume", "H2O"): True,
        ("speciesStandardVolume", "CO3-2"): True,
        ("speciesStandardVolume", "CO2"): True,
        ("speciesStandardVolume", "Ca+2"): True,
        ("speciesStandardVolume", "Cl-"): True,
        ("speciesStandardVolume", "HCO3-"): True,
        ("speciesStandardVolume", "Mg+2"): True,
        ("speciesStandardVolume", "MgCO3"): True,
        ("speciesStandardVolume", "MgOH+"): True,
        ("speciesStandardVolume", "Na+"): True,
        ("speciesStandardVolume", "OH-"): True,
        ("speciesStandardVolume", "Calcite"): True,
        ("speciesStandardVolumeT", "H+"): True,
        ("speciesStandardVolumeT", "H2O"): True,
        ("speciesStandardVolumeT", "CO3-2"): True,
        ("speciesStandardVolumeT", "CO2"): True,
        ("speciesStandardVolumeT", "Ca+2"): True,
        ("speciesStandardVolumeT", "Cl-"): True,
        ("speciesStandardVolumeT", "HCO3-"): True,
        ("speciesStandardVolumeT", "Mg+2"): True,
        ("speciesStandardVolumeT", "MgCO3"): True,
        ("speciesStandardVolumeT", "MgOH+"): True,
        ("speciesStandardVolumeT", "Na+"): True,
        ("speciesStandardVolumeT", "OH-"): True,
        ("speciesStandardVolumeT", "Calcite"): True,
        ("speciesStandardVolumeP", "H+"): True,
        ("speciesStandardVolumeP", "H2O"): True,
        ("speciesStandardVolumeP", "CO3-2"): True,
        ("speciesStandardVolumeP", "CO2"): True,
        ("speciesStandardVolumeP", "Ca+2"): True,
        ("speciesStandardVolumeP", "Cl-"): True,
        ("speciesStandardVolumeP", "HCO3-"): True,
        ("speciesStandardVolumeP", "Mg+2"): True,
        ("speciesStandardVolumeP", "MgCO3"): True,
        ("speciesStandardVolumeP", "MgOH+"): True,
        ("speciesStandardVolumeP", "Na+"): True,
        ("speciesStandardVolumeP", "OH-"): True,
        ("speciesStandardVolumeP", "Calcite"): True,
        ("speciesStandardHeatCapacityConstP", "H+"): True,
        ("speciesStandardHeatCapacityConstP", "H2O"): True,
        ("speciesStandardHeatCapacityConstP", "CO3-2"): True,
        ("speciesStandardHeatCapacityConstP", "CO2"): True,
        ("speciesStandardHeatCapacityConstP", "Ca+2"): True,
        ("speciesStandardHeatCapacityConstP", "Cl-"): True,
        ("speciesStandardHeatCapacityConstP", "HCO3-"): True,
        ("speciesStandardHeatCapacityConstP", "Mg+2"): True,
        ("speciesStandardHeatCapacityConstP", "MgCO3"): True,
        ("speciesStandardHeatCapacityConstP", "MgOH+"): True,
        ("speciesStandardHeatCapacityConstP", "Na+"): True,
        ("speciesStandardHeatCapacityConstP", "OH-"): True,
        ("speciesStandardHeatCapacityConstP", "Calcite"): True,
        ("correctiveMolarVolume", "AqueousPhase"): True,
        ("correctiveMolarVolume", "Calcite"): True,
        ("correctiveMolarVolumeT", "AqueousPhase"): True,
        ("correctiveMolarVolumeT", "Calcite"): True,
        ("correctiveMolarVolumeP", "AqueousPhase"): True,
        ("correctiveMolarVolumeP", "Calcite"): True,
        ("speciesCorrectiveMolarVolume", "H+"): True,
        ("speciesCorrectiveMolarVolume", "H2O"): True,
        ("speciesCorrectiveMolarVolume", "CO3-2"): True,
        ("speciesCorrectiveMolarVolume", "CO2"): True,
        ("speciesCorrectiveMolarVolume", "Ca+2"): True,
        ("speciesCorrectiveMolarVolume", "Cl-"): True,
        ("speciesCorrectiveMolarVolume", "HCO3-"): True,
        ("speciesCorrectiveMolarVolume", "Mg+2"): True,
        ("speciesCorrectiveMolarVolume", "MgCO3"): True,
        ("speciesCorrectiveMolarVolume", "MgOH+"): True,
        ("speciesCorrectiveMolarVolume", "Na+"): True,
        ("speciesCorrectiveMolarVolume", "OH-"): True,
        ("speciesCorrectiveMolarVolume", "Calcite"): True,
        ("correctiveMolarGibbsEnergy", "AqueousPhase"): True,
        ("correctiveMolarGibbsEnergy", "Calcite"): True,
        ("correctiveMolarEnthalpy", "AqueousPhase"): True,
        ("correctiveMolarEnthalpy", "Calcite"): True,
        ("correctiveMolarHeatCapacityConstP", "AqueousPhase"): True,
        ("correctiveMolarHeatCapacityConstP", "Calcite"): True,
        ("speciesActivityCoefficientLn", "H+"): True,
        ("speciesActivityCoefficientLn", "H2O"): True,
        ("speciesActivityCoefficientLn", "CO3-2"): True,
        ("speciesActivityCoefficientLn", "CO2"): True,
        ("speciesActivityCoefficientLn", "Ca+2"): True,
        ("speciesActivityCoefficientLn", "Cl-"): True,
        ("speciesActivityCoefficientLn", "HCO3-"): True,
        ("speciesActivityCoefficientLn", "Mg+2"): True,
        ("speciesActivityCoefficientLn", "MgCO3"): True,
        ("speciesActivityCoefficientLn", "MgOH+"): True,
        ("speciesActivityCoefficientLn", "Na+"): True,
        ("speciesActivityCoefficientLn", "OH-"): True,
        ("speciesActivityCoefficientLn", "Calcite"): True,
        ("speciesActivityLn", "H+"): True,
        ("speciesActivityLn", "H2O"): True,
        ("speciesActivityLn", "CO3-2"): True,
        ("speciesActivityLn", "CO2"): True,
        ("speciesActivityLn", "Ca+2"): True,
        ("speciesActivityLn", "Cl-"): True,
        ("speciesActivityLn", "HCO3-"): True,
        ("speciesActivityLn", "Mg+2"): True,
        ("speciesActivityLn", "MgCO3"): True,
        ("speciesActivityLn", "MgOH+"): True,
        ("speciesActivityLn", "Na+"): True,
        ("speciesActivityLn", "OH-"): True,
        ("speciesActivityLn", "Calcite"): True,
        ("speciesChemicalPotential", "H+"): True,
        ("speciesChemicalPotential", "H2O"): True,
        ("speciesChemicalPotential", "CO3-2"): True,
        ("speciesChemicalPotential", "CO2"): True,
        ("speciesChemicalPotential", "Ca+2"): True,
        ("speciesChemicalPotential", "Cl-"): True,
        ("speciesChemicalPotential", "HCO3-"): True,
        ("speciesChemicalPotential", "Mg+2"): True,
        ("speciesChemicalPotential", "MgCO3"): True,
        ("speciesChemicalPotential", "MgOH+"): True,
        ("speciesChemicalPotential", "Na+"): True,
        ("speciesChemicalPotential", "OH-"): True,
        ("speciesChemicalPotential", "Calcite"): True,
    }
    # jac_values = {}
    expected_values = {
        ("temperature", None): 293.15,
        ("pressure", None): 100000.0,
        ("speciesAmount", "H+"): 3.3073599548076276e-07,
        ("speciesAmount", "H2O"): 50.00283061005031,
        ("speciesAmount", "CO3-2"): 4.858515032822044e-07,
        ("speciesAmount", "CO2"): 0.0038308916484410043,
        ("speciesAmount", "Ca+2"): 0.007212090923707725,
        ("speciesAmount", "Cl-"): 0.5,
        ("speciesAmount", "HCO3-"): 0.004338167565255561,
        ("speciesAmount", "Mg+2"): 0.09995718896980808,
        ("speciesAmount", "MgCO3"): 4.254585850820247e-05,
        ("speciesAmount", "MgOH+"): 2.6517168390765324e-07,
        ("speciesAmount", "Na+"): 0.5,
        ("speciesAmount", "OH-"): 1.6426448993904945e-08,
        ("speciesAmount", "Calcite"): 0.002787909076292376,
        ("temperature", "AqueousPhase"): 293.15,
        ("temperature", "Calcite"): 293.15,
        ("pressure", "AqueousPhase"): 100000.0,
        ("pressure", "Calcite"): 100000.0,
        ("amount", "AqueousPhase"): 51.11821259320166,
        ("amount", "Calcite"): 0.002787909076292376,
        ("mass", "AqueousPhase"): 0.9332277420120733,
        ("mass", "Calcite"): 0.0002790448861460878,
        ("speciesMoleFraction", "H+"): 6.47002269255299e-09,
        ("speciesMoleFraction", "H2O"): 0.9781803406932955,
        ("speciesMoleFraction", "CO3-2"): 9.504469711188982e-09,
        ("speciesMoleFraction", "CO2"): 7.494181533550889e-05,
        ("speciesMoleFraction", "Ca+2"): 0.0001410865239186176,
        ("speciesMoleFraction", "Cl-"): 0.00978124966886061,
        ("speciesMoleFraction", "HCO3-"): 8.486540012223558e-05,
        ("speciesMoleFraction", "Mg+2"): 0.0019554124430223455,
        ("speciesMoleFraction", "MgCO3"): 8.323033288894915e-07,
        ("speciesMoleFraction", "MgOH+"): 5.187420890825887e-09,
        ("speciesMoleFraction", "Na+"): 0.00978124966886061,
        ("speciesMoleFraction", "OH-"): 3.2134239756437686e-10,
        ("speciesMoleFraction", "Calcite"): 1.0,
        ("speciesStandardGibbsEnergy", "H+"): 0.0,
        ("speciesStandardGibbsEnergy", "H2O"): 0.0,
        ("speciesStandardGibbsEnergy", "CO3-2"): 0.0,
        ("speciesStandardGibbsEnergy", "CO2"): -94047.6013718186,
        ("speciesStandardGibbsEnergy", "Ca+2"): 0.0,
        ("speciesStandardGibbsEnergy", "Cl-"): 0.0,
        ("speciesStandardGibbsEnergy", "HCO3-"): -58289.05139946059,
        ("speciesStandardGibbsEnergy", "Mg+2"): 0.0,
        ("speciesStandardGibbsEnergy", "MgCO3"): -16259.872017152507,
        ("speciesStandardGibbsEnergy", "MgOH+"): 67357.36289276853,
        ("speciesStandardGibbsEnergy", "Na+"): 0.0,
        ("speciesStandardGibbsEnergy", "OH-"): 79495.73931520512,
        ("speciesStandardGibbsEnergy", "Calcite"): -46945.31451388588,
        ("speciesStandardEnthalpy", "H+"): 0.0,
        ("speciesStandardEnthalpy", "H2O"): 0.0,
        ("speciesStandardEnthalpy", "CO3-2"): 0.0,
        ("speciesStandardEnthalpy", "CO2"): -27334.107086906843,
        ("speciesStandardEnthalpy", "Ca+2"): 0.0,
        ("speciesStandardEnthalpy", "Cl-"): 0.0,
        ("speciesStandardEnthalpy", "HCO3-"): -16358.697490147668,
        ("speciesStandardEnthalpy", "Mg+2"): 0.0,
        ("speciesStandardEnthalpy", "MgCO3"): 10079.550269621495,
        ("speciesStandardEnthalpy", "MgOH+"): 64513.091103668434,
        ("speciesStandardEnthalpy", "Na+"): 0.0,
        ("speciesStandardEnthalpy", "OH-"): 57312.73736577871,
        ("speciesStandardEnthalpy", "Calcite"): 12871.417396603083,
        ("speciesStandardVolume", "H+"): 0.0,
        ("speciesStandardVolume", "H2O"): 1.804844574465311e-05,
        ("speciesStandardVolume", "CO3-2"): -4.453260022294762e-06,
        ("speciesStandardVolume", "CO2"): 3.418595230492276e-05,
        ("speciesStandardVolume", "Ca+2"): -1.831312774455611e-05,
        ("speciesStandardVolume", "Cl-"): 1.7878459996233325e-05,
        ("speciesStandardVolume", "HCO3-"): 2.430528110419748e-05,
        ("speciesStandardVolume", "Mg+2"): -2.1743790327341462e-05,
        ("speciesStandardVolume", "MgCO3"): -1.7084061104585775e-05,
        ("speciesStandardVolume", "MgOH+"): 0.0,
        ("speciesStandardVolume", "Na+"): -1.806215053875956e-06,
        ("speciesStandardVolume", "OH-"): -4.358671824612883e-06,
        ("speciesStandardVolume", "Calcite"): 3.6899999999999996e-05,
        ("speciesStandardVolumeT", "H+"): 0.0,
        ("speciesStandardVolumeT", "H2O"): 0.0,
        ("speciesStandardVolumeT", "CO3-2"): 0.0,
        ("speciesStandardVolumeT", "CO2"): 0.0,
        ("speciesStandardVolumeT", "Ca+2"): 0.0,
        ("speciesStandardVolumeT", "Cl-"): 0.0,
        ("speciesStandardVolumeT", "HCO3-"): 0.0,
        ("speciesStandardVolumeT", "Mg+2"): 0.0,
        ("speciesStandardVolumeT", "MgCO3"): 0.0,
        ("speciesStandardVolumeT", "MgOH+"): 0.0,
        ("speciesStandardVolumeT", "Na+"): 0.0,
        ("speciesStandardVolumeT", "OH-"): 0.0,
        ("speciesStandardVolumeT", "Calcite"): 0.0,
        ("speciesStandardVolumeP", "H+"): 0.0,
        ("speciesStandardVolumeP", "H2O"): 0.0,
        ("speciesStandardVolumeP", "CO3-2"): 0.0,
        ("speciesStandardVolumeP", "CO2"): 0.0,
        ("speciesStandardVolumeP", "Ca+2"): 0.0,
        ("speciesStandardVolumeP", "Cl-"): 0.0,
        ("speciesStandardVolumeP", "HCO3-"): 0.0,
        ("speciesStandardVolumeP", "Mg+2"): 0.0,
        ("speciesStandardVolumeP", "MgCO3"): 0.0,
        ("speciesStandardVolumeP", "MgOH+"): 0.0,
        ("speciesStandardVolumeP", "Na+"): 0.0,
        ("speciesStandardVolumeP", "OH-"): 0.0,
        ("speciesStandardVolumeP", "Calcite"): 0.0,
        ("speciesStandardHeatCapacityConstP", "H+"): 0.0,
        ("speciesStandardHeatCapacityConstP", "H2O"): 0.0,
        ("speciesStandardHeatCapacityConstP", "CO3-2"): 0.0,
        ("speciesStandardHeatCapacityConstP", "CO2"): 672.5988861453108,
        ("speciesStandardHeatCapacityConstP", "Ca+2"): 0.0,
        ("speciesStandardHeatCapacityConstP", "Cl-"): 0.0,
        ("speciesStandardHeatCapacityConstP", "HCO3-"): 292.63763620645005,
        ("speciesStandardHeatCapacityConstP", "Mg+2"): 0.0,
        ("speciesStandardHeatCapacityConstP", "MgCO3"): 105.79595953289301,
        ("speciesStandardHeatCapacityConstP", "MgOH+"): 0.0,
        ("speciesStandardHeatCapacityConstP", "Na+"): 0.0,
        ("speciesStandardHeatCapacityConstP", "OH-"): -191.93637469177867,
        ("speciesStandardHeatCapacityConstP", "Calcite"): 366.4148868693987,
        ("correctiveMolarVolume", "AqueousPhase"): 0.0,
        ("correctiveMolarVolume", "Calcite"): 0.0,
        ("correctiveMolarVolumeT", "AqueousPhase"): 0.0,
        ("correctiveMolarVolumeT", "Calcite"): 0.0,
        ("correctiveMolarVolumeP", "AqueousPhase"): 0.0,
        ("correctiveMolarVolumeP", "Calcite"): 0.0,
        ("speciesCorrectiveMolarVolume", "H+"): 0.0,
        ("speciesCorrectiveMolarVolume", "H2O"): 0.0,
        ("speciesCorrectiveMolarVolume", "CO3-2"): 0.0,
        ("speciesCorrectiveMolarVolume", "CO2"): 0.0,
        ("speciesCorrectiveMolarVolume", "Ca+2"): 0.0,
        ("speciesCorrectiveMolarVolume", "Cl-"): 0.0,
        ("speciesCorrectiveMolarVolume", "HCO3-"): 0.0,
        ("speciesCorrectiveMolarVolume", "Mg+2"): 0.0,
        ("speciesCorrectiveMolarVolume", "MgCO3"): 0.0,
        ("speciesCorrectiveMolarVolume", "MgOH+"): 0.0,
        ("speciesCorrectiveMolarVolume", "Na+"): 0.0,
        ("speciesCorrectiveMolarVolume", "OH-"): 0.0,
        ("speciesCorrectiveMolarVolume", "Calcite"): 0.0,
        ("correctiveMolarGibbsEnergy", "AqueousPhase"): 0.0,
        ("correctiveMolarGibbsEnergy", "Calcite"): 0.0,
        ("correctiveMolarEnthalpy", "AqueousPhase"): 0.0,
        ("correctiveMolarEnthalpy", "Calcite"): 0.0,
        ("correctiveMolarHeatCapacityConstP", "AqueousPhase"): 0.0,
        ("correctiveMolarHeatCapacityConstP", "Calcite"): 0.0,
        ("speciesActivityCoefficientLn", "H+"): 0.0,
        ("speciesActivityCoefficientLn", "H2O"): 0.0,
        ("speciesActivityCoefficientLn", "CO3-2"): 0.0,
        ("speciesActivityCoefficientLn", "CO2"): 0.0,
        ("speciesActivityCoefficientLn", "Ca+2"): 0.0,
        ("speciesActivityCoefficientLn", "Cl-"): 0.0,
        ("speciesActivityCoefficientLn", "HCO3-"): 0.0,
        ("speciesActivityCoefficientLn", "Mg+2"): 0.0,
        ("speciesActivityCoefficientLn", "MgCO3"): 0.0,
        ("speciesActivityCoefficientLn", "MgOH+"): 0.0,
        ("speciesActivityCoefficientLn", "Na+"): 0.0,
        ("speciesActivityCoefficientLn", "OH-"): 0.0,
        ("speciesActivityCoefficientLn", "Calcite"): 0.0,
        ("speciesActivityLn", "H+"): -14.817529965555599,
        ("speciesActivityLn", "H2O"): -0.022306376849937075,
        ("speciesActivityLn", "CO3-2"): -14.432947397539397,
        ("speciesActivityLn", "CO2"): -5.460242285505204,
        ("speciesActivityLn", "Ca+2"): -4.8275809554418165,
        ("speciesActivityLn", "Cl-"): -0.5887317695618303,
        ("speciesActivityLn", "HCO3-"): -5.335887829018863,
        ("speciesActivityLn", "Mg+2"): -2.198597883963228,
        ("speciesActivityLn", "MgCO3"): -9.960512629041046,
        ("speciesActivityLn", "MgOH+"): -15.038472945932764,
        ("speciesActivityLn", "Na+"): -0.5887317695618303,
        ("speciesActivityLn", "OH-"): -17.819957646661674,
        ("speciesActivityLn", "Calcite"): 0.0,
        ("speciesChemicalPotential", "H+"): -36116.018823765444,
        ("speciesChemicalPotential", "H2O"): -54.36921862652199,
        ("speciesChemicalPotential", "CO3-2"): -35178.64320866276,
        ("speciesChemicalPotential", "CO2"): -107356.31163756711,
        ("speciesChemicalPotential", "Ca+2"): -11766.671305223217,
        ("speciesChemicalPotential", "Cl-"): -1434.96572782848,
        ("speciesChemicalPotential", "HCO3-"): -71294.66203242821,
        ("speciesChemicalPotential", "Mg+2"): -5358.828546167172,
        ("speciesChemicalPotential", "MgCO3"): -40537.47175482992,
        ("speciesChemicalPotential", "MgOH+"): 30702.82105897175,
        ("speciesChemicalPotential", "Na+"): -1434.96572782848,
        ("speciesChemicalPotential", "OH-"): 36061.64960513892,
        ("speciesChemicalPotential", "Calcite"): -46945.31451388588,
    }

    d = {}
    for key, available in av_dict.items():
        d[key] = rkt_jac.jac_rows.get_value(key)
    print(d)
    for key, available in av_dict.items():
        assert expected_dict[key] == available
        if available == False:
            assert rkt_jac.jac_rows.get_value(key) == 0
        else:
            assert (
                pytest.approx(rkt_jac.jac_rows.get_value(key), 1e-3)
                == expected_values[key]
            )
    types_jac = rkt_jac.display_jacobian_output_types()
    print(types_jac)
    expected_types = {
        ("speciesAmount", "H+"): "exact",
        ("speciesAmount", "H2O"): "exact",
        ("speciesAmount", "CO3-2"): "exact",
        ("speciesAmount", "CO2"): "exact",
        ("speciesAmount", "Ca+2"): "exact",
        ("speciesAmount", "Cl-"): "exact",
        ("speciesAmount", "HCO3-"): "exact",
        ("speciesAmount", "Mg+2"): "exact",
        ("speciesAmount", "MgCO3"): "exact",
        ("speciesAmount", "MgOH+"): "exact",
        ("speciesAmount", "Na+"): "exact",
        ("speciesAmount", "OH-"): "exact",
        ("speciesAmount", "Calcite"): "exact",
        ("pH", None): "calculated",
    }
    for key, t in types_jac.items():
        assert t == expected_types[key]
    assert len(types_jac) == len(expected_types)
    #     jac_values[key] = rkt_jac.jac_rows.get_value(prop, key)
    # print(jac_values)


def test_jacboian_output_types(build_standard_state):

    rkt_jac = build_standard_state

    # print(rkt_jac.output_specs.rkt_outputs.keys())
    assert (
        rkt_jac.output_specs.rkt_outputs[("speciesAmount", "H+")].jacobian_type
        == JacType.exact
    )
    assert (
        rkt_jac.output_specs.rkt_outputs[("pH", None)].jacobian_type
        == JacType.calculated
    )


def test_numeric_setup(build_standard_state):
    rkt_jac = build_standard_state
    for order in [2, 4, 6]:
        rkt_jac.configure_numerical_jacobian(jacobian_type="average", order=order)
        assert len(rkt_jac.numerical_steps) == order + 1
        assert len(rkt_jac.chem_prop_states) == order + 1
        assert len(rkt_jac.aqueous_prop_states) == order + 1
    for order in [2, 4, 6, 8, 10]:
        rkt_jac.configure_numerical_jacobian(
            jacobian_type="center_difference", order=order
        )
        assert len(rkt_jac.cdf_multipliers) == order
        assert pytest.approx(sum(rkt_jac.cdf_multipliers), 1e-3) == 0
        assert len(rkt_jac.chem_prop_states) == order
        assert len(rkt_jac.aqueous_prop_states) == order


def test_jacobian_matrix(build_standard_state):
    rkt_jac = build_standard_state

    dummyMatrix = np.ones(178)
    rkt_jac.update_jacobian_absolute_values()
    rkt_jac.partial_jac_vals = dummyMatrix
    jac_matrix = rkt_jac.process_jacobian_matrix(0.1, 0.01)
    assert pytest.approx(jac_matrix[0][0], 1e-5) == 2.93148000e02
    assert pytest.approx(jac_matrix[0][-1], 1e-5) == 2.93152000e02

    assert len(jac_matrix) == 178
