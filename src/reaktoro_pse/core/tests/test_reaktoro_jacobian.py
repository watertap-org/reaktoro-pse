import pytest
from reaktoro_pse.core.reaktoro_jacobian import (
    reaktoroJacobianSpec,
    jacType,
)
from reaktoro_pse.core.reaktoro_outputs import (
    reaktoroOutputSpec,
)
from reaktoro_pse.core.tests.test_reaktoro_state import (
    build_rkt_state_with_species,
)
import numpy as np

__author__ = "Alexander V. Dudchenko (SLAC)"


@pytest.fixture
def build_standard_state(build_rkt_state_with_species):
    m, rkt_state = build_rkt_state_with_species
    rkt_state.register_mineral_phases("Calcite")
    rkt_state.build_state()
    rkt_state.equilibrate_state()
    rkt_outputs = reaktoroOutputSpec(rkt_state)

    rkt_outputs.register_output("speciesAmount", get_all_indexes=True)
    rkt_outputs.register_output("pH")
    rkt_jacobian = reaktoroJacobianSpec(rkt_state, rkt_outputs)
    return rkt_jacobian


def test_available_and_eval_jacobian(build_standard_state):
    """testing setting out outputs"""
    rkt_jac = build_standard_state
    rkt_jac.update_jacobian_absolute_values()
    av_dict = rkt_jac.jacRows.display_available()
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
        ("correctiveMolarVolume", "AqueousPhase"): False,
        ("correctiveMolarVolume", "Calcite"): False,
        ("correctiveMolarVolumeT", "AqueousPhase"): False,
        ("correctiveMolarVolumeT", "Calcite"): False,
        ("correctiveMolarVolumeP", "AqueousPhase"): False,
        ("correctiveMolarVolumeP", "Calcite"): False,
        ("speciesCorrectiveMolarVolumeI", "H+"): False,
        ("speciesCorrectiveMolarVolumeI", "H2O"): False,
        ("speciesCorrectiveMolarVolumeI", "CO3-2"): False,
        ("speciesCorrectiveMolarVolumeI", "CO2"): False,
        ("speciesCorrectiveMolarVolumeI", "Ca+2"): False,
        ("speciesCorrectiveMolarVolumeI", "Cl-"): False,
        ("speciesCorrectiveMolarVolumeI", "HCO3-"): False,
        ("speciesCorrectiveMolarVolumeI", "Mg+2"): False,
        ("speciesCorrectiveMolarVolumeI", "MgCO3"): False,
        ("speciesCorrectiveMolarVolumeI", "MgOH+"): False,
        ("speciesCorrectiveMolarVolumeI", "Na+"): False,
        ("speciesCorrectiveMolarVolumeI", "OH-"): False,
        ("speciesCorrectiveMolarVolumeI", "Calcite"): False,
        ("correctiveMolarGibbsEnergy", "AqueousPhase"): False,
        ("correctiveMolarGibbsEnergy", "Calcite"): False,
        ("correctiveMolarEnthalpy", "AqueousPhase"): False,
        ("correctiveMolarEnthalpy", "Calcite"): False,
        ("correctiveMolarHeatCapacityConstP", "AqueousPhase"): False,
        ("correctiveMolarHeatCapacityConstP", "Calcite"): False,
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
        ("speciesAmount", "H+"): 2.2071333802987095e-08,
        ("speciesAmount", "H2O"): 49.99910142779299,
        ("speciesAmount", "CO3-2"): 3.0086828822209897e-06,
        ("speciesAmount", "CO2"): 0.00010563929464380397,
        ("speciesAmount", "Ca+2"): 0.0011644572460291775,
        ("speciesAmount", "Cl-"): 0.5,
        ("speciesAmount", "HCO3-"): 0.0017929108410269271,
        ("speciesAmount", "Mg+2"): 0.09973313624253165,
        ("speciesAmount", "MgCO3"): 0.000262898427477642,
        ("speciesAmount", "MgOH+"): 3.965329990911938e-06,
        ("speciesAmount", "Na+"): 0.5,
        ("speciesAmount", "OH-"): 2.4617165481522745e-07,
        ("speciesAmount", "Calcite"): 0.008835542753970915,
        ("temperature", "AqueousPhase"): 293.15,
        ("temperature", "Calcite"): 293.15,
        ("pressure", "AqueousPhase"): 100000.0,
        ("pressure", "Calcite"): 100000.0,
        ("amount", "AqueousPhase"): 51.10216771210056,
        ("amount", "Calcite"): 0.008835542753970915,
        ("mass", "AqueousPhase"): 0.9326123531906763,
        ("mass", "Calcite"): 0.0008843591933419782,
        ("speciesMoleFraction", "H+"): 4.319060186904907e-10,
        ("speciesMoleFraction", "H2O"): 0.9784144913279996,
        ("speciesMoleFraction", "CO3-2"): 5.887583671932099e-08,
        ("speciesMoleFraction", "CO2"): 2.0672174855468895e-06,
        ("speciesMoleFraction", "Ca+2"): 2.278684639347391e-05,
        ("speciesMoleFraction", "Cl-"): 0.009784320751653834,
        ("speciesMoleFraction", "HCO3-"): 3.5084829495449784e-05,
        ("speciesMoleFraction", "Mg+2"): 0.001951641989130643,
        ("speciesMoleFraction", "MgCO3"): 5.144565079093306e-06,
        ("speciesMoleFraction", "MgOH+"): 7.759612103446997e-08,
        ("speciesMoleFraction", "Na+"): 0.009784320751653834,
        ("speciesMoleFraction", "OH-"): 4.8172448613551885e-09,
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
        ("correctiveMolarVolume", "AqueousPhase"): 0,
        ("correctiveMolarVolume", "Calcite"): 0,
        ("correctiveMolarVolumeT", "AqueousPhase"): 0,
        ("correctiveMolarVolumeT", "Calcite"): 0,
        ("correctiveMolarVolumeP", "AqueousPhase"): 0,
        ("correctiveMolarVolumeP", "Calcite"): 0,
        ("speciesCorrectiveMolarVolumeI", "H+"): 0,
        ("speciesCorrectiveMolarVolumeI", "H2O"): 0,
        ("speciesCorrectiveMolarVolumeI", "CO3-2"): 0,
        ("speciesCorrectiveMolarVolumeI", "CO2"): 0,
        ("speciesCorrectiveMolarVolumeI", "Ca+2"): 0,
        ("speciesCorrectiveMolarVolumeI", "Cl-"): 0,
        ("speciesCorrectiveMolarVolumeI", "HCO3-"): 0,
        ("speciesCorrectiveMolarVolumeI", "Mg+2"): 0,
        ("speciesCorrectiveMolarVolumeI", "MgCO3"): 0,
        ("speciesCorrectiveMolarVolumeI", "MgOH+"): 0,
        ("speciesCorrectiveMolarVolumeI", "Na+"): 0,
        ("speciesCorrectiveMolarVolumeI", "OH-"): 0,
        ("speciesCorrectiveMolarVolumeI", "Calcite"): 0,
        ("correctiveMolarGibbsEnergy", "AqueousPhase"): 0,
        ("correctiveMolarGibbsEnergy", "Calcite"): 0,
        ("correctiveMolarEnthalpy", "AqueousPhase"): 0,
        ("correctiveMolarEnthalpy", "Calcite"): 0,
        ("correctiveMolarHeatCapacityConstP", "AqueousPhase"): 0,
        ("correctiveMolarHeatCapacityConstP", "Calcite"): 0,
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
        ("speciesActivityLn", "H+"): -17.524496190186596,
        ("speciesActivityLn", "H2O"): -0.022061722167158845,
        ("speciesActivityLn", "CO3-2"): -12.609518162424983,
        ("speciesActivityLn", "CO2"): -9.050990154335567,
        ("speciesActivityLn", "Ca+2"): -6.651010190556204,
        ("speciesActivityLn", "Cl-"): -0.5886571873576057,
        ("speciesActivityLn", "HCO3-"): -6.2194248185354475,
        ("speciesActivityLn", "Mg+2"): -2.200767304527362,
        ("speciesActivityLn", "MgCO3"): -8.139252814490765,
        ("speciesActivityLn", "MgOH+"): -12.33343148718312,
        ("speciesActivityLn", "Na+"): -0.5886571873576057,
        ("speciesActivityLn", "OH-"): -15.112746767347897,
        ("speciesActivityLn", "Calcite"): 0.0,
        ("speciesChemicalPotential", "H+"): -42713.93651661518,
        ("speciesChemicalPotential", "H2O"): -53.772901079057554,
        ("speciesChemicalPotential", "CO3-2"): -30734.24493632715,
        ("speciesChemicalPotential", "CO2"): -116108.34506847845,
        ("speciesChemicalPotential", "Ca+2"): -16211.069577558761,
        ("speciesChemicalPotential", "Cl-"): -1434.7839423151083,
        ("speciesChemicalPotential", "HCO3-"): -73448.18145294234,
        ("speciesChemicalPotential", "Mg+2"): -5364.116258364351,
        ("speciesChemicalPotential", "MgCO3"): -36098.3611946915,
        ("speciesChemicalPotential", "MgOH+"): 37296.047357171774,
        ("speciesChemicalPotential", "Na+"): -1434.7839423151083,
        ("speciesChemicalPotential", "OH-"): 42660.16361553613,
        ("speciesChemicalPotential", "Calcite"): -46945.31451388588,
    }

    d = {}
    for key, available in av_dict.items():
        d[key] = rkt_jac.jacRows.get_value(key)
    print(d)
    for key, available in av_dict.items():
        assert expected_dict[key] == available
        if available == False:
            assert rkt_jac.jacRows.get_value(key) == 0
        else:
            assert (
                pytest.approx(rkt_jac.jacRows.get_value(key), 1e-3)
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
        ("pH", None): "numeric",
    }
    for key, t in types_jac.items():
        assert t == expected_types[key]
    assert len(types_jac) == len(expected_types)
    #     jac_values[key] = rkt_jac.jacRows.get_value(prop, key)
    # print(jac_values)


def test_jacboian_output_types(build_standard_state):

    rkt_jac = build_standard_state

    # print(rkt_jac.rktOutputSpec.rktOutputs.keys())
    assert (
        rkt_jac.rktOutputSpec.rktOutputs[("speciesAmount", "H+")].jacobianType
        == jacType.exact
    )
    assert (
        rkt_jac.rktOutputSpec.rktOutputs[("pH", None)].jacobianType == jacType.numeric
    )


def test_numeric_setup(build_standard_state):
    rkt_jac = build_standard_state
    for order in [2, 4, 6]:
        rkt_jac.configure_numerical_jacobian(jacobian_type="average", order=order)
        assert len(rkt_jac.numericalSteps) == order + 1
        assert len(rkt_jac.chemPropStates) == order + 1
        assert len(rkt_jac.aqueousPropStates) == order + 1
    for order in [2, 4, 6, 8, 10]:
        rkt_jac.configure_numerical_jacobian(
            jacobian_type="center_difference", order=order
        )
        assert len(rkt_jac.cdfMultipliers) == order
        assert pytest.approx(sum(rkt_jac.cdfMultipliers), 1e-3) == 0
        assert len(rkt_jac.chemPropStates) == order
        assert len(rkt_jac.aqueousPropStates) == order


def test_jacobian_matrix(build_standard_state):
    rkt_jac = build_standard_state

    dummyMatrix = np.ones((len(rkt_jac.jacRows.keys), 2))
    dummyMatrix[:, 1] = 10
    rkt_jac.update_jacobian_absolute_values()
    jac_dict, jac_matrix = rkt_jac.process_jacobian_matrix(dummyMatrix, 0, 100)

    for key, value in jac_dict.items():
        assert value == 1

    assert pytest.approx(jac_matrix[0][0], 1e-3) == 293.14999
    assert pytest.approx(jac_matrix[0][0], 1e-3) == 293.15000293

    assert len(jac_matrix) == 178
