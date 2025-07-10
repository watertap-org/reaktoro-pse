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

from reaktoro_pse.reaktoro_block import ReaktoroBlock

from pyomo.environ import (
    ConcreteModel,
    Block,
    Var,
    Constraint,
    assert_optimal_termination,
    units as pyunits,
)

from pyomo.contrib.pynumero.interfaces.external_grey_box import (
    ExternalGreyBoxModel,
)
from reaktoro_pse.core.util_classes.cyipopt_solver import (
    get_cyipopt_watertap_solver,
)

from idaes.core.util.model_statistics import degrees_of_freedom


def build_comp(blk):
    blk.temp = Var(initialize=293.15, units=pyunits.K)
    blk.temp.fix()
    blk.pressure = Var(initialize=1e5, units=pyunits.Pa)
    blk.pressure.fix()
    blk.pH = Var(initialize=7, units=pyunits.dimensionless)
    blk.pH.fix()
    blk.composition = Var(
        ["H2O", "Mg", "Na", "Cl", "Ca", "HCO3", "SO4"],
        initialize=1,
        units=pyunits.mol / pyunits.s,
    )
    blk.composition["H2O"].fix(50)
    blk.composition["Mg"].fix(0.1)
    blk.composition["Na"].fix(0.5)
    blk.composition["Cl"].fix(0.5)
    blk.composition["Ca"].fix(0.01)
    blk.composition["HCO3"].fix(0.01)
    blk.composition["SO4"].fix(0.01)
    blk.outputs = Var(
        [("scalingTendency", "Calcite"), ("pH", None)],
        initialize=1,
    )


@pytest.fixture
def build_rkt_state_with_species():
    m = ConcreteModel()
    build_comp(m)
    return m


@pytest.fixture
def build_rkt_state_with_species_and_mixing():
    m = ConcreteModel()
    m.comp_a = Block()
    build_comp(m.comp_a)
    m.comp_b = Block()
    build_comp(m.comp_b)
    return m


@pytest.fixture
def build_rkt_state_with_indexed_species():
    m = ConcreteModel()
    m.temp = Var([0, 1], initialize=293.15, units=pyunits.K)
    m.temp.fix()
    m.pressure = Var([0, 1], initialize=1e5, units=pyunits.Pa)
    m.pressure.fix()
    m.pH = Var([0, 1], initialize=7, units=pyunits.dimensionless)
    m.pH.fix()
    m.composition = Var(
        [0, 1],
        ["H2O", "Mg", "Na", "Cl", "Ca", "HCO3"],
        initialize=1,
        units=pyunits.mol / pyunits.s,
    )
    for idx in [0, 1]:
        if idx == 0:
            m.composition[(idx, "H2O")].fix(50)
        else:
            m.composition[(idx, "H2O")].fix(20)
        m.composition[(idx, "Mg")].fix(0.1 * (1 + idx))
        m.composition[(idx, "Na")].fix(0.5 * (1 + idx))
        m.composition[(idx, "Cl")].fix(0.5 * (1 + idx))
        m.composition[(idx, "Ca")].fix(0.01 * (1 + idx))
        m.composition[(idx, "HCO3")].fix(0.01 * (1 + idx))

    m.outputs = Var(
        [0, 1], [("scalingTendency", "Calcite"), ("pH", None)], initialize=1
    )
    return m


def test_blockBuild(build_rkt_state_with_species):
    m = build_rkt_state_with_species
    m.outputs.display()
    m.property_block = ReaktoroBlock(
        aqueous_phase={
            "composition": m.composition,
            "convert_to_rkt_species": True,
        },
        system_state={
            "temperature": m.temp,
            "pressure": m.pressure,
            "pH": m.pH,
        },
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        outputs=m.outputs,
    )
    print("rkt block")
    m.property_block.reaktoro_model.display()
    print("rkt block")
    m.property_block.initialize()

    m.property_block.display_jacobian_scaling()
    cy_solver = get_cyipopt_watertap_solver()
    cy_solver.options["max_iter"] = 20
    m.pH.fix()
    m.composition["H2O"].unfix()
    m.composition["H2O"].setlb(30)
    m.outputs[("scalingTendency", "Calcite")].fix(5)
    m.property_block.output_constraints.pprint()
    print(degrees_of_freedom(m))
    assert degrees_of_freedom(m) == 0
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.display()
    assert pytest.approx(m.composition["H2O"].value, 1e-3) == 68.0601837


def test_activate_deactivate(build_rkt_state_with_species):
    m = build_rkt_state_with_species
    m.property_block = ReaktoroBlock(
        aqueous_phase={
            "composition": m.composition,
            "convert_to_rkt_species": True,
        },
        system_state={
            "temperature": m.temp,
            "pressure": m.pressure,
            "pH": m.pH,
        },
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        outputs=m.outputs,
    )
    m.property_block.initialize()
    m.property_block.deactivate()

    assert m.property_block.active == False
    for v in m.property_block.component_data_objects(Constraint):
        print(v.name)
        assert v.active == False
    for v in m.property_block.component_data_objects(ExternalGreyBoxModel):
        assert v.active == False

    cy_solver = get_cyipopt_watertap_solver()
    cy_solver.options["max_iter"] = 20
    m.pH.fix()
    m.composition["H2O"].unfix()
    m.composition["H2O"].setlb(30)
    m.outputs[("scalingTendency", "Calcite")].fix(5)

    # this should fail run solve, raising value error
    with pytest.raises(ValueError):
        result = cy_solver.solve(m, tee=True)

    assert pytest.approx(m.composition["H2O"].value, 1e-3) == 50
    m.property_block.activate()
    assert m.property_block.active == True
    for v in m.property_block.component_data_objects(Constraint):
        print(v.name)
        assert v.active == True
    for v in m.property_block.component_data_objects(ExternalGreyBoxModel):
        assert v.active == True

    # this solve should solve
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    assert pytest.approx(m.composition["H2O"].value, 1e-3) == 68.0601837


def test_blockBuild_solids_gas(build_rkt_state_with_species):
    m = build_rkt_state_with_species
    m.outputs.display()
    m.solid_gas_outputs = Var(
        [
            ("speciesAmount", "Calcite"),
            ("vaporPressure", "H2O(g)"),
            # ("speciesActivityLn", "H2O(g)"),
        ],
        initialize=0.5,
    )
    m.property_block = ReaktoroBlock(
        aqueous_phase={
            "composition": m.composition,
            "convert_to_rkt_species": True,
            "activity_model": "ActivityModelPitzer",
        },
        system_state={
            "temperature": m.temp,
            "pressure": m.pressure,
            "pH": m.pH,
        },
        mineral_phase={"phase_components": "Calcite"},
        gas_phase={
            "phase_components": ["H2O(g)"],
            "activity_model": "ActivityModelRedlichKwong",
        },
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        outputs=m.solid_gas_outputs,
    )
    m.display()
    m.property_block.initialize()
    cy_solver = get_cyipopt_watertap_solver()
    cy_solver.options["max_iter"] = 20
    m.temp.fix(273.15 + 50)
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.display()
    assert (
        pytest.approx(m.solid_gas_outputs[("vaporPressure", "H2O(g)")].value, 1e-1)
        == 49382.90
    )


@pytest.mark.parametrize(
    "coupling_type, ph_relax, water_relax",
    [
        (True, False, False),
        (False, True, False),
        (False, True, True),
    ],
)
def test_blockBuild_with_speciation_block(
    build_rkt_state_with_species, coupling_type, ph_relax, water_relax
):
    m = build_rkt_state_with_species
    m.CaO = Var(["CaO"], initialize=0.001, units=pyunits.mol / pyunits.s)
    m.CaO.fix()

    if ph_relax:
        # for pitzer we can't do pH relxacation and cahrge at as me time.
        charge_balance_pH = False
    else:
        charge_balance_pH = True
    m.property_block = ReaktoroBlock(
        aqueous_phase={
            "composition": m.composition,
            "convert_to_rkt_species": True,
        },
        system_state={
            "temperature": m.temp,
            "pressure": m.pressure,
            "pH": m.pH,
        },
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        chemistry_modifier=m.CaO,
        outputs=m.outputs,
        enable_pH_relaxation_on_property_block=ph_relax,
        enable_solvent_relaxation_on_property_block=water_relax,
        assert_charge_neutrality_on_property_block=charge_balance_pH,
        build_speciation_block=True,
        direct_speciation_to_property_block_coupling=coupling_type,
    )
    m.property_block.initialize()
    cy_solver = get_cyipopt_watertap_solver(limited_memory=False)
    cy_solver.options["max_iter"] = 50
    m.pH.unfix()
    new_scaling = m.property_block.display_jacobian_scaling()
    m.property_block.display_jacobian_scaling()
    m.outputs[("scalingTendency", "Calcite")].fix(5)
    print(degrees_of_freedom(m))
    assert degrees_of_freedom(m) == 0
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)

    assert pytest.approx(m.outputs[("pH", None)].value, 1e-2) == 6.7496301
    assert pytest.approx(m.pH.value, 1e-2) == 6.401
    m.property_block.update_block_scaling()
    m.property_block.update_jacobian_scaling()
    scaling_result = m.property_block.display_jacobian_scaling()
    print(scaling_result)
    if coupling_type:
        expected_scaling = {
            "property_block": {
                ("scalingTendency", "Calcite"): 4.843199291298157e-08,
                ("pH", None): 2.5891024643706047e-09,
            }
        }
        assert "speciation_block" not in scaling_result
        assert "property_block" in scaling_result
    else:
        expected_scaling = {
            "speciation_block": {
                ("speciesAmount", "H+"): 3.572346198193882e-17,
                ("speciesAmount", "H2O"): 4.999999999999997e-09,
                ("speciesAmount", "CO3-2"): 9.139929493679263e-15,
                ("speciesAmount", "CO2"): 2.8811688764488992e-11,
                ("speciesAmount", "Ca+2"): 1.0000000000000012e-12,
                ("speciesAmount", "Cl-"): 6.948233272771661e-11,
                ("speciesAmount", "HCO3-"): 2.8175490202179848e-11,
                ("speciesAmount", "SO4-2"): 9.99967082218832e-13,
                ("speciesAmount", "HSO4-"): 3.2917781167809654e-17,
                ("speciesAmount", "Mg+2"): 9.99535976796691e-12,
                ("speciesAmount", "MgCO3"): 6.27569290450919e-13,
                ("speciesAmount", "MgOH+"): 3.7201369273196507e-16,
                ("speciesAmount", "Na+"): 5.0000000000000015e-11,
                ("speciesAmount", "OH-"): 1.392175267414428e-16,
            }
        }
        if water_relax:
            expected_scaling["property_block"] = {
                ("scalingTendency", "Calcite"): 4.698386834294868e-08,
                ("pH", None): 6.749544872785636e-10,
                ("elementAmount", "H"): 1.0000508467563742e-08,
                ("elementAmount", "O"): 5.006613067418016e-09,
            }

        elif ph_relax:
            expected_scaling["property_block"] = {
                ("scalingTendency", "Calcite"): 4.698402870213914e-08,
                ("pH", None): 6.749544872785618e-10,
                ("elementAmount", "H"): 1.0000508467563742e-08,
            }

        else:
            expected_scaling["property_block"] = {
                ("scalingTendency", "Calcite"): 4.698387021360117e-08,
                ("pH", None): 4.860784942869509e-09,
            }
        assert "speciation_block" in scaling_result
        assert "property_block" in scaling_result
        new_scaling = {}
        for key in scaling_result["speciation_block"]:
            new_scaling[key] = 1
            assert (
                pytest.approx(scaling_result["speciation_block"][key], 1e-3)
                == expected_scaling["speciation_block"][key]
            )
        m.property_block.update_jacobian_scaling(new_scaling)
        scaling_result = m.property_block.display_jacobian_scaling()

        assert "speciation_block" in scaling_result
        for key in scaling_result["speciation_block"]:
            assert scaling_result["speciation_block"][key] == 1
    new_scaling = {}
    for key in scaling_result["property_block"]:
        new_scaling[key] = 1
        assert (
            pytest.approx(scaling_result["property_block"][key], 1e-3)
            == expected_scaling["property_block"][key]
        )
    m.property_block.update_jacobian_scaling(new_scaling)
    scaling_result = m.property_block.display_jacobian_scaling()
    m.property_block.display_reaktoro_state()
    assert "property_block" in scaling_result
    for key in scaling_result["property_block"]:
        assert scaling_result["property_block"][key] == 1


@pytest.mark.parametrize(
    "ph_relax, water_relax",
    [
        (True, False),
        (True, True),
    ],
)
def test_blockBuild_with_speciation_log_species_block(
    build_rkt_state_with_species, ph_relax, water_relax
):
    m = build_rkt_state_with_species
    m.CaO = Var(["CaO"], initialize=0.001, units=pyunits.mol / pyunits.s)
    m.CaO.fix()

    if ph_relax:
        # for pitzer we can't do pH relxacation and cahrge at as me time.
        charge_balance_pH = False
    else:
        charge_balance_pH = True
    m.property_block = ReaktoroBlock(
        aqueous_phase={
            "composition": m.composition,
            "convert_to_rkt_species": True,
        },
        system_state={
            "temperature": m.temp,
            "pressure": m.pressure,
            "pH": m.pH,
        },
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        chemistry_modifier=m.CaO,
        outputs=m.outputs,
        enable_pH_relaxation_on_property_block=ph_relax,
        enable_solvent_relaxation_on_property_block=water_relax,
        assert_charge_neutrality_on_property_block=charge_balance_pH,
        build_speciation_block=True,
        speciation_output_type="log10_species",
        direct_speciation_to_property_block_coupling=False,
    )
    m.property_block.initialize()
    cy_solver = get_cyipopt_watertap_solver(limited_memory=False)
    cy_solver.options["max_iter"] = 50
    m.pH.unfix()
    m.property_block.speciation_block.reaktoro_model.outputs.display()
    new_scaling = m.property_block.display_jacobian_scaling()
    m.property_block.display_jacobian_scaling()
    m.outputs[("scalingTendency", "Calcite")].fix(5)
    assert degrees_of_freedom(m) == 0
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)

    assert pytest.approx(m.outputs[("pH", None)].value, 1e-2) == 6.7496301
    assert pytest.approx(m.pH.value, 1e-2) == 6.401


@pytest.mark.parametrize(
    "ph_relax, water_relax",
    [
        (False, False),
        (True, False),
        (True, True),
    ],
)
def test_blockBuild_with_speciation_elements_block(
    build_rkt_state_with_species, ph_relax, water_relax
):
    m = build_rkt_state_with_species
    m.CaO = Var(["CaO"], initialize=0.001, units=pyunits.mol / pyunits.s)
    m.CaO.fix()

    if ph_relax:
        # for pitzer we can't do pH relxacation and cahrge at as me time.
        charge_balance_pH = False
    else:
        charge_balance_pH = True
    m.property_block = ReaktoroBlock(
        aqueous_phase={
            "composition": m.composition,
            "convert_to_rkt_species": True,
        },
        system_state={
            "temperature": m.temp,
            "pressure": m.pressure,
            "pH": m.pH,
        },
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        chemistry_modifier=m.CaO,
        outputs=m.outputs,
        enable_pH_relaxation_on_property_block=ph_relax,
        enable_solvent_relaxation_on_property_block=water_relax,
        assert_charge_neutrality_on_property_block=charge_balance_pH,
        build_speciation_block=True,
        speciation_output_type="elements",
        direct_speciation_to_property_block_coupling=False,
    )
    m.property_block.initialize()
    cy_solver = get_cyipopt_watertap_solver(limited_memory=False)
    cy_solver.options["max_iter"] = 50
    m.pH.unfix()
    m.property_block.speciation_block.reaktoro_model.outputs.display()
    new_scaling = m.property_block.display_jacobian_scaling()
    m.property_block.display_jacobian_scaling()
    m.outputs[("scalingTendency", "Calcite")].fix(5)
    assert degrees_of_freedom(m) == 0
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)

    assert pytest.approx(m.outputs[("pH", None)].value, 1e-2) == 6.7496301
    assert pytest.approx(m.pH.value, 1e-2) == 6.401


def test_blockBuild_with_speciation_and_mixing(
    build_rkt_state_with_species_and_mixing,
):
    m = build_rkt_state_with_species_and_mixing
    m.comp_a.CaO = Var(["CaO"], initialize=0.001, units=pyunits.mol / pyunits.s)
    m.comp_a.CaO.fix()
    m.comp_a.outputs.display()
    m.mixing_prop = ReaktoroBlock(
        aqueous_phase={
            "composition": m.comp_b.composition,
            "convert_to_rkt_species": True,
        },
        system_state={
            "temperature": m.comp_b.temp,
            "pressure": m.comp_b.pressure,
            "pH": m.comp_b.pH,
        },
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        outputs={"speciesAmount": True},
        build_graybox_model=False,
    )
    m.property_block = ReaktoroBlock(
        aqueous_phase={
            "composition": m.comp_a.composition,
            "convert_to_rkt_species": True,
        },
        system_state={
            "temperature": m.comp_a.temp,
            "pressure": m.comp_a.pressure,
            "pH": m.comp_a.pH,
        },
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        chemistry_modifier=m.comp_a.CaO,
        outputs=m.comp_a.outputs,
        build_speciation_block=True,
        external_speciation_reaktoro_blocks=[m.mixing_prop],
    )
    m.property_block.display()
    # iscale.set_scaling_factor(m.outputs[("scalingTendency", "Calcite")], 10)
    # iscale.set_scaling_factor(m.outputs[("pH", None)], 10)
    m.property_block.initialize()
    # m.property_block.display()
    # m.display()
    cy_solver = get_cyipopt_watertap_solver()
    cy_solver.options["max_iter"] = 100
    m.comp_a.pH.fix()
    m.comp_b.pH.unfix()
    m.comp_a.outputs[("scalingTendency", "Calcite")].fix(5)

    result = cy_solver.solve(m, tee=True)
    m.property_block.display_jacobian_scaling()
    m.display()
    assert_optimal_termination(result)

    # # m.display()
    # m.property_block.reaktoro_model.outputs.display()

    assert (
        pytest.approx(m.comp_a.outputs[("pH", None)].value, 1e-2) == 6.765353236871279
    )
    assert pytest.approx(m.comp_b.pH.value, 1e-2) == 6.222585017249812


def test_blockBuild_with_speciation_block_no_chem_addition(
    build_rkt_state_with_species,
):
    m = build_rkt_state_with_species
    m.outputs.display()
    m.property_block = ReaktoroBlock(
        aqueous_phase={
            "composition": m.composition,
            "convert_to_rkt_species": True,
        },
        system_state={
            "temperature": m.temp,
            "pressure": m.pressure,
            "pH": m.pH,
        },
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        outputs=m.outputs,
        build_speciation_block=True,
    )
    m.property_block.initialize()
    cy_solver = get_cyipopt_watertap_solver()
    cy_solver.options["max_iter"] = 20
    m.pH.unfix()
    m.outputs[("scalingTendency", "Calcite")].fix(5)
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.display()
    assert pytest.approx(m.outputs[("pH", None)].value, 1e-2) == m.pH.value


def test_blockBuild_with_temp_and_pressure_modification_in_speciation_block(
    build_rkt_state_with_species,
):
    m = build_rkt_state_with_species
    m.CaO = Var(["CaO"], initialize=0.001, units=pyunits.mol / pyunits.s)
    m.CaO.fix()
    m.outputs.display()
    m.temp_mod = Var(initialize=333.15, units=pyunits.K)
    m.temp_mod.fix()
    m.pressure_mod = Var(initialize=5e5, units=pyunits.Pa)
    m.pressure_mod.fix()

    m.outputs_mod = Var(
        [
            ("scalingTendency", "Calcite"),
            ("pH", None),
            ("temperature", None),
            ("pressure", None),
        ],
        initialize=1,
    )
    m.property_block = ReaktoroBlock(
        aqueous_phase={
            "composition": m.composition,
            "convert_to_rkt_species": True,
        },
        system_state={
            "temperature": m.temp,
            "pressure": m.pressure,
            "pH": m.pH,
        },
        system_state_modifier={
            "temperature": m.temp_mod,
            "pressure": m.pressure_mod,
        },
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        chemistry_modifier=m.CaO,
        outputs=m.outputs_mod,
        build_speciation_block=True,
    )
    m.property_block.initialize()
    cy_solver = get_cyipopt_watertap_solver()
    cy_solver.options["max_iter"] = 20
    m.pH.unfix()
    m.outputs_mod[("scalingTendency", "Calcite")].fix(5)
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.display()
    assert pytest.approx(m.outputs_mod[("pH", None)].value, 1e-2) == 6.250981308052
    assert pytest.approx(m.pH.value, 1e-2) == 5.995934005877454
    assert pytest.approx(m.outputs_mod[("temperature", None)].value, 1e-2) == 333.15
    assert pytest.approx(m.outputs_mod[("pressure", None)].value, 1e-2) == 5e5


@pytest.mark.parametrize(
    "direct_coupling, charge_relax",
    [
        (False, False),
        (False, True),
    ],
)
def test_blockBuild_with_speciation_block_no_chem_super_critical_db(
    build_rkt_state_with_species, direct_coupling, charge_relax
):
    translation_dict = {
        "H2O": "H2O(aq)",
        "Mg": "Mg+2",
        "Na": "Na+",
        "Cl": "Cl-",
        "SO4": "SO4-2",
        "Ca": "Ca+2",
        "HCO3": "HCO3-",
    }
    m = build_rkt_state_with_species
    m.outputs.display()
    m.CaO = Var(["CaO"], initialize=0.002, units=pyunits.mol / pyunits.s)
    m.CaO.fix()
    m.property_block = ReaktoroBlock(
        aqueous_phase={
            "composition": m.composition,
            "convert_to_rkt_species": True,
            "species_to_rkt_species_dict": translation_dict,
        },
        system_state={
            "temperature": m.temp,
            "pressure": m.pressure,
            "pH": m.pH,
        },
        chemistry_modifier=m.CaO,
        database="SupcrtDatabase",
        database_file="supcrtbl",
        outputs=m.outputs,
        build_speciation_block=True,
        enable_pH_relaxation_on_property_block=False,
        assert_charge_neutrality_on_property_block=True,
        enable_charge_relaxation_on_property_block=charge_relax,
        direct_speciation_to_property_block_coupling=direct_coupling,
    )

    m.property_block.initialize()
    m.property_block.display_jacobian_outputs()
    # sup critical does not like limited_memory option
    cy_solver = get_cyipopt_watertap_solver(limited_memory=False)

    cy_solver.options["max_iter"] = 50
    result = cy_solver.solve(m, tee=True)
    m.pH.unfix()
    # cy_solver = get_cyipopt_watertap_solver(limited_memory=True)
    m.outputs[("scalingTendency", "Calcite")].fix(5)
    # m.property_block.display()
    # m.property_block.speciation_block.outputs.display()
    # m.property_block.speciation_block.reaktoro_model.display()
    # m.property_block.reaktoro_model.display()
    result = cy_solver.solve(m, tee=True)

    m.display()
    assert_optimal_termination(result)
    m.property_block.display()
    # m.property_block.speciation_block.outputs.display()
    # m.property_block.speciation_block.reaktoro_model.display()
    # m.property_block.reaktoro_model.display()
    # m.property_block.display_reaktoro_state()
    m.property_block.display_jacobian_scaling()
    assert pytest.approx(m.outputs[("pH", None)].value, 1e-2) == 6.899783669305352
    assert pytest.approx(m.pH.value, 1e-2) == 6.330494020512032


def test_indexed_blockBuild(build_rkt_state_with_indexed_species):
    m = build_rkt_state_with_indexed_species
    m.outputs.display()
    m.property_block = ReaktoroBlock(
        [0, 1],
        aqueous_phase={
            "composition": m.composition,
            "convert_to_rkt_species": True,
        },
        system_state={
            "temperature": m.temp,
            "pressure": m.pressure,
            "pH": m.pH,
        },
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        outputs=m.outputs,
    )
    for blk in m.property_block:
        m.property_block[blk].initialize()
    m.property_block[0].reaktoro_model.display()
    cy_solver = get_cyipopt_watertap_solver()
    cy_solver.options["max_iter"] = 20
    m.pH.unfix()
    m.outputs[0, ("scalingTendency", "Calcite")].fix(5)
    m.outputs[1, ("scalingTendency", "Calcite")].fix(2.5)
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.display()
    assert pytest.approx(m.pH[0].value, 1e-3) == 6.78206
    assert pytest.approx(m.pH[1].value, 1e-3) == 5.719012533419923


def test_indexed_blockBuild_with_speciation_block(
    build_rkt_state_with_indexed_species,
):
    m = build_rkt_state_with_indexed_species
    m.CaO = Var([0, 1], ["CaO"], initialize=0.01, units=pyunits.mol / pyunits.s)
    m.CaO.fix()
    m.outputs.display()
    m.property_block = ReaktoroBlock(
        [0, 1],
        aqueous_phase={
            "composition": m.composition,
            "convert_to_rkt_species": True,
        },
        system_state={
            "temperature": m.temp,
            "pressure": m.pressure,
            "pH": m.pH,
        },
        chemistry_modifier=m.CaO,
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        outputs=m.outputs,
        build_speciation_block=True,
    )
    for blk in m.property_block:
        m.property_block[blk].initialize()
    m.property_block.display()
    cy_solver = get_cyipopt_watertap_solver()
    cy_solver.options["max_iter"] = 20
    m.CaO.unfix()
    m.outputs[(0, "pH", None)].fix(11.5)
    m.outputs[(1, "pH", None)].fix(10)
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.display()
    assert pytest.approx(m.CaO[(0, "CaO")].value, 1e-3) == 0.01732553618254949
    assert pytest.approx(m.CaO[(1, "CaO")].value, 1e-3) == 0.011351679127420139


def test_indexed_blockBuild_with_speciation_and_relaxation_block(
    build_rkt_state_with_indexed_species,
):
    m = build_rkt_state_with_indexed_species
    m.CaO = Var([0, 1], ["CaO"], initialize=0.01, units=pyunits.mol / pyunits.s)
    m.CaO.fix()
    m.outputs.display()
    m.property_block = ReaktoroBlock(
        [0, 1],
        aqueous_phase={
            "composition": m.composition,
            "convert_to_rkt_species": True,
        },
        system_state={
            "temperature": m.temp,
            "pressure": m.pressure,
            "pH": m.pH,
        },
        chemistry_modifier=m.CaO,
        database="PhreeqcDatabase",
        database_file="pitzer.dat",
        outputs=m.outputs,
        build_speciation_block=True,
        enable_pH_relaxation_on_property_block=True,
        enable_solvent_relaxation_on_property_block=True,
        assert_charge_neutrality_on_property_block=False,
        direct_speciation_to_property_block_coupling=False,
    )
    for blk in m.property_block:
        m.property_block[blk].initialize()
    m.property_block.display()
    cy_solver = get_cyipopt_watertap_solver()
    cy_solver.options["max_iter"] = 20
    m.CaO.unfix()
    m.outputs[(0, "pH", None)].fix(11.5)
    m.outputs[(1, "pH", None)].fix(10)
    result = cy_solver.solve(m, tee=True)
    assert_optimal_termination(result)
    m.display()
    assert pytest.approx(m.CaO[(0, "CaO")].value, 1e-3) == 0.01732553618254949
    assert pytest.approx(m.CaO[(1, "CaO")].value, 1e-3) == 0.011351679127420139
    assert (
        pytest.approx(m.property_block[0].relaxation_H2O.value, 1e-3)
        == 49.99237778905057
    )
    assert pytest.approx(m.property_block[0].relaxation_pH.value, 1e-3) == 11.5
    assert (
        pytest.approx(m.property_block[1].relaxation_H2O.value, 1e-3)
        == 20.006060396417787
    )
    assert pytest.approx(m.property_block[1].relaxation_pH.value, 1e-3) == 10
