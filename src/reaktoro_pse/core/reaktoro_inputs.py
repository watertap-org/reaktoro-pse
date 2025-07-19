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
import reaktoro as rkt
from reaktoro_pse.core.util_classes.rkt_inputs import RktInputs, RktInput, RktInputTypes
from reaktoro_pse.core.reaktoro_state import ReaktoroState
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)

__author__ = "Alexander V. Dudchenko"

""" class to setup input constraints, and specs for reaktoro solver class"""


class ReaktoroInputExport:
    def __init__(self):
        self.ignore_elements_for_constraints = []
        self.fixed_solvent_specie = {}
        self.fixed_solvent_speciation = {}
        self.fixed_solvent_type = {}

        self.fixed_species = {}
        self.rkt_chemical_inputs = None
        self.assert_charge_neutrality = None
        self.neutrality_ion = None
        self.dissolve_species_in_rkt = None
        self.exact_speciation = None

    def copy_chem_inputs(self, chem_inputs):
        self.rkt_chemical_inputs = RktInputs()
        for key, obj in chem_inputs.items():
            obj.update_values(True)
            self.rkt_chemical_inputs[key] = None
            self.rkt_chemical_inputs[key].time_unit = obj.time_unit
            self.rkt_chemical_inputs[key].main_unit = obj.main_unit
            self.rkt_chemical_inputs[key].conversion_unit = obj.conversion_unit
            self.rkt_chemical_inputs[key].conversion_value = obj.conversion_value
            self.rkt_chemical_inputs[key].required_unit = obj.required_unit
            self.rkt_chemical_inputs[key].lower_bound = obj.lower_bound
            self.rkt_chemical_inputs[key].input_type = obj.input_type
            self.rkt_chemical_inputs[key].io_type = obj.io_type
            self.rkt_chemical_inputs[key].value = obj.value
            self.rkt_chemical_inputs[key].converted_value = obj.converted_value

        self.rkt_chemical_inputs.registered_phases = chem_inputs.registered_phases
        self.rkt_chemical_inputs.all_species = chem_inputs.all_species
        self.rkt_chemical_inputs.species_list = chem_inputs.species_list
        self.rkt_chemical_inputs.convert_to_rkt_species = (
            chem_inputs.convert_to_rkt_species
        )
        self.rkt_chemical_inputs.composition_is_elements = (
            chem_inputs.composition_is_elements
        )
        self.rkt_chemical_inputs.conversion_method = chem_inputs.conversion_method
        self.rkt_chemical_inputs.rkt_input_list = chem_inputs.rkt_input_list


class ReaktoroInputSpec:
    def __init__(self, reaktor_state=None):
        # initialize parameters needed to build reaktor solver
        if reaktor_state is not None:
            self.state = reaktor_state
            if isinstance(self.state, ReaktoroState) == False:
                raise TypeError("Reator inputs require rektoroState class")
            self.user_inputs = reaktor_state.inputs  # user inputs provided to state
            self.rkt_inputs = RktInputs()  # inputs that will be fed to rkt spec
            self.rkt_chemical_inputs = RktInputs()
            # tracking which elements to not include in summation constraints
            self.ignore_elements_for_constraints = []
            self.fixed_solvent_specie = {}
            self.fixed_solvent_speciation = {}
            self.fixed_solvent_type = {}

            self.fixed_species = {}
            # execute default configuration options, user can update settings
            self.register_charge_neutrality()
            self.default_speciation()
            self.register_open_species()
            # register default for aqueous phase
            if RktInputTypes.aqueous_phase in self.state.inputs.registered_phases:
                self.register_fixed_solvent_specie(RktInputTypes.aqueous_phase, "H2O")

    def register_chemistry_modifiers(
        self, chemical_dict, index=None, log10_basis=False
    ):
        """registers chemistry modifiers being added to system
        chemistry_modifier -- chemicals to be added (pyo object should be mole flow of chemical that would enter a system with same species as in apparat_species_mol_flow
                        chemistry_modifier = {'HCl':m.fs.HCl_dose} example for HCl
        """
        for chemical, obj in chemical_dict.items():
            if index is None or index in chemical:
                if isinstance(chemical, tuple):
                    chemical = chemical[-1]
                self.register_chemistry_modifier(chemical, obj, log10_basis=log10_basis)

    def register_chemistry_modifier(self, chemical, pyomo_var, log10_basis=False):
        chemical = self.safe_modifier_name(chemical)
        if chemical not in self.chemical_to_elements:
            raise ValueError(
                f"{chemical} is not available in chemical_to_element dict, please add"
            )
        self.rkt_chemical_inputs[chemical] = RktInput(
            var_name=chemical, pyomo_var=pyomo_var
        )
        mw, mw_unit = self.get_modifier_mw(self.chemical_to_elements[chemical])
        self.state.verify_unit(self.rkt_chemical_inputs[chemical], mw, mw_unit)

    def register_open_species(self, specie=None):
        """registers species to open to optimization and write empty constraint for,
        this can help with solvability of some problems, but can
        lead to unexpected results depending on database, activity coefficients, and inputs chosen
        """
        self.empty_constraints = []
        if specie is not None:
            if isinstance(specie, str):
                self.empty_constraints = [specie]
            else:
                self.empty_constraints = specie
            for spc in specie:
                _log.warning(
                    f"Registered an empty constraint for {spc}, this can lead to unexpected results depending on reaktoro configuration, please use with caution"
                )

    def register_charge_neutrality(self, assert_neutrality=True, ion="Cl"):
        self.assert_charge_neutrality = assert_neutrality
        self.neutrality_ion = ion

    def register_fixed_solvent_specie(self, phase, specie, solvent_name=None):
        """defines aqueous species for system - used to set species when speciating - H/O change based on
        system speciation, so if we want to specify pH, we need to allow system to find eq. H/O and fix
        H2O"""
        self.fixed_solvent_specie[phase] = specie
        self.fixed_solvent_speciation[phase] = {}
        if solvent_name is not None:
            self.fixed_solvent_type[specie] = solvent_name

    def register_fixed_species(self, specie, input_name=None):
        if specie not in self.fixed_species:
            if input_name is not None:
                self.fixed_species[specie] = input_name
            else:
                self.fixed_species[specie] = specie

    def register_free_elements(self, elements):
        if elements is not None:
            if isinstance(elements, str):
                elements = [elements]
            for e in elements:
                if e not in self.ignore_elements_for_constraints:
                    self.ignore_elements_for_constraints.append(e)

    def configure_specs(
        self,
        dissolve_species_in_rkt=True,
        exact_speciation=False,
    ):
        """configures specification for the problem

        Keyword arguments:
        dissolve_species_in_rkt-- If true, species would be summed up to element amount in rkt, if false
        mode will contain conditions to build pyomo constraints via raktoro output class
        exact_speciation -- if True, will write exact element amount for all input species other wise
        will leave  H, and O open, while fixing aqueousSolvent to specified value (e.g. H2O)

        """
        self.dissolve_species_in_rkt = dissolve_species_in_rkt
        self.exact_speciation = exact_speciation

    def build_input_specs(self):
        """function to build all the input specs"""
        self.breakdown_species_to_elements()
        self.equilibrium_specs = rkt.EquilibriumSpecs(self.state.state.system())
        self.add_specs(
            self.equilibrium_specs,
            self.assert_charge_neutrality,
            self.dissolve_species_in_rkt,
        )
        # get input name order!
        self.register_jacobian_indexes()

    def register_jacobian_indexes(self):
        for idx, spec in enumerate(self.equilibrium_specs.namesInputs()):
            if spec == "T":
                spec_var_name = RktInputTypes.temperature
            elif spec == "P":
                spec_var_name = RktInputTypes.pressure
            elif spec == "H":
                spec_var_name = RktInputTypes.enthalpy
            elif "input" in spec:
                spec_var_name = spec.replace("input", "")
            else:
                spec_var_name = spec
            # only care for indexes that exists and were added to spec
            if self.rkt_inputs.get(spec_var_name) is not None:
                self.rkt_inputs[spec_var_name].set_jacobian_index(idx)
                # tracking inputs we are passing into our spec problem
                if spec_var_name not in self.rkt_inputs.rkt_input_list:
                    self.rkt_inputs.rkt_input_list.append(spec_var_name)

    def breakdown_species_to_elements(self):
        """this will take all species in rktstate and create a dictionary containing
        their elements and species amounts for use in summing
        eg. {'H2O':{'H':2,'O':1}}"""

        # TODO: probably want to make a class to track this
        self.specie_to_elements = {}
        self.element_to_species = {}
        for specie in self.state.state.system().species():
            self.specie_to_elements[specie.name()] = {}
            for i, el in enumerate(specie.elements().symbols()):
                if el not in self.element_to_species:
                    self.element_to_species[el] = []
                self.element_to_species[el].append(
                    (specie.elements().coefficients()[i], specie.name())
                )
                self.specie_to_elements[specie.name()][
                    el
                ] = specie.elements().coefficients()[i]
        self.state.specie_to_elements = self.specie_to_elements
        self.state.element_to_species = self.element_to_species
        self.chemical_to_elements.update(self.specie_to_elements)

    def add_specs(
        self, specs_object, assert_charge_neutrality, dissolve_species_in_rkt
    ):
        # ignore elements for constraints

        pressure_not_set = True
        temperature_not_set = True
        for input_name, _ in self.state.inputs.items():
            if input_name == RktInputTypes.temperature:
                specs_object.temperature()
                temperature_not_set = False
                self.rkt_inputs[RktInputTypes.temperature] = self.state.inputs[
                    RktInputTypes.temperature
                ]
                self.rkt_inputs[RktInputTypes.temperature].set_lower_bound(0)
            elif input_name == RktInputTypes.pressure:
                specs_object.pressure()
                pressure_not_set = False
                self.rkt_inputs[RktInputTypes.pressure] = self.state.inputs[
                    RktInputTypes.pressure
                ]
                self.rkt_inputs[RktInputTypes.pressure].set_lower_bound(0)
            elif input_name == RktInputTypes.enthalpy:
                specs_object.enthalpy()
                self.rkt_inputs[RktInputTypes.enthalpy] = self.state.inputs[
                    RktInputTypes.enthalpy
                ]
                self.rkt_inputs[RktInputTypes.enthalpy].set_lower_bound(None)
            elif input_name == RktInputTypes.pH:
                specs_object.pH()
                self.rkt_inputs[RktInputTypes.pH] = self.state.inputs[RktInputTypes.pH]
                self.rkt_inputs[RktInputTypes.pH].set_lower_bound(-1)
                self.rkt_inputs[RktInputTypes.pH].set_upper_bound(14)
            elif input_name == RktInputTypes.pOH:
                self.write_pOH_constraint(specs_object)
                self.rkt_inputs[RktInputTypes.pOH] = self.state.inputs[
                    RktInputTypes.pOH
                ]
                self.rkt_inputs[RktInputTypes.pOH].set_lower_bound(-1)
                self.rkt_inputs[RktInputTypes.pOH].set_upper_bound(14)
            else:
                pass
        if pressure_not_set:
            specs_object.unknownPressure()
            # self.write_empty_con(specs_object, "open_pressure")
        if temperature_not_set:
            specs_object.unknownTemperature()
            # self.write_empty_con(specs_object, "open_temperature")
        if assert_charge_neutrality:
            specs_object.charge()
            if self.neutrality_ion is not None:
                if self.neutrality_ion == RktInputTypes.pH:
                    specs_object.openTo("H+")
                elif self.neutrality_ion == RktInputTypes.pOH:
                    specs_object.openTo("OH-")
                else:
                    self.ignore_elements_for_constraints.append(self.neutrality_ion)

                    if self.neutrality_ion not in specs_object.namesInputs():
                        # needs to be a species!
                        specs_object.openTo(self.neutrality_ion)

        self._find_element_sums()
        # add/check if vars in rkt Inputs
        if dissolve_species_in_rkt:
            self.write_active_species(specs_object)
        else:
            for element in self.constraint_dict:
                if element not in self.rkt_inputs:
                    self.rkt_inputs[element] = RktInput(element)
                    self.rkt_inputs[element].set_rkt_input_name(f"input{element}")
                    self.rkt_inputs[element].set_lower_bound(0)
                    self.rkt_inputs[element].io_type = RktInputTypes.element

        # write reaktoro constraints to spec
        for element in self.constraint_dict:
            if dissolve_species_in_rkt:
                self.write_element_sum_constraint(specs_object, element)
            else:
                self.write_elementAmount_constraint(specs_object, element)
        for specie, input_name in self.fixed_species.items():
            self.write_speciesAmount_constraint(specs_object, specie, input_name)
            self.rkt_inputs[input_name] = self.state.inputs[input_name]
            self.rkt_inputs[input_name].set_rkt_input_name(input_name)
            self.rkt_inputs[input_name].set_lower_bound(0)
            self.rkt_inputs[input_name].io_type = RktInputTypes.specie
            self.rkt_inputs.rkt_input_list.append(input_name)
        if self.exact_speciation == False or self.fixed_solvent_type != {}:
            self.add_solvent_constraints(specs_object)

        self.write_empty_constraints(specs_object)

    def add_solvent_constraints(self, specs_object):
        """adds constraints for fixed solvent species
        args:
        specs_object -- reaktoro spec object
        auto_register -- if True, will add the fixed solvent species to the rkt_inputs
        """
        for phase in self.state.inputs.registered_phases:
            if phase in self.fixed_solvent_specie:
                specie = self.state.inputs.convert_rkt_species_fun(
                    self.fixed_solvent_specie[phase], phase
                )
                if self.fixed_solvent_type == {}:
                    spc_name = specie
                else:
                    spc_name = self.fixed_solvent_type[self.fixed_solvent_specie[phase]]
                self.write_speciesAmount_constraint(
                    specs_object, specie, input_name=spc_name
                )
                if "H2O_evaporation" in spc_name:
                    assert False
                if spc_name not in self.rkt_inputs and spc_name in self.state.inputs:
                    self.rkt_inputs[spc_name] = self.state.inputs[spc_name]
                    self.rkt_inputs[spc_name].set_rkt_input_name(spc_name)
                    self.rkt_inputs[spc_name].set_lower_bound(0)
                    self.rkt_inputs.rkt_input_list.append(spc_name)

                    self.rkt_inputs[spc_name].io_type = RktInputTypes.specie
        self.write_open_solvent_constraints(specs_object)

    def update_constraint_dict(self, element, specie, coeff):
        if (
            element
            not in self.ignore_elements_for_constraints
            # and specie not in self.fixed_species
        ):
            self.constraint_dict[element].append((coeff, specie))
            if specie not in self.active_species:
                self.active_species.append(specie)
        self.all_inclusive_constraint_dict[element].append((coeff, specie))

    def _find_element_sums(self):
        """
        Here in we will take all input species, elements, and chemicals and organize them such that
        amount of element == sum of all species inputs, except those we excluded due to
        specifying pH or for charge neutrality or otherwise, we will also track their
        respective coefficients"""
        # TODO: Should add some sort of way to check if all species were accounted for!
        self.constraint_dict = {}
        self.all_inclusive_constraint_dict = {}

        self.active_species = []
        rktState = self.state.state
        if self.exact_speciation == False:
            for phase in self.state.inputs.registered_phases:
                if phase in self.fixed_solvent_specie:
                    specie_elements = self.specie_to_elements[
                        self.state.inputs.convert_rkt_species_fun(
                            self.fixed_solvent_specie[phase], phase
                        )
                    ]

                    for element, coeff in specie_elements.items():
                        self.ignore_elements_for_constraints.append(element)
                        self.fixed_solvent_speciation[phase][element] = coeff
                        _log.info(
                            f"Exact speciation is not provided! Fixing aqueous solvent and, excluding {element}"
                        )
        self.rkt_elements = [specie.symbol() for specie in rktState.system().elements()]
        # loop over all elements in the rkt system
        for element in self.rkt_elements:
            # skip any we want to ignore
            # if element not in self.ignore_elements_for_constraints:
            self.constraint_dict[element] = []
            self.all_inclusive_constraint_dict[element] = []
            # if user provided true elements as inputs, we just add them as 1:1
            # constraints
            for phase in self.state.inputs.registered_phases:
                if self.state.inputs.composition_is_elements[phase]:
                    self.update_constraint_dict(element, element, 1)
                else:
                    # otherwise find break down of provided species to elements
                    for specie in self.state.inputs.species_list[phase]:
                        # check if specie is in list of rkt species
                        spc_dict = self.specie_to_elements.get(specie)
                        # might be empty as specie might not exist, thats okay
                        if spc_dict is not None:
                            coef = spc_dict.get(element)
                            # checks if element is in the actual species, might not exists and thats okay
                            # example C might not be in "H2O"
                            if coef is not None:
                                self.update_constraint_dict(element, specie, coef)
            # now lets also check if user provided chemical inputs and
            # add them to our elemental sum constraints (e.g. H = H(from H2O) + H (from HCL))
            for specie in self.rkt_chemical_inputs.keys():
                # if specie is a element add directly
                if specie == element:
                    self.update_constraint_dict(element, specie, 1)
                # if not then find species and add their coefficients
                elif specie in self.chemical_to_elements:
                    coef = self.chemical_to_elements[specie].get(element)
                    if coef is not None:
                        self.update_constraint_dict(element, specie, coef)
            # make sure we did not create empty element lists
            if len(self.constraint_dict[element]) == 0:
                del self.constraint_dict[element]
            if len(self.all_inclusive_constraint_dict[element]) == 0:
                del self.all_inclusive_constraint_dict[element]
        self.state.all_inclusive_constraint_dict = self.all_inclusive_constraint_dict

    def write_active_species(self, spec_object):
        # build inputs into rkt model, and track their indexes for writing rkt constraints
        for specie in self.active_species:
            self.register_input_species(spec_object, specie)

    def register_input_species(self, spec_object, specie):
        input_name = f"input{specie}"
        idx = spec_object.addInput(input_name)
        if specie in self.state.inputs:
            self.rkt_inputs[specie] = self.state.inputs[specie]
            self.rkt_inputs[specie].set_rkt_index(idx)
            self.rkt_inputs[specie].set_rkt_input_name(input_name)
            self.rkt_inputs[specie].set_lower_bound(0)

            self.rkt_inputs[specie].io_type = RktInputTypes.specie
        elif specie in self.rkt_chemical_inputs:
            self.rkt_inputs[specie] = self.rkt_chemical_inputs[specie]
            self.rkt_inputs[specie].set_rkt_index(idx)
            self.rkt_inputs[specie].set_rkt_input_name(input_name)

            self.rkt_inputs[specie].set_lower_bound(0)

            self.rkt_inputs[specie].io_type = RktInputTypes.chemical_specie
        # elif specie in self.rkt_inputs:
        #     self.rkt_inputs[specie].set_rkt_index(idx)
        #     self.rkt_inputs[specie].set_rkt_input_name(input_name)
        #     self.rkt_inputs[specie].set_lower_bound(0)
        else:
            raise KeyError(f"Specie is not found {specie}")

    def default_speciation(self):
        # TODO: probably want to make a class to track this stuff
        """defines species to element conversions"""
        self.chemical_to_elements = {
            "HCl": {"H": 1, "Cl": 1},
            "H2SO4": {"H": 2, "S": 1, "O": 4},
            "CaO": {"Ca": 1, "O": 1},
            "Ca(OH)2": {"Ca": 1, "O": 2, "H": 2},
            "Na2CO3": {"Na": 2, "C": 1, "O": 3},
            "CO2": {"C": 1, "O": 2},
            "NaOH": {"Na": 1, "O": 1, "H": 1},
            "H": {"H": 1},
            "OH": {"O": 1, "H": 1},
            "H2O_evaporation": {"O": -1, "H": -2},
        }
        self.ensure_safe_modifier_names()

    def ensure_safe_modifier_names(self):
        """we need to ensure any chemical has a safe subname added
        so it does not override exact species (e.g. HCl can exist in database)"""
        for key in list(self.chemical_to_elements.keys()):
            self.chemical_to_elements[self.safe_modifier_name(key)] = (
                self.chemical_to_elements.pop(key)
            )

    def safe_modifier_name(self, name):
        """ensures we use safe modifiers that do not replicate
        real species, e.g. exact species might contain HCl, but we might want
        to also add HCl to system, internally we want to track it as
        modifier_HCl it will be bound to provided pyomo var regardless"""
        if "modifier_" not in name:
            return f"modifier_{name}"
        else:
            return name

    def get_modifier_mw(self, elemental_composition):
        mw = 0
        for el, mol in elemental_composition.items():
            _mw, _unit = self.state.get_molar_mass_element(el)
            mw = mw + mol * _mw
        return mw, _unit

    def register_modifier(self, new_chemical):
        if new_chemical is not None:
            self.chemical_to_elements.update(new_chemical)
        self.ensure_safe_modifier_names()

    def write_element_sum_constraint(self, spec_object, element):
        """writes a sum of elements constraint for reaktoro"""
        # pull out all the input indexes  and their coefficients into a list
        # so we can write the constraints

        species_list = [
            (
                cv[0],
                self.rkt_inputs[cv[1]].get_rkt_index(),
            )
            for cv in self.constraint_dict[element]
        ]

        def _constraint_fn(props, w):
            """constraint function to sum up all species"""
            sum_species = []
            for mol, idx in species_list:
                sum_species.append(mol * w[idx])
            return props.elementAmount(element) - sum(sum_species)

        spec_object.openTo(element)
        constraint = rkt.EquationConstraint()
        constraint.id = f"{element}_constraint"
        constraint.fn = _constraint_fn
        spec_object.addConstraint(constraint)

    def write_elementAmount_constraint(self, spec_object, element, input_name=None):
        """writes a elements amount constraint for reaktoro"""
        spec_object.openTo(element)
        if input_name is None:
            idx = spec_object.addInput(f"input{element}")
        else:
            idx = spec_object.addInput(f"input{input_name}")
        constraint = rkt.EquationConstraint()
        constraint.id = f"{element}_constraint"
        constraint.fn = (
            lambda props, w: props.elementAmount(element) - w[idx]
        )  # - props.elementAmount(element)
        spec_object.addConstraint(constraint)

    def write_elementAmountInPhase_constraint(self, spec_object, element, phase):
        """writes a elements amount constraint for reaktoro"""
        spec_object.openTo(element)
        idx = spec_object.addInput(element)
        constraint = rkt.EquationConstraint()
        constraint.id = f"{element}_constraint"
        constraint.fn = (
            lambda props, w: props.elementAmountInPhase(element, phase) - w[idx]
        )

        spec_object.addConstraint(constraint)

    def write_speciesAmount_constraint(self, spec_object, species, input_name=None):
        """writes a elements amount constraint for reaktoro"""
        spec_object.openTo(species)
        if input_name is None:
            idx = spec_object.addInput(species)
        else:
            idx = spec_object.addInput(input_name)
        constraint = rkt.EquationConstraint()
        constraint.id = f"{species}_constraint"
        constraint.fn = lambda props, w: props.speciesAmount(species) - w[idx]
        spec_object.addConstraint(constraint)

    def write_pOH_constraint(self, spec_object):
        """writes a pOH constraint for reaktoro"""
        spec_object.openTo("OH-")
        idx = spec_object.addInput("pOH")
        constraint = rkt.EquationConstraint()
        constraint.id = f"pOH_constraint"
        constraint.fn = lambda props, w: w[idx] + props.speciesActivityLg("OH-")
        spec_object.addConstraint(constraint)

    def write_empty_con(self, spec_object, spc):
        constraint = rkt.EquationConstraint()
        constraint.id = f"{spc}_dummy_constraint"
        constraint.fn = lambda props, w: 0
        spec_object.addConstraint(constraint)

    def write_open_solvent_constraints(self, spec_object):
        """add redundant constraints for H2O"""
        opened_elements = []  # track elements we open, don't want duplicate calls
        for phase in self.fixed_solvent_speciation:
            for element, coeff in self.fixed_solvent_speciation[phase].items():
                if element not in opened_elements:
                    opened_elements.append(element)
                    spec_object.openTo(element)
                    self.write_empty_con(spec_object, element)

    def write_phase_volume_constraint(self, spec_object, phase, fixed_vlaue=1e-8):
        # spec_object.openTo(f"volume_{phase}")
        idx = spec_object.addInput(f"volume_{phase}")
        constraint = rkt.EquationConstraint()
        constraint.id = f"{phase}_volume_constraint"
        constraint.fn = lambda props, w: w[idx] - props.phaseProps(phase).volume()
        spec_object.addConstraint(constraint)

    def write_empty_constraints(self, spec_object):
        """add redundant constraints"""
        for specie in self.empty_constraints:
            spec_object.openTo(specie)
            self.write_empty_con(spec_object, specie)

    def export_config(self):
        export_object = ReaktoroInputExport()
        export_object.copy_chem_inputs(self.rkt_chemical_inputs)
        export_object.ignore_elements_for_constraints = (
            self.ignore_elements_for_constraints
        )
        export_object.fixed_solvent_specie = self.fixed_solvent_specie
        export_object.fixed_solvent_speciation = self.fixed_solvent_speciation
        export_object.assert_charge_neutrality = self.assert_charge_neutrality
        export_object.neutrality_ion = self.neutrality_ion
        export_object.dissolve_species_in_rkt = self.dissolve_species_in_rkt
        export_object.exact_speciation = self.exact_speciation
        export_object.chemical_to_elements = self.chemical_to_elements
        export_object.fixed_solvent_type = self.fixed_solvent_type
        export_object.empty_constraints = self.empty_constraints
        export_object.fixed_species = self.fixed_species
        return export_object

    def load_from_export_object(self, export_object):
        self.ignore_elements_for_constraints = (
            export_object.ignore_elements_for_constraints
        )
        self.fixed_solvent_specie = export_object.fixed_solvent_specie
        self.fixed_solvent_speciation = export_object.fixed_solvent_speciation
        self.fixed_solvent_type = export_object.fixed_solvent_type
        self.rkt_chemical_inputs = export_object.rkt_chemical_inputs
        self.assert_charge_neutrality = export_object.assert_charge_neutrality
        self.neutrality_ion = export_object.neutrality_ion
        self.dissolve_species_in_rkt = export_object.dissolve_species_in_rkt
        self.exact_speciation = export_object.exact_speciation
        self.chemical_to_elements = export_object.chemical_to_elements
        self.empty_constraints = export_object.empty_constraints
        self.fixed_species = export_object.fixed_species
