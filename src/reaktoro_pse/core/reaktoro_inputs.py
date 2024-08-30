###############################################################################
# #################################################################################
# # WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# # through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# # National Renewable Energy Laboratory, and National Energy Technology
# # Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# # of Energy). All rights reserved.
# #
# # Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# # information, respectively. These files are also available online at the URL
# # "https://github.com/watertap-org/reaktoro-pse/"
# #################################################################################
###############################################################################
import reaktoro as rkt
import reaktoro_pse.core.util_classes.rkt_inputs as RktInputs
from reaktoro_pse.core.reaktoro_state import ReaktoroState

import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)

__author__ = "Alexander Dudchenko"

""" class to setup input constraints, and specs for reaktoro solver class"""

# class RktChemModifier:
#     def __init__(self, elemental_makeup):
#         self.elemental_makeup=elemental_makeup

#     def get_mw(self):


class ReaktoroInputSpec:
    def __init__(self, reaktor_state):
        """initialize parameters needed to build reaktor solver"""
        self.state = reaktor_state
        if isinstance(self.state, ReaktoroState) == False:
            raise TypeError("Reator inputs require rektoroState class")
        self.user_inputs = reaktor_state.inputs  # user iputs provided to state
        self.rkt_inputs = RktInputs.RktInputs()  # inputs that will be fed to rkt spec
        self.rkt_chemical_inputs = RktInputs.RktInputs()

        """ execute default configuration options, user can update settings """
        self.register_charge_neutrality()
        self.default_speciation()
        self.register_aqueous_solvent()
        self.register_open_species()

    def register_chemistry_modifiers(self, chemical_dict, index=None):
        """registers chemistry modifiers being added to system
        chemistry_modifier -- chemicals to be added (pyo object should be mole flow of chemical that would enter a system with same species as in apparat_species_mol_flow
                        chemistry_modifier = {'HCl':m.fs.HCl_dose} example for HCl
        """
        for chemical, obj in chemical_dict.items():
            if index is None or index in chemical:
                if isinstance(chemical, tuple):
                    chemical = chemical[-1]
                self.register_chemistry_modifier(chemical, obj)

    def register_chemistry_modifier(self, chemical, pyomo_var):
        if chemical not in self.chemical_to_elements:
            raise ValueError(
                f"{chemical} is not avaialbe in chemical_to_element dict, please add"
            )
        self.rkt_chemical_inputs[chemical] = RktInputs.RktInput(
            var_name=chemical, pyomo_var=pyomo_var
        )

        mw, mw_unit = self.get_modifier_mw(self.chemical_to_elements[chemical])
        self.state.verify_unit(self.rkt_chemical_inputs[chemical], mw, mw_unit)

    def register_open_species(self, specie=None):
        """registers species to open to optimization and write empety constraint for,
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
                    f"Registered an empty constraint for {self.empty_constraints}, this can lead to unexpected results depending on reaktoro configuration, please use with caution"
                )

    def register_charge_neutrality(self, assert_neutrality=True, ion="Cl"):
        self.assert_charge_neutrality = assert_neutrality
        self.neutrality_ion = ion

    def register_aqueous_solvent(self, aqueous_solvent="H2O"):
        """defines aqueous species for system - used to set species when speciating - H/O change based on
        system specitation, so if we want to specify pH, we need to allow system to find eq. H/O and fix
        H2O"""
        self.aqueous_solvent = aqueous_solvent

        self.aqueous_solvent_speciation = {}

    def configure_specs(
        self,
        dissolve_species_in_rkt=True,
        exact_speciation=False,
    ):
        """configures specification for the problem

        Keyword arguments:
        dissolveSpeciesInRkt -- If true, species would be summed up to element amount in rkt, if false
        mode will contain conditions to build pyomo constraints via raktoroIO class
        exact_speciation -- if True, will write exact element amount for all input species other wise
        will leave  H, and O open, while fixing aqueousSolvent to specified value (e.g. H2O)

        """
        self.dissolve_species_in_rkt = dissolve_species_in_rkt
        self.exact_speciation = exact_speciation
        self.breakdown_species_to_elements()
        self.equilibrium_specs = rkt.EquilibriumSpecs(self.state.state.system())
        self.add_specs(
            self.equilibrium_specs,
            self.assert_charge_neutrality,
            dissolve_species_in_rkt,
        )

        """ get input name order!"""
        for idx, spec in enumerate(self.equilibrium_specs.namesInputs()):
            if spec == "T":
                spec_var_name = "temperature"
            elif spec == "P":
                spec_var_name = "pressure"
            elif "input" in spec:
                spec_var_name = spec.replace("input", "")
            else:
                spec_var_name = spec
            """ only care for indexes that exists and were added to spec"""

            if self.rkt_inputs.get(spec_var_name) is not None:
                self.rkt_inputs[spec_var_name].set_jacobian_index(idx)
                # tracking inputs we are passing into our spec problem
                if spec_var_name not in self.rkt_inputs.rkt_input_list:
                    self.rkt_inputs.rkt_input_list.append(spec_var_name)

    def breakdown_species_to_elements(self):
        """this will take all species in rktstate and create a dictionary containing
        thier elemnts and species amounts for use in summing
        eg. {'H2O':{'H':2,'O':1}}"""

        # TODO: probably want to make a class to track this
        self.specie_to_elements = {}
        for specie in self.state.state.system().species():
            self.specie_to_elements[specie.name()] = {}
            for i, el in enumerate(specie.elements().symbols()):
                self.specie_to_elements[specie.name()][
                    el
                ] = specie.elements().coefficients()[i]
        self.chemical_to_elements.update(self.specie_to_elements)

    def add_specs(self, specs_object, assert_charge_neutrality, dissolveSpeciesInRkt):
        # ignore elements for constraints
        self.ignore_elements_for_constraints = []
        # ignore elements for summation
        self.ignore_sum_elements = []
        pressure_not_set = True
        temperature_not_set = True
        for input_name, _ in self.state.inputs.items():
            if input_name == "temperature":
                specs_object.temperature()
                temperature_not_set = False
                self.rkt_inputs["temperature"] = self.state.inputs["temperature"]
                self.rkt_inputs["temperature"].set_lower_bound(0)
            elif input_name == "pressure":
                specs_object.pressure()
                pressure_not_set = False
                self.rkt_inputs["pressure"] = self.state.inputs["pressure"]
                self.rkt_inputs["pressure"].set_lower_bound(0)
            elif input_name == "pH":
                specs_object.pH()
                self.rkt_inputs["pH"] = self.state.inputs["pH"]
                self.rkt_inputs["pH"].set_lower_bound(0)
            else:
                pass
        if pressure_not_set:
            specs_object.unknownPressure()
            self.write_empty_con(specs_object, "open_pressure")

        if temperature_not_set:
            specs_object.unknownTemperature()
            self.write_empty_con(specs_object, "open_temperature")
        if assert_charge_neutrality:
            self.ignore_elements_for_constraints.append(self.neutrality_ion)
            specs_object.charge()
            if self.neutrality_ion not in specs_object.namesInputs():
                """needs to be a species!"""
                specs_object.openTo(self.neutrality_ion)

        self._find_element_sums()
        """ add/check if vars in rkt Inputs"""
        if dissolveSpeciesInRkt:
            self.write_active_species(specs_object)
        else:
            for element in self.constraint_dict:
                if element not in self.rkt_inputs:
                    self.rkt_inputs[element] = RktInputs.RktInput(element)
                    self.rkt_inputs[element].set_lower_bound(0)
        """ write reaktoro constraints to spec"""
        for element in self.constraint_dict:
            if dissolveSpeciesInRkt:
                self.write_element_sum_constraint(specs_object, element)
            else:
                self.write_elementAmount_constraint(specs_object, element)
        if self.exact_speciation == False:
            aq_specie = self.state.inputs.convert_rkt_species_fun(self.aqueous_solvent)
            self.write_speciesAmount_constraint(specs_object, aq_specie)
            if self.aqueous_solvent not in self.rkt_inputs:
                self.rkt_inputs[aq_specie] = self.state.inputs[self.aqueous_solvent]
                self.rkt_inputs[aq_specie].set_rkt_input_name(aq_specie)
                self.rkt_inputs[aq_specie].set_lower_bound(0)
            self.write_open_solvent_constraints(specs_object)
        self.write_empty_constraints(specs_object)
        """ legacy code """

    def _find_element_sums(self):
        """
        Here in we will take all input species, elements, and chemicals and organize them such that
        amount of element == sum of all species inputs, except those we excluded due to
        specifying pH or for charge neutrality or otherwise, we will also track their
        respective coefficients"""
        # TODO: Should add some sort of way to check if all species were accounted for!
        self.constraint_dict = {}
        self.active_species = []
        rktState = self.state.state
        if self.exact_speciation == False:
            # self.rktActiveSpecies.append(self.aqueousSolvent)
            aquous_phase_ions = self.specie_to_elements[
                self.state.inputs.convert_rkt_species_fun(self.aqueous_solvent)
            ]
            for element, coeff in aquous_phase_ions.items():
                self.ignore_elements_for_constraints.append(element)
                self.aqueous_solvent_speciation[element] = coeff
                _log.info(
                    f"Exact speciation is not provided! Fixing aqueous solvent and, excluding {element}"
                )
        self.rkt_elements = [specie.symbol() for specie in rktState.system().elements()]
        # loop over all elements in the rkt system
        for element in self.rkt_elements:
            # skip any we want to ignore
            if element not in self.ignore_elements_for_constraints:
                self.constraint_dict[element] = []
                # check if element is not in our inputs and not ignore list
                # if its not that means we need to find all species that are related to
                # that element
                if element not in list(self.state.inputs.keys()):
                    # loop over all input species
                    for specie in self.state.inputs.species_list:
                        # check if specie is in list of rkt species
                        spc_dict = self.specie_to_elements.get(specie)
                        # might be empty as specie might not exist, thats okay
                        if spc_dict is not None:
                            coef = spc_dict.get(element)
                            # checks if element is in the actual species, might not exists and thats okay
                            # example C might not be in "H2O"
                            if coef is not None:
                                self.constraint_dict[element].append((coef, specie))
                                if specie not in self.active_species:
                                    self.active_species.append(specie)
                # if element was in the list of inputs, this means user
                # provided excet elemental amounts (e.g. C instead of CO3-2)
                elif element in self.state.inputs:
                    self.constraint_dict[element].append((1, element))
                    if element not in self.active_species:
                        self.active_species.append(element)
                # now lets also check if user provided chemical inputs and
                # add them to our elemental sum constraints (e.g. H = H(from H2O) + H (from HCL))
                for specie in self.rkt_chemical_inputs.keys():
                    # if specie is a element add directly
                    if specie == element:
                        self.constraint_dict[element].append((1, specie))
                        if specie not in self.active_species:
                            self.active_species.append(specie)
                    # if not the nfind species and add their coefficients
                    elif specie in self.chemical_to_elements:
                        coef = self.chemical_to_elements[specie].get(element)
                        if coef is not None:
                            self.constraint_dict[element].append((coef, specie))
                            if specie not in self.active_species:
                                self.active_species.append(specie)
                # make sure we did not create empty element lists
                if len(self.constraint_dict[element]) == 0:
                    del self.constraint_dict[element]

    def write_active_species(self, spec_object):
        # build intputs into rkt model, and track thier indexes for writing rkt constraints
        for specie in self.active_species:
            input_name = f"input{specie}"
            idx = spec_object.addInput(input_name)
            if specie in self.state.inputs:
                self.rkt_inputs[specie] = self.state.inputs[specie]
                self.rkt_inputs[specie].set_rkt_index(idx)
                self.rkt_inputs[specie].set_rkt_input_name(input_name)
                self.rkt_inputs[specie].set_lower_bound(0)
            elif specie in self.rkt_chemical_inputs:
                self.rkt_inputs[specie] = self.rkt_chemical_inputs[specie]
                self.rkt_inputs[specie].set_rkt_index(idx)
                self.rkt_inputs[specie].set_rkt_input_name(input_name)
                self.rkt_inputs[specie].set_lower_bound(0)
            else:
                raise KeyError(f"Specie is not found {specie}")

    def default_speciation(self):
        # TODO: probably want to make a class to track this stuff
        """defines species to element conversions"""
        self.chemical_to_elements = {
            "HCl": {"H": 1, "Cl": 1},
            "CaO": {"Ca": 1, "O": 1},
            "Na2CO3": {"Na": 2, "C": 1, "O": 3},
            "CO2": {"C": 1, "O": 2},
            "NaOH": {"Na": 1, "O": 1, "H": 1},
            "H": {"H": 1},
            "OH": {"O": 1, "H": 1},
            "H2O_evaporation": {"O": 1, "H": 2},
        }

    def get_modifier_mw(self, elemental_composition):
        mw = 0
        for el, mol in elemental_composition.items():
            _mw, _unit = self.state.get_molar_mass_element(el)
            mw = mw + mol * _mw
        return mw, _unit

    def register_modifier(self, new_chemical):
        if new_chemical is not None:
            self.chemical_to_elements.update(new_chemical)

    def write_element_sum_constraint(self, spec_object, element):
        """writes a sum of elements constraint for reaktoro"""
        # pull out all the input indexes  and thier coefficents into a list
        # so we can write the constraints
        spec_object.openTo(element)
        species_list = [
            (cv[0], self.rkt_inputs[cv[1]].get_rkt_index())
            for cv in self.constraint_dict[element]
        ]
        constraint = rkt.EquationConstraint()
        constraint.id = f"{element}_constraint"
        constraint.fn = lambda props, w: sum(
            [mol * w[idx] for (mol, idx) in species_list]
        ) - props.elementAmount(element)
        spec_object.addConstraint(constraint)

    def write_elementAmount_constraint(self, spec_object, element):
        """writes a elements amount constraint for reaktoro"""
        spec_object.openTo(element)
        idx = spec_object.addInput(element)
        constraint = rkt.EquationConstraint()
        constraint.id = f"{element}_constraint"
        constraint.fn = lambda props, w: w[idx] - props.elementAmount(element)
        spec_object.addConstraint(constraint)

    def write_elementAmountInPhase_constraint(self, spec_object, element, phase):
        """writes a elements amount constraint for reaktoro"""
        spec_object.openTo(element)
        idx = spec_object.addInput(element)
        constraint = rkt.EquationConstraint()
        constraint.id = f"{element}_constraint"
        constraint.fn = lambda props, w: w[idx] - props.elementAmountInPhase(
            element, phase
        )
        spec_object.addConstraint(constraint)

    def write_speciesAmount_constraint(self, spec_object, species):
        """writes a elements amount constraint for reaktoro"""
        spec_object.openTo(species)
        idx = spec_object.addInput(species)
        constraint = rkt.EquationConstraint()
        constraint.id = f"{species}_constraint"
        constraint.fn = lambda props, w: w[idx] - props.speciesAmount(species)
        spec_object.addConstraint(constraint)

    def write_empty_con(self, spec_object, spc):
        constraint = rkt.EquationConstraint()
        constraint.id = f"{spc}_dummy_constraint"
        constraint.fn = lambda props, w: 0
        spec_object.addConstraint(constraint)

    def write_open_solvent_constraints(self, spec_object):
        """add redundant constraints for H2O"""
        for element, coeff in self.aqueous_solvent_speciation.items():
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

        # if "GaseousPhase" in [
        #     phase.name() for phase in self.state.state.system().phases()
        # ]:
        #     spec_object.openTo("H2O")
        #     spec_object.phaseVolume("GaseousPhase")
        # self.write_phase_volume_constraint(spec_object, "GaseousPhase")
        # _log.info("fixed volume")
        # spec_object.openTo("H2O(g)")
        # self.write_empty_con(spec_object, "H2O(g)")

    # def write_empty_constraints(
    #     self, rkt_chemical_inputs, spec_object, active_species=[]
    # ):
    #     ''' legacy code '''
    #     """this function will write in all species and open them to optimization of reaktoro"""
    #     def write_empty_con(spec_object, spc):
    #         constraint = rkt.EquationConstraint()
    #         constraint.id = f"{spc}_constraint"
    #         constraint.fn = lambda props, w: 0
    #         spec_object.addConstraint(constraint)

    # existing_constraints = spec_object.namesConstraints()
    # existing_variables = spec_object.namesControlVariables()
    # for chem in self.rktBase.rktState.system().species():
    #     if chem.name() not in str(spec_object.namesControlVariables()):
    #         spec_object.openTo(chem.name())
    #         write_empty_con(spec_object, chem.name())
    # if self.RktInputs.convert_to_rkt_species:
    #     aq_phase= RktInputs.specie_to_rkt_species(self.aqueousSolvent)
    # else:
    #     aq_phase=self.aqueousSolvent
    # if aq_phase not in str(spec_object.namesControlVariables()):
    #     spec_object.openTo(aq_phase)
    #     write_empty_con(spec_object, aq_phase)

    # print(spec_object.namesConstraints())
    # print(spec_object.namesControlVariables())
    # print(spec_object.namesControlVariablesP())
    # print(spec_object.namesControlVariablesQ())
    # print(spec_object.namesTitrants())
    # print(spec_object.namesTitrantsExplicit())
    # print(spec_object.namesTitrantsImplicit())
    # print(spec_object.namesInputs())
    # spec_object.openTo("H2O(g)")
    # write_empty_con(spec_object, "H2O(g)")
    # supported_species =[specie.name() for specie in self.rktBase.rktState.system().species()]
    # # print(self.rktBase.userInputs)
    # for specie in self.rktBase.RktInputs.species_list:
    #     print(specie)
    #     if specie not in str(spec_object.namesControlVariables()) and specie in supported_species:
    #         spec_object.openTo(specie)
    #         write_empty_con(spec_object, specie)
    #     existing_constraints = spec_object.namesConstraints()
    #     num_cons = len(existing_constraints)
    #     existing_variables = spec_object.namesInputs()
    #     num_vars = len(existing_variables)
    #     if num_cons < num_vars:
    #         for i in range(10):
    #             if f"dummy_{i}_constraint" not in existing_constraints:
    #                 write_empty_con(spec_object, f"dummy_{i}")
    #                 num_cons += 1
    #                 if num_cons == num_vars:
    #                     break#
