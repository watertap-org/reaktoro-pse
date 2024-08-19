import reaktoro as rkt
import reaktoro_pse.core.util_classes.rktInputs as rktInputs
from reaktoro_pse.core.reaktoro_state import reaktoroState
__author__ = "Alexander Dudchenko, Ben Knueven"

''' class to setup input constraints, and specs for reaktoro solver class'''

class reaktoroInputSpec:
    def __init__(self, reaktoroBase):
        """initialize parameters needed to build reaktor solver"""
        self.rktBase = reaktoroBase
        if isinstance(self.rktBase,reaktoroState) == False:
            raise TypeError('Reator inputs require rektoroState class')
        self.userInputs = reaktoroBase.rktInputs # user iputs provided to state
        self.rktInputs = rktInputs.rktInputs() # inputs that will be fed to rkt spec
        self.rktChemicalInputs=rktInputs.rktInputs()
        ''' execute default configuration options, user can update settings '''
        self.assert_charge_neutrality()
        self.default_chemical_speciation()



    def register_chemical_additions(self, chemical_dict, index=None):
        ''' registers chemicals being added to system
        chemical_addition -- chemicals to be added (pyo object should be mole flow of chemical that would enter a system with same species as in apparat_species_mol_flow
                        chemical_addition = {'HCl':m.fs.HCl_dose} example for HCl
        '''
        for chemical,obj in chemical_dict.items():
            if index is None or index in chemical:
                if isinstance(chemical, tuple):
                    chemical=chemical[-1]
                self.register_chemical_addition(chemical,obj)
                
    def register_chemical_addition(self, chemical, pyomo_var):
        self.rktChemicalInputs[chemical]=rktInputs.rktInput(var_name=chemical, pyomo_var=pyomo_var)

    def assert_charge_neutrality(self, assert_neutrality=True, ion="Cl"):
        self.assertChargeNeutrality = assert_neutrality
        self.neutralityIon = ion

    def configure_specs(self, dissolve_species_in_rkt=True):
        """configures specification for the problem

        Keyword arguments:
        dissolveSpeciesInRkt -- If true, species would be summed up to element amount in rkt, if false
        mode will contain conditions to build pyomo constraints via raktoroIO class"""
        self.dissolveSpeciesInRkt=dissolve_species_in_rkt
        self.breakdown_species_to_elements()
        self.rktEquilibriumSpecs = rkt.EquilibriumSpecs(self.rktBase.rktState.system())
        self.add_specs(self.rktEquilibriumSpecs, self.assertChargeNeutrality, dissolve_species_in_rkt)
        
        """ get input name order!"""
        for idx, spec in enumerate(self.rktEquilibriumSpecs.namesInputs()):
            if spec == "T":
                spec_var_name='temperature'
            elif spec == "P":
                spec_var_name='pressure'
            elif 'input' in spec:
                spec_var_name=spec.replace('input','')
            else:
                spec_var_name=spec
            ''' only care for indexes that exists and were added to spec'''

            if self.rktInputs.get(spec_var_name) is not None:
                self.rktInputs[spec_var_name].set_jacobian_index(idx)  
                # tracking inputs we are passing into our spec problem
                self.rktInputs.rktInputList.append(spec_var_name)

    def breakdown_species_to_elements(self):
        """this will take all species in rktstate and create a dictionary containing
        thier elemnts and species amounts for use in summing
        eg. {'H2O':{'H':2,'O':1}}"""

        #TODO: probably want to make a class to track this
        self.specieToElements = {}
        for specie in self.rktBase.rktState.system().species():
            self.specieToElements[specie.name()] = {}
            for i, el in enumerate(specie.elements().symbols()):
                self.specieToElements[specie.name()][el] = specie.elements().coefficients()[
                    i
                ]            
    def add_specs(self, specs_object, assert_charge_neutrality, dissolveSpeciesInRkt):
        # ignore elements for constraints 
        self.ignoreElementsForConstraints = []
        # ignore elements for summation 
        self.ignoreSumElements = []

        for input_name, _ in self.rktBase.rktInputs.items():
            if input_name == "temperature":
                specs_object.temperature()
                self.rktInputs["temperature"] = self.rktBase.rktInputs["temperature"]
            elif input_name == "pressure":
                specs_object.pressure()
                self.rktInputs["pressure"] = self.rktBase.rktInputs["pressure"]
            elif input_name == "pH":
                self.ignoreSumElements.append('H')
                specs_object.pH()
                self.rktInputs["pH"] = self.rktBase.rktInputs["pH"]
            else:
                pass

        if assert_charge_neutrality:
            self.ignoreElementsForConstraints.append(self.neutralityIon)
            specs_object.charge()
            if self.neutralityIon not in specs_object.namesInputs():
                ''' needs to be a species! '''
                specs_object.openTo(
                    rktInputs.specie_to_rkt_species(self.neutralityIon)
                )

        self._find_element_sums()
        ''' add/check if vars in rkt Inputs'''
        if dissolveSpeciesInRkt:            
            self.write_active_species(specs_object)
        else:
            for element in self.constraintDict:
                if element not in self.rktInputs:
                    self.rktInputs[element]=rktInputs.rktInput(element)
        ''' write reaktoro constraints to spec'''
        for element in self.constraintDict:
            if dissolveSpeciesInRkt:
                self.write_element_sum_constraint(
                    specs_object, element
                )
            else:
                self.write_elementAmount_constraint(specs_object, element)
        
        # ensure that all species in system are open, and have an empty constraint for 
        # to satisfy 0DOF 
        self.write_empty_constraints(
            self.rktChemicalInputs,
            specs_object,
            self.rktActiveSpecies,
        )


    def _find_element_sums(self):
        """
        Here in we will take all input species, elements, and chemicals and organize them such that 
        amount of element == sum of all species inputs, except those we excluded due to 
        specifying pH or for charge neutrality or otherwise, we will also track their 
        respective coefficients"""
        # TODO: Should add some sort of way to check if all species were accounted for!
        self.constraintDict = {}
        self.rktActiveSpecies = []
        rktState=self.rktBase.rktState
        self.rkt_elements=        [
            specie.symbol() for specie in rktState.system().elements()
        ]
        # loop over all elements in the rkt system
        for element in self.rkt_elements:
            # skip any we want to ignore
            if element not in self.ignoreElementsForConstraints:
                self.constraintDict[element] = []
                # check if element is not in our inputs and not ignore list 
                # if its not that means we need to find all species that are related to 
                # that element
                if element not in self.ignoreSumElements and element not in list(
                    self.rktBase.rktInputs.keys()
                ):
                    # loop over all input species
                    for specie in self.rktBase.rktInputs.speciesList:
                        # check if specie is in list of rkt species
                        spc_dict = self.specieToElements.get(specie)                       
                        # might be empty as specie might not exist, thats okay
                        if spc_dict is not None: 
                            coef = spc_dict.get(element)
                            # checks if element is in the actual species, might not exists and thats okay
                            # example C might not be in "H2O"
                            if coef is not None:
                                self.constraintDict[element].append((coef, specie))
                                if specie not in self.rktActiveSpecies:
                                    self.rktActiveSpecies.append(specie)
                # if element was in the list of inputs, this means user
                # provided excet elemental amounts (e.g. C instead of CO3-2)
                elif element in self.rktBase.rktInputs:
                    self.constraintDict[element].append((1, element))
                    if element not in self.rktActiveSpecies:
                        self.rktActiveSpecies.append(element)
                # now lets also check if user provided chemical inputs and 
                # add them to our elemental sum constraints (e.g. H = H(from H2O) + H (from HCL))
                for specie in self.rktChemicalInputs.keys():
                    # if specie is a element add directly
                    if specie == element:
                        self.constraintDict[element].append((1, specie))
                        if specie not in self.rktActiveSpecies:
                            self.rktActiveSpecies.append(specie)
                    # if not the nfind species and add their coefficients
                    elif specie in self.chemical_to_elements:
                        coef = self.chemical_to_elements[specie].get(element)
                        if coef is not None:
                            self.constraintDict[element].append((coef, specie))
                            if specie not in self.rktActiveSpecies:
                                self.rktActiveSpecies.append(specie)
                # make sure we did not create empty element lists
                if len(self.constraintDict[element]) == 0:
                    del self.constraintDict[element]

    def write_active_species(self, spec_object):
        # build intputs into rkt model, and track thier indexes for writing rkt constraints   
        for specie in self.rktActiveSpecies:
            input_name=f'input{specie}'
            idx = spec_object.addInput(input_name)
            if specie in self.rktBase.rktInputs:
                self.rktInputs[specie] = self.rktBase.rktInputs[specie]
                self.rktInputs[specie].set_rkt_index(idx)
                self.rktInputs[specie].set_rkt_input_name(input_name)
            elif specie in self.rktChemicalInputs:
                self.rktInputs[specie] = self.rktChemicalInputs[specie]
                self.rktInputs[specie].set_rkt_index(idx)
                self.rktInputs[specie].set_rkt_input_name(input_name)
            else:
                raise KeyError(f"Specie is not found {specie}")
        
    def default_chemical_speciation(self):
        #TODO: probably want to make a class to track this stuff
        """defines species to element conversions"""
        self.chemical_to_elements = {
            "HCl": {"H": 1, "Cl": 1},
            "CaO": {"Ca": 1, "O": 1},
            "Na2CO3": {"Na": 2, "C": 1, "O": 3},
            "CO2": {"C": 1, "O": 2},
            "NaOH": {"Na": 1, "O": 1, "H": 1},
            "H": {"H": 1},
            "OH": {"O": 1, "H": 1},
        }
    def register_chemical(self,new_chemical):
        if new_chemical is not None:
            self.chemical_to_elements.update(new_chemical)


    def write_element_sum_constraint(
        self, spec_object, element
    ):  
        ''' writes a sum of elements constraint for reaktoro'''
        # pull out all the input indexes  and thier coefficents into a list 
        # so we can write the constraints
        spec_object.openTo(element)
        species_list = [(cv[0], self.rktInputs[cv[1]].get_rkt_index()) for cv in self.constraintDict[element]]
        constraint = rkt.EquationConstraint()
        constraint.id = f"{element}_constraint"
        constraint.fn = lambda props, w: sum(
            [mol * w[idx] for (mol, idx) in species_list]
        ) - props.elementAmount(element)
        spec_object.addConstraint(constraint)

    def write_elementAmount_constraint(self, spec_object, element):
        ''' writes a elements amount constraint for reaktoro'''
        spec_object.openTo(element)
        idx = spec_object.addInput(element)
        constraint = rkt.EquationConstraint()
        constraint.id = f"{element}_constraint"
        constraint.fn = lambda props, w: (w[idx] - props.elementAmount(element))
        spec_object.addConstraint(constraint)

    def write_empty_constraints(
        self, rkt_chemical_inputs, spec_object, active_species=[]
    ):
        """this function will write in all species and open them to optimization of reaktoro"""
        def write_empty_con(spec_object, spc):
            constraint = rkt.EquationConstraint()
            constraint.id = f"{spc}_constraint"
            constraint.fn = lambda props, w: 0
            spec_object.addConstraint(constraint)

        # existing_constraints = spec_object.namesConstraints()
        # existing_variables = spec_object.namesControlVariables()
        for chem in self.rktBase.rktState.system().species():
            if chem.name() not in str(spec_object.namesControlVariables()):
                spec_object.openTo(chem.name())
                write_empty_con(spec_object, chem.name())
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
