# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.


from .mechanism import Mechanism
from .reaction import Reaction
from .species import Complex
from .component import Component
from .propensities import ProportionalHillPositive
from .mechanisms_binding import One_Step_Cooperative_Binding



class BasicIntegration(Mechanism):
    """Mechanism for the schema DNA1 + DNA2 --> DNA3 + DNA4."""
    def __init__(self, name: str="basic_integration", mechanism_type: str="integration",**keywords):
        """Initializes a BasicIntegration instance.

        :param name: name of the Mechanism, default: basic_integration
        :param mechanism_type: type of the Mechanism, default: integration
        :param keywords:
        """
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(self, DNA_inputs, DNA_outputs = None, **keywords):
        #this doesn't make any species because I use a Binding mechanism for that
        #maybe if we do the tetramerization mechanism then this would do something
        return []

    def update_reactions(self, DNA_inputs, DNA_outputs, component = None, part_id = None, kint = None, **keywords):
        if part_id is None and component is not None:
            part_id = component.name

        if kint is None and component is None:
            raise ValueError("Must pass in either a component or kint.")
        elif kint is None:
            kint = component.get_parameter("kint", part_id = part_id, mechanism = self)

        return [Reaction.from_massaction(inputs=DNA_inputs, outputs=DNA_outputs, k_forward=kint)]

class EnzymeIntegration(Mechanism):
    """Mechanism for the schema integrase+DNA1 + DNA2 --> integrase+DNA3 + DNA4."""
    def __init__(self, name: str="enzyme_integration", mechanism_type: str="integration",integrase="Int1",**keywords):
        """Initializes a BasicIntegration instance.

        :param name: name of the Mechanism, default: basic_integration
        :param mechanism_type: type of the Mechanism, default: integration
        :param keywords:
        """
        # TODO ZAT: remove unused keywords argument
        self.integrase = Component.set_species(integrase,material_type="protein")
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(self, DNA_inputs, DNA_outputs = None, **keywords):
        #this doesn't make any species because I use a Binding mechanism for that
        #maybe if we do the tetramerization mechanism then this would do something
        return []

    def update_reactions(self, DNA_inputs, DNA_outputs, component = None, part_id = None, kint = None, **keywords):
        
        
        if part_id is None and component is not None:
            part_id = component.name

        if kint is None and component is None:
            raise ValueError("Must pass in either a component or kint.")
        
        elif kint is None:
            kint = component.get_parameter("kint", part_id = part_id, mechanism = self)

        return [Reaction.from_massaction(inputs=[self.integrase]*4+DNA_inputs,\
                                        outputs=[self.integrase]*4+DNA_outputs, k_forward=kint)]


class EnzymeIntegrationBinding(One_Step_Cooperative_Binding,Mechanism):
    """Mechanism for the schema integrase+DNA1 + DNA2 --> integrase+DNA3 + DNA4."""
    def __init__(self, name: str="enzyme_integration", mechanism_type: str="integration",integrase="Int1",**keywords):
        """Initializes a BasicIntegration instance.

        :param name: name of the Mechanism, default: basic_integration
        :param mechanism_type: type of the Mechanism, default: integration
        :param keywords:
        """
        # TODO ZAT: remove unused keywords argument
        self.integrase = Component.set_species(integrase,material_type="protein")
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(self, binder, bindee, complex_species = None, cooperativity=None, component = None, part_id = None, **kwords):
        #this doesn't make any species because I use a Binding mechanism for that
        #maybe if we do the tetramerization mechanism then this would do something
        out_spec = One_Step_Cooperative_Binding.update_species(self,binder, bindee,\
                            complex_species = complex_species, cooperativity=cooperativity, \
                            component = component, part_id = part_id, **kwords)
        return out_spec

    def update_reactions(self, site1,site2, integrase_functions, DNA_outputs=None, \
                                        component = None, part_id = None, kint = None, **keywords):
        
        if part_id is None and component is not None:
            part_id = component.name

        if kint is None and component is None:
            raise ValueError("Must pass in either a component or kint.")
        
        elif kint is None:
            kint = component.get_parameter("kint", part_id = part_id, mechanism = self)
        
        out_rxns = One_Step_Cooperative_Binding.update_reactions(self,self.integrase,site1,\
                                                    component=component,part_id = self.integrase)
        if(self.integrase not in site2):
            #this means integrase is not bound to the other dna!
            return out_rxns
        
        int_cooperativity = component.get_parameter("cooperativity", part_id = part_id, mechanism = self, return_numerical = True)
        complexed_dna = Complex([site1]+[self.integrase]*int_cooperativity)
        complexed_parent = complexed_dna.parent


        


        if site1.parent == site1.parent:
            #this means this is an intramolecular reaction
            DNA_inputs = [complexed_parent]
            
        else:
            DNA_inputs = [complexed_parent,site2.parent]
        DNA_outputs = []
        for integrase_function in integrase_functions:
            DNA_outputs += [integrase_function.create_polymer(DNA_inputs)]
            
        return [Reaction.from_massaction(inputs=[self.integrase]*4+DNA_inputs,\
                                        outputs=[self.integrase]*4+DNA_outputs, k_forward=kint)]