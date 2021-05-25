# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.


from .mechanism import Mechanism
from .reaction import Reaction
from .species import Complex
from .component import Component
from .propensities import ProportionalHillPositive, GeneralPropensity



class BasicIntegration(Mechanism):
    """Mechanism for the schema DNA1 + DNA2 --> DNA3 + DNA4."""
    def __init__(self, name: str="basic_integration", mechanism_type: str="integration",**keywords):
        """Initializes a BasicIntegration instance.

        :param name: name of the Mechanism, default: basic_integration
        :param mechanism_type: type of the Mechanism, default: integration
        :param keywords:
        """
        # TODO ZAT: remove unused keywords argument
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
        #ku = component.get_parameter("ku",part_id=part_id,mechanism = self)

        #prop_hill = ProportionalHillPositive(k=kint, K=ku, s1=self.integrase, n=2, d=DNA_inputs[0])

        return [Reaction.from_massaction(inputs=[self.integrase]*4+DNA_inputs,\
                                        outputs=[self.integrase]*4+DNA_outputs, k_forward=kint)]

class RepressibleEnzymeIntegration(Mechanism):
    """Mechanism for the schema integrase+DNA1 + DNA2 --> integrase+DNA3 + DNA4. a 'repressor' will prevent this reaction from happening"""
    def __init__(self, name: str="repressible_enzyme_integration", mechanism_type: str="integration",integrase="Int1",repressor="rep1",**keywords):
        """Initializes a BasicIntegration instance.

        :param name: name of the Mechanism, default: basic_integration
        :param mechanism_type: type of the Mechanism, default: integration
        :param keywords:
        """
        #TODO currently this reaction will not work because reactions are generated from EITHER attB or attP. Thus if this reaction has a specific 
        #repressor in mind, the attB reaction cannot also be specific (since there is only one of those).
        #not sure how to solve this
        # TODO ZAT: remove unused keywords argument
        self.integrase = Component.set_species(integrase,material_type="protein")
        self.repressor = Component.set_species(repressor,material_type="protein")
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(self, DNA_inputs, DNA_outputs = None, **keywords):
        #this doesn't make any species because I use a Binding mechanism for that
        #maybe if we do the tetramerization mechanism then this would do something

        return [self.repressor]

    def update_reactions(self, DNA_inputs, DNA_outputs, component = None, part_id = None, kint = None, kb_rep = None,nr = None, **keywords):
        
        
        if part_id is None and component is not None:
            part_id = component.name

        if kint is None and component is None:
            raise ValueError("Must pass in either a component or kint.")
        elif kint is None:
            kint = component.get_parameter("kint", part_id = part_id, mechanism = self)
        if kb_rep is None:
            kb_rep = component.get_parameter("kb",part_id = self.repressor.name,mechanism=self)
        if nr is None:
            nr = component.get_parameter("cooperativity",part_id = self.repressor.name,mechanism=self)
        #ku = component.get_parameter("ku",part_id=part_id,mechanism = self)
        prop = GeneralPropensity(f"kint * {str(self.integrase)}^4 * ({' * '.join([str(a) for a in DNA_inputs])})* 1/(kb + {str(self.repressor)}^cooperativity)",propensity_species=[self.integrase,self.repressor]+DNA_inputs,propensity_parameters=[kb_rep,kint,nr])
        #prop_hill = ProportionalHillPositive(k=kint, K=ku, s1=self.integrase, n=2, d=DNA_inputs[0])

        return [Reaction(inputs=[self.integrase]*4+DNA_inputs+[self.repressor]*nr.value,\
                                        outputs=[self.integrase]*4+DNA_outputs+[self.repressor]*nr.value, propensity_type=prop)]