#  Copyright (c) 2019, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.


from unittest import TestCase

#This file tests ComplexSpecies, OrderedComplexSpecies, and Multimers which are all subclasses of species.

class TestComplexSpecies(TestCase):

    def test_species_initialization(self):
        from biocrnpyler import Species
        from biocrnpyler import ComplexSpecies
        from biocrnpyler import OrderedComplexSpecies
        from biocrnpyler import Multimer
        
        
        s1 = Species(name='s1')
        s2 = Species(name = 's2', material_type = "m2")
        
        
        
        #Check invalidity of ComplexSpecies with fewer than 1 component
        with self.assertRaises(ValueError):
            c = ComplexSpecies([s1])
        
        #Check invalidity of OrderedComplexSpecies with fewer than 1 component
        with self.assertRaises(ValueError):
            oc = OrderedComplexSpecies([s1])
        
        #Check invalidity of multimers with fewer than 1 component
        with self.assertRaises(ValueError):
            m = Multimer(s1, 1)
        
        

        
        
        #Check the naming conventions
        oc1 = OrderedComplexSpecies([s2, s1])
        c1 = ComplexSpecies([s2, s1])
        m1 = Multimer(s1, 2)
        c3 = ComplexSpecies([s1, s1])

        #Check validity of ComplexSpecies, Multimers and OrderedComplexSpecies with strings instead of species
        self.assertEqual(OrderedComplexSpecies([s2, "s1"]), oc1)
        self.assertEqual(ComplexSpecies([s2, "s1"]), c1)
        self.assertEqual(Multimer("s1", 2), m1)
        
        #ComplexSpecies should sort the species added alphabetically by representation
        self.assertEqual(repr(c1), "complex_"+repr(s2)+"_"+repr(s1))
        
        #OrderedComplexSpecies do not sort their internal species
        self.assertEqual(repr(oc1), "ordered_complex_"+repr(s2)+"_"+repr(s1))
        
        #Multimers are just complexes with multiplicity
        self.assertEqual(repr(m1), "complex_2x_"+repr(s1))
        self.assertEqual(repr(c3), repr(m1))
        
        
        
        
    def test_species_equality(self):
        from biocrnpyler import Species
        from biocrnpyler import ComplexSpecies
        from biocrnpyler import OrderedComplexSpecies
        
        s1 = Species(name='s1', material_type = "m1")
        s2 = Species(name = 's2', material_type = "m2")
        
        c1 = ComplexSpecies([s1, s2])
        c2 = ComplexSpecies([s2, s1])
        #check equality of differently ordered complexes
        self.assertEqual(c1, c2)
                         
        oc1 = OrderedComplexSpecies([s2, s1])
        oc2 = OrderedComplexSpecies([s1, s2])
        self.assertFalse(oc1==oc2)


        