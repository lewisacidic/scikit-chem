#! /usr/bin/env python
#
# Copyright (C) 2007-2009 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

from unittest import TestCase
from skchem.tests.utils import resource
#required for running test
import pandas as pd
import skchem as skc

# required for test
# fix for python2 and 3 compatability
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO


single_molecule_props = {'PUBCHEM_IUPAC_INCHIKEY', 'PUBCHEM_COMPOUND_CANONICALIZED', 'PUBCHEM_IUPAC_INCHI', 'PUBCHEM_COMPOUND_CID', 'PUBCHEM_OPENEYE_ISO_SMILES', 'PUBCHEM_ATOM_UDEF_STEREO_COUNT', 'PUBCHEM_MOLECULAR_FORMULA', 'PUBCHEM_ISOTOPIC_ATOM_COUNT', 'PUBCHEM_CACTVS_COMPLEXITY', 'PUBCHEM_COORDINATE_TYPE', 'PUBCHEM_BOND_DEF_STEREO_COUNT', 'PUBCHEM_CACTVS_HBOND_DONOR', 'PUBCHEM_IUPAC_OPENEYE_NAME', 'PUBCHEM_EXACT_MASS', 'PUBCHEM_ATOM_DEF_STEREO_COUNT', 'PUBCHEM_IUPAC_TRADITIONAL_NAME', 'PUBCHEM_OPENEYE_CAN_SMILES', 'PUBCHEM_IUPAC_NAME', 'PUBCHEM_MOLECULAR_WEIGHT', 'PUBCHEM_CACTVS_TAUTO_COUNT', 'PUBCHEM_CACTVS_HBOND_ACCEPTOR', 'PUBCHEM_CACTVS_ROTATABLE_BOND', 'PUBCHEM_TOTAL_CHARGE', 'PUBCHEM_IUPAC_CAS_NAME', 'PUBCHEM_MONOISOTOPIC_WEIGHT', 'PUBCHEM_HEAVY_ATOM_COUNT', 'PUBCHEM_BOND_UDEF_STEREO_COUNT', 'PUBCHEM_CACTVS_SUBSKEYS', 'PUBCHEM_IUPAC_SYSTEMATIC_NAME', 'PUBCHEM_CACTVS_TPSA', 'PUBCHEM_XLOGP3_AA', 'PUBCHEM_COMPONENT_COUNT'}
single_molecule_name = '297'
single_molecule_num_atoms = 1
single_molecule_num_atoms_w_Hs = 5

multi_molecule_num_molecules = 3
multi_molecule_names = ['297', '6324', '6334']
multi_molecule_props = single_molecule_props.union({'DUMMY_PROPERTY_A', 'DUMMY_PROPERTY_B', 'DUMMY_PROPERTY_C'})

class TestSDF(TestCase):

    def test_opening_with_file(self):
        
        """ Can an sdf file be opened with a file-like object? """
        
        with open(resource('sdf', 'single_molecule-simple.sdf'), 'rb') as f:
            df = skc.read_sdf(f)
            self.assertTrue(len(df) == 1) 

    def test_with_file_correct_structure(self):
        
        """ When opened with a file-like object, is the structure correct?
        Done by checking atom number (should be one, as rdkit ignores Hs by default """
        
        with open(resource('sdf', 'single_molecule-simple.sdf'), 'rb') as f:
            df = skc.read_sdf(f)
            self.assertTrue(df.structure[single_molecule_name].GetNumAtoms() == 1)

    def test_opening_with_path(self):

        """ Do we find a molecule in example file? """

        df = skc.read_sdf(resource('sdf', 'single_molecule-simple.sdf'))   
        self.assertTrue(len(df) == 1)

    def test_with_path_correct_structure(self):

        """ When opened with a path, is the structure correct? """

        single_molecule_df = skc.read_sdf(resource('sdf', 'single_molecule-simple.sdf'))  
        self.assertTrue(single_molecule_df.structure[single_molecule_name].GetNumAtoms() == single_molecule_num_atoms)

    def test_arg_forwarding(self):

        """ Check that kwargs can still be parsed to the rdkit object """

        single_molecule_df = skc.read_sdf(resource('sdf', 'single_molecule-simple.sdf'), removeHs=False)  
        self.assertTrue(single_molecule_df.structure[single_molecule_name].GetNumAtoms() == single_molecule_num_atoms_w_Hs)

    def test_single_index_detected(self):

        """ Does molecule have a name set to index? """

        single_molecule_df = skc.read_sdf(resource('sdf', 'single_molecule-simple.sdf'))  
        self.assertFalse((single_molecule_df.index == pd.DataFrame(['dummy']).index).all())

    def test_single_index_correct(self):

        """ is name correct? """

        single_molecule_df = skc.read_sdf(resource('sdf', 'single_molecule-simple.sdf'))  
        self.assertTrue(single_molecule_df.index[0] == single_molecule_name) 

    def test_single_properties_detected(self):

        """ Does the dataframe have properties? """

        df = skc.read_sdf(resource('sdf', 'single_molecule-properties.sdf'))
        test = set(df.columns)
        test.remove('structure')
        self.assertTrue(len(test) > 1)

    def test_single_properties_correct(self):
        '''are they the right properties?'''

        single_molecule_df = skc.read_sdf(resource('sdf', 'single_molecule-properties.sdf'))  
        props = set(single_molecule_df.columns)
        props.remove('structure')
        self.assertTrue(props == single_molecule_props)

    def test_multi_parsed(self):
        '''do we find right number of molecules?'''

        multi_molecule_df = skc.read_sdf(resource('sdf', 'multi_molecule-simple.sdf'))
        self.assertTrue(multi_molecule_df.shape[0] == multi_molecule_num_molecules)

    def test_multi_index_detected(self):
        '''is index set?'''

        multi_molecule_df = skc.read_sdf(resource('sdf', 'multi_molecule-simple.sdf'))
        self.assertFalse((multi_molecule_df.index == pd.DataFrame(['dummy'] * multi_molecule_num_molecules).index).all())

    def test_multi_index_correct(self):
        '''is it the right index?'''

        multi_molecule_df = skc.read_sdf(resource('sdf', 'multi_molecule-simple.sdf'))
        self.assertTrue((multi_molecule_df.index == multi_molecule_names).all())

    def test_multi_properties_different_correct(self):
        '''if there are properties not common for all, are they all detected?'''

        multi_molecule_df = skc.read_sdf(resource('sdf', 'multi_molecule-properties.sdf'))
        props = set(multi_molecule_df.columns)
        props.remove('structure')
        self.assertTrue(props == multi_molecule_props)

    def test_bad_structure(self):
        """ Does it throw an error if bad structures are given? """

        with self.assertRaises(ValueError):
            multi_molecule_df = skc.read_sdf(resource('sdf', 'multi_molecule-bad_structure.sdf'))