from unittest import TestCase

#required for running test
import pandas as pd
import skchem as skc

# required for test

# need to fix for python2 and 3 compatability
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
    def test_single_parsed(self):
        '''do we find a molecule in example file?'''

        single_molecule_sdf = open('skchem/tests/test_resources/methane.sdf')
        single_molecule_df = skc.read_sdf(single_molecule_sdf)   
        self.assertTrue(single_molecule_df.shape[0] == 1)

    def test_single_index_detected(self):
        '''does molecule have a name set to index?'''

        single_molecule_sdf = open('skchem/tests/test_resources/methane.sdf')
        single_molecule_df = skc.read_sdf(single_molecule_sdf)  
        self.assertFalse((single_molecule_df.index == pd.DataFrame(['dummy']).index).all())

    def test_single_index_correct(self):
        '''is name correct?'''

        single_molecule_sdf = open('skchem/tests/test_resources/methane.sdf')
        single_molecule_df = skc.read_sdf(single_molecule_sdf)  
        self.assertTrue(single_molecule_df.index[0] == single_molecule_name) 


    def test_single_properties_detected(self):
        '''does dataframe have properties?'''

        single_molecule_sdf = open('skchem/tests/test_resources/methane.sdf')
        single_molecule_df = skc.read_sdf(single_molecule_sdf)  
        self.assertTrue(single_molecule_df.columns.shape[0] > 1)

    def test_single_properties_correct(self):
        '''are they the right properties?'''

        single_molecule_sdf = open('skchem/tests/test_resources/methane.sdf')
        single_molecule_df = skc.read_sdf(single_molecule_sdf)  
        props = set(single_molecule_df.columns)
        props.remove('structure')
        self.assertTrue(props == single_molecule_props)

    def test_single_correct_structure(self):
        '''check the structure correct. Done by checking atom number (should be one, as rdkit ignores Hs by default'''

        single_molecule_sdf = open('skchem/tests/test_resources/methane.sdf')
        single_molecule_df = skc.read_sdf(single_molecule_sdf)  
        self.assertTrue(single_molecule_df.structure[single_molecule_name].GetNumAtoms() == single_molecule_num_atoms)

    def test_arg_forwarding(self):
        '''check that kwargs can still be parsed to the rdkit object'''

        single_molecule_sdf = open('skchem/tests/test_resources/methane.sdf')
        single_molecule_df = skc.read_sdf(single_molecule_sdf, removeHs=False)  
        self.assertTrue(single_molecule_df.structure[single_molecule_name].GetNumAtoms() == single_molecule_num_atoms_w_Hs)

    def test_multi_parsed(self):
        '''do we find right number of molecules?'''

        multi_molecule_sdf = open('skchem/tests/test_resources/hydrocarbons.sdf')
        multi_molecule_df = skc.read_sdf(multi_molecule_sdf)
        self.assertTrue(multi_molecule_df.shape[0] == multi_molecule_num_molecules)

    def test_multi_index_detected(self):
        '''is index set?'''

        multi_molecule_sdf = open('skchem/tests/test_resources/hydrocarbons.sdf')
        multi_molecule_df = skc.read_sdf(multi_molecule_sdf)
        self.assertFalse((multi_molecule_df.index == pd.DataFrame(['dummy'] * multi_molecule_num_molecules).index).all())

    def test_multi_index_correct(self):
        '''is it the right index?'''

        multi_molecule_sdf = open('skchem/tests/test_resources/hydrocarbons.sdf')
        multi_molecule_df = skc.read_sdf(multi_molecule_sdf)
        self.assertTrue((multi_molecule_df.index == multi_molecule_names).all())

    def test_multi_properties_different_correct(self):
        '''if there are properties not common for all, are they all detected?'''

        multi_molecule_sdf = open('skchem/tests/test_resources/hydrocarbons.sdf')
        multi_molecule_df = skc.read_sdf(multi_molecule_sdf)
        props = set(multi_molecule_df.columns)
        props.remove('structure')
        self.assertTrue(props == multi_molecule_props)

