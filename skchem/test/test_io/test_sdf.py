#! /usr/bin/env python
#
# Copyright (C) 2015-2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

""" Tests for sdf io functionality """

import pandas as pd
import pytest

from ...resource import resource
from ...io import read_sdf

SINGLE_MOLECULE_PROPS = {
    'PUBCHEM_IUPAC_INCHIKEY', 'PUBCHEM_COMPOUND_CANONICALIZED',
    'PUBCHEM_IUPAC_INCHI', 'PUBCHEM_COMPOUND_CID', 'PUBCHEM_OPENEYE_ISO_SMILES',
    'PUBCHEM_ATOM_UDEF_STEREO_COUNT', 'PUBCHEM_MOLECULAR_FORMULA',
    'PUBCHEM_ISOTOPIC_ATOM_COUNT', 'PUBCHEM_CACTVS_COMPLEXITY',
    'PUBCHEM_COORDINATE_TYPE', 'PUBCHEM_BOND_DEF_STEREO_COUNT',
    'PUBCHEM_CACTVS_HBOND_DONOR', 'PUBCHEM_IUPAC_OPENEYE_NAME',
    'PUBCHEM_EXACT_MASS', 'PUBCHEM_ATOM_DEF_STEREO_COUNT',
    'PUBCHEM_IUPAC_TRADITIONAL_NAME', 'PUBCHEM_OPENEYE_CAN_SMILES',
    'PUBCHEM_IUPAC_NAME', 'PUBCHEM_MOLECULAR_WEIGHT',
    'PUBCHEM_CACTVS_TAUTO_COUNT', 'PUBCHEM_CACTVS_HBOND_ACCEPTOR',
    'PUBCHEM_CACTVS_ROTATABLE_BOND', 'PUBCHEM_TOTAL_CHARGE',
    'PUBCHEM_IUPAC_CAS_NAME', 'PUBCHEM_MONOISOTOPIC_WEIGHT',
    'PUBCHEM_HEAVY_ATOM_COUNT', 'PUBCHEM_BOND_UDEF_STEREO_COUNT',
    'PUBCHEM_CACTVS_SUBSKEYS', 'PUBCHEM_IUPAC_SYSTEMATIC_NAME',
    'PUBCHEM_CACTVS_TPSA', 'PUBCHEM_XLOGP3_AA', 'PUBCHEM_COMPONENT_COUNT'
}

NON_SHARED_PROPS = {'DUMMY_PROPERTY_A', 'DUMMY_PROPERTY_B', 'DUMMY_PROPERTY_C'}
SINGLE_MOLECULE_NAME = '297'
SINGLE_MOLECULE_NUM_ATOMS = 1
SINGLE_MOLECULE_NUM_ATOMS_W_HS = 5

MULTI_MOLECULE_NUM_MOLECULES = 3
MULTI_MOLECULE_NAMES = ['297', '6324', '6334']
MULTI_MOLECULE_PROPS = SINGLE_MOLECULE_PROPS.union(NON_SHARED_PROPS)

class TestSDF(object):

    """ Test class for sdf file parser """

    def test_opening_with_file(self):

        """ Can an sdf file be opened with a file-like object? """

        with open(resource('test_sdf', 'single_molecule-simple.sdf'), 'rb') as f:
            df = read_sdf(f)
            assert len(df) == 1

    def test_file_correct_structure(self):

        """ When opened with a file-like object, is the structure correct?
        Done by checking atom number (should be one, as rdkit ignores Hs by default """

        with open(resource('test_sdf', 'single_molecule-simple.sdf'), 'rb') as f:
            df = read_sdf(f)
            assert df[SINGLE_MOLECULE_NAME].GetNumAtoms() == 1

    def test_opening_with_path(self):

        """ Do we find a molecule in example file? """

        df = read_sdf(resource('test_sdf', 'single_molecule-simple.sdf'))
        assert len(df) == 1

    def test_path_correct_structure(self):

        """ When opened with a path, is the structure correct? """

        df = read_sdf(resource('test_sdf', 'single_molecule-simple.sdf'))
        assert df[SINGLE_MOLECULE_NAME].GetNumAtoms() \
            == SINGLE_MOLECULE_NUM_ATOMS

    def test_arg_forwarding(self):

        """ Check that kwargs can still be parsed to the rdkit object """

        df = read_sdf(resource('test_sdf', 'single_molecule-simple.sdf'), removeHs=False)
        assert df[SINGLE_MOLECULE_NAME].GetNumAtoms() \
            == SINGLE_MOLECULE_NUM_ATOMS_W_HS

    def test_single_index_detected(self):

        """ Does molecule have a name set to index? """

        df = read_sdf(resource('test_sdf', 'single_molecule-simple.sdf'))
        assert not (df.index == pd.DataFrame(['dummy']).index).all()

    def test_single_index_correct(self):

        """ is name correct? """

        single_molecule_df = read_sdf(resource('test_sdf', 'single_molecule-simple.sdf'))
        assert single_molecule_df.index[0] == SINGLE_MOLECULE_NAME

    def test_single_properties_detected(self):

        """ Does the dataframe have properties? """

        df = read_sdf(resource('test_sdf', 'single_molecule-properties.sdf'))
        test = set(df.columns)
        test.remove('structure')
        assert len(test) > 1

    def test_single_properties_correct(self):

        """ Are they the right properties? """

        df = read_sdf(resource('test_sdf', 'single_molecule-properties.sdf'))
        props = set(df.columns)
        props.remove('structure')
        assert props == SINGLE_MOLECULE_PROPS

    def test_multi_parsed(self):

        """ Do we find right number of molecules?"""

        df = read_sdf(resource('test_sdf', 'multi_molecule-simple.sdf'))
        assert df.shape[0] == MULTI_MOLECULE_NUM_MOLECULES

    def test_multi_index_detected(self):

        """ Is index set? """

        df = read_sdf(resource('test_sdf', 'multi_molecule-simple.sdf'))
        dummy_df = pd.DataFrame(['dummy'] * MULTI_MOLECULE_NUM_MOLECULES)
        assert not (df.index == dummy_df.index).all()

    def test_multi_index_correct(self):
        '''is it the right index?'''

        df = read_sdf(resource('test_sdf', 'multi_molecule-simple.sdf'))
        assert (df.index == MULTI_MOLECULE_NAMES).all()

    def test_multi_diff_properties(self):
        '''if there are properties not common for all, are they all detected?'''

        df = read_sdf(resource('test_sdf', 'multi_molecule-properties.sdf'))
        props = set(df.columns)
        props.remove('structure')
        assert props == MULTI_MOLECULE_PROPS

    def test_bad_structure(self):
        """ Does it throw an error if bad structures are given? """

        with pytest.raises(ValueError):
            read_sdf(resource('test_sdf', 'multi_molecule-bad_structure.sdf'), error_bad_mol=True)
