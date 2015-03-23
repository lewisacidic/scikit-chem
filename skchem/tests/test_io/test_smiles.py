#! /usr/bin/env python
#
# Copyright (C) 2007-2009 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD


#test libraries
import pytest
import os

#required for test
from skchem import read_smiles, Mol

from skchem.tests.utils import resource 

class TestSmiles(object):

    def test_single_parsed(self):

        """ Do we find a molecule in a single smiles file"""

        df = read_smiles(resource('smiles', 'single_molecule.smiles'))

        assert len(df) == 1
        assert df.iloc[0, 0].to_smiles() == 'C'

    def test_multiple_parsed(self):

        """ Do we find the exact number of molecules expected in a multi molecule smiles file?"""
        
        df = read_smiles(resource('smiles', 'multi_molecule-simple.smiles'))
        assert len(df) == 3

    def test_title_line(self):

        """ Test parsing a smiles file with a header. """
        df = read_smiles(resource('smiles', 'header_set.smiles'), title_line=True)
        assert len(df) == 3

    def test_header_correct(self):

        """ Is the header line correctly set? """
        df = read_smiles(resource('smiles', 'header_set.smiles'), title_line=True)
        assert list(df.columns) == ['structure', 'name']

    def test_configure_header(self):

        """ Can you pass header directly through to pandas? """

        df = read_smiles(resource('smiles', 'custom_header.smiles'), header=2)
        assert len(df) == 3
        assert (df.structure.apply(type) == Mol).all()

    def test_change_smiles_column(self):

        """ Does it work with smiles at different positions """
        df = read_smiles(resource('smiles', 'smiles_col_changed.smiles'), smiles_column=1)
        assert list(map(lambda m: m.to_smiles(), df.structure)) == ['C', 'CC', 'CCC']

    def test_name_column(self):
        
        """ Can it set the index? """

        df = read_smiles(resource('smiles', 'name_set.smiles'), name_column=1)
        assert list(df.index) == ['methane', 'ethane', 'propane']

    def test_properties(self):
    
        """ Can we read other properties? """

        multi_molecule_props = {'TEST_PROPERTY_A', 'TEST_PROPERTY_B', 'TEST_PROPERTY_C'}
        df = read_smiles(resource('smiles', 'multi_molecule-properties.smiles'), title_line=True)
        test = set(df.columns)
        test.remove('structure')
        assert test ==  multi_molecule_props

    def test_bad_smiles(self):

        """ Does it throw an error for an improper smiles code?"""

        with pytest.raises(ValueError):
            df = read_smiles(resource('smiles', 'multi_molecule-bad_smiles.smiles'), name_column=1, title_line=False, force=False)

    def test_bad_chemistry(self):

        """ Does it throw an error without force?"""

        with pytest.raises(ValueError):
            df = read_smiles(resource('smiles', 'multi_molecule-bad_chemistry.smiles'), name_column=1, title_line=False, force=False)

    def test_bad_chemistry_force(self):

        """ Can we force the parse? """

        df = read_smiles(resource('smiles', 'multi_molecule-bad_chemistry.smiles'), name_column=1, title_line=False, force=True)
        assert len(df) == 5


