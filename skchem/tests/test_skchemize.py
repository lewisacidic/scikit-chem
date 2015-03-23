#! /usr/bin/env python
#
# Copyright (C) 2007-2009 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

from unittest import TestCase
from skchem.tests.utils import resource

from rdkit import Chem

from skchem.io import read_sdf
from skchem.core import Mol
from skchem.descriptors import skchemize
import pandas as pd

class TestSKChemize(TestCase):
    def test_modifies_rdkit_fps(self):

        '''Check to see if the fingerprint now gives a Series object'''

        m = Mol.from_smiles('C')
        f = skchemize(Chem.RDKFingerprint)
        self.assertTrue(type(f(m)) == pd.Series)

    def test_fingerprint_is_same(self):

        '''Check if the fingerprint is the same'''

        m = Mol.from_smiles('c1ccccc1')
        f = skchemize(Chem.RDKFingerprint)
        gold = list(Chem.RDKFingerprint(m))
        self.assertTrue(f(m).tolist() == gold)
        
    def test_gives_dataframe(self):

        '''check to see if a dataframe is obtained'''

        f = skchemize(Chem.RDKFingerprint)
        df = read_sdf(resource('sdf', 'multi_molecule-simple.sdf'))
        self.assertTrue(type(df.structure.apply(f)) == pd.DataFrame)
