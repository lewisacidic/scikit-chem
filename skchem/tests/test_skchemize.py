#! /usr/bin/env python
#
# Copyright (C) 2007-2009 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

""" Tests for scikit chemize """

from unittest import TestCase
from skchem.data import resource

from rdkit import Chem

from skchem.io import read_sdf
from skchem.core import Mol
from skchem.descriptors import skchemize
import pandas as pd

class TestSKChemize(TestCase):

    """ test class for skchemize """

    def test_modifies_rdkit_fps(self):

        '''Check to see if the fingerprint now gives a Series object'''

        m = Mol.from_smiles('C')
        func = skchemize(Chem.RDKFingerprint)
        self.assertTrue(isinstance(func(m), pd.Series))

    def test_fingerprint_is_same(self):

        '''Check if the fingerprint is the same'''

        m = Mol.from_smiles('c1ccccc1')
        f = skchemize(Chem.RDKFingerprint)
        gold = list(Chem.RDKFingerprint(m))
        self.assertTrue(f(m).tolist() == gold)

    def test_gives_dataframe(self):

        '''check to see if a dataframe is obtained'''

        f = skchemize(Chem.RDKFingerprint)
        df = read_sdf(resource('test_sdf', 'multi_molecule-simple.sdf'))
        self.assertTrue(isinstance(df.structure.apply(f), pd.DataFrame))
