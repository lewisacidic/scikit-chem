from unittest import TestCase
from rdkit import Chem
import skchem as skc
import pandas as pd

class TestSKChemize(TestCase):
    def test_modifies_rdkit_fps(self):

        '''Check to see if the fingerprint now gives a Series object'''

        m = Chem.MolFromSmiles('C')
        f = skc.skchemize(Chem.RDKFingerprint)
        self.assertTrue(type(f(m)) == pd.Series)

    def test_fingerprint_is_same(self):

        '''Check if the fingerprint is the same'''

        m = Chem.MolFromSmiles('c1ccccc1')
        f = skc.skchemize(Chem.RDKFingerprint)
        gold = list(Chem.RDKFingerprint(m))
        self.assertTrue(f(m).tolist() == gold)
        
    def test_gives_dataframe(self):

        '''check to see if a dataframe is obtained'''

        f = skc.skchemize(Chem.RDKFingerprint)
        df = skc.read_sdf('skchem/tests/test_resources/hydrocarbons.sdf')
        self.assertTrue(type(df.structure.apply(f)) == pd.DataFrame)
