#! /usr/bin/env python
#
# Copyright (C) 2007-2009 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
skchem.descriptors.physicochemical

Physicochemical descriptors and associated functions are defined.

"""

from rdkit.Chem import Descriptors
import pandas as pd

from .fingerprints import Fingerprinter
from ..utils import camel_to_snail

def molecular_weight(m):
    return Descriptors.MolWt(m)

DESCRIPTORS = [(camel_to_snail(s), f) for (s, f) in Descriptors.descList]

class PhysicochemicalFingerprinter(Fingerprinter):

    NAME = 'physchem'
    
    def __init__(self, descriptors='all'):
        if descriptors == 'all':
            self.descriptors = DESCRIPTORS
        else:
            self.descriptors = descriptors
        self.descriptor_names, _ = zip(*self.descriptors)

    def _transform(self, mol):

        return pd.Series({n: f(mol) for (n, f) in self.descriptors}, name=mol.name)
