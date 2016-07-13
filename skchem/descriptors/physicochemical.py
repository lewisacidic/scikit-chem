#! /usr/bin/env python
#
# Copyright (C) 2015-2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
## skchem.descriptors.physicochemical

Physicochemical descriptors and associated functions are defined.

"""

from rdkit.Chem import Descriptors
import pandas as pd
import numpy as np

from .fingerprints import Fingerprinter
from ..utils import camel_to_snail

DESCRIPTORS = [(camel_to_snail(s), f) for (s, f) in Descriptors.descList]

class PhysicochemicalFingerprinter(Fingerprinter):

    """ Physicochemical descriptor generator using RDKit descriptor """

    NAME = 'physchem'
    sparse = False

    def __init__(self, descriptors='all'):

        """ Create a physicochemical descriptor generator.

        Args:
            descriptors (list<(str, func)> or 'all'):
                Descriptors to calculate, or if 'all', use all descriptors."""

        if descriptors == 'all':
            self.descriptors = DESCRIPTORS
        else:
            self.descriptors = descriptors
        self.descriptor_names, _ = zip(*self.descriptors)

    @property
    def index(self):
        return self.descriptor_names

    def _transform(self, mol):
        res = []
        for (n, f) in self.descriptors:
            try:
                res.append(f(mol))
            except ValueError:
                print(mol)
                return res.append(np.NaN)

        return np.array(res)
