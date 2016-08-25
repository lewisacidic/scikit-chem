#! /usr/bin/env python
#
# Copyright (C) 2015-Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
## skchem.descriptors.physicochemical

Physicochemical descriptors and associated functions are defined.

"""

from rdkit.Chem import Descriptors
import pandas as pd
import numpy as np

from ..base import Transformer, Featurizer
from ..utils import camel_to_snail

DESCRIPTORS = {camel_to_snail(s): f for (s, f) in Descriptors.descList}


class PhysicochemicalFeaturizer(Transformer, Featurizer):

    """ Physicochemical descriptor generator using RDKit descriptor """

    def __init__(self, features='all', **kwargs):

        """ Create a physicochemical descriptor generator.

        Args:
            descriptors (list<(str, func)> or 'all'):
                Descriptors to calculate, or if 'all', use all descriptors."""

        super(PhysicochemicalFeaturizer, self).__init__(**kwargs)

        self.features = features

    @property
    def features(self):
        return self._features

    @features.setter
    def features(self, features):
        if features == 'all':
            features = DESCRIPTORS
        elif isinstance(features, str):
            features = {features: DESCRIPTORS[features]}
        elif isinstance(features, list):
            features = {feature: DESCRIPTORS[feature] for feature in features}
        elif isinstance(features, (dict, pd.Series)):
            features = features
        else:
            raise NotImplementedError('Cannot use features {}'.format(features))

        self._features = pd.Series(features)
        self._features.index.name = 'physicochemical_features'

    @property
    def name(self):
        return 'physchem'

    @property
    def columns(self):
        return self.features.index

    def _transform_mol(self, mol):
        res = []
        for (n, f) in self.features.iteritems():
            try:
                res.append(f(mol))
            except ValueError:
                return res.append(np.NaN)

        return np.array(res)
