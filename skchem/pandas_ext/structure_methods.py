#! /usr/bin/env python
#
# Copyright (C) 2015-2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

""" # skchem.pandas.structure_methods

 Tools for adding a default attribute to pandas objects."""

from sklearn.manifold import TSNE, MDS
from sklearn.decomposition import PCA

from matplotlib import pyplot as plt
import pandas as pd

from pandas.core.base import NoNewAttributesMixin, AccessorProperty
from pandas.core.series import Series
from pandas.core.index import Index

from .. import core
from .. import descriptors

DIM_RED = {
    'tsne': TSNE,
    'pca': PCA,
    'mds': MDS
}


class StructureMethods(NoNewAttributesMixin):

    """ Accessor for calling chemical methods on series of molecules. """

    def __init__(self, data):
        self._data = data

    def add_hs(self, **kwargs):
        return self._data.apply(lambda m: m.add_hs(**kwargs))

    def remove_hs(self, **kwargs):
        return self._data.apply(lambda m: m.remove_hs(**kwargs))

    def visualize(self, fper='morgan', dim_red='tsne', dim_red_kw={}, **kwargs):

        if isinstance(dim_red, str):
            dim_red = DIM_RED.get(dim_red.lower())(**dim_red_kw)

        fper = descriptors.get(fper)
        fper.verbose = False
        feats = fper.transform(self._data)
        feats = feats.fillna(feats.mean())
        twod = pd.DataFrame(dim_red.fit_transform(feats))
        return twod.plot.scatter(x=0, y=1, **kwargs)

    @property
    def atoms(self):
        return self._data.apply(lambda m: m.atoms)


def only_contains_mols(ser):
    return ser.apply(lambda s: isinstance(s, core.Mol)).all()


class StructureAccessorMixin(object):

    """ Mixin to bind chemical methods to objects. """

    def _make_structure_accessor(self):
        if isinstance(self, Index):
            raise AttributeError('Can only use .mol accessor with molecules,'
                                 'which use np.object_ in scikit-chem.')
        if not only_contains_mols(self):
            raise AttributeError('Can only use .mol accessor with '
                                 'Series that only contain mols.')

        return StructureMethods(self)
    mol = AccessorProperty(StructureMethods, _make_structure_accessor)

Series.__bases__ += StructureAccessorMixin,
