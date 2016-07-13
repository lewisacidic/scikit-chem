#! /usr/bin/env python
#
# Copyright (C) 2015-2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

""" # skchem.pandas.structure_methods

 Tools for adding a default attribute to pandas objects."""

from pandas.core.base import NoNewAttributesMixin, AccessorProperty
from pandas.core.series import Series
from pandas.core.index import Index

from .. import core

class StructureMethods(NoNewAttributesMixin):
    def __init__(self, data):
        self._data = data

    def add_hs(self, **kwargs):
        return self._data.apply(lambda m: m.add_hs(**kwargs))

    def remove_hs(self, **kwargs):
        return self._data.apply(lambda m: m.remove_hs(**kwargs))

    @property
    def atoms(self):
        return self._data.apply(lambda m: m.atoms)

def only_contains_mols(ser):
    return ser.apply(lambda s: isinstance(s, core.Mol)).all()

class StructureAccessorMixin(object):

    def _make_structure_accessor(self):
        if isinstance(self, Index):
            raise AttributeError('Can only use .mol accessor with molecules,'
                                 'which use np.object_ in scikit-chem.')
        if not only_contains_mols(self):
            raise AttributeError('Can only use .mol accessor with '
                                 'Series that only contain mols.')

        return StructureMethods(self)
    ms = AccessorProperty(StructureMethods, _make_structure_accessor)


Series.__bases__ += StructureAccessorMixin,
