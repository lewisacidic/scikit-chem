#! /usr/bin/env python
#
# Copyright (C) 2015-2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
## skchem.core.base

Define base classes for scikit chem objects
"""

from abc import ABCMeta, abstractmethod
import warnings
import numpy as np
import pandas as pd


class ChemicalObject(object):

    """ A mixin for each chemical object in scikit-chem """

    @classmethod
    def from_super(cls, obj):

        """A method that converts the class of an object of parent class to that of the child. """

        obj.__class__ = cls
        return obj


class AtomView(object):

    """ Atom interface wrapper """

    def __init__(self, owner):
        self.owner = owner
        self.props = AtomPropertyView(self)

    def __getitem__(self, index):
        from .atom import Atom
        return Atom.from_super(self.owner.GetAtomWithIdx(index))

    def __len__(self):
        return self.owner.GetNumAtoms()

    def __iter__(self):
        return AtomIterator(self.owner)

    def __str__(self):
        return str(list(str(atom) for atom in self))

    @property
    def elements(self):
        return pd.Series((atom.element for atom in self), index=self.index)

    @property
    def atomic_number(self):
        return pd.Series((atom.atomic_number for atom in self), index=self.index)

    @property
    def atomic_mass(self):
        return pd.Series((atom.mass for atom in self), index=self.index)

    @property
    def index(self):
        return pd.RangeIndex(len(self), name='atom_idx')

    def __repr__(self):
        return '<{class_} values="{values}" at {address}>'.format(
            class_=self.__class__.__name__,
            values=str(self),
            address=hex(id(self)))


class AtomIterator(AtomView):

    """ Atom iterator """

    def __init__(self, owner):
        super(AtomIterator, self).__init__(owner)
        self._current = 0
        self._high = self.owner.GetNumAtoms()

    def __next__(self):
        if self._current >= self._high:
            raise StopIteration
        else:
            self._current += 1
            return self[self._current - 1]

    # py2 compat
    next = __next__


class View(object):

    """ View wrapper interface """
    __metaclass__ = ABCMeta

    @abstractmethod
    def keys(self):
        return []

    def get(self, index, default=None):
        if index in self.keys():
            return self[index]
        else:
            return default

    def pop(self, index, default=None):
        if default:
            val = self.get(index, default)
        else:
            val = self[index]
        self.remove(index)
        return val

    def clear(self):
        for idx in self.keys():
            self.remove(idx)

    def items(self):
        return list((k, self[k]) for k in self.keys())

    def remove(self, key):
        self.__delitem__(key)

    def __getitem__(self, key):
        raise NotImplemented

    def __setitem__(self, key, value):
        raise NotImplemented

    def __delitem__(self, key):
        raise NotImplemented

    def __iter__(self):
        return iter(self.keys())

    def __str__(self):
        return str(dict(self))

    def __len__(self):
        return len(self.keys())

    def __repr__(self):
        return '<{klass} values="{values}" at {address}>'.format(
            klass=self.__class__.__name__,
            values=str(self),
            address=hex(id(self)))


class PropertyView(View):

    """ Property object wrapper """

    def __init__(self, owner):
        self._owner = owner

    def keys(self):
        return list(k for k in self._owner.GetPropNames() if k[:1] != '_')

    def __getitem__(self, key):

        # we manually work out if it was a float that was stored, as GetProp
        # returns floats and ints set by SetDoubleProp and SetIntProp as strings
        value = self._owner.GetProp(str(key))
        try:
            return int(value)
        except ValueError:
            try:
                return float(value)
            except ValueError:
                return value

    def __setitem__(self, key, value):

        if not isinstance(key, str):
            warnings.warn("RDKit property keys can only be of type `str`.  Using `{key}` as a `str`.".format(key=key))
            key = str(key)

        if key[0] == '_':
            warnings.warn("`{value}` is a private RDKit property key. "
                          "Using this may have unintended consequences.".format(value=value))

        if isinstance(value, str):
            self._owner.SetProp(key, value)
        elif isinstance(value, (int, np.int64, np.int32)):
            self._owner.SetIntProp(key, value)
        elif isinstance(value, (float, np.float64, np.float32)):
            self._owner.SetDoubleProp(key, value)
        else:
            warnings.warn("RDKit property keys can only be `str`, `int` or `float`."
                          "Using `{value}` as a `str`.".format(value=value))
            self._owner.SetProp(key, str(value))

    def __delitem__(self, index):
        self._owner.ClearProp(index)


class AtomPropertyView(View):

    """ Atom property wrapper """

    def __init__(self, atom_view):
        self._atom_view = atom_view

    def keys(self):
        res = set()
        for atom in self._atom_view:
            res = res.union(set(atom.props.keys()))
        return list(res)

    def get(self, key, default=None):
        return [a.props.get(key, default) for a in self._atom_view]

    def __getitem__(self, key):
        if key not in self.keys():
            raise KeyError('No atoms have the property set.')
        return self.get(key, None)

    def __setitem__(self, key, value):
        assert len(self._atom_view) == len(value), "Must pass same number of values as atoms."
        for atom, val in zip(self._atom_view, value):
            atom.props[key] = val

    def __delitem__(self, key):
        for atom in self._atom_view:
            atom.props.remove(key)
