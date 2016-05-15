#! /usr/bin/env python
#
# Copyright (C) 2015-2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
## skchem.core.base

Define base classes for scikit chem objects
"""

import warnings

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
        self._current = 0
        self._high = self.owner.GetNumAtoms()

    def __getitem__(self, index):
        from .atom import Atom
        return Atom.from_super(self.owner.GetAtomWithIdx(index))

    def __iter__(self):
        return AtomIterator(self.owner)

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


class PropertyView(object):

    """ Property object wrapper """

    def __init__(self, owner):
        self.owner = owner

    def __getitem__(self, index):
        return self.owner.GetProp(str(index))

    def __setitem__(self, key, value):
        if not isinstance(value, str):
            warnings.warn("""RDKit property keys and values can only be
                        `str`.  Using `{value}` as a `str`.""".format(value=value))
        self.owner.SetProp(str(key), str(value))

    def __iter__(self):
        return iter(self.keys())

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

    def remove(self, index):
        self.owner.ClearProp(index)

    def clear(self):
        for idx in self.keys():
            self.remove(idx)

    def keys(self):
        return list(k for k in self.owner.GetPropNames() if k[:1] != '_')

    def items(self):
        return list((k, self[k]) for k in self.keys())

    def __str__(self):
        return str(dict(self))

    def __repr__(self):
        return '<{klass} values="{values}" at {address}>'.format(
            klass=self.__class__.__name__,
            values=str(self),
            address=hex(id(self)))
