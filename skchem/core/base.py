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

    """ A mixin for each chemical object in scikit-chem. """

    @classmethod
    def from_super(cls, obj):

        """A method that converts the class of an object of parent class to that of the child. """

        obj.__class__ = cls
        return obj


class ChemicalObjectView(object):

    """ Abstract iterable view of chemical objects.

    Concrete classes inheriting from it should implement `__getitem__` and `__len__`.
    """

    __metaclass__ = ABCMeta

    def __init__(self, owner):
        """ Return a view """
        self.owner = owner

    @abstractmethod
    def __getitem__(self, index):
        # subclass call this, then implement the actual get if None
        if isinstance(index, slice):
            return self.to_list()[index]
        if isinstance(index, list) \
                and all(isinstance(i, bool) for i in index) \
                and len(index) == len(self):
            return [self[i] for i, ix in enumerate(index) if ix]
        elif isinstance(index, (list, tuple)):
            return [self[ix] for ix in index]
        else:
            return None

    @abstractmethod
    def __len__(self):
        pass

    def __iter__(self):
        return ChemicalObjectIterator(self)

    def __str__(self):
        return str(list(str(obj) for obj in self))

    @property
    def props(self):
        """ Return a property view of the objects in the view. """

        return MolPropertyView(self)

    def to_list(self):
        """ Return a list of objects in the view. """

        return list(self)

    def __repr__(self):
        return '<{class_} values="{values}" at {address}>'.format(
            class_=self.__class__.__name__,
            values=str(self),
            address=hex(id(self)))


class ChemicalObjectIterator(object):

    """  Iterator for chemical object views.  """

    def __init__(self, view):
        """ Create an iterator from a chemical object view. """
        self.view = view
        self._current = 0
        self._high = len(self.view)

    def __next__(self):
        if self._current >= self._high:
            raise StopIteration
        else:
            self._current += 1
            return self.view[self._current - 1]

    # py2 compat
    next = __next__


class View(object):

    """ View wrapper interface.  Conforms to the dictionary interface.

    Objects inheriting from this should implement the `keys`, `getitem`,
    `setitem` and `delitem` methods."""

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
        """ Remove all properties from the object. """
        for idx in self.keys():
            self.remove(idx)

    def items(self):
        """ Return an iterable of key, value pairs. """

        return list((k, self[k]) for k in self.keys())

    def remove(self, key):
        """ Remove a property from the object. """
        self.__delitem__(key)

    @abstractmethod
    def __getitem__(self, key):
        pass

    @abstractmethod
    def __setitem__(self, key, value):
        pass

    @abstractmethod
    def __delitem__(self, key):
        pass

    def __iter__(self):
        return iter(self.keys())

    def __str__(self):
        return str(self.to_dict())

    def __len__(self):
        return len(self.keys())

    def __repr__(self):
        return '<{klass} values="{values}" at {address}>'.format(
            klass=self.__class__.__name__,
            values=str(self),
            address=hex(id(self)))

    def to_dict(self):
        """ Return a dict of the properties on the object. """
        return dict(self)

    def to_series(self):
        """ Return a pd.Series of the properties on the object. """
        return pd.Series(self.to_dict())


class PropertyView(View):
    """ Property object wrapper.

     This provides properties for rdkit objects. """

    def __init__(self, owner):
        """ Initialize a PropertyView.

        Args:
            owner(skchem.ChemicalObject):
                A chemical object with properties, such as `skchem.Atom`."""

        self._owner = owner

    def keys(self):
        """ The available property keys on the object. """

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


class MolPropertyView(View):

    """ Mol property wrapper.

    This provides properties for the atom and bond views. """

    def __init__(self, obj_view):
        self._obj_view = obj_view

    def keys(self):
        """ The available property keys on the object. """

        res = set()
        for atom in self._obj_view:
            res = res.union(set(atom.props.keys()))
        return list(res)

    def get(self, key, default=None):
        return pd.Series((a.props.get(key, default) for a in self._obj_view), index=self._obj_view.index)

    def __getitem__(self, key):
        if key not in self.keys():
            raise KeyError('No atoms have the property set.')
        return self.get(key, None)

    def __setitem__(self, key, value):
        if isinstance(value, (pd.Series, dict)):
            for idx, val in pd.compat.iteritems(value):
                self._obj_view[int(idx)].props[key] = val
        else:
            assert len(self._obj_view) == len(value), "Must pass same number of values as atoms."
            for atom, val in zip(self._obj_view, value):
                atom.props[key] = val

    def __delitem__(self, key):
        for atom in self._obj_view:
            atom.props.remove(key)

    def to_frame(self):

        """ Return a DataFrame of the properties of the objects of the molecular view. """

        return pd.DataFrame(dict(self), index=self._obj_view.index)

    def to_dict(self):

        """ Return a dict of the properties of the objectos fo the molecular view. """

        return {k: v.tolist() for k, v in self.items()}

    def __str__(self):
        return str(self.to_dict())

