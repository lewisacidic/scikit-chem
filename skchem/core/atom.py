#! /usr/bin/env python
#
# Copyright (C) 2015-2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD


"""
## skchem.core.atom

Defining atoms in scikit-chem.
"""

import pandas as pd

from rdkit import Chem

from .base import ChemicalObject, PropertyView, ChemicalObjectView


class Atom(Chem.rdchem.Atom, ChemicalObject):

    """ Object representing an Atom in scikit-chem. """

    @property
    def element(self):

        """ str: the element symbol of the atom. """

        return self.GetSymbol()

    @property
    def atomic_number(self):

        """ int: the atomic number of the atom. """

        return self.GetAtomicNum()

    @property
    def mass(self):

        """ float: the mass of the atom.

        Usually relative atomic mass unless explicitly set. """

        return self.GetMass()

    @property
    def props(self):

        """ PropertyView: rdkit properties of the atom. """

        if not hasattr(self, '_props'):
            self._props = PropertyView(self)
        return PropertyView(self)

    def __repr__(self):

        return '<{klass} element="{element}" at {address}>'.format(
            klass=self.__class__.__name__,
            element=self.element,
            address=hex(id(self))
            )

    def __str__(self):

        return self.element


class AtomView(ChemicalObjectView):

    def __getitem__(self, index):
        res = super(AtomView, self).__getitem__(index)
        if res is None:
            if index >= len(self):
                msg = 'Index {} out of range for molecule with' \
                    '{} atoms.'.format(index, len(self))
                raise IndexError(msg)
            else:
                return Atom.from_super(self.owner.GetAtomWithIdx(index))
        else:
            return res

    def __len__(self):
        return self.owner.GetNumAtoms()

    @property
    def element(self):
        """ A `pd.Series` of the element of the atoms in `AtomView`. """

        return pd.Series((atom.element for atom in self), index=self.index)

    @property
    def atomic_number(self):
        """ A `pd.Series` of the atomic number of the atoms in `AtomView`. """

        return pd.Series((atom.atomic_number for atom in self),
                         index=self.index)

    @property
    def atomic_mass(self):
        """ A `pd.Series` of the atomic mass of the atoms in `AtomView`. """

        return pd.Series((atom.mass for atom in self), index=self.index)

    @property
    def index(self):
        """ A `pd.Index` of the atoms in the `AtomView`. """

        return pd.RangeIndex(len(self), name='atom_idx')
