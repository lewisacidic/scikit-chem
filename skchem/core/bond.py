#! /usr/bin/env python
#
# Copyright (C) 2015-2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
## skchem.core.bond

Defining chemical bonds in scikit-chem.
"""

import pandas as pd
import rdkit.Chem

from . import Atom
from .base import ChemicalObject, PropertyView, ChemicalObjectView


class Bond(rdkit.Chem.rdchem.Bond, ChemicalObject):

    """
    Class representing a chemical bond in scikit-chem.

    """

    @property
    def props(self):

        """ PropertyView: rdkit properties of the atom. """

        if not hasattr(self, '_props'):
            self._props = PropertyView(self)
        return PropertyView(self)

    @property
    def order(self):

        """ int: the order of the bond. """

        return self.GetBondTypeAsDouble()

    @property
    def atoms(self):

        """ list[Atom]: list of atoms involved in the bond. """

        return [Atom.from_super(self.GetBeginAtom()),
                Atom.from_super(self.GetEndAtom())]

    def draw(self):

        """ str: Draw the bond in ascii. """

        return '{}{}{}'.format(self.atoms[0].element,
                               '-' if self.order == 1 else self.GetSmarts(),
                               self.atoms[1].element)

    def to_dict(self):

        """ dict: Convert to a dictionary representation. """

        return {"b": self.GetBeginAtomIdx(),
                "e": self.GetEndAtomIdx(),
                "o": self.order}

    def __repr__(self):
        return '<{klass} type="{bond}" at {address}>'.format(
            klass=self.__class__.__name__,
            bond=self.draw(),
            address=hex(id(self)))

    def __str__(self):
        return self.draw()


class BondView(ChemicalObjectView):

    """ Bond interface wrapper """
    def __getitem__(self, index):
        res = super(BondView, self).__getitem__(index)
        if res is None:
            if index >= len(self):
                raise IndexError('Index {} out of range for molecule with '
                                 '{} bonds.'.format(index, len(self)))
            else:
                return Bond.from_super(self.owner.GetBondWithIdx(index))
        else:
            return res

    def __len__(self):
        return self.owner.GetNumBonds()

    @property
    def order(self):
        """ A `pd.Series` of the bond orders of the view's bonds. """

        return pd.Series((bond.order for bond in self), index=self.index)

    @property
    def index(self):
        """ A `pd.Index` of the bonds in the `BondView`. """
        return pd.RangeIndex(len(self), name='bond_idx')
