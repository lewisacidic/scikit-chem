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
import numpy as np

from .atom import Atom
from .base import ChemicalObject, PropertyView, ChemicalObjectView


class Bond(rdkit.Chem.rdchem.Bond, ChemicalObject):

    """
    Class representing a chemical bond in scikit-chem.

    """

    @property
    def atoms(self):

        """ tuple[Atom]: list of atoms involved in the bond. """

        return (Atom.from_super(self.GetBeginAtom()),
                Atom.from_super(self.GetEndAtom()))

    @property
    def atom_idxs(self):

        """ tuple[int]: list of atom indexes involved in the bond. """

        return (self.GetBeginAtomIdx(), self.GetEndAtomIdx())

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
    def is_aromatic(self):

        """ bool: whether the bond is aromatic. """

        return self.GetIsAromatic()

    @property
    def is_conjugated(self):

        """ bool: whether the bond is conjugated. """

        return self.GetIsConjugated()

    @property
    def owner(self):

        """ skchem.Mol: the molecule this bond is a part of. """

        from .mol import Mol
        return Mol.from_super(self.GetOwningMol())

    @property
    def stereo_symbol(self):

        """ str: the stereo label of the bond ('Z', 'E', 'ANY', 'NONE') """

        return self.GetStereo().name.lstrip('STEREO')

    @property
    def is_in_ring(self):

        """ bool: whether the bond is in a ring. """

        return self.IsInRing()

    def draw(self):

        """ str: Draw the bond in ascii. """

        return '{}{}{}'.format(self.atoms[0].symbol,
                               '-' if self.order == 1 else self.GetSmarts(),
                               self.atoms[1].symbol)

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
            if abs(index) >= len(self):
                raise IndexError('Index {} out of range for molecule with '
                                 '{} bonds.'.format(index, len(self)))

            # index is negative, so adding gives desired indexing from back
            if index < 0:
                index += len(self)

            return Bond.from_super(self.owner.GetBondWithIdx(index))

        else:
            return res

    def __len__(self):
        return self.owner.GetNumBonds()

    @property
    def atom_idxs(self):

        """ The atom indices for the bonds in the view. """

        return np.array([atom.atom_idxs for atom in self])

    @property
    def order(self):

        """ np.array<int> the bond orders of the bonds in the view. """

        return np.array([bond.order for bond in self])

    @property
    def is_aromatic(self):

        """ np.array<bool> whether each of the bonds in the view are
        aromatic. """

        return np.array([bond.is_aromatic for bond in self])

    @property
    def is_conjugated(self):

        """ np.array<bool> whether each of the bonds in the view are c
        onjugated. """

        return np.array([bond.is_conjugated for bond in self])

    @property
    def is_in_ring(self):

        """ np.array<bool> whether each of the bonds in the view are in a
        ring. """

        return np.array([bond.is_in_ring for bond in self])

    @property
    def stereo_symbol(self):

        """ np.array<str> the stereo symbol of the bonds in the view. """

        return np.array([bond.stereo_symbol for bond in self])

    @property
    def index(self):

        """ A `pd.Index` of the bonds in the `BondView`. """

        return pd.RangeIndex(len(self), name='bond_idx')

__all__ = ['Atom', 'AtomView']
