#! /usr/bin/env python
#
# Copyright (C) 2015-2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
## skchem.core.bond

Defining chemical bonds in scikit-chem.
"""

import rdkit.Chem
from . import Atom
from .base import ChemicalObject

class Bond(rdkit.Chem.rdchem.Bond, ChemicalObject):

    """
    Class representing a chemical bond in scikit-chem.

    """

    @property
    def order(self):

        """ int: the order of the bond. """

        return self.GetBondTypeAsDouble()

    @property
    def atoms(self):

        """ list[Atom]: list of atoms involved in the bond. """

        return [Atom.from_super(self.GetBeginAtom()), Atom.from_super(self.GetEndAtom())]

    def draw(self):

        """ str: Draw the bond in ascii. """

        return '{}{}{}'.format(self.atoms[0].element, \
            '-' if self.order == 1 else self.GetSmarts(), \
            self.atoms[1].element)

    def to_dict(self):

        """ dict: Convert to a dictionary representation. """

        return {"b": self.GetBeginAtomIdx(), "e": self.GetEndAtomIdx(), "o": self.order}

    def __repr__(self):
        return '<{klass} type="{bond}" at {address}>'.format(klass=self.__class__.__name__, \
            bond=self.draw(), \
            address=hex(id(self)))

    def __str__(self):
        return self.draw()
