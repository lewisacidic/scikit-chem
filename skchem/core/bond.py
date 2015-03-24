#! /usr/bin/env python
#
# Copyright (C) 2007-2009 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
skchem.core.bond

Defining chemical bonds in scikit-chem.
"""

import rdkit.Chem
from skchem.core import Atom
from skchem.core import ChemicalObject

class Bond(rdkit.Chem.rdchem.Bond, ChemicalObject):

    """
    Class representing a chemical bond in scikit-chem.

    """

    @property
    def order(self):

        """
        The order of the bond.

        Parameters
        ----------

        None

        Returns
        -------

        bond :int

        """

        return self.GetBondTypeAsDouble()

    @order.setter
    def order(self, value):

        """ Set the order of the bond.  Not implemented."""

        raise NotImplementedError

    @property
    def atoms(self):

        """ Return an iterable of the atoms involved in the bond. """

        return [Atom.from_super(self.GetBeginAtom()), Atom.from_super(self.GetEndAtom())]

    @atoms.setter
    def atoms(self, value):

        """ Set the atoms of the molecule.  Not implemented. """

        raise NotImplementedError

    def draw(self):

        """ Draw the bond inline. """

        return '{}{}{}'.format(self.atoms[0].element, \
            '-' if self.order == 1 else self.GetSmarts(), \
            self.atoms[0].element)

    def to_dict(self):

        """ Convert to a dictionary representation. """

        return {"b": self.GetBeginAtomIdx(), "e": self.GetEndAtomIdx(), "o": self.order}

    def __repr__(self):
        return '<{klass} type="{bond}" at {address}>'.format(klass=self.__class__.__name__, \
            bond=self.draw(), \
            address=hex(id(self)))

    def __str__(self):
        return self.draw()
