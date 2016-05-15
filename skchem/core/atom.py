#! /usr/bin/env python
#
# Copyright (C) 2015-2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD


"""
## skchem.core.atom

Defining atoms in scikit-chem.
"""

from rdkit import Chem
from .base import ChemicalObject, PropertyView

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
