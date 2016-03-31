#! /usr/bin/env python
#
# Copyright (C) 2007-2009 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD


"""
skchem.core.atom

Defining atoms in scikit-chem.
"""

from rdkit import Chem
from skchem.core import ChemicalObject

class Atom(Chem.rdchem.Atom, ChemicalObject):

    """ Object representing an Atom in scikit-chem. """

    @property
    def element(self):

        """ Get the element of the atom as a string. """

        return self.GetSymbol()

    @property
    def atomic_number(self):

        """ Get the atomic number of the atom as a float. """

        return self.GetAtomicNum()

    @property
    def mass(self):

        """ Get the mass of the atom as a float. """

        return self.GetMass()

    @property
    def atomic_mass(self):

        """ Get the mass of the atom as a float. """

        return self.mass

    @property
    def props(self):

        """ Return a dictionary of properties of the atom. """

        # Some atom properties are inaccessible, but still give values.
        #

        props = {}

        for prop in self.GetPropNames():
            try:
                props[prop] = self.GetProp(prop)
            except RuntimeError:
                pass

        return props

    def __repr__(self):

        return '<{klass} element="{element}" at {address}>'.format(
            klass=self.__class__.__name__,
            element=self.element,
            address=hex(id(self))
            )

    def __str__(self):
        return self.element
