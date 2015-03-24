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

    @element.setter
    def element(self, value):

        """ Set the element of the atom.  Not implemented. """

        raise NotImplementedError

    @property
    def props(self):

        """ Return a dictionary of properties of the atom. """

        return {i: self.GetProp() for i in self.GetProps()}

    @props.setter
    def props(self, value):

        """ Set the properties of the Atom.  Not implemented. """

        raise NotImplementedError

    def __repr__(self):

        return '<{klass} element="{element}" at {address}>'.format(
            klass=self.__class__.__name__,
            element=self.element,
            address=hex(id(self))
            )

    def __str__(self):
        return self.element
