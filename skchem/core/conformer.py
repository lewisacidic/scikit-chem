#! /usr/bin/env python
#
# Copyright (C) 2015-2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
## skchem.core.conformer

Defining conformers in scikit-chem.
"""

import rdkit.Chem
from .point import Point3D
from .base import ChemicalObject

class Conformer(rdkit.Chem.rdchem.Conformer, ChemicalObject):

    """ Class representing a Conformer in scikit-chem. """
    #should use a view, list for now
    @property
    def atom_positions(self):

        """ Return the atom positions in the conformer for the atoms in the molecule. """

        return [Point3D.from_super(self.GetAtomPosition(i)) for i in range(self.GetNumAtoms())]

    @atom_positions.setter
    def atom_positions(self, value):

        """ Set the atom positions in the conformer.  Not implemented. """

        raise NotImplementedError

    @property
    def is_three_d(self):

        """ Return whether the conformer is three dimensional. """

        return self.is3D()

    @is_three_d.setter
    def is_three_d(self, value):

        """ Set whether the conformer is three dimensional. """

        self.set3D(value)

    def __repr__(self):
        return '<{klass} id="{id}" at {address}>'.format(klass=self.__class__.__name__, \
            id=self.GetId(), address=hex(id(self)))
