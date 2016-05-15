#! /usr/bin/env python
#
# Copyright (C) 2015-2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
## skchem.core.point

Defining points in scikit-chem.
"""

# the functionality here perhaps should be replaced by numpy arrays to get closer
# integration with the scientific python stack.

import rdkit.Geometry.rdGeometry
from .base import ChemicalObject

class Point3D(rdkit.Geometry.rdGeometry.Point3D, ChemicalObject):

    """ Class representing a point in scikit-chem """

    def to_dict(self, two_d=True):

        """ Dictionary representation of the point.

        Args:
            two_d (bool):
                Whether the point is in two dimensions or three.

        Returns:
            dict[str: float]: dictionary of coordinates to values.
        """

        if two_d:
            return {"x": round(self.x), "y": round(self.y)}
        else:
            return {"x": round(self.x), "y": round(self.y), "z": round(self.z)}

    def __repr__(self):
        return '<{klass} coords="({x:.2f}, {y:.2f}, {z:.2f})" at {address}>'.format(
            klass=self.__class__.__name__, \
            x=self.x, \
            y=self.y, \
            z=self.z, \
            address=hex(id(self)))

    def __str__(self):
        return '({x:.2f}, {y:.2f}, {z:.2f})'.format(x=self.x, y=self.y, z=self.z)
