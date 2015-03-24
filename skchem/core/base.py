#! /usr/bin/env python
#
# Copyright (C) 2007-2009 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
skchem.core.base

Define base classes for scikit chem objects
"""

class ChemicalObject(object):

    """ A mixin for each chemical object in scikit-chem """

    @classmethod
    def from_super(cls, obj):

        """A method that converts the class of an object of parent class to that of the child. """

        obj.__class__ = cls
        return obj
        