#! /usr/bin/env python
#
# Copyright (C) 2015 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
## skchem.forcefields.uff

Module specifying the universal force field.
"""
import warnings

from .base import ForceField
from rdkit.Chem.rdForceFieldHelpers import UFFOptimizeMolecule


class UFF(ForceField):

    def __init__(self, **kwargs):
        super(UFF, self).__init__(**kwargs)

    def _optimize(self, mol):
        return UFFOptimizeMolecule(mol)
