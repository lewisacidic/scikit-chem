#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
## skchem.forcefields.uff

Module specifying the universal force field.
"""
from rdkit.Chem.rdForceFieldHelpers import UFFOptimizeMolecule

from .base import ForceField

class UFF(ForceField):

    def __init__(self, **kwargs):
        super(UFF, self).__init__(**kwargs)

    def _optimize(self, mol):
        try:
            return UFFOptimizeMolecule(mol)
        except RuntimeError:
            return None