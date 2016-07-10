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

    def __init__(self):
        pass

    def optimize(self, mol):
        res = UFFOptimizeMolecule(mol)

        if res == -1:
            msg = 'Failed to optimize molecule \'{}\' using MMFF'.format(mol.name)
            if self.error_on_fail:
                raise RuntimeError(msg)
            elif self.warn_on_fail:
                warnings.warn(msg)
            else:
                pass

        return mol