#! /usr/bin/env python
#
# Copyright (C) 2015 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
## skchem.forcefields.mmff

Module specifying the Merck Molecular Force Field.
"""
import warnings

from rdkit.Chem.rdForceFieldHelpers import MMFFOptimizeMolecule
from .base import ForceField


class MMFF(ForceField):
    def optimize(self, mol):
        res = MMFFOptimizeMolecule(mol)

        if res == -1:
            msg = 'Failed to optimize molecule \'{}\' using MMFF'.format(mol.name)
            if self.error_on_fail:
                raise RuntimeError(msg)
            elif self.warn_on_fail:
                warnings.warn(msg)
            else:
                pass

        return mol