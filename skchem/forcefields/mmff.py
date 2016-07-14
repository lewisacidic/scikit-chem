#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
## skchem.forcefields.mmff

Module specifying the Merck Molecular Force Field.
"""
from rdkit.Chem.rdForceFieldHelpers import MMFFOptimizeMolecule

from .base import ForceField


class MMFF(ForceField):

    def __init__(self, **kwargs):
        super(MMFF, self).__init__(**kwargs)

    def _optimize(self, mol):

        return MMFFOptimizeMolecule(mol)