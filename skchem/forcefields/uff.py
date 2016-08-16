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

    """ Universal Force Field transformer. """

    def __init__(self, preembed=True, warn_on_fail=True,
                 error_on_fail=False, add_hs=True, verbose=True):

        """ Initialize a UFF object.

        Args:
            preembed (bool):
                Whether to embed before optimizing.
            warn_on_fail (bool):
                Whether to warn if a molecule fails to optimise.
            error_on_fail (bool):
                Whether to raise an error if a molecule fails to optimise.
            add_hs (bool):
                Whether to automatically add hydrogens.
        """

        super(UFF, self).__init__(preembed=preembed, warn_on_fail=warn_on_fail,
                                  error_on_fail=error_on_fail, add_hs=add_hs,
                                  verbose=verbose)

    def _optimize(self, mol):
        try:
            return UFFOptimizeMolecule(mol)
        except RuntimeError:
            return None