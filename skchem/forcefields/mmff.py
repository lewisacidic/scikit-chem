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

    """ Merck Molecular Force Field transformer. """

    def __init__(self, preembed=True, warn_on_fail=True, error_on_fail=False,
                 add_hs=True, n_jobs=1, verbose=True):

        """ Initialize a MMFF object.

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
        super(MMFF, self).__init__(preembed=preembed,
                                   warn_on_fail=warn_on_fail,
                                   error_on_fail=error_on_fail, add_hs=add_hs,
                                   verbose=verbose, n_jobs=n_jobs)

    def _optimize(self, mol):

        return MMFFOptimizeMolecule(mol)