#! /usr/bin/env python
#
# Copyright (C) 2016 Rich Lewis <rl403@cam.ac.uk>
# License: 3-clause BSD

"""
## skchem.forcefields.base

Module specifying base class for forcefields.
"""
import warnings
from abc import ABCMeta, abstractmethod

import pandas as pd
from rdkit.Chem.rdDistGeom import EmbedMolecule

from ..utils import Suppressor
from ..base import Transformer
from ..filters.base import TransformFilter


class ForceField(Transformer, TransformFilter):
    # TODO: Multiple conformer generation handling.

    """ Base forcefield class.

    Filter drops those that fail to be optimized.

     """

    __metaclass__ = ABCMeta

    def __init__(self, preembed=True, warn_on_fail=True,
                 error_on_fail=False, add_hs=True, verbose=True):

        self.add_hs = add_hs
        self.warn_on_fail = warn_on_fail
        self.error_on_fail = error_on_fail
        self.preembed = preembed
        super(ForceField, self).__init__(verbose=verbose)

    @property
    def columns(self):
        return pd.Index(['structure'])

    def embed(self, mol):

        success = EmbedMolecule(mol)
        if success == -1:
            msg = 'Failed to Embed Molecule {}'.format(mol.name)
            if self.error_on_fail:
                raise RuntimeError(msg)
            elif self.warn_on_fail:
                warnings.warn(msg)
            return None

        if self.add_hs:
            return mol.add_hs(add_coords=True)
        else:
            return mol

    def _transform_mol(self, mol):

        mol = mol.copy()

        with Suppressor():
            if self.preembed:
                mol = self.embed(mol)

            if mol is None:  # embedding failed
                return None

            res = self._optimize(mol)

        if res == -1:
            msg = 'Failed to optimize molecule \'{}\' using {}'.format(
                mol.name, self.__class__)
            if self.error_on_fail:
                raise RuntimeError(msg)
            elif self.warn_on_fail:
                warnings.warn(msg)
            return None

        return mol

    @abstractmethod
    def _optimize(self, mol):
        pass


class RoughEmbedding(ForceField):
    def _optimize(self, mol):
        return mol
